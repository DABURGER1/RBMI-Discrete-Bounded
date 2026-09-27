# BtheB reference-based multiple imputation with discrete Kumaraswamy margins.
# Gaussian copula with continuous-time AR(1) correlation at months 2, 3, 5, 8.
# The median anchor is continuous; its ceiling is the discrete lower median.
# Run from the application directory.
# Required packages: HSAUR2, posterior, and DHARMa.
# Installation: install.packages(c("HSAUR2", "posterior", "DHARMa"))
# Sourcing this file runs imputation, posterior summaries, and predictive checks.
# Treatment analysis and Rubin pooling are performed in a separate program.
# Exported LaTeX tables require booktabs, tabularx, siunitx, and bm.

DefaultBthebDkSettings <- function() {
  list(
    OutputRoot = "BtheB_DK_RBMI",
    ScoreMaximum = 63,
    Tau = 0.5,
    VisitVariables = c("bdi.2m", "bdi.3m", "bdi.5m", "bdi.8m"),
    VisitMonths = c(2, 3, 5, 8),
    ExcludeNoPostBaseline = TRUE,
    BaselineScale = 10,
    # The imputation model includes a separate baseline coefficient at each visit.
    Rules = c(MAR = 0, Half = 0.5, J2R = 1),
    NChains = 4,
    Warmup = 5000,
    Sample = 10000,
    Seed = 20260921,
    UseParallel = TRUE,
    CoreFraction = 0.8,
    AdaptEvery = 100,
    PitSweeps = 2,
    PrintEvery = 1000,
    PrintProgress = TRUE,
    ForceRefit = FALSE,
    # Normal priors on regression coefficients and log shape; uniform prior on rho.
    InterceptPriorSd = 5,
    TreatmentPriorSd = 5,
    BaselinePriorSd = 5,
    LogShapePriorMean = 0,
    LogShapePriorSd = 2,
    RhoPrior = "Uniform",
    RunPredictiveChecks = TRUE,
    SavePredictiveDraws = FALSE,
    MakePlots = TRUE,
    # Every post-warmup draw supplies one imputation under every rule.
    # Compact RDS files are the default; optional CSV export also uses all draws.
    ExportCompletedCsv = FALSE,
    RunDistributionChecks = TRUE
  )
}

Settings <- DefaultBthebDkSettings()

Softplus <- function(Value) {
  pmax(Value, 0) + log1p(exp(-abs(Value)))
}

LogOneMinusExp <- function(LogValue) {
  Result <- rep(NA_real_, length(LogValue))
  Low <- !is.na(LogValue) & LogValue < -log(2)
  High <- !is.na(LogValue) & !Low & LogValue <= 0
  Result[Low] <- log1p(-exp(LogValue[Low]))
  Result[High] <- log(-expm1(LogValue[High]))
  Result
}

LogAddExp <- function(First, Second) {
  Maximum <- pmax(First, Second)
  Result <- Maximum + log1p(exp(-abs(First - Second)))
  BothZero <- is.infinite(Maximum) & Maximum < 0
  Result[BothZero] <- -Inf
  Result
}

LogDiffExp <- function(First, Second) {
  Result <- First + LogOneMinusExp(pmin(0, Second - First))
  Equal <- !is.na(First) & !is.na(Second) & First == Second
  Result[Equal] <- -Inf
  Result
}

LogNegativeLogOneMinusExp <- function(LogValue) {
  # log(-log(1 - exp(LogValue))). The tail branch avoids underflow.
  Result <- LogValue
  Interior <- !is.na(LogValue) & LogValue > -35
  Result[Interior] <- log(-LogOneMinusExp(LogValue[Interior]))
  Result
}

DkLogShapes <- function(Eta, LogShape, Tau = 0.5) {
  Shape <- exp(LogShape)
  list(
    LogA = LogShape - log(Softplus(-Eta)),
    LogB = log(-log1p(-Tau)) - LogNegativeLogOneMinusExp(-Shape)
  )
}

DkLogSurvival <- function(Y, Eta, LogShape, ScoreMaximum = 63, Tau = 0.5) {
  Length <- max(length(Y), length(Eta), length(LogShape))
  Y <- rep_len(Y, Length)
  Eta <- rep_len(Eta, Length)
  LogShape <- rep_len(LogShape, Length)
  Result <- rep(NA_real_, Length)
  Result[!is.na(Y) & Y < 0] <- 0
  Result[!is.na(Y) & Y >= ScoreMaximum] <- -Inf
  Interior <- which(!is.na(Y) & Y >= 0 & Y < ScoreMaximum)
  if (length(Interior) > 0) {
    Shapes <- DkLogShapes(Eta[Interior], LogShape[Interior], Tau)
    LogFraction <- log((floor(Y[Interior]) + 1)/(ScoreMaximum + 1))
    LogPower <- -exp(Shapes$LogA + log(-LogFraction))
    Result[Interior] <- -exp(Shapes$LogB + LogNegativeLogOneMinusExp(LogPower))
  }
  Result
}

DkCdf <- function(Y, Eta, LogShape, ScoreMaximum = 63, Tau = 0.5,
                  LowerTail = TRUE, LogProbability = FALSE) {
  LogSurvival <- DkLogSurvival(Y, Eta, LogShape, ScoreMaximum, Tau)
  Result <- if (LowerTail) LogOneMinusExp(LogSurvival) else LogSurvival
  if (LogProbability) Result else exp(Result)
}

DkInterval <- function(Y, Eta, LogShape, ScoreMaximum = 63, Tau = 0.5) {
  LogSLower <- DkLogSurvival(Y - 1, Eta, LogShape, ScoreMaximum, Tau)
  LogSUpper <- DkLogSurvival(Y, Eta, LogShape, ScoreMaximum, Tau)
  LogFLower <- LogOneMinusExp(LogSLower)
  LogFUpper <- LogOneMinusExp(LogSUpper)
  UseLower <- LogFUpper < -log(2)
  LogMass <- LogDiffExp(LogSLower, LogSUpper)
  LogMass[UseLower] <- LogDiffExp(LogFUpper[UseLower], LogFLower[UseLower])
  Invalid <- is.na(Y) | Y < 0 | Y > ScoreMaximum | Y != floor(Y)
  LogMass[Invalid] <- -Inf
  list(LogFLower = LogFLower, LogSUpper = LogSUpper, LogMass = LogMass)
}

DkPmf <- function(Y, Eta, LogShape, ScoreMaximum = 63, Tau = 0.5,
                  LogProbability = FALSE) {
  Result <- DkInterval(Y, Eta, LogShape, ScoreMaximum, Tau)$LogMass
  if (LogProbability) Result else exp(Result)
}

DkQuantileFromLogSurvival <- function(LogSurvivalProbability, Eta, LogShape,
                                      ScoreMaximum = 63, Tau = 0.5) {
  # The input is log(1 - U), computed directly from the Gaussian upper tail.
  # This avoids rounding U to 1 for a large positive latent Gaussian value.
  Shapes <- DkLogShapes(Eta, LogShape, Tau)
  LogExponent <- log(-LogSurvivalProbability) - Shapes$LogB
  LogOneMinusPower <- -exp(LogExponent)
  LogInside <- LogOneMinusExp(LogOneMinusPower)
  LogRoot <- -exp(log(-LogInside) - Shapes$LogA)
  Anchor <- (ScoreMaximum + 1)*exp(LogRoot) - 1
  # Remove only floating-point roundoff at an integer threshold before ceiling.
  NearestInteger <- round(Anchor)
  NearInteger <- is.finite(Anchor) & abs(Anchor - NearestInteger) <=
    8*.Machine$double.eps*pmax(1, abs(Anchor))
  Anchor[NearInteger] <- NearestInteger[NearInteger]
  Result <- pmin(ScoreMaximum, pmax(0, ceiling(Anchor)))
  Result[LogSurvivalProbability == 0] <- 0
  Result[LogSurvivalProbability == -Inf] <- ScoreMaximum
  Result
}

DkQuantile <- function(Probability, Eta, LogShape, ScoreMaximum = 63, Tau = 0.5) {
  if (anyNA(Probability) || any(Probability < 0 | Probability > 1)) {
    stop("Quantile probabilities must be between 0 and 1.")
  }
  DkQuantileFromLogSurvival(log1p(-Probability), Eta, LogShape, ScoreMaximum, Tau)
}

PitToGaussian <- function(PitEta, LogFLower, LogSUpper, LogMass) {
  LogV <- -Softplus(-PitEta)
  LogOneMinusV <- -Softplus(PitEta)
  LogU <- LogAddExp(LogFLower, LogV + LogMass)
  LogOneMinusU <- LogAddExp(LogSUpper, LogOneMinusV + LogMass)
  Result <- rep(NA_real_, length(LogU))
  Lower <- LogU <= -log(2)
  Result[Lower] <- qnorm(LogU[Lower], log.p = TRUE)
  Result[!Lower] <- qnorm(LogOneMinusU[!Lower], lower.tail = FALSE, log.p = TRUE)
  Result
}

CheckDkDistribution <- function() {
  Grid <- expand.grid(Eta = c(-6, -2, 0, 3), Shape = c(0.15, 0.5, 1, 4, 10),
                      Maximum = c(2, 15, 63))
  Probabilities <- c(1e-10, 0.001, 0.1, 0.5, 0.9, 0.999, 1 - 1e-10)
  for (Index in seq_len(nrow(Grid))) {
    Row <- Grid[Index, ]
    Mass <- DkPmf(0:Row$Maximum, Row$Eta, log(Row$Shape), Row$Maximum)
    if (any(!is.finite(Mass)) || any(Mass < 0) || abs(sum(Mass) - 1) > 1e-10) {
      stop("DK probability normalization check failed.")
    }
    Quantiles <- DkQuantile(Probabilities, Row$Eta, log(Row$Shape), Row$Maximum)
    Lower <- DkCdf(Quantiles - 1, Row$Eta, log(Row$Shape), Row$Maximum)
    Upper <- DkCdf(Quantiles, Row$Eta, log(Row$Shape), Row$Maximum)
    if (any(Lower > Probabilities + 1e-12) || any(Upper < Probabilities - 1e-12)) {
      stop("DK quantile inversion check failed.")
    }
    ExpectedMedian <- ceiling((Row$Maximum + 1)*plogis(Row$Eta) - 1)
    if (DkQuantile(0.5, Row$Eta, log(Row$Shape), Row$Maximum) != ExpectedMedian) {
      stop("DK median anchor check failed.")
    }
  }
  invisible(TRUE)
}

CheckSettings <- function(Settings) {
  if (!identical(Settings$Tau, 0.5)) stop("This application uses the median, Tau = 0.5.")
  if (length(Settings$VisitMonths) != length(Settings$VisitVariables) ||
      any(!is.finite(Settings$VisitMonths)) || any(diff(Settings$VisitMonths) <= 0)) {
    stop("Visit months must be finite and strictly increasing.")
  }
  for (Name in c("NChains", "Warmup", "Sample", "AdaptEvery", "PitSweeps", "PrintEvery")) {
    Value <- Settings[[Name]]
    if (length(Value) != 1 || !is.finite(Value) || Value < 1 || Value != floor(Value)) {
      stop(Name, " must be a positive integer.")
    }
  }
  if (Settings$NChains < 2 || Settings$Sample < 4) stop("Use at least two chains and four draws per chain.")
  if (Settings$CoreFraction <= 0 || Settings$CoreFraction > 1) stop("CoreFraction must be in (0, 1].")
  if (Settings$BaselineScale <= 0) stop("BaselineScale must be positive.")
  if (is.null(names(Settings$Rules)) || anyDuplicated(names(Settings$Rules)) ||
      any(!is.finite(Settings$Rules)) || any(Settings$Rules < 0 | Settings$Rules > 1) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9]*$", names(Settings$Rules)))) {
    stop("Rules need unique simple names and borrowing weights between 0 and 1.")
  }
  for (Name in c("InterceptPriorSd", "TreatmentPriorSd", "BaselinePriorSd", "LogShapePriorSd")) {
    Value <- Settings[[Name]]
    if (length(Value) != 1 || !is.finite(Value) || Value <= 0) {
      stop(Name, " must be a positive finite number.")
    }
  }
  if (length(Settings$LogShapePriorMean) != 1 || !is.finite(Settings$LogShapePriorMean)) {
    stop("LogShapePriorMean must be a finite number.")
  }
  if (!identical(Settings$RhoPrior, "Uniform")) {
    stop("RhoPrior must be 'Uniform'.")
  }
  invisible(TRUE)
}

PrepareBthebData <- function(Settings, DataWide = NULL) {
  if (is.null(DataWide)) {
    if (!requireNamespace("HSAUR2", quietly = TRUE)) stop("Install HSAUR2 before running this program.")
    DataEnvironment <- new.env(parent = baseenv())
    utils::data("BtheB", package = "HSAUR2", envir = DataEnvironment)
    DataWide <- DataEnvironment$BtheB
  }
  Required <- c("treatment", "bdi.pre", Settings$VisitVariables)
  Missing <- setdiff(Required, names(DataWide))
  if (length(Missing) > 0) stop("Missing BtheB variables: ", paste(Missing, collapse = ", "))
  DataWide <- as.data.frame(DataWide)
  if (!"OriginalSubjectId" %in% names(DataWide)) DataWide$OriginalSubjectId <- seq_len(nrow(DataWide))
  if (anyDuplicated(DataWide$OriginalSubjectId)) stop("Subject identifiers must be unique.")
  YAll <- as.matrix(DataWide[, Settings$VisitVariables, drop = FALSE])
  if (!is.numeric(YAll)) stop("BDI scores must be numeric.")
  if (any(!is.na(YAll) & (!is.finite(YAll) | YAll < 0 | YAll > Settings$ScoreMaximum |
                          YAll != floor(YAll)))) stop("BDI scores must be integers between 0 and 63.")
  Keep <- if (Settings$ExcludeNoPostBaseline) rowSums(!is.na(YAll)) > 0 else rep(TRUE, nrow(YAll))
  Dropped <- DataWide[!Keep, , drop = FALSE]
  DataWide <- DataWide[Keep, , drop = FALSE]
  rownames(DataWide) <- NULL
  if (nrow(DataWide) == 0) stop("No participants remain in the analysis population.")
  Groups <- as.character(DataWide$treatment)
  if (anyNA(Groups) || !setequal(unique(Groups), c("TAU", "BtheB"))) {
    stop("Treatment must contain TAU and BtheB.")
  }
  Baseline <- DataWide$bdi.pre
  if (anyNA(Baseline) || any(!is.finite(Baseline)) || any(Baseline < 0 | Baseline > Settings$ScoreMaximum |
                                                          Baseline != floor(Baseline))) {
    stop("A valid baseline BDI score is required for every included participant.")
  }
  DataWide$SubjectId <- seq_len(nrow(DataWide))
  DataWide$Treatment <- as.numeric(Groups == "BtheB")
  DataWide$treatment <- factor(Groups, levels = c("TAU", "BtheB"))
  BaselineCenter <- mean(Baseline)
  DataWide$BaselineCentered <- Baseline - BaselineCenter
  DataWide$BaselineScaled <- DataWide$BaselineCentered/Settings$BaselineScale
  Y <- as.matrix(DataWide[, Settings$VisitVariables, drop = FALSE])
  storage.mode(Y) <- "double"
  NSubjects <- nrow(Y)
  NVisits <- ncol(Y)
  Observed <- !is.na(Y)
  Trigger <- rep(NVisits + 1, NSubjects)
  PostTrigger <- matrix(FALSE, NSubjects, NVisits)
  Monotone <- rep(TRUE, NSubjects)
  for (Subject in seq_len(NSubjects)) {
    Visits <- which(Observed[Subject, ])
    LastObserved <- if (length(Visits) == 0) 0 else max(Visits)
    Trigger[Subject] <- LastObserved + 1
    PostTrigger[Subject, ] <- !Observed[Subject, ] & seq_len(NVisits) >= Trigger[Subject]
    Gaps <- which(!Observed[Subject, ])
    if (length(Gaps) > 0) Monotone[Subject] <- !any(Observed[Subject, seq.int(min(Gaps), NVisits)])
  }
  # Clinical event dates are unavailable, so the trigger follows the last
  # recorded visit. All recorded scores are retained; intermittent gaps use
  # continuation under the assigned treatment.
  PatternKeys <- apply(Observed, 1, function(Row) paste(as.integer(Row), collapse = ""))
  Patterns <- lapply(unique(PatternKeys), function(Key) {
    Subjects <- which(PatternKeys == Key)
    list(Subjects = Subjects, Visits = which(Observed[Subjects[1], ]))
  })
  ObservedCells <- which(Observed, arr.ind = TRUE)
  MissingCells <- which(!Observed, arr.ind = TRUE)
  Index <- list(
    Intercept = seq_len(NVisits),
    Treatment = NVisits + seq_len(NVisits),
    Baseline = 2*NVisits + seq_len(NVisits),
    LogShape = 3*NVisits + seq_len(NVisits),
    LogitRho = 4*NVisits + 1
  )
  ParameterNames <- c(paste0("Intercept[", seq_len(NVisits), "]"),
                      paste0("Treatment[", seq_len(NVisits), "]"),
                      paste0("Baseline[", seq_len(NVisits), "]"),
                      paste0("LogShape[", seq_len(NVisits), "]"), "LogitRho")
  TriggerTable <- data.frame(SubjectId = DataWide$SubjectId, OriginalSubjectId = DataWide$OriginalSubjectId,
                             Treatment = DataWide$treatment, TriggerVisit = Trigger,
                             TriggerMonth = c(Settings$VisitMonths, NA_real_)[Trigger],
                             Monotone = Monotone, NObserved = rowSums(Observed),
                             NIntermittentMissing = rowSums(!Observed & !PostTrigger),
                             NPostTriggerMissing = rowSums(PostTrigger))
  list(Data = DataWide, Dropped = Dropped, Y = Y, Observed = Observed,
       ObservedCells = ObservedCells, MissingCells = MissingCells, Patterns = Patterns,
       Treatment = DataWide$Treatment, Baseline = DataWide$BaselineScaled,
       BaselineCenter = BaselineCenter, BaselineScale = Settings$BaselineScale,
       Trigger = Trigger, PostTrigger = PostTrigger, TriggerTable = TriggerTable,
       NSubjects = NSubjects, NVisits = NVisits, Index = Index,
       ParameterNames = ParameterNames, VisitMonths = Settings$VisitMonths,
       VisitVariables = Settings$VisitVariables, ScoreMaximum = Settings$ScoreMaximum, Tau = Settings$Tau)
}

MakeContinuousAr1 <- function(VisitMonths, Rho) {
  if (length(Rho) != 1 || !is.finite(Rho) || Rho < 0 || Rho >= 1) {
    stop("The correlation per month must be in [0, 1).")
  }
  Distance <- abs(outer(VisitMonths, VisitMonths, "-"))
  if (Rho == 0) return(diag(length(VisitMonths)))
  exp(log(Rho)*Distance)
}

BuildCopulaGeometry <- function(RMatrix, Patterns) {
  lapply(Patterns, function(Pattern) {
    Dimension <- length(Pattern$Visits)
    if (Dimension <= 1) return(list(PrecisionMinusIdentity = matrix(0, Dimension, Dimension), LogDet = 0))
    Factor <- chol(RMatrix[Pattern$Visits, Pattern$Visits, drop = FALSE])
    list(PrecisionMinusIdentity = chol2inv(Factor) - diag(Dimension),
         LogDet = 2*sum(log(diag(Factor))))
  })
}

CopulaLogBySubject <- function(Z, Patterns, Geometry) {
  Result <- numeric(nrow(Z))
  for (PatternIndex in seq_along(Patterns)) {
    Pattern <- Patterns[[PatternIndex]]
    if (length(Pattern$Visits) <= 1) next
    Values <- Z[Pattern$Subjects, Pattern$Visits, drop = FALSE]
    Current <- Geometry[[PatternIndex]]
    Result[Pattern$Subjects] <- -0.5*Current$LogDet -
      0.5*rowSums((Values%*%Current$PrecisionMinusIdentity)*Values)
  }
  Result
}

ComputeDkPrior <- function(Parameter, Model, Settings) {
  Index <- Model$Index
  LogitRho <- Parameter[Index$LogitRho]
  # Uniform density on rho; change of variables to the sampled logit coordinate.
  LogRhoJacobian <- -Softplus(-LogitRho) - Softplus(LogitRho)
  sum(dnorm(Parameter[Index$Intercept], 0, Settings$InterceptPriorSd, log = TRUE)) +
    sum(dnorm(Parameter[Index$Treatment], 0, Settings$TreatmentPriorSd, log = TRUE)) +
    sum(dnorm(Parameter[Index$Baseline], 0, Settings$BaselinePriorSd, log = TRUE)) +
    sum(dnorm(Parameter[Index$LogShape], Settings$LogShapePriorMean, Settings$LogShapePriorSd, log = TRUE)) +
    LogRhoJacobian
}

EvaluateDkState <- function(Parameter, PitEta, Model, Settings, Current = NULL,
                            ChangedVisit = NULL, ChangeCorrelation = FALSE) {
  Invalid <- list(LogTarget = -Inf)
  if (any(!is.finite(Parameter))) return(Invalid)
  LogPrior <- ComputeDkPrior(Parameter, Model, Settings)
  if (!is.finite(LogPrior)) return(Invalid)
  Index <- Model$Index
  if (is.null(Current)) {
    State <- list(LogFLower = matrix(NA_real_, Model$NSubjects, Model$NVisits),
                  LogSUpper = matrix(NA_real_, Model$NSubjects, Model$NVisits),
                  LogMass = matrix(NA_real_, Model$NSubjects, Model$NVisits),
                  Z = matrix(NA_real_, Model$NSubjects, Model$NVisits))
    Visits <- seq_len(Model$NVisits)
    ChangeCorrelation <- TRUE
  } else {
    State <- Current
    Visits <- if (is.null(ChangedVisit)) integer(0) else ChangedVisit
  }
  if (ChangeCorrelation) {
    Rho <- plogis(Parameter[Index$LogitRho])
    if (!is.finite(Rho) || Rho <= 0 || Rho >= 1) return(Invalid)
    State$Correlation <- MakeContinuousAr1(Model$VisitMonths, Rho)
    State$Geometry <- tryCatch(BuildCopulaGeometry(State$Correlation, Model$Patterns),
                               error = function(Error) NULL)
    if (is.null(State$Geometry)) return(Invalid)
  }
  for (Visit in Visits) {
    Rows <- which(Model$Observed[, Visit])
    if (length(Rows) == 0) next
    Eta <- Parameter[Index$Intercept[Visit]] + Parameter[Index$Treatment[Visit]]*Model$Treatment[Rows] +
      Parameter[Index$Baseline[Visit]]*Model$Baseline[Rows]
    Interval <- DkInterval(Model$Y[Rows, Visit], Eta, Parameter[Index$LogShape[Visit]],
                           Model$ScoreMaximum, Model$Tau)
    if (any(!is.finite(Interval$LogMass))) return(Invalid)
    State$LogFLower[Rows, Visit] <- Interval$LogFLower
    State$LogSUpper[Rows, Visit] <- Interval$LogSUpper
    State$LogMass[Rows, Visit] <- Interval$LogMass
    State$Z[Rows, Visit] <- PitToGaussian(PitEta[Rows, Visit], Interval$LogFLower,
                                          Interval$LogSUpper, Interval$LogMass)
    if (any(!is.finite(State$Z[Rows, Visit]))) return(Invalid)
  }
  State$LogCopulaBySubject <- CopulaLogBySubject(State$Z, Model$Patterns, State$Geometry)
  State$LogPrior <- LogPrior
  State$LogMarginal <- sum(State$LogMass[Model$Observed])
  State$LogJacobianBySubject <- rowSums(-Softplus(-PitEta) - Softplus(PitEta), na.rm = TRUE)
  State$LogTarget <- LogPrior + State$LogMarginal +
    sum(State$LogCopulaBySubject) + sum(State$LogJacobianBySubject)
  if (!is.finite(State$LogTarget)) return(Invalid)
  State
}

UpdateDkPit <- function(PitEta, State, Model, ProposalSd) {
  # Subjects are independent conditional on model parameters. Their MH proposals
  # are evaluated together, with a separate acceptance decision per participant.
  CandidateEta <- PitEta
  CandidateEta[Model$Observed] <- PitEta[Model$Observed] + rnorm(sum(Model$Observed), sd = ProposalSd)
  CandidateZ <- State$Z
  CandidateZ[Model$Observed] <- PitToGaussian(CandidateEta[Model$Observed],
                                              State$LogFLower[Model$Observed],
                                              State$LogSUpper[Model$Observed],
                                              State$LogMass[Model$Observed])
  CandidateCopula <- CopulaLogBySubject(CandidateZ, Model$Patterns, State$Geometry)
  CandidateJacobian <- rowSums(-Softplus(-CandidateEta) - Softplus(CandidateEta), na.rm = TRUE)
  LogRatio <- CandidateCopula + CandidateJacobian - State$LogCopulaBySubject - State$LogJacobianBySubject
  Eligible <- rowSums(Model$Observed) > 0
  Accepted <- Eligible & is.finite(LogRatio) & log(runif(Model$NSubjects)) < LogRatio
  if (any(Accepted)) {
    PitEta[Accepted, ] <- CandidateEta[Accepted, ]
    State$Z[Accepted, ] <- CandidateZ[Accepted, ]
    State$LogCopulaBySubject[Accepted] <- CandidateCopula[Accepted]
    State$LogJacobianBySubject[Accepted] <- CandidateJacobian[Accepted]
    State$LogTarget <- State$LogPrior + State$LogMarginal +
      sum(State$LogCopulaBySubject) + sum(State$LogJacobianBySubject)
  }
  list(PitEta = PitEta, State = State, Accepted = sum(Accepted), Proposed = sum(Eligible))
}

RegularizeProposalCovariance <- function(Covariance, Minimum = 1e-6, Maximum = 4) {
  Covariance <- (Covariance + t(Covariance))/2
  Decomposition <- eigen(Covariance, symmetric = TRUE)
  Values <- pmin(Maximum, pmax(Minimum, Decomposition$values))
  Decomposition$vectors%*%diag(Values, nrow = length(Values))%*%t(Decomposition$vectors)
}

InitialDkParameters <- function(Model, Settings) {
  Index <- Model$Index
  Parameter <- setNames(numeric(length(Model$ParameterNames)), Model$ParameterNames)
  Covariances <- vector("list", Model$NVisits + 1)
  for (Visit in seq_len(Model$NVisits)) {
    Rows <- which(Model$Observed[, Visit])
    if (length(Rows) < 5 || length(unique(Model$Treatment[Rows])) < 2) {
      stop("Too few observed scores in one treatment group at visit ", Visit, ".")
    }
    X <- cbind(1, Model$Treatment[Rows], Model$Baseline[Rows])
    Transformed <- qlogis((Model$Y[Rows, Visit] + 0.5)/(Model$ScoreMaximum + 1))
    Start <- as.numeric(qr.solve(X, Transformed))
    Start <- c(Start, 0)
    PriorSd <- c(Settings$InterceptPriorSd, Settings$TreatmentPriorSd,
                 Settings$BaselinePriorSd, Settings$LogShapePriorSd)
    PriorMean <- c(0, 0, 0, Settings$LogShapePriorMean)
    Objective <- function(Values) {
      LogMass <- DkPmf(Model$Y[Rows, Visit], as.numeric(X%*%Values[1:3]), Values[4],
                       Model$ScoreMaximum, Model$Tau, LogProbability = TRUE)
      Value <- -sum(LogMass) - sum(dnorm(Values, PriorMean, PriorSd, log = TRUE))
      if (is.finite(Value)) Value else 1e100
    }
    Fit <- optim(Start, Objective, method = "BFGS", hessian = TRUE,
                 control = list(maxit = 1000, reltol = 1e-9))
    Values <- if (is.finite(Fit$value) && Fit$value < 1e90) Fit$par else Start
    BlockIndex <- c(Index$Intercept[Visit], Index$Treatment[Visit],
                    Index$Baseline[Visit], Index$LogShape[Visit])
    Parameter[BlockIndex] <- Values
    Hessian <- (Fit$hessian + t(Fit$hessian))/2
    Covariance <- tryCatch(solve(Hessian), error = function(Error) diag(0.02, 4))
    if (any(!is.finite(Covariance))) Covariance <- diag(0.02, 4)
    Covariances[[Visit]] <- RegularizeProposalCovariance(Covariance)
  }
  Parameter[Index$LogitRho] <- 0
  Covariances[[Model$NVisits + 1]] <- matrix(0.04, 1, 1)
  list(Parameter = Parameter, Covariances = Covariances)
}

BuildDkBlocks <- function(Model, Initial) {
  Index <- Model$Index
  lapply(seq_len(Model$NVisits + 1), function(Block) {
    IsCorrelation <- Block > Model$NVisits
    Positions <- if (IsCorrelation) Index$LogitRho else
      c(Index$Intercept[Block], Index$Treatment[Block], Index$Baseline[Block], Index$LogShape[Block])
    Dimension <- length(Positions)
    list(Name = if (IsCorrelation) "LogitRho" else paste0("Visit", Block),
         Index = Positions, Visit = if (IsCorrelation) NULL else Block,
         IsCorrelation = IsCorrelation, Dimension = Dimension,
         Covariance = Initial$Covariances[[Block]],
         Cholesky = chol(Initial$Covariances[[Block]]), Scale = 1,
         Count = 0, Mean = numeric(Dimension), M2 = matrix(0, Dimension, Dimension),
         AcceptedWindow = 0, ProposedWindow = 0,
         AcceptedSampling = 0, ProposedSampling = 0)
  })
}

UpdateBlockMoments <- function(Block, Values) {
  Block$Count <- Block$Count + 1
  Difference <- Values - Block$Mean
  Block$Mean <- Block$Mean + Difference/Block$Count
  Block$M2 <- Block$M2 + outer(Difference, Values - Block$Mean)
  Block
}

AdaptDkBlock <- function(Block, AdaptationNumber) {
  Rate <- Block$AcceptedWindow/max(1, Block$ProposedWindow)
  Target <- if (Block$Dimension == 1) 0.44 else 0.30
  Gain <- min(0.5, 1/sqrt(AdaptationNumber))
  Block$Scale <- min(20, max(0.02, Block$Scale*exp(2*Gain*(Rate - Target))))
  if (Block$Count > 2*Block$Dimension && AdaptationNumber >= 5) {
    Empirical <- Block$M2/(Block$Count - 1)
    Covariance <- 0.95*Empirical + 0.05*diag(diag(Block$Covariance), Block$Dimension)
    Block$Cholesky <- chol(RegularizeProposalCovariance(Covariance))
  }
  Block$AcceptedWindow <- 0
  Block$ProposedWindow <- 0
  Block
}

SaveRdsSafely <- function(Object, Path) {
  Temporary <- tempfile(pattern = "DkSave", tmpdir = dirname(Path))
  on.exit(unlink(Temporary), add = TRUE)
  saveRDS(Object, Temporary, compress = "gzip")
  if (!file.copy(Temporary, Path, overwrite = TRUE)) stop("Could not save: ", Path)
  invisible(Path)
}

HashDkObject <- function(Object) {
  Path <- tempfile(pattern = "DkSignature")
  on.exit(unlink(Path), add = TRUE)
  saveRDS(Object, Path, compress = FALSE, version = 2)
  unname(tools::md5sum(Path))
}

RunDkChain <- function(Chain, Model, Settings, Initial, Signature, CacheDirectory) {
  CacheFile <- file.path(CacheDirectory, paste0("Chain", Chain, ".rds"))
  if (!Settings$ForceRefit && file.exists(CacheFile)) {
    Cached <- tryCatch(readRDS(CacheFile), error = function(Error) NULL)
    if (!is.null(Cached) && identical(Cached$Signature, Signature)) {
      if (Settings$PrintProgress) cat("Loaded completed chain", Chain, "\n")
      return(Cached)
    }
  }
  set.seed(Settings$Seed + 1000*Chain)
  StartTime <- proc.time()[[3]]
  Parameter <- Initial$Parameter
  Blocks <- BuildDkBlocks(Model, Initial)
  PitEta <- matrix(NA_real_, Model$NSubjects, Model$NVisits)
  PitEta[Model$Observed] <- qlogis(runif(sum(Model$Observed), 0.1, 0.9))
  for (Attempt in seq_len(20)) {
    Candidate <- Initial$Parameter
    for (Block in Blocks) {
      Candidate[Block$Index] <- Candidate[Block$Index] +
        as.numeric(rnorm(Block$Dimension)%*%Block$Cholesky)*0.5
    }
    State <- EvaluateDkState(Candidate, PitEta, Model, Settings)
    if (is.finite(State$LogTarget)) {
      Parameter <- Candidate
      break
    }
  }
  if (!is.finite(State$LogTarget)) stop("Could not initialize chain ", Chain, ".")
  Draws <- matrix(NA_real_, Settings$Sample, length(Parameter),
                  dimnames = list(NULL, Model$ParameterNames))
  ZObserved <- matrix(NA_real_, Settings$Sample, sum(Model$Observed))
  LogLikelihood <- numeric(Settings$Sample)
  PitSd <- 0.8
  PitAcceptedWindow <- 0
  PitProposedWindow <- 0
  PitAcceptedSampling <- 0
  PitProposedSampling <- 0
  TotalIterations <- Settings$Warmup + Settings$Sample
  for (Iteration in seq_len(TotalIterations)) {
    for (BlockIndex in seq_along(Blocks)) {
      Block <- Blocks[[BlockIndex]]
      Candidate <- Parameter
      Candidate[Block$Index] <- Candidate[Block$Index] +
        as.numeric(rnorm(Block$Dimension)%*%Block$Cholesky)*Block$Scale
      CandidateState <- EvaluateDkState(Candidate, PitEta, Model, Settings, Current = State,
                                        ChangedVisit = Block$Visit, ChangeCorrelation = Block$IsCorrelation)
      LogRatio <- CandidateState$LogTarget - State$LogTarget
      Accepted <- is.finite(LogRatio) && log(runif(1)) < LogRatio
      if (Accepted) {
        Parameter <- Candidate
        State <- CandidateState
      }
      Block$ProposedWindow <- Block$ProposedWindow + 1
      Block$AcceptedWindow <- Block$AcceptedWindow + as.numeric(Accepted)
      if (Iteration > Settings$Warmup) {
        Block$ProposedSampling <- Block$ProposedSampling + 1
        Block$AcceptedSampling <- Block$AcceptedSampling + as.numeric(Accepted)
      } else {
        Block <- UpdateBlockMoments(Block, Parameter[Block$Index])
      }
      Blocks[[BlockIndex]] <- Block
    }
    for (Sweep in seq_len(Settings$PitSweeps)) {
      Update <- UpdateDkPit(PitEta, State, Model, PitSd)
      PitEta <- Update$PitEta
      State <- Update$State
      PitAcceptedWindow <- PitAcceptedWindow + Update$Accepted
      PitProposedWindow <- PitProposedWindow + Update$Proposed
      if (Iteration > Settings$Warmup) {
        PitAcceptedSampling <- PitAcceptedSampling + Update$Accepted
        PitProposedSampling <- PitProposedSampling + Update$Proposed
      }
    }
    if (Iteration <= Settings$Warmup && Iteration %% Settings$AdaptEvery == 0) {
      AdaptationNumber <- Iteration/Settings$AdaptEvery
      Blocks <- lapply(Blocks, AdaptDkBlock, AdaptationNumber = AdaptationNumber)
      Gain <- min(0.5, 1/sqrt(AdaptationNumber))
      Rate <- PitAcceptedWindow/max(1, PitProposedWindow)
      PitSd <- min(10, max(0.02, PitSd*exp(2*Gain*(Rate - 0.4))))
      PitAcceptedWindow <- 0
      PitProposedWindow <- 0
    }
    if (Iteration > Settings$Warmup) {
      Row <- Iteration - Settings$Warmup
      Draws[Row, ] <- Parameter
      ZObserved[Row, ] <- State$Z[Model$Observed]
      LogLikelihood[Row] <- State$LogMarginal + sum(State$LogCopulaBySubject)
    }
    if (Settings$PrintProgress && Chain == 1 && Iteration %% Settings$PrintEvery == 0) {
      cat("Chain 1:", Iteration, "of", TotalIterations, "iterations; rho per month =",
          round(plogis(Parameter[Model$Index$LogitRho]), 3), "\n")
    }
  }
  Acceptance <- do.call(rbind, lapply(Blocks, function(Block) {
    data.frame(Chain = Chain, Block = Block$Name,
               Acceptance = Block$AcceptedSampling/max(1, Block$ProposedSampling),
               ProposalScale = Block$Scale)
  }))
  Acceptance <- rbind(Acceptance, data.frame(Chain = Chain, Block = "PIT",
                                             Acceptance = PitAcceptedSampling/max(1, PitProposedSampling),
                                             ProposalScale = PitSd))
  Result <- list(Chain = Chain, Signature = Signature, Draws = Draws, ZObserved = ZObserved,
                 LogAugmentedLikelihood = LogLikelihood, Acceptance = Acceptance,
                 RuntimeSeconds = proc.time()[[3]] - StartTime)
  SaveRdsSafely(Result, CacheFile)
  if (Settings$PrintProgress) cat("Completed chain", Chain, "in", round(Result$RuntimeSeconds, 1), "seconds.\n")
  Result
}

RunDkChainTask <- function(Task) {
  RunDkChain(Task$Chain, Task$Model, Task$Settings, Task$Initial, Task$Signature, Task$CacheDirectory)
}

DkSamplerFunctionNames <- function() {
  c("Softplus", "LogOneMinusExp", "LogAddExp", "LogDiffExp", "LogNegativeLogOneMinusExp",
    "DkLogShapes", "DkLogSurvival", "DkCdf", "DkInterval", "DkPmf", "DkQuantileFromLogSurvival",
    "DkQuantile", "PitToGaussian", "MakeContinuousAr1", "BuildCopulaGeometry", "CopulaLogBySubject",
    "ComputeDkPrior", "EvaluateDkState", "UpdateDkPit", "RegularizeProposalCovariance",
    "BuildDkBlocks", "UpdateBlockMoments", "AdaptDkBlock", "SaveRdsSafely", "RunDkChain", "RunDkChainTask")
}

FitDkCopula <- function(Model, Settings, ModelDirectory) {
  Initial <- InitialDkParameters(Model, Settings)
  FitSettingNames <- c(
    "Tau", "ScoreMaximum", "NChains", "Warmup", "Sample", "Seed",
    "AdaptEvery", "PitSweeps", "InterceptPriorSd", "TreatmentPriorSd",
    "BaselinePriorSd", "LogShapePriorMean", "LogShapePriorSd", "RhoPrior"
  )
  FunctionEnvironment <- environment(RunDkChain)
  FunctionBodies <- lapply(c(DkSamplerFunctionNames(), "InitialDkParameters"), function(Name) {
    Function <- get(Name, envir = FunctionEnvironment)
    # Text excludes source-reference environments and runtime compilation state,
    # so executing an unchanged function does not invalidate its cache signature.
    list(Name = Name, Formals = deparse(formals(Function), width.cutoff = 500),
         Body = deparse(body(Function), width.cutoff = 500))
  })
  Signature <- HashDkObject(list(Version = "DKMedianGaussianCAR1v1", Model = Model,
                                 Settings = Settings[FitSettingNames], Functions = FunctionBodies))
  Tasks <- lapply(seq_len(Settings$NChains), function(Chain) {
    list(Chain = Chain, Model = Model, Settings = Settings, Initial = Initial,
         Signature = Signature, CacheDirectory = ModelDirectory)
  })
  Detected <- parallel::detectCores(logical = TRUE)
  if (!is.finite(Detected)) Detected <- 1
  NCores <- if (Settings$UseParallel) min(Settings$NChains, max(1, floor(Settings$CoreFraction*Detected))) else 1
  if (Settings$PrintProgress) cat("Fitting Gaussian copula with DK median margins using", NCores, "worker(s).\n")
  StartTime <- proc.time()[[3]]
  if (NCores > 1) {
    Cluster <- parallel::makeCluster(NCores, type = "PSOCK", outfile = "")
    on.exit(parallel::stopCluster(Cluster), add = TRUE)
    parallel::clusterExport(Cluster, DkSamplerFunctionNames(), envir = FunctionEnvironment)
    Chains <- parallel::parLapply(Cluster, Tasks, RunDkChainTask)
  } else {
    Chains <- lapply(Tasks, RunDkChainTask)
  }
  list(Chains = Chains, Signature = Signature, RuntimeSeconds = proc.time()[[3]] - StartTime,
       TotalDraws = Settings$NChains*Settings$Sample, Initial = Initial)
}

DescribeDkDraws <- function(Values) {
  Quantiles <- as.numeric(quantile(Values, c(0.025, 0.5, 0.975), names = FALSE, type = 8))
  c(Mean = mean(Values), Sd = sd(Values), Q025 = Quantiles[1], Median = Quantiles[2], Q975 = Quantiles[3])
}

DerivedDkParameters <- function(Draws, Model) {
  Index <- Model$Index
  Rho <- plogis(Draws[, Index$LogitRho])
  Derived <- matrix(Rho, ncol = 1, dimnames = list(NULL, "RhoPerMonth"))
  for (Visit in seq_len(Model$NVisits)) {
    ReferenceAnchor <- (Model$ScoreMaximum + 1)*plogis(Draws[, Index$Intercept[Visit]]) - 1
    ActiveAnchor <- (Model$ScoreMaximum + 1)*plogis(Draws[, Index$Intercept[Visit]] +
                                                      Draws[, Index$Treatment[Visit]]) - 1
    Values <- cbind(exp(Draws[, Index$LogShape[Visit]]), ReferenceAnchor, ActiveAnchor)
    colnames(Values) <- paste0(c("Shape[", "TauMedianAnchor[", "BtheBMedianAnchor["), Visit, "]")
    Derived <- cbind(Derived, Values)
  }
  cbind(Draws, Derived)
}

SummariseDkFit <- function(Fit, Model, SummaryDirectory, PlotDirectory, Settings) {
  if (!requireNamespace("posterior", quietly = TRUE)) stop("Install posterior before running this program.")
  Derived <- lapply(Fit$Chains, function(Chain) DerivedDkParameters(Chain$Draws, Model))
  Array <- array(NA_real_, c(Settings$Sample, Settings$NChains, ncol(Derived[[1]])),
                 dimnames = list(iteration = seq_len(Settings$Sample), chain = seq_len(Settings$NChains),
                                 variable = colnames(Derived[[1]])))
  for (Chain in seq_len(Settings$NChains)) Array[, Chain, ] <- Derived[[Chain]]
  Diagnostics <- posterior::summarise_draws(posterior::as_draws_array(Array),
                                            "rhat", "ess_bulk", "ess_tail", "mcse_mean")
  Combined <- do.call(rbind, Derived)
  Summary <- as.data.frame(t(apply(Combined, 2, DescribeDkDraws)))
  Summary$Parameter <- rownames(Summary)
  rownames(Summary) <- NULL
  Match <- match(Summary$Parameter, Diagnostics$variable)
  Summary$Rhat <- Diagnostics$rhat[Match]
  Summary$EssBulk <- Diagnostics$ess_bulk[Match]
  Summary$EssTail <- Diagnostics$ess_tail[Match]
  Summary$McseMean <- Diagnostics$mcse_mean[Match]
  Summary <- Summary[, c("Parameter", "Mean", "Sd", "Q025", "Median", "Q975",
                         "Rhat", "EssBulk", "EssTail", "McseMean")]
  Summary$NDraws <- Fit$TotalDraws
  write.csv(Summary, file.path(SummaryDirectory, "ParameterPosteriorSummary.csv"), row.names = FALSE)
  Acceptance <- do.call(rbind, lapply(Fit$Chains, function(Chain) Chain$Acceptance))
  write.csv(Acceptance, file.path(SummaryDirectory, "SamplerAcceptance.csv"), row.names = FALSE)
  CorrelationRows <- list()
  Rho <- Combined[, "RhoPerMonth"]
  Counter <- 0
  for (First in seq_len(Model$NVisits - 1)) {
    for (Second in seq.int(First + 1, Model$NVisits)) {
      Counter <- Counter + 1
      CorrelationRows[[Counter]] <- data.frame(FirstMonth = Model$VisitMonths[First],
                                               SecondMonth = Model$VisitMonths[Second],
                                               GapMonths = Model$VisitMonths[Second] - Model$VisitMonths[First],
                                               t(DescribeDkDraws(Rho^(Model$VisitMonths[Second] - Model$VisitMonths[First]))))
    }
  }
  write.csv(do.call(rbind, CorrelationRows), file.path(SummaryDirectory, "LatentCorrelations.csv"), row.names = FALSE)
  Core <- Summary[Summary$Parameter %in% Model$ParameterNames, ]
  Flag <- !is.finite(Core$Rhat) | Core$Rhat > 1.01 | !is.finite(Core$EssBulk) |
    Core$EssBulk < 400 | !is.finite(Core$EssTail) | Core$EssTail < 400
  if (any(Flag)) {
    warning("Some sampling diagnostics need review; see ParameterPosteriorSummary.csv. All draws are retained.", call. = FALSE)
  }
  if (Settings$MakePlots) PlotDkTraces(Derived, Model, PlotDirectory)
  list(Parameters = Summary, Acceptance = Acceptance, DiagnosticFlag = any(Flag))
}

PlotDkTraces <- function(Chains, Model, PlotDirectory) {
  Colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#555555", "#E69F00")
  pdf(file.path(PlotDirectory, "ParameterTracePlots.pdf"), width = 10, height = 8, family = "Helvetica")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(3, 2), mar = c(3, 3.5, 2, 0.5), mgp = c(1.8, 0.6, 0))
  for (Name in c(Model$ParameterNames, "RhoPerMonth")) {
    Limits <- range(vapply(Chains, function(Chain) range(Chain[, Name]), numeric(2)))
    plot(seq_len(nrow(Chains[[1]])), Chains[[1]][, Name], type = "n", ylim = Limits,
         xlab = "Post-warmup iteration", ylab = "Value", main = Name)
    for (Chain in seq_along(Chains)) {
      lines(Chains[[Chain]][, Name], col = Colors[(Chain - 1) %% length(Colors) + 1], lwd = 0.6)
    }
    if (Name == Model$ParameterNames[1]) {
      legend("topright", legend = paste("Chain", seq_along(Chains)),
             col = rep_len(Colors, length(Chains)), lty = 1, bty = "n", cex = 0.7)
    }
  }
}

DkTableNumberFormat <- function(Values, Digits) {
  Precision <- max(Digits)
  Rounded <- round(Values, Digits)
  Rounded <- Rounded[is.finite(Rounded)]
  if (length(Rounded) == 0) return(paste0("1.", Precision))
  IntegerDigits <- max(nchar(format(floor(abs(Rounded)), trim = TRUE, scientific = FALSE)))
  paste0(if (any(Rounded < 0)) "-" else "", IntegerDigits, ".", Precision)
}

ExportDkObservedTable <- function(ObservedSummary, Model, SummaryDirectory) {
  # Format the same observed summaries as the CSV; no sampling is performed.
  # Separate numeric columns align the count, percentage, mean, and SD independently.
  FormatNumber <- function(Value) {
    if (is.finite(Value)) sprintf("%.1f", Value) else "{--}"
  }
  ParenthesizedNumber <- function(Value) {
    if (is.finite(Value)) paste0("{(}", FormatNumber(Value), "{)}") else "{(--)}"
  }
  Formats <- c(
    Observed = DkTableNumberFormat(ObservedSummary$NObserved, 0),
    Missing = DkTableNumberFormat(ObservedSummary$NMissing, 0),
    Percent = DkTableNumberFormat(100*ObservedSummary$NMissing/ObservedSummary$NParticipants, 1),
    Mean = DkTableNumberFormat(ObservedSummary$Mean, 1),
    Sd = DkTableNumberFormat(ObservedSummary$Sd, 1),
    Median = DkTableNumberFormat(ObservedSummary$Median, 1)
  )
  Treatments <- unique(as.character(ObservedSummary$Treatment))
  Lines <- c(
    "% Requires \\usepackage{booktabs,siunitx,bm}.",
    "% Generated from the recorded scores in the analysis population.",
    "\\begin{table}[tbp]",
    paste0("\\caption{Observed post-baseline BDI-II scores in the ", Model$NSubjects,
           "-participant analysis population.}"),
    "\\label{tab:observed}",
    "\\centering\\small",
    "\\sisetup{table-align-text-pre = false, table-align-text-post = false,",
    "  group-digits = false, output-decimal-marker = {.}, round-mode = none}",
    "\\begin{tabular}{",
    "l c",
    paste0("S[table-format = ", Formats["Observed"], "]"),
    paste0("S[table-format = ", Formats["Missing"], "] @{\\;}"),
    paste0("S[table-format = ", Formats["Percent"],
           ", table-space-text-pre = {(}, table-space-text-post = {)}]"),
    paste0("S[table-format = ", Formats["Mean"], "] @{\\;}"),
    paste0("S[table-format = ", Formats["Sd"],
           ", table-space-text-pre = {(}, table-space-text-post = {)}]"),
    paste0("S[table-format = ", Formats["Median"], "]"),
    "}",
    "\\toprule",
    "\\textbf{Treatment} & {\\textbf{Month}} & {\\textbf{Observed $\\bm n$}} &",
    "\\multicolumn{2}{c}{\\textbf{Missing $\\bm n$ (\\%)}} &",
    "\\multicolumn{2}{c}{\\textbf{Mean (SD)}} & {\\textbf{Median}}\\\\",
    "\\midrule"
  )
  Denominators <- character(length(Treatments))
  for (Group in seq_along(Treatments)) {
    GroupRows <- ObservedSummary[ObservedSummary$Treatment == Treatments[Group], , drop = FALSE]
    Denominators[Group] <- paste(GroupRows$NParticipants[1], Treatments[Group])
    if (Group > 1) Lines <- c(Lines, "\\addlinespace")
    for (Visit in seq_len(nrow(GroupRows))) {
      Row <- GroupRows[Visit, , drop = FALSE]
      MissingPercent <- 100*Row$NMissing/Row$NParticipants
      Cells <- c(
        if (Visit == 1) Treatments[Group] else "",
        sprintf("%g", Row$Month), as.character(Row$NObserved), as.character(Row$NMissing),
        ParenthesizedNumber(MissingPercent), FormatNumber(Row$Mean),
        ParenthesizedNumber(Row$Sd), FormatNumber(Row$Median)
      )
      Lines <- c(Lines, paste0(paste(Cells, collapse = " & "), "\\\\"))
    }
  }
  Lines <- c(
    Lines, "\\bottomrule", "\\end{tabular}", "\\par\\smallskip\\raggedright\\footnotesize",
    paste0("Missing percentages use ", paste(Denominators, collapse = " and "),
           " participants as denominators. Summary statistics use recorded scores only; ",
           "the contributing participants change over follow-up."),
    "\\end{table}"
  )
  File <- file.path(SummaryDirectory, "ObservedScoreSummary.tex")
  writeLines(Lines, File)
  invisible(File)
}

ExportDkInputs <- function(Model, Settings, SummaryDirectory) {
  write.csv(Model$Data, file.path(SummaryDirectory, "AnalysisPopulation.csv"), row.names = FALSE)
  write.csv(Model$Dropped, file.path(SummaryDirectory, "ExcludedNoPostBaseline.csv"), row.names = FALSE)
  write.csv(Model$TriggerTable, file.path(SummaryDirectory, "MissingnessTriggers.csv"), row.names = FALSE)
  Rows <- list()
  Counter <- 0
  for (Treatment in 0:1) {
    Subjects <- which(Model$Treatment == Treatment)
    for (Visit in seq_len(Model$NVisits)) {
      Values <- Model$Y[Subjects, Visit]
      Recorded <- Values[!is.na(Values)]
      Counter <- Counter + 1
      Rows[[Counter]] <- data.frame(Treatment = c("TAU", "BtheB")[Treatment + 1],
                                    Month = Model$VisitMonths[Visit], NParticipants = length(Subjects),
                                    NObserved = length(Recorded), NMissing = sum(is.na(Values)),
                                    Mean = mean(Recorded), Sd = sd(Recorded), Median = median(Recorded),
                                    PZero = mean(Recorded == 0), PCeiling = mean(Recorded == Model$ScoreMaximum))
    }
  }
  ObservedSummary <- do.call(rbind, Rows)
  write.csv(ObservedSummary, file.path(SummaryDirectory, "ObservedScoreSummary.csv"), row.names = FALSE)
  ExportDkObservedTable(ObservedSummary, Model, SummaryDirectory)
  Population <- data.frame(NOriginal = nrow(Model$Data) + nrow(Model$Dropped), NIncluded = nrow(Model$Data),
                           NExcludedNoPostBaseline = nrow(Model$Dropped),
                           AllMonotone = all(Model$TriggerTable$Monotone), BaselineCenter = Model$BaselineCenter,
                           BaselineScale = Model$BaselineScale)
  write.csv(Population, file.path(SummaryDirectory, "PopulationSummary.csv"), row.names = FALSE)
}

SampleDkMissingGaussians <- function(Model, Correlation, ZObserved) {
  Z <- matrix(NA_real_, Model$NSubjects, Model$NVisits)
  Z[Model$Observed] <- ZObserved
  for (Pattern in Model$Patterns) {
    Observed <- Pattern$Visits
    Missing <- setdiff(seq_len(Model$NVisits), Observed)
    if (length(Missing) == 0) next
    Subjects <- Pattern$Subjects
    if (length(Observed) == 0) {
      Mean <- matrix(0, length(Subjects), length(Missing))
      Variance <- Correlation[Missing, Missing, drop = FALSE]
    } else {
      Regression <- t(solve(Correlation[Observed, Observed, drop = FALSE],
                            Correlation[Observed, Missing, drop = FALSE]))
      Mean <- Z[Subjects, Observed, drop = FALSE]%*%t(Regression)
      Variance <- Correlation[Missing, Missing, drop = FALSE] -
        Regression%*%Correlation[Observed, Missing, drop = FALSE]
    }
    Variance <- (Variance + t(Variance))/2
    Factor <- chol(Variance)
    Noise <- matrix(rnorm(length(Subjects)*length(Missing)), length(Subjects), length(Missing))
    Z[Subjects, Missing] <- Mean + Noise%*%Factor
  }
  Z[Model$MissingCells]
}

MapDkMissing <- function(Parameter, MissingLatents, Model, Omega) {
  Cells <- Model$MissingCells
  if (nrow(Cells) == 0) return(integer(0))
  Subjects <- Cells[, 1]
  Visits <- Cells[, 2]
  Index <- Model$Index
  ReferenceEta <- Parameter[Index$Intercept[Visits]] + Parameter[Index$Baseline[Visits]]*Model$Baseline[Subjects]
  AssignedDifference <- Parameter[Index$Treatment[Visits]]*Model$Treatment[Subjects]
  Weight <- Omega*Model$Treatment[Subjects]*Model$PostTrigger[Cells]
  RuleEta <- ReferenceEta + (1 - Weight)*AssignedDifference
  Values <- DkQuantileFromLogSurvival(pnorm(MissingLatents, lower.tail = FALSE, log.p = TRUE),
                                      RuleEta, Parameter[Index$LogShape[Visits]], Model$ScoreMaximum, Model$Tau)
  if (any(!is.finite(Values)) || any(Values < 0 | Values > Model$ScoreMaximum) ||
      any(Values != floor(Values))) stop("An imputed score is outside the integer support.")
  as.integer(Values)
}

BuildDkDrawIndex <- function(Fit) {
  do.call(rbind, lapply(Fit$Chains, function(Chain) {
    data.frame(Chain = Chain$Chain, PostWarmupIteration = seq_len(nrow(Chain$Draws)))
  }))
}

ImputeDkOutcomes <- function(Fit, Model, Settings, CompletedDirectory, SummaryDirectory, PlotDirectory) {
  set.seed(Settings$Seed + 100000)
  DrawIndex <- BuildDkDrawIndex(Fit)
  DrawIndex$Imputation <- seq_len(nrow(DrawIndex))
  DrawIndex <- DrawIndex[, c("Imputation", "Chain", "PostWarmupIteration")]
  Cells <- Model$MissingCells
  NMissing <- nrow(Cells)
  RuleNames <- names(Settings$Rules)
  Imputed <- setNames(lapply(RuleNames, function(Rule) matrix(NA_integer_, NMissing, Fit$TotalDraws)), RuleNames)
  Draw <- 0
  if (Settings$PrintProgress) cat("Generating", Fit$TotalDraws, "completed datasets per rule using every post-warmup draw.\n")
  for (Chain in Fit$Chains) {
    for (Iteration in seq_len(nrow(Chain$Draws))) {
      Draw <- Draw + 1
      Parameter <- Chain$Draws[Iteration, ]
      Correlation <- MakeContinuousAr1(Model$VisitMonths, plogis(Parameter[Model$Index$LogitRho]))
      Latents <- SampleDkMissingGaussians(Model, Correlation, Chain$ZObserved[Iteration, ])
      for (Rule in RuleNames) Imputed[[Rule]][, Draw] <- MapDkMissing(Parameter, Latents, Model, Settings$Rules[[Rule]])
      if (Settings$PrintProgress && Draw %% Settings$PrintEvery == 0) cat("Imputed", Draw, "of", Fit$TotalDraws, "draws.\n")
    }
  }
  CellTable <- data.frame(MissingCell = seq_len(NMissing), SubjectId = Model$Data$SubjectId[Cells[, 1]],
                          OriginalSubjectId = Model$Data$OriginalSubjectId[Cells[, 1]],
                          Treatment = Model$Data$treatment[Cells[, 1]], Visit = Cells[, 2],
                          Month = Model$VisitMonths[Cells[, 2]],
                          MissingType = ifelse(Model$PostTrigger[Cells], "PostTrigger", "Intermittent"),
                          BaselineBdi = Model$Data$bdi.pre[Cells[, 1]])
  Common <- which(Model$Treatment[Cells[, 1]] == 0 | !Model$PostTrigger[Cells])
  for (Rule in RuleNames) {
    if (!identical(Imputed[[Rule]][Common, , drop = FALSE], Imputed[[1]][Common, , drop = FALSE])) {
      stop("Reference-arm or intermittent imputations differ across rules.")
    }
  }
  Bundle <- list(Version = "DKMedianGaussianCAR1v1", FitSignature = Fit$Signature,
                 Data = Model$Data, ObservedScores = Model$Y, MissingCells = Cells,
                 MissingCellTable = CellTable, ImputedValues = Imputed, DrawIndex = DrawIndex,
                 VisitVariables = Model$VisitVariables, VisitMonths = Model$VisitMonths,
                 ScoreMaximum = Model$ScoreMaximum, Tau = Model$Tau, Rules = Settings$Rules,
                 BaselineCenter = Model$BaselineCenter, BaselineScale = Model$BaselineScale,
                 TriggerTable = Model$TriggerTable,
                 Note = "Each column is one joint posterior draw. Observed scores are retained exactly.")
  File <- file.path(CompletedDirectory, "BtheB_DK_CompletedDatasets.rds")
  SaveRdsSafely(Bundle, File)
  write.csv(DrawIndex, file.path(CompletedDirectory, "ImputationDrawIndex.csv"), row.names = FALSE)
  write.csv(CellTable, file.path(SummaryDirectory, "ImputedCells.csv"), row.names = FALSE)
  SummaryRows <- list()
  Counter <- 0
  if (NMissing > 0) {
    for (Rule in RuleNames) {
      Summary <- as.data.frame(t(apply(Imputed[[Rule]], 1, DescribeDkDraws)))
      Counter <- Counter + 1
      SummaryRows[[Counter]] <- cbind(Rule = Rule, CellTable, Summary,
                                      PZero = rowMeans(Imputed[[Rule]] == 0),
                                      PCeiling = rowMeans(Imputed[[Rule]] == Model$ScoreMaximum),
                                      NDraws = Fit$TotalDraws)
    }
    write.csv(do.call(rbind, SummaryRows), file.path(SummaryDirectory, "SubjectImputationSummary.csv"), row.names = FALSE)
  }
  GroupSummary <- SummariseDkCompletions(Bundle, SummaryDirectory)
  if (Settings$MakePlots) PlotDkCompletions(GroupSummary, Model, PlotDirectory)
  if (Settings$ExportCompletedCsv) ExportDkCompletedCsv(Bundle, CompletedDirectory)
  list(File = normalizePath(File, winslash = "/"), Summary = GroupSummary,
       NImputationsPerRule = Fit$TotalDraws, NMissingCells = NMissing)
}

GetCompletedDkData <- function(Bundle, Rule = "MAR", Imputation = 1, Long = FALSE) {
  # Bundle can be the saved RDS filename or the result of readRDS().
  if (is.character(Bundle) && length(Bundle) == 1) Bundle <- readRDS(Bundle)
  if (!Rule %in% names(Bundle$ImputedValues)) stop("Unknown imputation rule: ", Rule)
  if (length(Imputation) != 1 || !is.finite(Imputation) || Imputation != floor(Imputation) ||
      Imputation < 1 || Imputation > nrow(Bundle$DrawIndex)) stop("Invalid imputation index.")
  Scores <- Bundle$ObservedScores
  Scores[Bundle$MissingCells] <- Bundle$ImputedValues[[Rule]][, Imputation]
  Data <- Bundle$Data
  for (Visit in seq_along(Bundle$VisitVariables)) Data[[Bundle$VisitVariables[Visit]]] <- Scores[, Visit]
  Data$Rule <- Rule
  Data$Imputation <- Imputation
  Data$Chain <- Bundle$DrawIndex$Chain[Imputation]
  Data$PostWarmupIteration <- Bundle$DrawIndex$PostWarmupIteration[Imputation]
  if (!Long) return(Data)
  NSubjects <- nrow(Data)
  NVisits <- ncol(Scores)
  Subjects <- rep(seq_len(NSubjects), times = NVisits)
  Visits <- rep(seq_len(NVisits), each = NSubjects)
  Result <- Data[Subjects, setdiff(names(Data), Bundle$VisitVariables), drop = FALSE]
  Result$Visit <- Visits
  Result$Month <- Bundle$VisitMonths[Visits]
  Result$Bdi <- as.vector(Scores)
  Result$Imputed <- as.vector(is.na(Bundle$ObservedScores))
  Result$PostTrigger <- Visits >= Bundle$TriggerTable$TriggerVisit[Subjects]
  rownames(Result) <- NULL
  Result
}

ExportDkCompletedCsv <- function(Bundle, Directory) {
  for (Rule in names(Bundle$ImputedValues)) {
    Path <- file.path(Directory, paste0("CompletedDatasets_", Rule, "_Wide.csv"))
    Connection <- file(Path, open = "wt")
    tryCatch({
      for (Draw in seq_len(nrow(Bundle$DrawIndex))) {
        write.table(GetCompletedDkData(Bundle, Rule, Draw), Connection, sep = ",", row.names = FALSE,
                    col.names = Draw == 1, quote = TRUE)
      }
    }, finally = close(Connection))
  }
  invisible(NULL)
}

SummariseDkCompletions <- function(Bundle, SummaryDirectory) {
  Rows <- list()
  Counter <- 0
  for (Rule in names(Bundle$ImputedValues)) {
    for (Treatment in c("TAU", "BtheB")) {
      Subjects <- which(as.character(Bundle$Data$treatment) == Treatment)
      for (Visit in seq_along(Bundle$VisitMonths)) {
        Observed <- Bundle$ObservedScores[Subjects, Visit]
        CellRows <- which(Bundle$MissingCells[, 1] %in% Subjects & Bundle$MissingCells[, 2] == Visit)
        Values <- Bundle$ImputedValues[[Rule]][CellRows, , drop = FALSE]
        CompletedMean <- (sum(Observed, na.rm = TRUE) + colSums(Values))/length(Subjects)
        CompletedZero <- (sum(Observed == 0, na.rm = TRUE) + colSums(Values == 0))/length(Subjects)
        CompletedCeiling <- (sum(Observed == Bundle$ScoreMaximum, na.rm = TRUE) +
                               colSums(Values == Bundle$ScoreMaximum))/length(Subjects)
        for (Quantity in c("CompletedMean", "CompletedZero", "CompletedCeiling")) {
          Counter <- Counter + 1
          Rows[[Counter]] <- data.frame(Rule = Rule, Treatment = Treatment, Month = Bundle$VisitMonths[Visit],
                                        NParticipants = length(Subjects), NImputedParticipants = length(CellRows),
                                        Quantity = Quantity, t(DescribeDkDraws(get(Quantity))))
        }
      }
    }
  }
  Summary <- do.call(rbind, Rows)
  write.csv(Summary, file.path(SummaryDirectory, "CompletedScoreProfiles.csv"), row.names = FALSE)
  Summary
}

PlotDkCompletions <- function(Summary, Model, PlotDirectory) {
  RuleNames <- unique(Summary$Rule)
  Colors <- rep_len(c("#0072B2", "#E69F00", "#D55E00", "#009E73"), length(RuleNames))
  Symbols <- rep_len(c(16, 17, 15, 18), length(RuleNames))
  pdf(file.path(PlotDirectory, "CompletedScoreProfiles.pdf"), width = 10, height = 8, family = "Helvetica")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(3, 2), mar = c(3, 3.5, 2, 0.5), mgp = c(1.8, 0.6, 0))
  Labels <- c(CompletedMean = "Completed mean BDI-II", CompletedZero = "Proportion at zero",
              CompletedCeiling = "Proportion at 63")
  for (Quantity in names(Labels)) {
    for (Treatment in c("TAU", "BtheB")) {
      Data <- Summary[Summary$Treatment == Treatment & Summary$Quantity == Quantity, ]
      Limits <- range(c(Data$Q025, Data$Q975))
      if (diff(Limits) == 0) Limits <- Limits + c(-0.01, 0.01)
      plot(range(Model$VisitMonths), Limits, type = "n", xlab = "Month", ylab = Labels[[Quantity]],
           main = Treatment, xaxt = "n")
      axis(1, at = Model$VisitMonths)
      for (Index in seq_along(RuleNames)) {
        Values <- Data[Data$Rule == RuleNames[Index], ]
        Offset <- if (Treatment == "TAU") 0 else (Index - (length(RuleNames) + 1)/2)*0.07
        X <- Values$Month + Offset
        lines(X, Values$Mean, col = Colors[Index], lty = Index)
        segments(X, Values$Q025, X, Values$Q975, col = Colors[Index])
        points(X, Values$Mean, col = Colors[Index], pch = Symbols[Index])
      }
      if (Quantity == names(Labels)[1] && Treatment == "BtheB") {
        legend("topright", legend = RuleNames, col = Colors, pch = Symbols,
               lty = seq_along(RuleNames), bty = "n", cex = 0.7)
      }
    }
  }
}

PredictDkScores <- function(Parameter, Model) {
  Correlation <- MakeContinuousAr1(Model$VisitMonths, plogis(Parameter[Model$Index$LogitRho]))
  Z <- matrix(rnorm(Model$NSubjects*Model$NVisits), Model$NSubjects, Model$NVisits)%*%chol(Correlation)
  Scores <- matrix(NA_integer_, Model$NSubjects, Model$NVisits)
  for (Visit in seq_len(Model$NVisits)) {
    Eta <- Parameter[Model$Index$Intercept[Visit]] +
      Parameter[Model$Index$Treatment[Visit]]*Model$Treatment +
      Parameter[Model$Index$Baseline[Visit]]*Model$Baseline
    Scores[, Visit] <- as.integer(DkQuantileFromLogSurvival(pnorm(Z[, Visit], lower.tail = FALSE, log.p = TRUE),
                                                            Eta, Parameter[Model$Index$LogShape[Visit]],
                                                            Model$ScoreMaximum, Model$Tau))
  }
  Scores
}

RunDkPredictiveChecks <- function(Fit, Model, Settings, SummaryDirectory, PlotDirectory, ModelDirectory) {
  set.seed(Settings$Seed + 200000)
  NObserved <- sum(Model$Observed)
  Simulated <- matrix(NA_integer_, NObserved, Fit$TotalDraws)
  Draw <- 0
  for (Chain in Fit$Chains) {
    for (Iteration in seq_len(nrow(Chain$Draws))) {
      Draw <- Draw + 1
      Scores <- PredictDkScores(Chain$Draws[Iteration, ], Model)
      Simulated[, Draw] <- Scores[Model$Observed]
    }
  }
  ObservedValues <- Model$Y[Model$Observed]
  Cells <- Model$ObservedCells
  Groups <- list(list(Treatment = "All", Month = NA_real_, Rows = seq_len(NObserved)))
  for (Treatment in 0:1) {
    for (Visit in seq_len(Model$NVisits)) {
      Rows <- which(Model$Treatment[Cells[, 1]] == Treatment & Cells[, 2] == Visit)
      if (length(Rows) == 0) next
      Groups[[length(Groups) + 1]] <- list(Treatment = c("TAU", "BtheB")[Treatment + 1],
                                           Month = Model$VisitMonths[Visit], Rows = Rows)
    }
  }
  FrequencyRows <- list()
  StatisticRows <- list()
  for (Group in Groups) {
    Values <- Simulated[Group$Rows, , drop = FALSE]
    Recorded <- ObservedValues[Group$Rows]
    Counts <- vapply(seq_len(Fit$TotalDraws), function(Draw) {
      tabulate(Values[, Draw] + 1, nbins = Model$ScoreMaximum + 1)
    }, integer(Model$ScoreMaximum + 1))
    Frequencies <- Counts/length(Group$Rows)
    FrequencySummary <- as.data.frame(t(apply(Frequencies, 1, DescribeDkDraws)))
    FrequencyRows[[length(FrequencyRows) + 1]] <- cbind(
      Treatment = Group$Treatment, Month = Group$Month, Score = 0:Model$ScoreMaximum,
      Observed = tabulate(Recorded + 1, nbins = Model$ScoreMaximum + 1)/length(Recorded),
      FrequencySummary, NObserved = length(Recorded), NDraws = Fit$TotalDraws)
    Predictions <- list(Mean = colMeans(Values), Variance = apply(Values, 2, var),
                        PZero = colMeans(Values == 0), PCeiling = colMeans(Values == Model$ScoreMaximum))
    ObservedStatistics <- c(Mean = mean(Recorded), Variance = var(Recorded),
                            PZero = mean(Recorded == 0), PCeiling = mean(Recorded == Model$ScoreMaximum))
    for (Name in names(Predictions)) {
      StatisticRows[[length(StatisticRows) + 1]] <- data.frame(
        Treatment = Group$Treatment, Month = Group$Month, Quantity = Name,
        Observed = ObservedStatistics[[Name]], t(DescribeDkDraws(Predictions[[Name]])),
        NObserved = length(Recorded), NDraws = Fit$TotalDraws)
    }
  }
  Frequencies <- do.call(rbind, FrequencyRows)
  Statistics <- do.call(rbind, StatisticRows)
  write.csv(Frequencies, file.path(SummaryDirectory, "PredictiveScoreFrequencies.csv"), row.names = FALSE)
  write.csv(Statistics, file.path(SummaryDirectory, "PredictiveSummary.csv"), row.names = FALSE)
  # These randomized predictive ranks are descriptive. No independent-uniform
  # goodness-of-fit p-value is claimed for dependent posterior predictive ranks.
  RankLower <- rowSums(Simulated < ObservedValues)/(Fit$TotalDraws + 1)
  RankWidth <- (rowSums(Simulated == ObservedValues) + 1)/(Fit$TotalDraws + 1)
  Ranks <- RankLower + runif(NObserved)*RankWidth
  Residuals <- data.frame(SubjectId = Model$Data$SubjectId[Cells[, 1]],
                          OriginalSubjectId = Model$Data$OriginalSubjectId[Cells[, 1]],
                          Treatment = Model$Data$treatment[Cells[, 1]], Month = Model$VisitMonths[Cells[, 2]],
                          ObservedScore = ObservedValues, PredictiveRank = Ranks,
                          PredictiveMean = rowMeans(Simulated),
                          RankLower = RankLower, RankUpper = RankLower + RankWidth)
  write.csv(Residuals, file.path(SummaryDirectory, "PredictiveRanks.csv"), row.names = FALSE)
  if (Settings$SavePredictiveDraws) {
    SaveRdsSafely(list(ObservedCells = Cells, SimulatedScores = Simulated, DrawIndex = BuildDkDrawIndex(Fit)),
                  file.path(ModelDirectory, "PosteriorPredictiveScores.rds"))
  }
  if (Settings$MakePlots) PlotDkPredictiveChecks(Frequencies, Residuals, Model, PlotDirectory)
  list(Statistics = Statistics, Frequencies = Frequencies, Ranks = Residuals, NDraws = Fit$TotalDraws)
}

PlotDkPredictiveChecks <- function(Frequencies, Residuals, Model, PlotDirectory) {
  pdf(file.path(PlotDirectory, "PredictiveScoreFrequencies.pdf"), width = 10, height = 7, family = "Helvetica")
  tryCatch({
    par(mfrow = c(2, 2), mar = c(3.5, 4, 2.5, 0.8), mgp = c(2.2, 0.7, 0),
        cex = 1.1, cex.axis = 1.05, cex.lab = 1.1, cex.main = 1.1)
    Groups <- unique(Frequencies[, c("Treatment", "Month")])
    for (Index in seq_len(nrow(Groups))) {
      Match <- Frequencies$Treatment == Groups$Treatment[Index] &
        (if (is.na(Groups$Month[Index])) is.na(Frequencies$Month) else
          !is.na(Frequencies$Month) & Frequencies$Month == Groups$Month[Index])
      Data <- Frequencies[Match, ]
      Title <- if (is.na(Groups$Month[Index])) "All observed scores" else
        paste(Groups$Treatment[Index], "- month", Groups$Month[Index])
      plot(Data$Score, Data$Observed, type = "h", col = "#333333", lwd = 1,
           xlab = "BDI-II score", ylab = "Proportion", main = Title,
           ylim = c(0, max(c(Data$Observed, Data$Q975))*1.1 + 0.001))
      segments(Data$Score + 0.2, Data$Q025, Data$Score + 0.2, Data$Q975, col = "#0072B2")
      points(Data$Score + 0.2, Data$Median, col = "#0072B2", pch = 16, cex = 0.65)
      if (Index == 1) legend("topright", legend = c("Observed", "Predictive median and 95% interval"),
                             col = c("#333333", "#0072B2"), pch = c(NA, 16), lty = c(1, 1),
                             bty = "n", cex = 0.85)
    }
  }, finally = dev.off())
  pdf(file.path(PlotDirectory, "PredictiveRanks.pdf"), width = 9, height = 4.8, family = "Helvetica")
  tryCatch({
    layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.16))
    par(mar = c(4.5, 4.5, 2.5, 1), mgp = c(2.7, 0.8, 0),
        cex = 1.1, cex.axis = 1.05, cex.lab = 1.1, cex.main = 1.1)
    hist(Residuals$PredictiveRank, breaks = seq(0, 1, length.out = 11),
         col = "#0072B2", border = "white", main = "Randomized predictive ranks",
         xlab = "Rank", axes = FALSE)
    axis(1)
    axis(2)
    box(bty = "l")
    plot(Residuals$PredictiveMean, Residuals$PredictiveRank,
         col = ifelse(Residuals$Treatment == "TAU", "#0072B2", "#D55E00"),
         pch = ifelse(Residuals$Treatment == "TAU", 16, 17), cex = 0.9,
         xlab = "Predictive mean BDI-II", ylab = "Predictive rank",
         ylim = c(0, 1), axes = FALSE)
    axis(1)
    axis(2)
    box(bty = "l")
    abline(h = 0.5, lty = 2, col = "grey50")
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = c("TAU", "BtheB"), col = c("#0072B2", "#D55E00"),
           pch = c(16, 17), bty = "n", horiz = TRUE, cex = 1.05)
  }, finally = dev.off())
}

WriteDkRunNotes <- function(Model, Settings, Fit, Diagnostics, Directory) {
  Notes <- c(
    "BtheB DK reference-based imputation",
    paste("Included participants:", Model$NSubjects),
    paste("Excluded with no post-baseline scores:", nrow(Model$Dropped)),
    paste("Observed visits are monotone for every participant:", all(Model$TriggerTable$Monotone)),
    paste("Post-baseline months:", paste(Model$VisitMonths, collapse = ", ")),
    "Margins: discrete Kumaraswamy on the integers 0,...,63, parameterized at Tau = 0.5.",
    "The relaxed median is 64*expit(Eta) - 1; the integer median is its ceiling.",
    "Eta has visit-specific intercept, treatment and baseline coefficients.",
    paste("Baseline is centered at", signif(Model$BaselineCenter, 8), "and divided by", Model$BaselineScale, "."),
    "Shape varies by visit and is shared by treatment groups.",
    "The copula is Gaussian: R[j,k] = RhoPerMonth^abs(Month[j] - Month[k]).",
    "RhoPerMonth is positive and less than one; all margins have standard normal latent scores.",
    "The dependence matrix is shared by arms and retained across sensitivity rules within each draw.",
    "The trigger is the first visit after the last observed post-baseline score.",
    "Under monotone missingness this is exactly the first missing visit.",
    "Intermittent gaps before the trigger use assigned-arm continuation under every rule.",
    "Active-arm post-trigger margins use ReferenceEta + (1 - Omega)*TreatmentCoefficient.",
    paste("Rules:", paste(paste(names(Settings$Rules), Settings$Rules, sep = "="), collapse = ", ")),
    "Half means half of the treatment coefficient retained on the logit scaled-median-anchor scale.",
    "It does not mean half of the observed mean difference or half of the integer median difference.",
    "Reference-arm imputations are identical across rules, using common latent random draws.",
    "Shape is unchanged across rules within a posterior draw, so J2R uses the reference marginal distribution.",
    "Conditioning on the retained history still affects individual imputation distributions.",
    paste("Prior SDs (intercept, treatment, scaled baseline):", paste(c(Settings$InterceptPriorSd,
                                                                        Settings$TreatmentPriorSd, Settings$BaselinePriorSd), collapse = ", ")),
    paste("Log-shape prior: Normal mean", Settings$LogShapePriorMean, "SD", Settings$LogShapePriorSd),
    "Correlation prior: rho ~ Uniform(0, 1), where rho is the one-month latent correlation.",
    "The logit-rho target includes the transformation Jacobian; the shape prior is specified directly on log shape.",
    "Metropolis-within-Gibbs with randomized PIT augmentation; proposal adaptation stops after warmup.",
    "PIT interval probabilities use log CDF/survival calculations without a positive probability floor.",
    "The auxiliary PIT variables have a uniform base measure; the copula weights their posterior distribution.",
    "The logit-PIT proposal includes its Jacobian in the MH target.",
    paste("Chains:", Settings$NChains, "Warmup per chain:", Settings$Warmup,
          "Post-warmup draws per chain:", Settings$Sample),
    paste("All", Fit$TotalDraws, "post-warmup draws are used; no thinning or draw selection is performed."),
    paste("Sampling diagnostic review flag:", Diagnostics$DiagnosticFlag),
    "Diagnostic flags do not remove chains or draws and do not automatically prevent output creation.",
    "PPC replicates use the original covariates and recorded-score observation mask, with Gaussian copula dependence.",
    "Randomized predictive ranks are descriptive; no independent-uniform test is claimed.",
    "CompletedScoreProfiles contains descriptive summaries of completed cohorts, not treatment-effect analyses.",
    "Completed datasets are stored compactly as original data plus a matrix of missing scores for each rule.",
    "Columns correspond exactly to ImputationDrawIndex.csv, retaining chain and iteration identifiers.",
    "Use GetCompletedDkData() to reconstruct any completed dataset in wide or long form.",
    "No completed-data model fitting or Rubin pooling is performed.",
    paste("Fit signature:", Fit$Signature)
  )
  writeLines(Notes, file.path(Directory, "RunNotes.txt"))
  writeLines(capture.output(sessionInfo()), file.path(Directory, "SessionInfo.txt"))
}

ExportDkPosteriorTable <- function(
    Result = file.path("BtheB_DK_RBMI", "ModelOutputs",
                       "BtheB_DK_ApplicationResult.rds"),
    OutputDirectory = NULL, ChainDirectory = NULL,
    Digits = 2, AnchorDigits = 1, RhoDigits = 3) {
  
  ResultPath <- NULL
  if (is.character(Result) && length(Result) == 1 && !is.na(Result)) {
    if (!file.exists(Result)) stop("The saved application result was not found: ", Result)
    ResultPath <- normalizePath(Result, winslash = "/", mustWork = TRUE)
    Result <- readRDS(ResultPath)
  }
  Required <- c("Settings", "Data", "ChainFiles", "FitSignature", "TotalPostWarmupDraws")
  if (!is.list(Result) || !all(Required %in% names(Result))) {
    stop("Supply BthebDkResult or the path to BtheB_DK_ApplicationResult.rds.")
  }
  Settings <- Result$Settings
  Months <- Settings$VisitMonths
  NVisits <- length(Months)
  if (!is.numeric(Months) || NVisits < 1 || any(!is.finite(Months)) ||
      any(diff(Months) <= 0)) stop("Invalid visit times in the saved settings.")
  if (!isTRUE(Settings$Tau == 0.5)) stop("This table requires the median parameterization.")
  if (length(Settings$ScoreMaximum) != 1 || !is.finite(Settings$ScoreMaximum) ||
      Settings$ScoreMaximum < 2 || Settings$ScoreMaximum != floor(Settings$ScoreMaximum)) {
    stop("Invalid score maximum in the saved settings.")
  }
  if (length(Settings$BaselineScale) != 1 || !is.finite(Settings$BaselineScale) ||
      Settings$BaselineScale <= 0) stop("Invalid baseline scale in the saved settings.")
  if (is.null(Result$FitSignature)) stop("The fit signature is missing.")
  for (Precision in list(Digits, AnchorDigits, RhoDigits)) {
    if (!is.numeric(Precision) || length(Precision) != 1 ||
        !is.finite(Precision) || Precision < 0 || Precision > 6 ||
        Precision != floor(Precision)) stop("Digits must be integers from 0 to 6.")
  }
  Baseline <- Result$Data$bdi.pre
  if (!is.numeric(Baseline) || length(Baseline) == 0 ||
      any(!is.finite(Baseline))) stop("The included participants' baseline scores are missing.")
  BaselineCenter <- mean(Baseline)
  
  # Prefer chains beside a result loaded from disk, so a moved output folder works.
  SavedPaths <- Result$ChainFiles
  if (!is.character(SavedPaths) || length(SavedPaths) != Settings$NChains ||
      anyNA(SavedPaths)) stop("The saved chain list does not match the number of chains.")
  if (!is.null(ChainDirectory)) {
    ChainFiles <- file.path(ChainDirectory, basename(SavedPaths))
  } else {
    Candidates <- list()
    if (!is.null(ResultPath)) {
      Candidates <- c(Candidates, list(file.path(dirname(ResultPath), basename(SavedPaths))))
    }
    Candidates <- c(Candidates, list(
      SavedPaths,
      file.path(Settings$OutputRoot, "ModelOutputs", basename(SavedPaths))
    ))
    Available <- vapply(Candidates, function(Paths) all(file.exists(Paths)), logical(1))
    if (!any(Available)) {
      stop("The saved chains were not found. Set ChainDirectory to their folder.")
    }
    ChainFiles <- Candidates[[which(Available)[1]]]
  }
  if (!all(file.exists(ChainFiles))) stop("Some chain files are missing from ChainDirectory.")
  ChainFiles <- normalizePath(ChainFiles, winslash = "/", mustWork = TRUE)
  if (anyDuplicated(ChainFiles)) stop("The chain list contains a duplicate file.")
  
  ParameterNames <- c(
    paste0("Intercept[", seq_len(NVisits), "]"),
    paste0("Treatment[", seq_len(NVisits), "]"),
    paste0("Baseline[", seq_len(NVisits), "]"),
    paste0("LogShape[", seq_len(NVisits), "]"),
    "LogitRho"
  )
  ChainDraws <- vector("list", length(ChainFiles))
  ChainIds <- numeric(length(ChainFiles))
  for (Index in seq_along(ChainFiles)) {
    Chain <- readRDS(ChainFiles[Index])
    if (!identical(Chain$Signature, Result$FitSignature)) {
      stop("A chain belongs to a different fit: ", ChainFiles[Index])
    }
    if (!is.matrix(Chain$Draws) || !is.numeric(Chain$Draws) ||
        !identical(colnames(Chain$Draws), ParameterNames)) {
      stop("Unexpected parameter columns in: ", ChainFiles[Index])
    }
    if (nrow(Chain$Draws) != Settings$Sample || any(!is.finite(Chain$Draws))) {
      stop("A chain has an unexpected draw count or nonfinite values: ", ChainFiles[Index])
    }
    if (length(Chain$Chain) != 1 || !is.numeric(Chain$Chain) ||
        !is.finite(Chain$Chain)) stop("A chain identifier is missing or invalid.")
    ChainIds[Index] <- Chain$Chain
    ChainDraws[[Index]] <- Chain$Draws
  }
  if (anyDuplicated(ChainIds) || !setequal(ChainIds, seq_len(Settings$NChains))) {
    stop("The saved files do not contain exactly the expected chains.")
  }
  Draws <- do.call(rbind, ChainDraws)
  NDraws <- nrow(Draws)
  if (NDraws < 2 || NDraws != Result$TotalPostWarmupDraws) {
    stop("The total draw count does not match the saved application result.")
  }
  rm(Chain, ChainDraws)
  
  Summarise <- function(Values, Panel, Quantity, Parameter, Description,
                        Visit = NA_real_, Arm = "Shared", Precision = Digits) {
    if (any(!is.finite(Values))) stop("Nonfinite derived values for ", Quantity, ".")
    Bounds <- as.numeric(quantile(Values, c(0.025, 0.975), names = FALSE, type = 8))
    data.frame(
      Panel = Panel, Quantity = Quantity, Parameter = Parameter,
      Description = Description, Visit = Visit,
      Month = if (is.na(Visit)) NA_real_ else Months[Visit], Arm = Arm,
      Mean = mean(Values), Lower95 = Bounds[1], Upper95 = Bounds[2],
      NDraws = length(Values), Digits = Precision, stringsAsFactors = FALSE
    )
  }
  MonthLabel <- function(Visit) format(Months[Visit], trim = TRUE, scientific = FALSE)
  Rows <- list()
  AddRow <- function(Row) Rows[[length(Rows) + 1]] <<- Row
  
  # Transform within each draw; in particular, preserve the alpha/beta covariance.
  for (Quantity in c("Intercept", "MedianAnchor")) {
    for (Visit in seq_len(NVisits)) {
      Alpha <- Draws[, paste0("Intercept[", Visit, "]")]
      Beta <- Draws[, paste0("Treatment[", Visit, "]")]
      for (Arm in c("TAU", "BtheB")) {
        Eta <- Alpha + (if (Arm == "BtheB") Beta else 0)
        if (Quantity == "Intercept") {
          Values <- Eta
          Symbol <- sprintf("$\\alpha_{%s} + \\beta_{%s} A_i$", Visit, Visit)
          Description <- paste0("Intercept, month ", MonthLabel(Visit))
          Precision <- Digits
        } else {
          Values <- (Settings$ScoreMaximum + 1)*plogis(Eta) - 1
          Symbol <- sprintf("$\\xi_{i%s}$", Visit)
          Description <- paste0("Median anchor, month ", MonthLabel(Visit))
          Precision <- AnchorDigits
        }
        AddRow(Summarise(Values, "ArmSpecific", Quantity, Symbol, Description,
                         Visit, Arm, Precision))
      }
    }
  }
  for (Quantity in c("Treatment", "Baseline", "Shape")) {
    for (Visit in seq_len(NVisits)) {
      if (Quantity == "Treatment") {
        Values <- Draws[, paste0("Treatment[", Visit, "]")]
        Symbol <- sprintf("$\\beta_{%s}$", Visit)
        Description <- paste0("Treatment coefficient, month ", MonthLabel(Visit))
        Arm <- "BtheB minus TAU"
      } else if (Quantity == "Baseline") {
        Values <- Draws[, paste0("Baseline[", Visit, "]")]
        Symbol <- sprintf("$\\gamma_{%s1}$", Visit)
        Description <- paste0("Baseline coefficient, month ", MonthLabel(Visit))
        Arm <- "Shared"
      } else {
        Values <- exp(Draws[, paste0("LogShape[", Visit, "]")])
        Symbol <- sprintf("$\\phi_{%s}$", Visit)
        Description <- paste0("Shape, month ", MonthLabel(Visit))
        Arm <- "Shared"
      }
      AddRow(Summarise(Values, "Common", Quantity, Symbol, Description, Visit, Arm))
    }
  }
  AddRow(Summarise(plogis(Draws[, "LogitRho"]), "Common", "Correlation",
                   "$\\rho$", "One-month latent correlation", Precision = RhoDigits))
  Summary <- do.call(rbind, Rows)
  rownames(Summary) <- NULL
  
  Number <- function(Value, Precision) {
    # Avoid displaying "-0.00" after rounding.
    if (round(Value, Precision) == 0) Value <- 0
    formatC(Value, format = "f", digits = Precision, decimal.mark = ".")
  }
  ArmSummary <- Summary[Summary$Panel == "ArmSpecific", ]
  CommonRows <- Summary[Summary$Panel == "Common", ]
  ArmFormats <- vapply(c("Mean", "Lower95", "Upper95"), function(Name) {
    DkTableNumberFormat(ArmSummary[[Name]], ArmSummary$Digits)
  }, character(1))
  CommonFormats <- vapply(c("Mean", "Lower95", "Upper95"), function(Name) {
    DkTableNumberFormat(CommonRows[[Name]], CommonRows$Digits)
  }, character(1))
  Estimate <- function(Row) Number(Row$Mean, Row$Digits)
  Lower <- function(Row) paste0("{[}", Number(Row$Lower95, Row$Digits), "{,}")
  Upper <- function(Row) paste0(Number(Row$Upper95, Row$Digits), "{]}")
  # Fixed-width numeric boxes preserve decimal alignment within spanning cells.
  CommonEstimate <- function(Row) paste0(
    "\\tablenum[table-format = ", CommonFormats["Mean"], "]{", Estimate(Row), "}"
  )
  CommonInterval <- function(Row) paste0(
    "\\tablenum[table-format = ", CommonFormats["Lower95"],
    ", table-space-text-pre = {[}, table-space-text-post = {,}]{", Lower(Row), "}",
    "\\;\\tablenum[table-format = ", CommonFormats["Upper95"],
    ", table-space-text-post = {]}]{", Upper(Row), "}"
  )
  ArmColumns <- c(
    paste0("S[table-format = ", ArmFormats["Mean"], "]"),
    paste0("S[table-format = ", ArmFormats["Lower95"],
           ", table-space-text-pre = {[}, table-space-text-post = {,}] @{\\;}"),
    paste0("S[table-format = ", ArmFormats["Upper95"], ", table-space-text-post = {]}]")
  )
  Lines <- c(
    "% Generated from all saved post-warmup draws; do not edit numerical cells by hand.",
    "% Requires \\usepackage{booktabs,tabularx,siunitx}.",
    "\\begin{table}[tbp]",
    "\\caption{Posterior summaries of the DK Gaussian copula imputation model.}",
    "\\label{tab:dk-parameters}",
    "\\centering",
    "\\begingroup",
    "\\small",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\sisetup{table-align-text-pre = false, table-align-text-post = false,",
    "  group-digits = false, output-decimal-marker = {.}, round-mode = none}",
    "\\begin{tabularx}{\\linewidth}{@{}l>{\\raggedright\\arraybackslash}X",
    ArmColumns, ArmColumns, "@{}}",
    "\\toprule",
    "\\textbf{Parameter} & \\textbf{Description} & \\multicolumn{3}{c}{\\textbf{TAU}} & \\multicolumn{3}{c}{\\textbf{BtheB}} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(l){6-8}",
    "& & {\\textbf{Mean}} & \\multicolumn{2}{c}{\\textbf{95\\% CrI}} & {\\textbf{Mean}} & \\multicolumn{2}{c}{\\textbf{95\\% CrI}} \\\\",
    "\\midrule",
    "\\multicolumn{8}{@{}l}{\\textbf{\\textit{Arm-specific summaries at the centered baseline value}}} \\\\",
    "\\addlinespace[2pt]"
  )
  ArmRows <- ArmSummary[ArmSummary$Arm == "TAU", ]
  for (Index in seq_len(nrow(ArmRows))) {
    Reference <- ArmRows[Index, ]
    Active <- Summary[Summary$Panel == "ArmSpecific" & Summary$Arm == "BtheB" &
                        Summary$Quantity == Reference$Quantity &
                        !is.na(Summary$Visit) & Summary$Visit == Reference$Visit, ]
    if (Index > 1 && Reference$Quantity != ArmRows$Quantity[Index - 1]) {
      Lines <- c(Lines, "\\addlinespace[3pt]")
    }
    Lines <- c(Lines, paste0(
      Reference$Parameter, " & ", Reference$Description, " & ",
      Estimate(Reference), " & ", Lower(Reference), " & ", Upper(Reference), " & ",
      Estimate(Active), " & ", Lower(Active), " & ", Upper(Active), " \\\\"
    ))
  }
  Lines <- c(
    Lines, "\\midrule",
    "\\multicolumn{8}{@{}l}{\\textbf{\\textit{Treatment contrasts and shared parameters}}} \\\\",
    "\\addlinespace[2pt]",
    "\\textbf{Parameter} & \\textbf{Description} & \\multicolumn{3}{c}{\\textbf{Mean}} & \\multicolumn{3}{c}{\\textbf{95\\% CrI}} \\\\",
    "\\cmidrule(l){1-8}"
  )
  for (Index in seq_len(nrow(CommonRows))) {
    Row <- CommonRows[Index, ]
    if (Index > 1 && Row$Quantity != CommonRows$Quantity[Index - 1]) {
      Lines <- c(Lines, "\\addlinespace[3pt]")
    }
    Lines <- c(Lines, paste0(
      Row$Parameter, " & ", Row$Description, " & \\multicolumn{3}{c}{",
      CommonEstimate(Row), "} & \\multicolumn{3}{c}{", CommonInterval(Row), "} \\\\"
    ))
  }
  DrawCount <- format(NDraws, big.mark = ",", scientific = FALSE, trim = TRUE)
  Notes <- paste0(
    "Entries are posterior means and 95\\% equal-tailed credible intervals (CrI), ",
    "using all ", DrawCount, " post-warmup draws. ",
    "Visit indices $j=1,\\ldots,", NVisits, "$ correspond to months ",
    paste(vapply(seq_len(NVisits), MonthLabel, character(1)), collapse = ", "), ". ",
    "Arm-specific summaries use $B_i=0$, corresponding to a baseline BDI-II score of ",
    Number(BaselineCenter, 4), ". ",
    "The intercept is $\\alpha_j$ for TAU and $\\alpha_j+\\beta_j$ for BtheB. ",
    "Median anchors $\\xi_{ij}$ are in score units, before conditioning on retained ",
    "follow-up outcomes; they are not the discrete medians or covariate-standardized means. ",
    "Regression coefficients refer to the logit scale of the normalized median anchor. ",
    "The baseline coefficient $\\gamma_{j1}$ refers to a ",
    format(Settings$BaselineScale, trim = TRUE, scientific = FALSE),
    "-point increase in baseline score. ",
    "Baseline coefficients, shapes $\\phi_j$, and $\\rho$ are shared across arms. ",
    "All imputation rules use this same fitted posterior."
  )
  Lines <- c(
    Lines, "\\bottomrule", "\\end{tabularx}", "\\endgroup",
    "\\par\\smallskip",
    "\\begin{minipage}{\\linewidth}", "\\footnotesize\\raggedright",
    Notes, "\\end{minipage}", "\\end{table}"
  )
  
  if (is.null(OutputDirectory)) {
    OutputRoot <- if (is.null(ResultPath)) Settings$OutputRoot else dirname(dirname(ResultPath))
    OutputDirectory <- file.path(OutputRoot, "Tables")
  }
  dir.create(OutputDirectory, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OutputDirectory)) stop("The table output directory could not be created.")
  OutputDirectory <- normalizePath(OutputDirectory, winslash = "/", mustWork = TRUE)
  TableFile <- file.path(OutputDirectory, "DkPosteriorParameterTable.tex")
  StandaloneFile <- file.path(OutputDirectory, "DkPosteriorParameterTableStandalone.tex")
  CsvFile <- file.path(OutputDirectory, "DkPosteriorParameterTable.csv")
  writeLines(Lines, TableFile, useBytes = TRUE)
  writeLines(c(
    "\\documentclass[12pt]{article}",
    "\\usepackage[margin=0.9in]{geometry}",
    "\\usepackage[T1]{fontenc}",
    "\\usepackage{helvet}",
    "\\renewcommand{\\familydefault}{\\sfdefault}",
    "\\usepackage{booktabs,tabularx,siunitx}",
    "\\usepackage{caption}",
    "\\captionsetup{justification=raggedright,singlelinecheck=false}",
    "\\captionsetup[table]{name=TABLE,labelfont=bf,labelsep=period}",
    "\\pagestyle{empty}",
    "\\begin{document}", Lines, "\\end{document}"
  ), StandaloneFile, useBytes = TRUE)
  write.csv(Summary, CsvFile, row.names = FALSE, na = "")
  message("Used all ", DrawCount, " post-warmup draws from ", length(ChainFiles), " chains.")
  message("LaTeX table: ", TableFile)
  message("Standalone LaTeX: ", StandaloneFile)
  message("Unrounded summaries: ", CsvFile)
  invisible(list(Table = TableFile, Standalone = StandaloneFile, Csv = CsvFile,
                 Summary = Summary, NDraws = NDraws))
}


RunBthebDkRbmi <- function(Settings = DefaultBthebDkSettings(), DataWide = NULL) {
  CheckSettings(Settings)
  if (!requireNamespace("posterior", quietly = TRUE)) stop("Install posterior before running this program.")
  if (Settings$RunDistributionChecks) CheckDkDistribution()
  Model <- PrepareBthebData(Settings, DataWide)
  Directories <- list(Model = file.path(Settings$OutputRoot, "ModelOutputs"),
                      Summary = file.path(Settings$OutputRoot, "Summaries"),
                      Completed = file.path(Settings$OutputRoot, "CompletedDatasets"),
                      Plot = file.path(Settings$OutputRoot, "Plots"))
  for (Directory in Directories) {
    dir.create(Directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(Directory)) stop("Could not create output directory: ", Directory)
  }
  ExportDkInputs(Model, Settings, Directories$Summary)
  if (Settings$PrintProgress) {
    cat("BtheB population:", Model$NSubjects, "participants;", nrow(Model$Dropped), "excluded.\n")
    cat("Monotone missingness:", all(Model$TriggerTable$Monotone), "\n")
  }
  Fit <- FitDkCopula(Model, Settings, Directories$Model)
  Diagnostics <- SummariseDkFit(Fit, Model, Directories$Summary, Directories$Plot, Settings)
  Imputation <- ImputeDkOutcomes(Fit, Model, Settings, Directories$Completed, Directories$Summary, Directories$Plot)
  Predictive <- if (Settings$RunPredictiveChecks) {
    if (Settings$PrintProgress) cat("Running predictive checks with every post-warmup draw.\n")
    RunDkPredictiveChecks(Fit, Model, Settings, Directories$Summary, Directories$Plot, Directories$Model)
  } else NULL
  WriteDkRunNotes(Model, Settings, Fit, Diagnostics, Directories$Summary)
  Result <- list(Settings = Settings, FitSignature = Fit$Signature, Data = Model$Data,
                 Dropped = Model$Dropped, TriggerTable = Model$TriggerTable,
                 ParameterSummary = Diagnostics$Parameters, SamplerAcceptance = Diagnostics$Acceptance,
                 DiagnosticFlag = Diagnostics$DiagnosticFlag, Imputation = Imputation,
                 Predictive = Predictive, TotalPostWarmupDraws = Fit$TotalDraws,
                 FitRuntimeSeconds = Fit$RuntimeSeconds,
                 ChainFiles = file.path(Directories$Model, paste0("Chain", seq_len(Settings$NChains), ".rds")))
  SaveRdsSafely(Result, file.path(Directories$Model, "BtheB_DK_ApplicationResult.rds"))
  ExportDkPosteriorTable(Result)
  if (Settings$PrintProgress) {
    cat("\nFinished. Completed datasets:\n", Imputation$File, "\n", sep = "")
    cat("Review ParameterPosteriorSummary.csv and the predictive checks before interpreting results.\n")
  }
  invisible(Result)
}

if (isTRUE(getOption("BthebDkRun", TRUE))) {
  BthebDkResult <- RunBthebDkRbmi(Settings)
}

# DHARMa posterior predictive diagnostics.
# Uses saved chains and all post-warmup draws without refitting the model.

PrepareDkDharmaSimulations <- function(Result) {
  Settings <- Result$Settings
  Model <- PrepareBthebData(Settings, Result$Data)
  ChainFiles <- Result$ChainFiles
  if (!all(file.exists(ChainFiles))) {
    ChainFiles <- file.path(Settings$OutputRoot, "ModelOutputs", basename(ChainFiles))
  }
  if (length(ChainFiles) == 0 || !all(file.exists(ChainFiles))) {
    stop("The saved chains could not be found.")
  }
  Chains <- lapply(ChainFiles, function(Path) {
    Chain <- readRDS(Path)
    if (!identical(Chain$Signature, Result$FitSignature)) {
      stop("A saved chain does not match this fit.")
    }
    if (!identical(colnames(Chain$Draws), Model$ParameterNames)) {
      stop("The saved parameter names do not match this model.")
    }
    list(Chain = Chain$Chain, Draws = Chain$Draws)
  })
  DrawIndex <- BuildDkDrawIndex(list(Chains = Chains))
  TotalDraws <- nrow(DrawIndex)
  if (TotalDraws != Result$TotalPostWarmupDraws || anyDuplicated(DrawIndex)) {
    stop("The saved draw counts do not match this fit.")
  }
  DrawIndex$Simulation <- seq_len(TotalDraws)
  DrawIndex <- DrawIndex[, c("Simulation", "Chain", "PostWarmupIteration")]
  Cells <- Model$ObservedCells
  ObservedData <- data.frame(
    ObservationIndex = seq_len(nrow(Cells)),
    SubjectId = Model$Data$SubjectId[Cells[, 1]],
    OriginalSubjectId = Model$Data$OriginalSubjectId[Cells[, 1]],
    Treatment = Model$Data$treatment[Cells[, 1]],
    Visit = Cells[, 2],
    Month = Model$VisitMonths[Cells[, 2]],
    Baseline = Model$Data$bdi.pre[Cells[, 1]],
    ObservedScore = as.numeric(Model$Y[Model$Observed])
  )
  SimulatedResponse <- matrix(NA_integer_, nrow(Cells), TotalDraws)
  PredictiveSeed <- Settings$Seed + 200000
  set.seed(PredictiveSeed)
  Draw <- 0
  ProgressEvery <- max(1, floor(TotalDraws/10))
  for (Chain in Chains) {
    for (Iteration in seq_len(nrow(Chain$Draws))) {
      Draw <- Draw + 1
      # Fresh joint Gaussian trajectories retain continuous-time dependence.
      # Only the originally observed cells enter the residual diagnostics.
      Scores <- PredictDkScores(Chain$Draws[Iteration, ], Model)
      SimulatedResponse[, Draw] <- Scores[Model$Observed]
      if (Settings$PrintProgress && (Draw %% ProgressEvery == 0 || Draw == TotalDraws)) {
        cat("DHARMa predictive simulation", Draw, "of", TotalDraws, "\n")
      }
    }
  }
  list(ObservedData = ObservedData, SimulatedResponse = SimulatedResponse,
       DrawIndex = DrawIndex, PredictiveSeed = PredictiveSeed)
}

SaveDkDharmaPlot <- function(Directory, Name, PlotFunction, Width = 8, Height = 6) {
  pdf(file.path(Directory, paste0(Name, ".pdf")), width = Width,
      height = Height, family = "Helvetica")
  tryCatch(PlotFunction(), finally = dev.off())
  png(file.path(Directory, paste0(Name, ".png")), width = Width,
      height = Height, units = "in", res = 300)
  tryCatch(PlotFunction(), finally = dev.off())
}

RunBthebDkDharma <- function(Result, ResidualSeed = 123, NBoot = 100, MakePlots = TRUE) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    stop('Install DHARMa first: install.packages("DHARMa").')
  }
  if (length(NBoot) != 1 || !is.finite(NBoot) || NBoot < 1 || NBoot != floor(NBoot)) {
    stop("NBoot must be a positive integer.")
  }
  HadRandomSeed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (HadRandomSeed) PreviousRandomSeed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (HadRandomSeed) {
      assign(".Random.seed", PreviousRandomSeed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  Settings <- Result$Settings
  ModelDirectory <- file.path(Settings$OutputRoot, "ModelOutputs")
  SummaryDirectory <- file.path(Settings$OutputRoot, "Summaries")
  PlotDirectory <- file.path(Settings$OutputRoot, "Plots")
  for (Directory in c(ModelDirectory, SummaryDirectory, PlotDirectory)) {
    dir.create(Directory, recursive = TRUE, showWarnings = FALSE)
  }
  Inputs <- PrepareDkDharmaSimulations(Result)
  SimulatedResponse <- Inputs$SimulatedResponse
  ResidualData <- Inputs$ObservedData
  FittedMedian <- apply(SimulatedResponse, 1, median)
  DharmaResidualObject <- DHARMa::createDHARMa(
    simulatedResponse = SimulatedResponse,
    observedResponse = ResidualData$ObservedScore,
    fittedPredictedResponse = FittedMedian,
    integerResponse = TRUE,
    seed = ResidualSeed,
    method = "PIT"
  )
  ResidualData$FittedMedian <- FittedMedian
  ResidualData$ScaledResidual <- DharmaResidualObject$scaledResiduals
  # The integer outlier test uses DHARMa's bootstrap, with NBoot = 100
  # bootstrap replicates by default.
  # NBoot controls this diagnostic only; all posterior draws remain in use.
  # Treat the p-values as exploratory, particularly the uniformity test:
  # repeated observations can retain dependence in these marginal residuals.
  set.seed(ResidualSeed + 1)
  EffectiveNBoot <- min(NBoot, ncol(SimulatedResponse))
  if (Settings$PrintProgress) {
    cat("Running DHARMa tests with", EffectiveNBoot, "outlier bootstrap replicates...\n")
  }
  Tests <- list(
    Uniformity = DHARMa::testUniformity(DharmaResidualObject, plot = FALSE),
    Dispersion = DHARMa::testDispersion(DharmaResidualObject, plot = FALSE),
    ZeroInflation = DHARMa::testZeroInflation(DharmaResidualObject, plot = FALSE),
    Outliers = DHARMa::testOutliers(DharmaResidualObject, type = "bootstrap",
                                    nBoot = EffectiveNBoot, plot = FALSE)
  )
  TestSummary <- data.frame(
    Test = c("Uniformity", "Dispersion", "Zero inflation", "Outliers"),
    Method = vapply(Tests, function(Test) Test$method, character(1)),
    Statistic = vapply(Tests, function(Test) as.numeric(Test$statistic)[1], numeric(1)),
    PValue = vapply(Tests, function(Test) Test$p.value, numeric(1)),
    row.names = NULL
  )
  PredictiveRows <- list()
  if (Settings$PrintProgress) cat("Summarising posterior predictive scores by treatment and visit...\n")
  Groups <- unique(ResidualData[, c("Treatment", "Visit", "Month")])
  Groups <- Groups[order(Groups$Treatment, Groups$Visit), ]
  Statistics <- list(
    Mean = mean, Variance = var,
    PZero = function(Values) mean(Values == 0),
    PCeiling = function(Values) mean(Values == Settings$ScoreMaximum),
    Q90 = function(Values) as.numeric(quantile(Values, 0.9, names = FALSE))
  )
  for (Group in seq_len(nrow(Groups))) {
    Rows <- which(ResidualData$Treatment == Groups$Treatment[Group] &
                    ResidualData$Visit == Groups$Visit[Group])
    for (Name in names(Statistics)) {
      Statistic <- Statistics[[Name]]
      Predicted <- apply(SimulatedResponse[Rows, , drop = FALSE], 2, Statistic)
      PredictiveRows[[length(PredictiveRows) + 1]] <- data.frame(
        Treatment = Groups$Treatment[Group], Month = Groups$Month[Group],
        Quantity = Name, NObserved = length(Rows),
        Observed = Statistic(ResidualData$ObservedScore[Rows]),
        t(DescribeDkDraws(Predicted)), NDraws = ncol(SimulatedResponse)
      )
    }
  }
  PredictiveSummary <- do.call(rbind, PredictiveRows)
  SimulationSummary <- data.frame(
    TotalPosteriorDraws = Result$TotalPostWarmupDraws,
    NDharmaSimulations = ncol(SimulatedResponse),
    UsedAllRetainedPosteriorDraws = TRUE,
    NObservedCells = nrow(SimulatedResponse),
    PredictiveSeed = Inputs$PredictiveSeed, ResidualSeed = ResidualSeed,
    OutlierBootstrapSeed = ResidualSeed + 1, NBootstrapReplicates = EffectiveNBoot,
    DharmaVersion = as.character(utils::packageVersion("DHARMa")),
    FitSignature = Result$FitSignature
  )
  write.csv(TestSummary, file.path(SummaryDirectory, "DHARMaResidualTestSummary.csv"), row.names = FALSE)
  write.csv(ResidualData, file.path(SummaryDirectory, "DHARMaResidualData.csv"), row.names = FALSE)
  write.csv(PredictiveSummary, file.path(SummaryDirectory, "DHARMaPredictiveCheckSummary.csv"), row.names = FALSE)
  write.csv(SimulationSummary, file.path(SummaryDirectory, "DHARMaSimulationSummary.csv"), row.names = FALSE)
  write.csv(Inputs$DrawIndex, file.path(SummaryDirectory, "DHARMaDrawIndex.csv"), row.names = FALSE)
  SaveRdsSafely(DharmaResidualObject, file.path(ModelDirectory, "DHARMaResidualObject.rds"))
  SaveRdsSafely(Tests, file.path(ModelDirectory, "DHARMaResidualTests.rds"))
  
  if (MakePlots) {
    Colors <- c(TAU = "#0072B2", BtheB = "#D55E00")
    Symbols <- c(TAU = 16, BtheB = 17)
    PlotQq <- function(ShowTests = FALSE) {
      Expected <- (seq_len(nrow(ResidualData)) - 0.5)/nrow(ResidualData)
      plot(Expected, sort(ResidualData$ScaledResidual), xlim = c(0, 1), ylim = c(0, 1),
           xaxs = "i", yaxs = "i", asp = 1, pch = 16, cex = 0.85,
           col = adjustcolor(Colors[["TAU"]], alpha.f = 0.75),
           xlab = "Expected quantiles under Uniform(0,1)",
           ylab = "Empirical posterior predictive PIT quantiles", bty = "l")
      abline(0, 1, lwd = 1)
      if (ShowTests) {
        Labels <- paste0(TestSummary$Test, ": p = ", format.pval(TestSummary$PValue, digits = 3))
        legend("bottomright", legend = Labels, title = "Exploratory checks", bty = "n", cex = 0.9)
      }
    }
    PlotFitted <- function(Rows, Title = "", ShowLegend = TRUE) {
      Data <- ResidualData[Rows, ]
      Treatment <- as.character(Data$Treatment)
      plot(Data$FittedMedian, Data$ScaledResidual, ylim = c(0, 1),
           pch = Symbols[Treatment], col = adjustcolor(Colors[Treatment], alpha.f = 0.5),
           cex = 1, xlab = "Posterior predictive median BDI-II score",
           ylab = "Scaled residual", main = Title, bty = "l")
      abline(h = c(0.25, 0.5, 0.75), lty = c(3, 1, 3), col = "grey60")
      for (Group in names(Colors)) {
        Keep <- which(Treatment == Group)
        if (length(Keep) >= 5 && length(unique(Data$FittedMedian[Keep])) >= 4) {
          Smooth <- lowess(Data$FittedMedian[Keep], Data$ScaledResidual[Keep])
          lines(Smooth, col = Colors[[Group]], lty = match(Group, names(Colors)), lwd = 1.5)
        }
      }
      if (ShowLegend) {
        legend("topright", legend = names(Colors), col = Colors, pch = Symbols,
               lty = 1:2, bty = "n", cex = 0.95)
      }
    }
    SaveDkDharmaPlot(PlotDirectory, "DHARMa_QQ_Residuals", function() {
      par(mar = c(5, 5, 2, 1), mgp = c(3, 0.9, 0),
          cex = 1.1, cex.axis = 1.05, cex.lab = 1.1)
      PlotQq(TRUE)
    }, Width = 6, Height = 6)
    SaveDkDharmaPlot(PlotDirectory, "DHARMa_QQ_Residuals_Clean", function() {
      par(mar = c(5, 5, 2, 1), mgp = c(3, 0.9, 0),
          cex = 1.1, cex.axis = 1.05, cex.lab = 1.1)
      PlotQq(FALSE)
    }, Width = 6, Height = 6)
    SaveDkDharmaPlot(PlotDirectory, "DHARMa_Residual_Diagnostics", function() {
      par(mfrow = c(1, 2), mar = c(4.8, 4.8, 2, 1), mgp = c(2.8, 0.8, 0),
          cex = 1.05, cex.axis = 1, cex.lab = 1.05, cex.main = 1.05)
      PlotQq(TRUE)
      PlotFitted(seq_len(nrow(ResidualData)))
    }, Width = 11, Height = 5.5)
    SaveDkDharmaPlot(PlotDirectory, "DHARMa_Residuals_By_Fitted_Visit", function() {
      layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
             heights = c(1, 1, 0.18))
      par(mar = c(4.5, 4.5, 2.5, 1), mgp = c(2.7, 0.8, 0),
          cex = 1.05, cex.axis = 1, cex.lab = 1.05, cex.main = 1.05)
      for (Month in Settings$VisitMonths) {
        PlotFitted(which(ResidualData$Month == Month), paste("Month", Month), ShowLegend = FALSE)
      }
      par(mar = c(0, 0, 0, 0))
      plot.new()
      legend("center", legend = names(Colors), col = Colors, pch = Symbols,
             lty = 1:2, bty = "n", horiz = TRUE, cex = 1.05)
    }, Width = 10, Height = 8.8)
    SaveDkDharmaPlot(PlotDirectory, "DHARMa_Residuals_By_Visit_Treatment", function() {
      par(mar = c(8, 5, 2, 1), mgp = c(3, 0.9, 0),
          cex = 1.1, cex.axis = 1.05, cex.lab = 1.1)
      Values <- list()
      for (Month in Settings$VisitMonths) {
        for (Group in names(Colors)) {
          Rows <- ResidualData$Month == Month & ResidualData$Treatment == Group
          Values[[paste(Month, Group)]] <- ResidualData$ScaledResidual[Rows]
        }
      }
      Positions <- rep(seq_along(Settings$VisitMonths), each = 2) + rep(c(-0.18, 0.18), length(Settings$VisitMonths))
      boxplot(Values, at = Positions, boxwex = 0.28, axes = FALSE, ylim = c(0, 1),
              xlim = c(0.5, length(Settings$VisitMonths) + 0.5),
              col = rep(adjustcolor(Colors, alpha.f = 0.45), length(Settings$VisitMonths)),
              border = rep(Colors, length(Settings$VisitMonths)),
              xlab = "Month", ylab = "Scaled residual")
      axis(1, at = seq_along(Settings$VisitMonths), labels = Settings$VisitMonths)
      axis(2)
      box(bty = "l")
      abline(h = c(0.25, 0.5, 0.75), lty = c(3, 1, 3), col = "grey60")
      par(xpd = NA)
      legend("bottom", inset = c(0, -0.36), legend = names(Colors),
             fill = adjustcolor(Colors, alpha.f = 0.45), border = Colors,
             bty = "n", horiz = TRUE, cex = 1.05)
    }, Width = 8.5, Height = 6.2)
  }
  cat("\nDHARMa checks completed using all", ncol(SimulatedResponse), "post-warmup draws.\n")
  print(TestSummary[, c("Test", "Statistic", "PValue")], row.names = FALSE)
  cat("\nTables:", SummaryDirectory, "\nResidual object:", ModelDirectory, "\n")
  if (MakePlots) cat("Plots:", PlotDirectory, "\n")
  invisible(list(ResidualObject = DharmaResidualObject, Tests = Tests,
                 TestSummary = TestSummary, ResidualData = ResidualData,
                 PredictiveSummary = PredictiveSummary, SimulationSummary = SimulationSummary,
                 DrawIndex = Inputs$DrawIndex))
}

DharmaDkResult <- RunBthebDkDharma(BthebDkResult)