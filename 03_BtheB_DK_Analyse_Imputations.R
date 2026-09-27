# Beta-binomial analysis of completed BtheB datasets with Rubin pooling.
# Run from the same application directory as the imputation program:
# source("03_BtheB_DK_Analyse_Imputations.R")
# Sourcing this file selects imputations and runs the analysis.
#
# With the default settings and four chains, selection uses 25 evenly spaced
# draws per chain and the same 100 draw indices for every rule. The analysis
# fits 300 models across MAR, halfway borrowing, and J2R.
# Results: BtheB_DK_RBMI/BbAnalysis100Auto
# Sourcing the file again resumes the analysis and reuses successful fits
# from matching checkpoints. The original input bundle is preserved.
# The LaTeX table requires booktabs and siunitx for decimal alignment.
# If glmmTMB is missing, it is installed automatically from CRAN.
#
# Model:
#   cbind(Score, K - Score) ~ 0 + Visit + Visit:Treatment +
#                           BaselineScaled + (1 | Subject)
# Treatment is numeric: TAU = 0, BtheB = 1.
# Visit-specific intercepts and treatment coefficients; common baseline effect,
# BB precision and participant random intercept variance; maximum likelihood.
#
# Target: the empirical baseline distribution of all included participants,
# held fixed and used for both arms. Random intercepts are integrated over
# N(0, sigma^2), not replaced by zero or fitted participant effects.
# Within-imputation variances use a joint delta method for fixed effects and
# log(sigma), retaining their covariance. Dispersion is estimated jointly.
# Rubin intervals use the large-sample complete-data convention. With correlated
# posterior imputations, the finite-M term is an approximation, not an ESS
# correction. Overlapping batch means within chains assess Monte Carlo error.
#
# Outputs: pooled CSV and LaTeX tables, chain summaries, all fit diagnostics,
# compact estimates/variances in RDS, and resumable batch checkpoints.
# All 100 selected imputations enter pooling; failed fits are not excluded.
# With 25 selected draws per chain, batch-means MCSE is reported as NA.
#
# glmmTMB references:
# https://glmmtmb.github.io/glmmTMB/reference/glmmTMB.html
# https://glmmtmb.github.io/glmmTMB/reference/vcov.glmmTMB.html

DefaultBthebBbSettings <- function() {
  list(
    InputFile = file.path("BtheB_DK_RBMI", "CompletedDatasets",
                          "BtheB_DK_CompletedDatasets.rds"),
    OutputRoot = file.path("BtheB_DK_RBMI", "BbAnalysis100Auto"),
    NImputations = 100,
    UseParallel = TRUE,
    CoreFraction = 0.8,
    Workers = NULL,
    BatchSize = 10,
    MaxIterations = 2000,
    MaxEvaluations = 3000,
    MaxGradient = 0.01,
    QuadratureNodes = c(30, 60, 120, 240),
    QuadratureTolerance = 1e-6,
    ConfidenceLevel = 0.95,
    RetryFailed = TRUE,
    PrintProgress = TRUE
  )
}

BbSaveRds <- function(Object, File) {
  Temporary <- tempfile("BbWrite", tmpdir = dirname(File))
  Backup <- paste0(File, ".previous")
  on.exit(unlink(Temporary), add = TRUE)
  saveRDS(Object, Temporary)
  if (file.exists(File)) {
    if (file.exists(Backup)) unlink(Backup)
    if (!file.rename(File, Backup)) stop("Could not back up: ", File)
  }
  if (!file.rename(Temporary, File)) {
    if (file.exists(Backup)) file.rename(Backup, File)
    stop("Could not save: ", File)
  }
  if (file.exists(Backup)) unlink(Backup)
  invisible(File)
}

BbReadCheckpoint <- function(File) {
  if (!file.exists(File) && file.exists(paste0(File, ".previous"))) {
    File <- paste0(File, ".previous")
  }
  if (!file.exists(File)) return(NULL)
  tryCatch(readRDS(File), error = function(Error) {
    stop("Could not read checkpoint: ", File, ". ", conditionMessage(Error))
  })
}

BbHash <- function(Object) {
  File <- tempfile("BbHash")
  on.exit(unlink(File), add = TRUE)
  saveRDS(Object, File, compress = FALSE, version = 2)
  unname(tools::md5sum(File))
}

BbNormalQuadrature <- function(Nodes) {
  Jacobi <- matrix(0, Nodes, Nodes)
  Indices <- seq_len(Nodes - 1)
  Jacobi[cbind(Indices, Indices + 1)] <- sqrt(Indices)
  Jacobi[cbind(Indices + 1, Indices)] <- sqrt(Indices)
  Eigen <- eigen(Jacobi, symmetric = TRUE)
  Order <- order(Eigen$values)
  list(Nodes = Eigen$values[Order], Weights = Eigen$vectors[1, Order]^2)
}

BbSelectImputations <- function(Bundle, NImputations) {
  DrawIndex <- Bundle$DrawIndex
  Required <- c("Imputation", "Chain", "PostWarmupIteration")
  if (!is.data.frame(DrawIndex) || !all(Required %in% names(DrawIndex)) ||
      anyNA(DrawIndex[, Required]) || !is.list(Bundle$ImputedValues) ||
      length(Bundle$ImputedValues) == 0) stop("Invalid completed-data bundle.")
  Available <- nrow(DrawIndex)
  if (NImputations > Available) stop("The bundle contains fewer than ", NImputations, " imputations.")
  if (!identical(as.numeric(DrawIndex$Imputation), as.numeric(seq_len(Available))) ||
      anyDuplicated(DrawIndex[, c("Chain", "PostWarmupIteration")])) {
    stop("Invalid source draw identifiers.")
  }
  for (Values in Bundle$ImputedValues) {
    if (!is.matrix(Values) || ncol(Values) != Available) {
      stop("Imputed values and draw identifiers do not agree.")
    }
  }
  # Preserve original draw identifiers when selecting from a previously subsetted bundle.
  if (!"OriginalImputation" %in% names(DrawIndex)) {
    DrawIndex$OriginalImputation <- DrawIndex$Imputation
  }
  if (!"OriginalPostWarmupIteration" %in% names(DrawIndex)) {
    DrawIndex$OriginalPostWarmupIteration <- DrawIndex$PostWarmupIteration
  }
  DrawIndex$PostWarmupIteration <- DrawIndex$OriginalPostWarmupIteration
  ChainIndices <- split(seq_len(Available), DrawIndex$Chain)
  NChains <- length(ChainIndices)
  if (NImputations %% NChains != 0) {
    stop("NImputations must be divisible by the number of chains.")
  }
  PerChain <- NImputations/NChains
  if (any(lengths(ChainIndices) < PerChain)) stop("A chain has too few available imputations.")
  Selected <- unlist(lapply(ChainIndices, function(Indices) {
    Indices <- Indices[order(DrawIndex$PostWarmupIteration[Indices])]
    Positions <- if (PerChain == 1) ceiling(length(Indices)/2) else {
      round(seq(1, length(Indices), length.out = PerChain))
    }
    Indices[Positions]
  }), use.names = FALSE)
  if (length(Selected) != NImputations || anyDuplicated(Selected)) {
    stop("The selected imputation indices are not valid.")
  }
  Bundle$DrawIndex <- DrawIndex[Selected, , drop = FALSE]
  Bundle$DrawIndex$Imputation <- seq_along(Selected)
  rownames(Bundle$DrawIndex) <- NULL
  Bundle$ImputedValues <- lapply(Bundle$ImputedValues, function(Values) {
    Values[, Selected, drop = FALSE]
  })
  Bundle$AnalysisSelection <- list(
    Method = "Evenly spaced within each chain; identical indices across rules",
    AvailablePerRule = Available, SelectedPerRule = NImputations,
    SelectedPerChain = PerChain, SourceColumns = Selected
  )
  Bundle
}

BbPrepareContext <- function(Bundle, Settings) {
  Required <- c("Data", "ObservedScores", "MissingCells", "ImputedValues",
                "DrawIndex", "VisitMonths", "ScoreMaximum", "BaselineCenter",
                "BaselineScale", "FitSignature")
  if (!all(Required %in% names(Bundle))) stop("This is not a DK completed-data bundle.")
  Scores <- Bundle$ObservedScores
  if (!is.matrix(Scores) || !is.numeric(Scores)) stop("Invalid recorded score matrix.")
  NSubjects <- nrow(Scores)
  NVisits <- ncol(Scores)
  NImputations <- nrow(Bundle$DrawIndex)
  K <- Bundle$ScoreMaximum
  if (length(K) != 1 || !is.finite(K) || K < 2 || K != floor(K)) stop("Invalid score bound.")
  if (nrow(Bundle$Data) != NSubjects || NVisits != length(Bundle$VisitMonths) ||
      NImputations < 2) stop("The bundle dimensions do not agree.")
  if (any(!is.finite(Bundle$VisitMonths)) || any(diff(Bundle$VisitMonths) <= 0)) {
    stop("Invalid visit times.")
  }
  Recorded <- Scores[!is.na(Scores)]
  if (any(!is.finite(Recorded) | Recorded < 0 | Recorded > K |
          Recorded != floor(Recorded))) stop("Invalid recorded scores.")
  Cells <- Bundle$MissingCells
  if (!is.matrix(Cells) || ncol(Cells) != 2 || anyNA(Cells) ||
      any(Cells != floor(Cells)) || any(Cells[, 1] < 1 | Cells[, 1] > NSubjects) ||
      any(Cells[, 2] < 1 | Cells[, 2] > NVisits)) stop("Invalid missing-cell indices.")
  Positions <- (Cells[, 2] - 1)*NSubjects + Cells[, 1]
  if (anyDuplicated(Positions) || !setequal(Positions, which(is.na(Scores)))) {
    stop("Missing-cell indices do not match the recorded score matrix.")
  }
  Rules <- names(Bundle$ImputedValues)
  if (length(Rules) == 0 || anyDuplicated(Rules) ||
      any(!Rules %in% c("MAR", "Half", "J2R"))) stop("Unexpected imputation rules.")
  for (Rule in Rules) {
    Values <- Bundle$ImputedValues[[Rule]]
    if (!is.matrix(Values) || !is.numeric(Values) ||
        !identical(dim(Values), as.integer(c(nrow(Cells), NImputations))) ||
        any(!is.finite(Values) | Values < 0 | Values > K | Values != floor(Values))) {
      stop("Invalid imputations for ", Rule, ".")
    }
  }
  DrawIndex <- Bundle$DrawIndex
  RequiredIndex <- c("Imputation", "Chain", "PostWarmupIteration")
  if (!all(RequiredIndex %in% names(DrawIndex)) ||
      anyNA(DrawIndex[, RequiredIndex]) ||
      !identical(as.numeric(DrawIndex$Imputation), as.numeric(seq_len(NImputations))) ||
      anyDuplicated(DrawIndex[, c("Chain", "PostWarmupIteration")])) {
    stop("Invalid imputation draw index.")
  }
  for (Chain in unique(DrawIndex$Chain)) {
    Iterations <- DrawIndex$PostWarmupIteration[DrawIndex$Chain == Chain]
    if (!is.numeric(Iterations) || any(!is.finite(Iterations)) ||
        any(Iterations < 1 | Iterations != floor(Iterations)) || anyDuplicated(Iterations)) {
      stop("Invalid original post-warmup iteration numbers.")
    }
  }
  Baseline <- Bundle$Data$bdi.pre
  Arm <- as.character(Bundle$Data$treatment)
  if (!is.numeric(Baseline) || length(Baseline) != NSubjects ||
      any(!is.finite(Baseline)) || !setequal(Arm, c("TAU", "BtheB")) || anyNA(Arm)) {
    stop("Baseline scores or treatment labels are missing.")
  }
  if (length(Bundle$BaselineScale) != 1 || !is.finite(Bundle$BaselineScale) ||
      Bundle$BaselineScale <= 0 || length(Bundle$BaselineCenter) != 1 ||
      !is.finite(Bundle$BaselineCenter)) stop("Invalid baseline centering or scale.")
  BaselineScaled <- (Baseline - Bundle$BaselineCenter)/Bundle$BaselineScale
  LongData <- data.frame(
    Score = as.vector(Scores), Failures = K - as.vector(Scores),
    Subject = factor(rep(seq_len(NSubjects), NVisits)),
    Visit = factor(rep(seq_len(NVisits), each = NSubjects), levels = seq_len(NVisits)),
    Treatment = rep(as.numeric(Arm == "BtheB"), NVisits),
    BaselineScaled = rep(BaselineScaled, NVisits)
  )
  FixedFormula <- ~0 + Visit + Visit:Treatment + BaselineScaled
  Design <- model.matrix(FixedFormula, LongData)
  if (qr(Design)$rank != ncol(Design)) stop("The analysis design matrix is not full rank.")
  PredictionDesign <- list()
  for (Visit in seq_len(NVisits)) {
    for (Treatment in 0:1) {
      NewData <- data.frame(
        Visit = factor(rep(Visit, NSubjects), levels = seq_len(NVisits)),
        Treatment = Treatment, BaselineScaled = BaselineScaled
      )
      PredictionDesign[[length(PredictionDesign) + 1]] <- model.matrix(FixedFormula, NewData)
    }
  }
  Transformation <- matrix(0, 3*NVisits, 2*NVisits)
  Targets <- data.frame(
    Visit = rep(seq_len(NVisits), each = 3),
    Month = rep(Bundle$VisitMonths, each = 3),
    Summary = rep(c("TAU", "BtheB", "BtheB minus TAU"), NVisits)
  )
  for (Visit in seq_len(NVisits)) {
    Rows <- 3*(Visit - 1) + 1:3
    Columns <- 2*(Visit - 1) + 1:2
    Transformation[Rows, Columns] <- rbind(c(1, 0), c(0, 1), c(-1, 1))
  }
  list(
    Bundle = Bundle, LongData = LongData, PredictionDesign = PredictionDesign,
    FixedNames = colnames(Design), Transformation = Transformation, Targets = Targets,
    Rules = Rules, NImputations = NImputations, NSubjects = NSubjects, NVisits = NVisits,
    Quadrature = lapply(Settings$QuadratureNodes, BbNormalQuadrature)
  )
}

BbMarginalMoments <- function(Beta, LogSd, Context, Quadrature) {
  Sd <- exp(LogSd)
  Shift <- Sd*Quadrature$Nodes
  NMeans <- length(Context$PredictionDesign)
  Means <- numeric(NMeans)
  Gradient <- matrix(0, NMeans, length(Beta) + 1)
  for (Index in seq_len(NMeans)) {
    Design <- Context$PredictionDesign[[Index]]
    Probability <- plogis(outer(as.vector(Design%*%Beta), Shift, "+"))
    Derivative <- Probability*(1 - Probability)
    Means[Index] <- Context$Bundle$ScoreMaximum*mean(Probability%*%Quadrature$Weights)
    Gradient[Index, seq_along(Beta)] <- Context$Bundle$ScoreMaximum*
      colMeans(Design*as.vector(Derivative%*%Quadrature$Weights))
    Gradient[Index, length(Beta) + 1] <- Context$Bundle$ScoreMaximum*
      mean(Derivative%*%(Quadrature$Weights*Shift))
  }
  list(Estimate = as.vector(Context$Transformation%*%Means),
       Gradient = Context$Transformation%*%Gradient)
}

BbStandardize <- function(Fit, Context, Settings) {
  Beta <- glmmTMB::fixef(Fit)$cond
  if (!identical(names(Beta), Context$FixedNames)) stop("Unexpected fixed-effect ordering.")
  Parameter <- Fit$fit$par
  ParameterNames <- names(Parameter)
  BetaIndex <- which(ParameterNames == "beta")
  ThetaIndex <- which(ParameterNames == "theta")
  if (length(BetaIndex) != length(Beta) || length(ThetaIndex) != 1 ||
      !isTRUE(all.equal(unname(Parameter[BetaIndex]), unname(Beta)))) {
    stop("Unexpected glmmTMB parameter layout.")
  }
  # TMB's covariance uses the same unconstrained coordinates as par.fixed.
  # Verify that they match the fitted optimizer vector before subsetting.
  if (!isTRUE(all.equal(as.numeric(Fit$sdr$par.fixed), as.numeric(Parameter),
                        tolerance = 1e-8))) stop("The parameter covariance is not aligned.")
  FullCovariance <- as.matrix(Fit$sdr$cov.fixed)
  if (!identical(dim(FullCovariance), as.integer(rep(length(Parameter), 2))) ||
      any(!is.finite(FullCovariance))) stop("The parameter covariance is unavailable.")
  PublicCovariance <- as.matrix(vcov(Fit)$cond)
  if (!isTRUE(all.equal(unname(FullCovariance[BetaIndex, BetaIndex, drop = FALSE]),
                        unname(PublicCovariance), tolerance = 1e-7))) {
    stop("Fixed-effect covariance verification failed.")
  }
  Keep <- c(BetaIndex, ThetaIndex)
  Covariance <- FullCovariance[Keep, Keep, drop = FALSE]
  if (inherits(try(chol(Covariance), silent = TRUE), "try-error")) {
    stop("The parameter covariance is not positive definite.")
  }
  Previous <- NULL
  for (Index in seq_along(Context$Quadrature)) {
    Current <- BbMarginalMoments(Beta, Parameter[ThetaIndex], Context,
                                 Context$Quadrature[[Index]])
    Current$Variance <- rowSums((Current$Gradient%*%Covariance)*Current$Gradient)
    if (any(!is.finite(c(Current$Estimate, Current$Variance))) ||
        any(Current$Variance < 0)) stop("Invalid standardized estimates or variances.")
    if (!is.null(Previous)) {
      MeanError <- max(abs(Current$Estimate - Previous$Estimate))
      VarianceError <- max(abs(Current$Variance - Previous$Variance)/
                             (1 + abs(Current$Variance)))
      if (MeanError <= Settings$QuadratureTolerance &&
          VarianceError <= Settings$QuadratureTolerance) {
        return(list(Estimate = Current$Estimate, Variance = Current$Variance,
                    Nodes = Settings$QuadratureNodes[Index],
                    QuadratureMeanError = MeanError, QuadratureVarianceError = VarianceError))
      }
    }
    Previous <- Current
  }
  stop("Random intercept integration did not reach the requested accuracy.")
}

BbAnalyseOne <- function(Imputation, Rule, Context, Settings, Start = NULL) {
  Started <- proc.time()[[3]]
  Data <- Context$LongData
  Scores <- Context$Bundle$ObservedScores
  Scores[Context$Bundle$MissingCells] <- Context$Bundle$ImputedValues[[Rule]][, Imputation]
  Data$Score <- as.vector(Scores)
  Data$Failures <- Context$Bundle$ScoreMaximum - Data$Score
  Errors <- character()
  Warnings <- character()
  # A fresh BFGS fit is the fallback, so a poor warm start is not reused.
  for (Attempt in 1:2) {
    AttemptWarnings <- character()
    Candidate <- tryCatch(withCallingHandlers({
      Control <- if (Attempt == 1) {
        glmmTMB::glmmTMBControl(
          optCtrl = list(iter.max = Settings$MaxIterations, eval.max = Settings$MaxEvaluations),
          parallel = 1
        )
      } else {
        glmmTMB::glmmTMBControl(
          optimizer = optim, optArgs = list(method = "BFGS"),
          optCtrl = list(maxit = Settings$MaxIterations, reltol = 1e-10), parallel = 1
        )
      }
      Fit <- glmmTMB::glmmTMB(
        cbind(Score, Failures) ~ 0 + Visit + Visit:Treatment +
          BaselineScaled + (1 | Subject),
        family = glmmTMB::betabinomial(link = "logit"), data = Data,
        dispformula = ~1, ziformula = ~0, REML = FALSE, na.action = na.fail,
        control = Control, start = if (Attempt == 1) Start else NULL
      )
      if (!isTRUE(Fit$fit$convergence == 0)) stop("The optimizer did not converge.")
      if (!isTRUE(Fit$sdr$pdHess)) stop("The Hessian is not positive definite.")
      Gradient <- max(abs(Fit$obj$gr(Fit$fit$par)))
      if (!is.finite(Gradient) || Gradient > Settings$MaxGradient) {
        stop("The likelihood gradient needs review.")
      }
      Marginal <- BbStandardize(Fit, Context, Settings)
      Parameter <- Fit$fit$par
      DispersionIndex <- which(names(Parameter) %in% c("betad", "betadisp"))
      if (length(DispersionIndex) != 1) stop("Unexpected BB dispersion model.")
      ThetaIndex <- which(names(Parameter) == "theta")
      Precision <- exp(Parameter[DispersionIndex])
      RandomSd <- exp(Parameter[ThetaIndex])
      if (any(!is.finite(c(Precision, RandomSd)))) stop("Invalid scale estimates.")
      list(
        Ok = TRUE, Estimate = Marginal$Estimate, Variance = Marginal$Variance,
        Start = lapply(split(unname(Parameter), names(Parameter)), unname),
        Coefficients = c(glmmTMB::fixef(Fit)$cond, LogRandomSd = unname(Parameter[ThetaIndex]),
                         LogPrecision = unname(Parameter[DispersionIndex])),
        Attempts = Attempt, Optimizer = c("nlminb", "BFGS")[Attempt],
        Convergence = Fit$fit$convergence, PdHess = Fit$sdr$pdHess,
        MaxGradient = Gradient, Nodes = Marginal$Nodes,
        QuadratureMeanError = Marginal$QuadratureMeanError,
        QuadratureVarianceError = Marginal$QuadratureVarianceError,
        RandomSd = unname(RandomSd), Precision = unname(Precision)
      )
    }, warning = function(Warning) {
      AttemptWarnings <<- c(AttemptWarnings, conditionMessage(Warning))
      invokeRestart("muffleWarning")
    }), error = function(Error) Error)
    Warnings <- c(Warnings, AttemptWarnings)
    if (!inherits(Candidate, "error")) {
      Candidate$Warnings <- paste(unique(Warnings), collapse = " | ")
      Candidate$EarlierErrors <- paste(Errors, collapse = " | ")
      Candidate$Error <- ""
      Candidate$RuntimeSeconds <- proc.time()[[3]] - Started
      return(Candidate)
    }
    Errors <- c(Errors, conditionMessage(Candidate))
  }
  list(Ok = FALSE, Attempts = 2, Optimizer = "nlminb; BFGS",
       Warnings = paste(unique(Warnings), collapse = " | "),
       Error = paste(unique(Errors), collapse = " | "), EarlierErrors = "",
       RuntimeSeconds = proc.time()[[3]] - Started)
}

BbRunTask <- function(Task, Context, Settings, Signature) {
  Cached <- BbReadCheckpoint(Task$File)
  if (!is.null(Cached) && (!identical(Cached$Signature, Signature) ||
                           !identical(Cached$Indices, Task$Indices) || !identical(Cached$Rule, Task$Rule))) {
    stop("A checkpoint does not match this analysis: ", Task$File)
  }
  Results <- if (is.null(Cached)) vector("list", length(Task$Indices)) else Cached$Results
  if (length(Results) != length(Task$Indices)) stop("An incomplete checkpoint needs review.")
  Start <- NULL
  for (Index in seq_along(Task$Indices)) {
    Existing <- Results[[Index]]
    if (is.null(Existing) || (!isTRUE(Existing$Ok) && Settings$RetryFailed)) {
      Results[[Index]] <- BbAnalyseOne(Task$Indices[Index], Task$Rule, Context, Settings, Start)
    }
    if (isTRUE(Results[[Index]]$Ok)) Start <- Results[[Index]]$Start
  }
  BbSaveRds(list(Signature = Signature, Rule = Task$Rule,
                 Indices = Task$Indices, Results = Results), Task$File)
  list(N = length(Results), Failed = sum(!vapply(Results, function(Value) isTRUE(Value$Ok), logical(1))))
}

BbWorkerTask <- function(Task) {
  BbRunTask(Task, BbWorkerContext, BbWorkerSettings, BbWorkerSignature)
}

BbMcse <- function(Values, DrawIndex) {
  # Overlapping batch means use every iteration, including chain ends.
  Total <- length(Values)
  Variance <- 0
  for (Chain in unique(DrawIndex$Chain)) {
    Positions <- which(DrawIndex$Chain == Chain)
    Positions <- Positions[order(DrawIndex$PostWarmupIteration[Positions])]
    ChainValues <- Values[Positions]
    N <- length(ChainValues)
    if (N < 100) return(NA_real_)
    BatchSize <- floor(sqrt(N))
    CumSum <- c(0, cumsum(ChainValues))
    Starts <- seq_len(N - BatchSize + 1)
    BatchMeans <- (CumSum[Starts + BatchSize] - CumSum[Starts])/BatchSize
    LongRunVariance <- N*BatchSize/((N - BatchSize)*(N - BatchSize + 1))*
      sum((BatchMeans - mean(ChainValues))^2)
    Variance <- Variance + N*LongRunVariance/Total^2
  }
  sqrt(Variance)
}

BbPool <- function(Estimates, Variances, DrawIndex, Targets, Rule, Level) {
  M <- nrow(Estimates)
  Mean <- colMeans(Estimates)
  Within <- colMeans(Variances)
  Between <- apply(Estimates, 2, var)
  Total <- Within + (1 + 1/M)*Between
  Fraction <- (1 + 1/M)*Between/Total
  Degrees <- ifelse(Between == 0, Inf, (M - 1)/Fraction^2)
  Critical <- qt((1 + Level)/2, Degrees)
  Mcse <- apply(Estimates, 2, BbMcse, DrawIndex = DrawIndex)
  data.frame(
    Rule = Rule, Targets, M = M, Estimate = Mean, WithinVariance = Within,
    BetweenVariance = Between, TotalVariance = Total, SE = sqrt(Total),
    DegreesOfFreedom = Degrees, Lower = Mean - Critical*sqrt(Total),
    Upper = Mean + Critical*sqrt(Total), MCSE = Mcse,
    MCSEOverSE = Mcse/sqrt(Total), ConfidenceLevel = Level, row.names = NULL
  )
}

# Derive numeric column widths from displayed values for each treatment summary.
BbTableNumberFormat <- function(Values, Digits = 2) {
  Text <- formatC(Values, format = "f", digits = Digits, decimal.mark = ".")
  Integer <- sub("[.].*$", "", Text)
  paste0(if (any(startsWith(Integer, "-"))) "-" else "",
         max(nchar(sub("^-", "", Integer))), ".", Digits)
}

BbWriteTable <- function(Pooled, Context, Settings) {
  Level <- format(100*Settings$ConfidenceLevel, trim = TRUE)
  RuleLabels <- c(MAR = "MAR", Half = "Halfway", J2R = "J2R")
  Targets <- c("TAU", "BtheB", "BtheB minus TAU")
  Columns <- unlist(lapply(seq_along(Targets), function(Index) {
    Rows <- Pooled[Pooled$Summary == Targets[Index], ]
    Formats <- vapply(c("Estimate", "Lower", "Upper"), function(Name) {
      BbTableNumberFormat(Rows[[Name]])
    }, character(1))
    c(
      paste0("S[table-format = ", Formats["Estimate"], "]"),
      "@{\\hspace{\\fontdimen2\\font}\\extracolsep{0pt}}",
      paste0("S[table-format = ", Formats["Lower"],
             ", table-space-text-pre = {(}, table-space-text-post = {,}]"),
      "@{\\hspace{\\fontdimen2\\font}}",
      paste0("S[table-format = ", Formats["Upper"], ", table-space-text-post = {)}]"),
      if (Index < length(Targets)) "@{\\hspace{1em}\\extracolsep{\\fill}}" else "@{}"
    )
  }), use.names = FALSE)
  Lines <- c(
    "% Requires booktabs and siunitx. Numerical results are unchanged by this export.",
    "\\begin{table}[tbp]",
    "\\caption{Standardized mean BDI-II scores and treatment differences after reference-based imputation.}",
    "\\label{tab:bb-analysis}", "\\centering", "\\begingroup\\small",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\sisetup{table-align-text-pre = false, table-align-text-post = false,",
    "  group-digits = false, output-decimal-marker = {.}, round-mode = none}",
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lc",
    Columns, "}",
    "\\toprule",
    "& & \\multicolumn{6}{c}{\\textbf{Standardized mean}} & \\multicolumn{3}{c}{\\textbf{Mean difference}} \\\\",
    "\\cmidrule(lr){3-8}\\cmidrule(l){9-11}",
    "\\textbf{Rule} & \\textbf{Month} & \\multicolumn{3}{c}{\\textbf{TAU}} & \\multicolumn{3}{c}{\\textbf{BtheB}} & \\multicolumn{3}{c}{\\textbf{BtheB minus TAU}} \\\\",
    "\\midrule"
  )
  Cell <- function(Row) paste0(
    formatC(Row$Estimate, format = "f", digits = 2, decimal.mark = "."), " & {(}",
    formatC(Row$Lower, format = "f", digits = 2, decimal.mark = "."), "{,} & ",
    formatC(Row$Upper, format = "f", digits = 2, decimal.mark = "."), "{)}"
  )
  for (Rule in Context$Rules) {
    for (Visit in seq_len(Context$NVisits)) {
      Rows <- Pooled[Pooled$Rule == Rule & Pooled$Visit == Visit, ]
      Cells <- vapply(Targets,
                      function(Target) Cell(Rows[Rows$Summary == Target, ]), character(1))
      Lines <- c(Lines, paste0(
        if (Visit == 1) RuleLabels[[Rule]] else "", " & ",
        Context$Bundle$VisitMonths[Visit], " & ", paste(Cells, collapse = " & "), " \\\\"
      ))
    }
    if (Rule != tail(Context$Rules, 1)) Lines <- c(Lines, "\\addlinespace")
  }
  Lines <- c(
    Lines, "\\bottomrule", "\\end{tabular*}", "\\endgroup", "\\par\\smallskip",
    "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright",
    paste0(
      "Predictions from the BB random intercept model are integrated over the fitted ",
      "random intercept distribution and averaged over the baseline scores of all ",
      Context$NSubjects, " included participants, using the same distribution for both arms. ",
      "Estimates and delta-method variances are pooled separately within each rule ",
      "using ", format(Context$NImputations, big.mark = ",", scientific = FALSE, trim = TRUE),
      " completed datasets selected at evenly spaced iterations within each chain. ",
      "The same draw indices are used across rules. ",
      "Intervals are ", Level, "\\% confidence intervals, using Rubin's large-sample complete-data approximation. ",
      "The primary comparison is the mean difference at month ",
      tail(Context$Bundle$VisitMonths, 1), "; negative differences favor BtheB."
    ), "\\end{minipage}", "\\end{table}"
  )
  writeLines(Lines, file.path(Settings$OutputRoot, "BbAnalysisTable.tex"))
  writeLines(c(
    "\\documentclass[12pt]{article}", "\\usepackage[margin=0.9in]{geometry}",
    "\\usepackage[T1]{fontenc}", "\\usepackage{helvet}",
    "\\renewcommand{\\familydefault}{\\sfdefault}", "\\usepackage{booktabs,siunitx}",
    "\\usepackage{caption}",
    "\\captionsetup{justification=raggedright,singlelinecheck=false}",
    "\\captionsetup[table]{name=TABLE,labelfont=bf,labelsep=period}",
    "\\pagestyle{empty}", "\\begin{document}", Lines, "\\end{document}"
  ), file.path(Settings$OutputRoot, "BbAnalysisTableStandalone.tex"))
}

BbCollect <- function(Tasks, Context, Settings, Signature) {
  M <- Context$NImputations
  Q <- nrow(Context$Targets)
  Estimates <- Variances <- setNames(lapply(Context$Rules, function(Rule) matrix(NA_real_, M, Q)),
                                     Context$Rules)
  Coefficients <- setNames(lapply(Context$Rules, function(Rule) {
    matrix(NA_real_, M, length(Context$FixedNames) + 2,
           dimnames = list(NULL, c(Context$FixedNames, "LogRandomSd", "LogPrecision")))
  }), Context$Rules)
  StatusRows <- vector("list", M*length(Context$Rules))
  Counter <- 0
  for (Task in Tasks) {
    Checkpoint <- BbReadCheckpoint(Task$File)
    if (is.null(Checkpoint) || !identical(Checkpoint$Signature, Signature) ||
        !identical(Checkpoint$Indices, Task$Indices) ||
        !identical(Checkpoint$Rule, Task$Rule)) stop("A required checkpoint is missing or mismatched.")
    for (Index in seq_along(Task$Indices)) {
      Imputation <- Task$Indices[Index]
      Fit <- Checkpoint$Results[[Index]]
      Ok <- isTRUE(Fit$Ok)
      Value <- function(Name, Default = NA_real_) if (is.null(Fit[[Name]])) Default else Fit[[Name]]
      Counter <- Counter + 1
      StatusRows[[Counter]] <- data.frame(
        Rule = Task$Rule, Context$Bundle$DrawIndex[Imputation, ], Success = Ok,
        Attempts = Value("Attempts"), Optimizer = Value("Optimizer", ""),
        Convergence = Value("Convergence"), PdHess = Value("PdHess", FALSE),
        MaxGradient = Value("MaxGradient"), QuadratureNodes = Value("Nodes"),
        QuadratureMeanError = Value("QuadratureMeanError"),
        QuadratureVarianceError = Value("QuadratureVarianceError"),
        RandomInterceptSd = Value("RandomSd"), BbPrecision = Value("Precision"),
        RuntimeSeconds = Value("RuntimeSeconds"), Warnings = Value("Warnings", ""),
        EarlierErrors = Value("EarlierErrors", ""), Error = Value("Error", ""),
        row.names = NULL
      )
      if (Ok) {
        Estimates[[Task$Rule]][Imputation, ] <- Fit$Estimate
        Variances[[Task$Rule]][Imputation, ] <- Fit$Variance
        Coefficients[[Task$Rule]][Imputation, ] <- Fit$Coefficients
      }
    }
  }
  Status <- do.call(rbind, StatusRows)
  write.csv(Status, file.path(Settings$OutputRoot, "BbFitDiagnostics.csv"), row.names = FALSE)
  Result <- list(
    Complete = all(Status$Success), Signature = Signature, Settings = Settings,
    FitSignature = Context$Bundle$FitSignature, DrawIndex = Context$Bundle$DrawIndex,
    Selection = Context$Bundle$AnalysisSelection,
    Targets = Context$Targets, Estimates = Estimates, Variances = Variances,
    Coefficients = Coefficients, Diagnostics = Status
  )
  if (!Result$Complete) {
    BbSaveRds(Result, file.path(Settings$OutputRoot, "BbAnalysisIncomplete.rds"))
    stop(sum(!Status$Success), " fits need review; see BbFitDiagnostics.csv. ",
         "No selected imputations were excluded and no pooled results were produced.")
  }
  Pooled <- ChainRows <- list()
  for (Rule in Context$Rules) {
    Pooled[[Rule]] <- BbPool(Estimates[[Rule]], Variances[[Rule]], Context$Bundle$DrawIndex,
                             Context$Targets, Rule, Settings$ConfidenceLevel)
    for (Chain in unique(Context$Bundle$DrawIndex$Chain)) {
      Positions <- which(Context$Bundle$DrawIndex$Chain == Chain)
      if (length(Positions) < 2) next
      Row <- BbPool(Estimates[[Rule]][Positions, , drop = FALSE],
                    Variances[[Rule]][Positions, , drop = FALSE],
                    Context$Bundle$DrawIndex[Positions, ], Context$Targets,
                    Rule, Settings$ConfidenceLevel)
      Row$Chain <- Chain
      ChainRows[[length(ChainRows) + 1]] <- Row
    }
  }
  Result$Pooled <- do.call(rbind, Pooled)
  Result$ByChain <- do.call(rbind, ChainRows)
  write.csv(Result$Pooled, file.path(Settings$OutputRoot, "BbPooledResults.csv"), row.names = FALSE)
  write.csv(Result$Pooled[Result$Pooled$Visit == Context$NVisits &
                            Result$Pooled$Summary == "BtheB minus TAU", ],
            file.path(Settings$OutputRoot, "BbPrimaryComparison.csv"), row.names = FALSE)
  write.csv(Result$ByChain, file.path(Settings$OutputRoot, "BbResultsByChain.csv"), row.names = FALSE)
  BbWriteTable(Result$Pooled, Context, Settings)
  BbSaveRds(Result, file.path(Settings$OutputRoot, "BbAnalysisResult.rds"))
  capture.output(sessionInfo(), file = file.path(Settings$OutputRoot, "SessionInfo.txt"))
  writeLines(c(
    "All selected imputations were analyzed; none were excluded because of fit diagnostics.",
    paste("Available imputations per rule:", Context$Bundle$AnalysisSelection$AvailablePerRule),
    paste("Imputations per rule:", M),
    paste("Selected imputations per chain:", Context$Bundle$AnalysisSelection$SelectedPerChain),
    "Selection uses evenly spaced iterations within each chain and identical indices across rules.",
    "Original draw identifiers are in SelectedImputationIndex.csv.",
    "The full input bundle and original posterior summaries are unchanged.",
    paste("Standardization population:", Context$NSubjects, "included participants."),
    "The empirical covariate distribution is held fixed for variance estimation.",
    "The BB model uses a common baseline coefficient, precision and random intercept variance.",
    "Random intercept integration and its variance contribution are included in each fit.",
    "Nominal Rubin intervals are not clipped to the score bounds.",
    "The finite-M factor uses the actual imputation count, not an effective sample size.",
    "MCSE uses overlapping batch means within each chain; fewer than 100 draws per chain gives NA.",
    "Review MCSE/SE and agreement across chains before reporting results.",
    "Rubin intervals are a pooling convention; coverage under DK imputation and BB analysis requires simulation.",
    paste("Maximum gradient:", max(Status$MaxGradient)),
    paste("Fits requiring the fallback optimizer:", sum(Status$Attempts > 1)),
    paste("Fits recording warnings:", sum(nzchar(Status$Warnings)))
  ), file.path(Settings$OutputRoot, "AnalysisNotes.txt"))
  Result
}

RunBthebBbAnalysis <- function(Settings = DefaultBthebBbSettings()) {
  if (!file.exists(Settings$InputFile)) stop("The completed-data bundle was not found: ", Settings$InputFile)
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    message("Installing glmmTMB from CRAN.")
    install.packages("glmmTMB", repos = "https://cloud.r-project.org")
    if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("glmmTMB could not be installed.")
  }
  for (Name in c("NImputations", "BatchSize", "MaxIterations", "MaxEvaluations")) {
    Value <- Settings[[Name]]
    if (length(Value) != 1 || !is.finite(Value) || Value < 1 || Value != floor(Value)) {
      stop(Name, " must be a positive integer.")
    }
  }
  if (length(Settings$ConfidenceLevel) != 1 || !is.finite(Settings$ConfidenceLevel) ||
      Settings$ConfidenceLevel <= 0 || Settings$ConfidenceLevel >= 1) stop("Invalid confidence level.")
  if (length(Settings$CoreFraction) != 1 || !is.finite(Settings$CoreFraction) ||
      Settings$CoreFraction <= 0 || Settings$CoreFraction > 1) stop("Invalid core fraction.")
  for (Name in c("QuadratureTolerance", "MaxGradient")) {
    if (length(Settings[[Name]]) != 1 || !is.finite(Settings[[Name]]) ||
        Settings[[Name]] <= 0) stop(Name, " must be positive.")
  }
  Nodes <- Settings$QuadratureNodes
  if (length(Nodes) < 2 || any(!is.finite(Nodes)) || any(Nodes < 2 | Nodes != floor(Nodes)) ||
      any(diff(Nodes) <= 0)) stop("QuadratureNodes must contain increasing integers.")
  dir.create(Settings$OutputRoot, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(Settings$OutputRoot)) stop("The output folder could not be created.")
  Settings$OutputRoot <- normalizePath(Settings$OutputRoot, winslash = "/", mustWork = TRUE)
  Settings$InputFile <- normalizePath(Settings$InputFile, winslash = "/", mustWork = TRUE)
  CheckpointDirectory <- file.path(Settings$OutputRoot, "Checkpoints")
  dir.create(CheckpointDirectory, showWarnings = FALSE)
  Bundle <- readRDS(Settings$InputFile)
  Bundle <- BbSelectImputations(Bundle, Settings$NImputations)
  Context <- BbPrepareContext(Bundle, Settings)
  NumericalSettings <- Settings[setdiff(names(Settings), c(
    "InputFile", "OutputRoot", "UseParallel", "CoreFraction", "Workers", "RetryFailed", "PrintProgress"
  ))]
  Manifest <- list(
    Version = "BtheBDkBbSelected2", InputMd5 = unname(tools::md5sum(Settings$InputFile)),
    FitSignature = Bundle$FitSignature, Settings = NumericalSettings,
    RVersion = as.character(getRversion()), GlmmTmbVersion = as.character(packageVersion("glmmTMB")),
    TmbVersion = as.character(packageVersion("TMB")), MatrixVersion = as.character(packageVersion("Matrix"))
  )
  Signature <- BbHash(Manifest)
  ManifestFile <- file.path(Settings$OutputRoot, "AnalysisManifest.rds")
  ExistingManifest <- BbReadCheckpoint(ManifestFile)
  if (!is.null(ExistingManifest) && !identical(ExistingManifest$Signature, Signature)) {
    stop("This output folder belongs to a different analysis or software version. Use a new OutputRoot.")
  }
  if (is.null(ExistingManifest)) BbSaveRds(list(Signature = Signature, Manifest = Manifest), ManifestFile)
  write.csv(Bundle$DrawIndex, file.path(Settings$OutputRoot, "SelectedImputationIndex.csv"),
            row.names = FALSE)
  BbSaveRds(Bundle, file.path(Settings$OutputRoot, "SelectedCompletedDatasets.rds"))
  Tasks <- list()
  Batches <- split(seq_len(Context$NImputations),
                   ceiling(seq_len(Context$NImputations)/Settings$BatchSize))
  for (Rule in Context$Rules) {
    for (Batch in seq_along(Batches)) {
      Tasks[[length(Tasks) + 1]] <- list(
        Rule = Rule, Indices = Batches[[Batch]],
        File = file.path(CheckpointDirectory, sprintf("%sBatch%05d.rds", Rule, Batch))
      )
    }
  }
  Pending <- logical(length(Tasks))
  for (Index in seq_along(Tasks)) {
    Task <- Tasks[[Index]]
    Cached <- BbReadCheckpoint(Task$File)
    if (is.null(Cached)) {
      Pending[Index] <- TRUE
    } else {
      if (!identical(Cached$Signature, Signature) || !identical(Cached$Indices, Task$Indices) ||
          !identical(Cached$Rule, Task$Rule) || length(Cached$Results) != length(Task$Indices)) {
        stop("A checkpoint does not match this analysis: ", Task$File)
      }
      Pending[Index] <- any(vapply(Cached$Results, function(Value) {
        is.null(Value) || (!isTRUE(Value$Ok) && Settings$RetryFailed)
      }, logical(1)))
    }
  }
  AvailableCores <- parallel::detectCores(logical = TRUE)
  if (is.na(AvailableCores)) AvailableCores <- 1
  Workers <- if (is.null(Settings$Workers)) max(1, floor(Settings$CoreFraction*AvailableCores)) else Settings$Workers
  if (length(Workers) != 1 || !is.finite(Workers) || Workers < 1 || Workers != floor(Workers)) {
    stop("Workers must be a positive integer.")
  }
  if (!Settings$UseParallel) Workers <- 1
  Workers <- min(Workers, max(1, sum(Pending)))
  if (Settings$PrintProgress) {
    cat("Selected", Context$NImputations, "of", Bundle$AnalysisSelection$AvailablePerRule,
        "imputations per rule;", Bundle$AnalysisSelection$SelectedPerChain, "from each chain.\n")
    cat("Analyzing", Context$NImputations*length(Context$Rules), "completed datasets across",
        length(Context$Rules), "rules.\n")
    cat(sum(Pending), "batches to process using", Workers, "workers.\n")
  }
  Cluster <- NULL
  if (Workers > 1 && any(Pending)) {
    Cluster <- parallel::makePSOCKcluster(Workers, outfile = file.path(Settings$OutputRoot, "WorkerLog.txt"))
    on.exit(parallel::stopCluster(Cluster), add = TRUE)
    FunctionNames <- c("BbSaveRds", "BbReadCheckpoint", "BbMarginalMoments",
                       "BbStandardize", "BbAnalyseOne", "BbRunTask", "BbWorkerTask")
    parallel::clusterExport(Cluster, FunctionNames, envir = environment(RunBthebBbAnalysis))
    BbWorkerContext <- Context
    BbWorkerSettings <- Settings
    BbWorkerSignature <- Signature
    parallel::clusterExport(Cluster, c("BbWorkerContext", "BbWorkerSettings", "BbWorkerSignature"),
                            envir = environment())
    parallel::clusterEvalQ(Cluster, {
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("glmmTMB is missing on a worker.")
      NULL
    })
  }
  Remaining <- Tasks[Pending]
  Started <- proc.time()[[3]]
  CompletedFits <- 0
  if (length(Remaining) > 0) {
    Waves <- split(seq_along(Remaining), ceiling(seq_along(Remaining)/Workers))
    for (Wave in Waves) {
      Reports <- if (is.null(Cluster)) {
        lapply(Remaining[Wave], BbRunTask, Context = Context, Settings = Settings, Signature = Signature)
      } else {
        parallel::parLapplyLB(Cluster, Remaining[Wave], BbWorkerTask)
      }
      CompletedFits <- CompletedFits + sum(vapply(Reports, function(Value) Value$N, numeric(1)))
      if (sum(vapply(Reports, function(Value) Value$Failed, numeric(1))) ==
          sum(vapply(Reports, function(Value) Value$N, numeric(1)))) {
        FailureRows <- list()
        for (Task in Remaining[Wave]) {
          Cached <- BbReadCheckpoint(Task$File)
          for (Index in seq_along(Task$Indices)) {
            FailureRows[[length(FailureRows) + 1]] <- data.frame(
              Rule = Task$Rule, Imputation = Task$Indices[Index],
              Error = Cached$Results[[Index]]$Error,
              Warnings = Cached$Results[[Index]]$Warnings
            )
          }
        }
        write.csv(do.call(rbind, FailureRows),
                  file.path(Settings$OutputRoot, "BbFitFailures.csv"), row.names = FALSE)
        stop("All fits in the latest batch group failed. See BbFitFailures.csv before resuming.")
      }
      if (Settings$PrintProgress) {
        cat("Checkpointed", CompletedFits, "fits in this run;",
            sum(vapply(Reports, function(Value) Value$Failed, numeric(1))),
            "unresolved in the latest batch group; elapsed",
            round((proc.time()[[3]] - Started)/60, 1), "minutes.\n")
        flush.console()
      }
    }
  }
  Result <- BbCollect(Tasks, Context, Settings, Signature)
  if (Settings$PrintProgress) {
    cat("Analysis complete. Results:", Settings$OutputRoot, "\n")
    print(Result$Pooled[Result$Pooled$Visit == Context$NVisits &
                          Result$Pooled$Summary == "BtheB minus TAU",
                        c("Rule", "Estimate", "SE", "Lower", "Upper", "MCSE", "MCSEOverSE")],
          row.names = FALSE)
  }
  invisible(Result)
}

# Source this file to run the complete selection and analysis.
BbSettings <- DefaultBthebBbSettings()
BbResult <- RunBthebBbAnalysis(BbSettings)