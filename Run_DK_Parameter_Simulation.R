# Simulation study 1: complete-data DK Gaussian copula parameter recovery.
#
# How to run
#   Source this file. The default is the full study: 300 replicates in each
#   of four scenarios (1,200 datasets). Set Stage to "Pilot" for five per
#   scenario. Compatible pilot fits are reused when Stage becomes "Full".
#
#   Workers = max(1, floor(CoreFraction*available cores)), capped by pending tasks.
#   The default CoreFraction is 0.60.
#   Available cores are detected at run time, including CPU restrictions.
#   Replicates run in parallel; chains within each replicate are sequential.
#   Workers take the next dataset immediately when they finish. Cumulative
#   reports refresh after each worker-count group of completed attempts.
#   A console heartbeat is printed every 30 seconds while jobs are running.
#   All post-warmup parameter draws are saved. There is no thinning.
#
# Outputs
#   The console prints the report folder. ParameterPerformance.tex is ready
#   for \input{} in a manuscript using booktabs, amsmath, pdflscape, and siunitx.
#   Numeric columns and coverage MCSEs are aligned at their decimal points.
#   CSV files contain performance measures, Monte Carlo errors, diagnostics,
#   and the failure log.
#   Replicate folders contain the simulated data and all parameter draws.
#   An interrupted run resumes from completed chains when sourced again.
#   Changing model, prior, sampler, seed, or software settings automatically
#   creates a separate study folder. Changing Stage, replication count, or
#   available cores does not invalidate completed compatible fits.
#   Compatible results from statistical specification 1.1.0 are reused.
#   Keep the same OutputRoot and statistical settings when resuming.
#   With all scheduled fits completed successfully, sourcing this file only
#   regenerates the reports. Table formatting does not change the study identity.
#
# Model and scope
#   Y in {0,...,63}; visits at months 2, 3, 5, 8; complete data only.
#   eta[i,j] = alpha[j] + beta[j]*A[i] + gamma[j,1]*B[i].
#   Phi[j] and Rho are summarized after transforming each posterior draw.
#   No imputation, BB analysis, median-anchor recovery, or treatment-mean
#   comparison is performed in this first simulation study.
#
# Priors
#   Independent Normal(0, 5^2) priors for alpha, beta, and gamma.
#   log(phi) ~ Normal(0, 2^2); rho ~ Uniform(0, 1).
#   Normal second arguments above are variances; settings below use SDs.
#   The logit(rho) target includes the Uniform-prior transformation Jacobian.
#   A changed prior creates a new study folder and requires new fits.
#
# The DK calculations and sampler follow the BtheB imputation implementation,
# with the priors specified above. This script is self-contained.
# Required CRAN packages are installed automatically if missing.
# To load functions without running: options(DkSimulationRun = FALSE).
# Version: 1.1.1, 2026-09-23 (statistical specification 1.1.0).

DefaultDkSimulationSettings <- function() {
  list(
    Stage = "Full", # Use "Pilot" for five replicates per scenario.
    OutputRoot = "DkParameterSimulation",
    PilotReplicates = 5, # Per scenario: 20 pilot replicates in total.
    FullReplicates = 300, # Per scenario, not across all scenarios.
    PilotScenarios = "All", # Or a vector such as c("ApplicationN100").
    SampleSizes = c(100, 200),
    VisitMonths = c(2, 3, 5, 8),
    ScoreMaximum = 63,
    Tau = 0.5,
    BaselineScale = 10,
    InterceptTruth = c(-0.95, -1.11, -1.48, -1.74),
    TreatmentTruth = c(-0.38, -0.43, -0.28, -0.16),
    BaselineTruth = c(0.44, 0.60, 0.63, 0.52),
    ShapeTruth = c(2.32, 1.99, 2.03, 2.15),
    RhoTruth = 0.78,
    FloorInterceptShift = -2,
    NChains = 4,
    Warmup = 5000,
    Sample = 10000,
    AdaptEvery = 100,
    PitSweeps = 2,
    Seed = 20260924,
    DesignSeed = 20260923,
    UseParallel = TRUE,
    CoreFraction = 0.60,
    RetryFailures = TRUE, # Same dataset and seeds; no replacement.
    InstallMissingPackages = TRUE,
    PrintProgress = TRUE,
    ProgressEvery = 1000,
    RhatThreshold = 1.01,
    EssThreshold = 400,
    InterceptPriorSd = 5,
    TreatmentPriorSd = 5,
    BaselinePriorSd = 5,
    LogShapePriorMean = 0,
    LogShapePriorSd = 2,
    RhoPrior = "Uniform" # Uniform on rho in (0, 1).
  )
}

DkCheckSimulationSettings <- function(Settings) {
  if (!Settings$Stage %in% c("Pilot", "Full")) stop("Stage must be 'Pilot' or 'Full'.")
  PositiveIntegers <- c("PilotReplicates", "FullReplicates", "NChains", "Warmup",
                        "Sample", "AdaptEvery", "PitSweeps", "ProgressEvery")
  for (Name in PositiveIntegers) {
    Value <- Settings[[Name]]
    if (length(Value) != 1 || !is.finite(Value) || Value < 1 || Value != floor(Value)) {
      stop(Name, " must be a positive integer.")
    }
  }
  if (Settings$NChains < 2 || Settings$Sample < 4) stop("Use at least two chains and four saved draws per chain.")
  if (length(Settings$SampleSizes) == 0 || anyDuplicated(Settings$SampleSizes) ||
      any(!is.finite(Settings$SampleSizes)) || any(Settings$SampleSizes < 10) ||
      any(Settings$SampleSizes %% 2 != 0)) stop("SampleSizes must be distinct even integers, at least 10.")
  if (length(Settings$VisitMonths) < 2 || any(!is.finite(Settings$VisitMonths)) ||
      any(diff(Settings$VisitMonths) <= 0)) stop("VisitMonths must be strictly increasing.")
  for (Name in c("InterceptTruth", "TreatmentTruth", "BaselineTruth", "ShapeTruth")) {
    if (length(Settings[[Name]]) != length(Settings$VisitMonths) ||
        any(!is.finite(Settings[[Name]]))) stop(Name, " needs one finite value per visit.")
  }
  if (any(Settings$ShapeTruth <= 0) || length(Settings$RhoTruth) != 1 ||
      !is.finite(Settings$RhoTruth) || Settings$RhoTruth <= 0 || Settings$RhoTruth >= 1) {
    stop("Shape truths must be positive and RhoTruth must lie strictly between zero and one.")
  }
  if (length(Settings$ScoreMaximum) != 1 || !is.finite(Settings$ScoreMaximum) ||
      Settings$ScoreMaximum < 2 || Settings$ScoreMaximum != floor(Settings$ScoreMaximum)) {
    stop("ScoreMaximum must be an integer of at least two.")
  }
  if (!identical(Settings$Tau, 0.5)) stop("This study uses the median parameterization, Tau = 0.5.")
  for (Name in c("BaselineScale", "InterceptPriorSd", "TreatmentPriorSd", "BaselinePriorSd",
                 "LogShapePriorSd", "EssThreshold")) {
    if (length(Settings[[Name]]) != 1 || !is.finite(Settings[[Name]]) || Settings[[Name]] <= 0) {
      stop(Name, " must be finite and positive.")
    }
  }
  for (Name in c("LogShapePriorMean", "FloorInterceptShift", "RhatThreshold")) {
    if (length(Settings[[Name]]) != 1 || !is.finite(Settings[[Name]])) stop(Name, " must be finite.")
  }
  if (Settings$FloorInterceptShift >= 0) stop("FloorInterceptShift must be negative.")
  if (!identical(Settings$RhoPrior, "Uniform")) stop("RhoPrior must be 'Uniform'.")
  if (Settings$RhatThreshold <= 1) stop("RhatThreshold must exceed one.")
  if (length(Settings$CoreFraction) != 1 || !is.finite(Settings$CoreFraction) ||
      Settings$CoreFraction <= 0 || Settings$CoreFraction > 0.75) {
    stop("CoreFraction must be in (0, 0.75].")
  }
  for (Name in c("Seed", "DesignSeed")) {
    Value <- Settings[[Name]]
    if (length(Value) != 1 || !is.finite(Value) || Value < 1 ||
        Value > .Machine$integer.max || Value != floor(Value)) stop(Name, " is not a valid positive R seed.")
  }
  for (Name in c("UseParallel", "RetryFailures", "InstallMissingPackages", "PrintProgress")) {
    if (!is.logical(Settings[[Name]]) || length(Settings[[Name]]) != 1 || is.na(Settings[[Name]])) {
      stop(Name, " must be TRUE or FALSE.")
    }
  }
  if (!is.character(Settings$OutputRoot) || length(Settings$OutputRoot) != 1 ||
      !nzchar(Settings$OutputRoot)) stop("OutputRoot must be a nonempty path.")
  invisible(TRUE)
}

DkEnsurePackages <- function(Settings) {
  Required <- c("HSAUR2", "posterior", "parallelly", "filelock")
  Missing <- Required[!vapply(Required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(Missing) > 0 && Settings$InstallMissingPackages) {
    message("Installing missing packages: ", paste(Missing, collapse = ", "))
    utils::install.packages(Missing, repos = "https://cloud.r-project.org")
  }
  Missing <- Required[!vapply(Required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(Missing) > 0) stop("Install these packages and rerun: ", paste(Missing, collapse = ", "))
  invisible(Required)
}

DkWithSeed <- function(Seed, Code) {
  OldKind <- RNGkind()
  HadSeed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (HadSeed) OldSeed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(OldKind))
    if (HadSeed) assign(".Random.seed", OldSeed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  if (length(Seed) == 1) set.seed(Seed) else assign(".Random.seed", Seed, envir = .GlobalEnv)
  force(Code)
}

DkHashObject <- function(Object) {
  Path <- tempfile("DkHash")
  on.exit(unlink(Path), add = TRUE)
  saveRDS(Object, Path, compress = FALSE, version = 2)
  unname(tools::md5sum(Path))
}

DkCommitFile <- function(Temporary, Path) {
  Backup <- paste0(Path, ".previous")
  if (file.exists(Path)) {
    if (file.exists(Backup)) unlink(Backup)
    if (!file.rename(Path, Backup)) stop("Cannot replace output file: ", Path)
  }
  if (!file.rename(Temporary, Path)) {
    if (file.exists(Backup)) file.rename(Backup, Path)
    stop("Cannot finish saving: ", Path)
  }
  if (file.exists(Backup)) unlink(Backup)
  invisible(Path)
}

DkSaveRds <- function(Object, Path) {
  dir.create(dirname(Path), recursive = TRUE, showWarnings = FALSE)
  Temporary <- tempfile("DkSave", tmpdir = dirname(Path))
  on.exit(unlink(Temporary), add = TRUE)
  saveRDS(Object, Temporary, compress = "gzip", version = 2)
  DkCommitFile(Temporary, Path)
}

DkReadRds <- function(Path) {
  if (!file.exists(Path) && file.exists(paste0(Path, ".previous"))) {
    file.rename(paste0(Path, ".previous"), Path)
  }
  if (!file.exists(Path)) return(NULL)
  tryCatch(readRDS(Path), error = function(Error) {
    warning("Unreadable checkpoint will be rebuilt from the same seeds: ", Path, call. = FALSE)
    NULL
  })
}

DkWriteCsv <- function(Object, Path) {
  Temporary <- tempfile("DkCsv", tmpdir = dirname(Path))
  on.exit(unlink(Temporary), add = TRUE)
  utils::write.csv(Object, Temporary, row.names = FALSE, na = "")
  DkCommitFile(Temporary, Path)
}

DkWriteLines <- function(Lines, Path) {
  Temporary <- tempfile("DkText", tmpdir = dirname(Path))
  on.exit(unlink(Temporary), add = TRUE)
  writeLines(Lines, Temporary, useBytes = TRUE)
  DkCommitFile(Temporary, Path)
}

DkParameterMetadata <- function(Settings) {
  J <- length(Settings$VisitMonths)
  data.frame(
    Parameter = c(paste0("Alpha", seq_len(J)), paste0("Beta", seq_len(J)),
                  paste0("Gamma", seq_len(J)), paste0("Phi", seq_len(J)), "Rho"),
    Latex = c(paste0("\\alpha_{", seq_len(J), "}"), paste0("\\beta_{", seq_len(J), "}"),
              paste0("\\gamma_{", seq_len(J), "1}"), paste0("\\phi_{", seq_len(J), "}"), "\\rho"),
    Description = c(paste("Intercept, month", Settings$VisitMonths),
                    paste("Treatment coefficient, month", Settings$VisitMonths),
                    paste("Baseline coefficient, month", Settings$VisitMonths),
                    paste("Shape, month", Settings$VisitMonths), "One-month latent correlation"),
    Family = c(rep(c("Intercept", "Treatment", "Baseline", "Shape"), each = J), "Correlation"),
    Visit = c(rep(seq_len(J), 4), NA_real_),
    Month = c(rep(Settings$VisitMonths, 4), NA_real_),
    stringsAsFactors = FALSE
  )
}

DkPrepareStudy <- function(Settings) {
  DataEnvironment <- new.env(parent = baseenv())
  utils::data("BtheB", package = "HSAUR2", envir = DataEnvironment)
  Data <- DataEnvironment$BtheB
  Followup <- c("bdi.2m", "bdi.3m", "bdi.5m", "bdi.8m")
  Included <- rowSums(!is.na(Data[, Followup])) > 0
  BaselinePool <- Data$bdi.pre[Included]
  if (any(!is.finite(BaselinePool))) stop("The BtheB baseline pool is not fully observed.")
  BaselineCenter <- mean(BaselinePool)
  Designs <- DkWithSeed(Settings$DesignSeed, {
    lapply(Settings$SampleSizes, function(N) {
      Baseline <- sample(BaselinePool, N, replace = TRUE)
      Treatment <- sample(rep(c(0, 1), each = N/2))
      Design <- data.frame(Participant = seq_len(N), Treatment = Treatment,
                           BaselineScore = Baseline,
                           Baseline = (Baseline - BaselineCenter)/Settings$BaselineScale)
      if (qr(cbind(1, Design$Treatment, Design$Baseline))$rank < 3) stop("The fixed design is rank deficient.")
      Design
    })
  })
  names(Designs) <- paste0("N", Settings$SampleSizes)
  Scenarios <- list()
  for (N in Settings$SampleSizes) {
    for (Setting in c("Application", "Floor")) {
      Shift <- if (Setting == "Floor") Settings$FloorInterceptShift else 0
      Name <- paste0(Setting, "N", N)
      Truth <- c(Settings$InterceptTruth + Shift, Settings$TreatmentTruth,
                 Settings$BaselineTruth, Settings$ShapeTruth, Settings$RhoTruth)
      names(Truth) <- DkParameterMetadata(Settings)$Parameter
      Scenarios[[Name]] <- list(Name = Name, N = N, Setting = Setting,
                                InterceptShift = Shift, Truth = Truth,
                                Design = Designs[[paste0("N", N)]])
    }
  }
  if (Settings$Stage == "Pilot" && !identical(Settings$PilotScenarios, "All") &&
      (length(Settings$PilotScenarios) == 0 || anyDuplicated(Settings$PilotScenarios) ||
       any(!Settings$PilotScenarios %in% names(Scenarios)))) {
    stop("PilotScenarios must name scenarios such as ApplicationN100 or FloorN100.")
  }
  list(Scenarios = Scenarios, Designs = Designs, BaselineCenter = BaselineCenter,
       BaselinePool = BaselinePool, Metadata = DkParameterMetadata(Settings))
}

DkMakeCompleteModel <- function(Y, Design, Settings) {
  J <- length(Settings$VisitMonths)
  N <- nrow(Design)
  if (!identical(dim(Y), c(as.integer(N), as.integer(J))) || anyNA(Y) ||
      any(!is.finite(Y)) || any(Y < 0 | Y > Settings$ScoreMaximum | Y != floor(Y))) {
    stop("The simulated outcome matrix must be complete and within the integer support.")
  }
  Index <- list(Intercept = seq_len(J), Treatment = J + seq_len(J),
                Baseline = 2*J + seq_len(J), LogShape = 3*J + seq_len(J), LogitRho = 4*J + 1)
  Names <- c(paste0("Intercept[", seq_len(J), "]"), paste0("Treatment[", seq_len(J), "]"),
             paste0("Baseline[", seq_len(J), "]"), paste0("LogShape[", seq_len(J), "]"), "LogitRho")
  list(Y = Y, Observed = matrix(TRUE, N, J),
       Patterns = list(list(Subjects = seq_len(N), Visits = seq_len(J))),
       Treatment = Design$Treatment, Baseline = Design$Baseline,
       NSubjects = N, NVisits = J, Index = Index, ParameterNames = Names,
       VisitMonths = Settings$VisitMonths, ScoreMaximum = Settings$ScoreMaximum, Tau = Settings$Tau)
}

DkGenerateOutcomes <- function(Scenario, Settings) {
  J <- length(Settings$VisitMonths)
  Truth <- Scenario$Truth
  Correlation <- MakeContinuousAr1(Settings$VisitMonths, Truth[4*J + 1])
  Z <- matrix(rnorm(Scenario$N*J), Scenario$N, J)%*%chol(Correlation)
  Y <- matrix(0, Scenario$N, J, dimnames = list(NULL, paste0("Month", Settings$VisitMonths)))
  for (Visit in seq_len(J)) {
    Eta <- Truth[Visit] + Truth[J + Visit]*Scenario$Design$Treatment +
      Truth[2*J + Visit]*Scenario$Design$Baseline
    # Upper-tail log probabilities avoid rounding Phi(Z) to one.
    Y[, Visit] <- DkQuantileFromLogSurvival(pnorm(Z[, Visit], lower.tail = FALSE, log.p = TRUE),
                                            Eta, log(Truth[3*J + Visit]), Settings$ScoreMaximum, Settings$Tau)
  }
  Y
}

DkScenarioProbabilities <- function(Study, Settings) {
  Rows <- list()
  J <- length(Settings$VisitMonths)
  for (Scenario in Study$Scenarios) {
    Truth <- Scenario$Truth
    for (Arm in c(0, 1)) {
      # Actual fixed baseline distribution in this assigned arm.
      B <- Scenario$Design$Baseline[Scenario$Design$Treatment == Arm]
      for (Visit in seq_len(J)) {
        Eta <- Truth[Visit] + Truth[J + Visit]*Arm + Truth[2*J + Visit]*B
        LogShape <- log(Truth[3*J + Visit])
        Rows[[length(Rows) + 1]] <- data.frame(
          Scenario = Scenario$Name, N = Scenario$N, Setting = Scenario$Setting,
          Treatment = Arm, Visit = Visit, Month = Settings$VisitMonths[Visit],
          ProbabilityZero = mean(DkPmf(0, Eta, LogShape, Settings$ScoreMaximum)),
          ProbabilityMaximum = mean(DkPmf(Settings$ScoreMaximum, Eta, LogShape, Settings$ScoreMaximum)),
          ProbabilityAtMostFive = mean(DkCdf(min(5, Settings$ScoreMaximum), Eta, LogShape, Settings$ScoreMaximum))
        )
      }
    }
  }
  do.call(rbind, Rows)
}

# ---- Stable DK calculations and sampler -----------------------------------

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
  # Parameter proposals do not change PitEta. UpdateDkPit maintains this
  # cached term after accepted auxiliary proposals, including partial rows.
  if (is.null(Current)) {
    State$LogJacobianBySubject <- rowSums(-Softplus(-PitEta) - Softplus(PitEta), na.rm = TRUE)
  }
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
  # Start at rho = 0.5; chains are subsequently jittered independently.
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

# ---- Chain checkpoints and parameter summaries ----------------------------

RunDkSimulationChain <- function(Chain, Model, Settings, Initial, Signature,
                                 Directory, Seed, TaskLabel) {
  CacheFile <- file.path(Directory, paste0("Chain", Chain, ".rds"))
  Cached <- DkReadRds(CacheFile)
  if (!is.null(Cached)) {
    if (!identical(Cached$Signature, Signature) || !identical(Cached$Seed, Seed)) {
      stop("Chain checkpoint does not match this dataset and RNG stream: ", CacheFile)
    }
    if (!identical(dim(Cached$Draws), c(as.integer(Settings$Sample), length(Model$ParameterNames))) ||
        any(!is.finite(Cached$Draws))) stop("Invalid draws in checkpoint: ", CacheFile)
    return(Cached)
  }
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", Seed, envir = .GlobalEnv)
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
    }
    if (Iteration %% Settings$ProgressEvery == 0 || Iteration == TotalIterations) {
      DkWriteLines(paste(format(Sys.time(), tz = "UTC", usetz = TRUE), TaskLabel,
                         "chain", Chain, "iteration", Iteration, "of", TotalIterations),
                   file.path(Directory, "Progress.txt"))
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
  Result <- list(Chain = Chain, Signature = Signature, Seed = Seed, Draws = Draws,
                 Acceptance = Acceptance,
                 RuntimeSeconds = proc.time()[[3]] - StartTime)
  DkSaveRds(Result, CacheFile)
  if (Settings$PrintProgress) cat(TaskLabel, "chain", Chain, "completed in",
                                  round(Result$RuntimeSeconds, 1), "seconds.\n")
  Result
}

DkApplyDiagnosticFlags <- function(Summary, Settings) {
  Summary$DiagnosticFlag <- !is.finite(Summary$Rhat) | Summary$Rhat > Settings$RhatThreshold |
    !is.finite(Summary$EssBulk) | Summary$EssBulk < Settings$EssThreshold |
    !is.finite(Summary$EssTail) | Summary$EssTail < Settings$EssThreshold |
    !is.finite(Summary$McseMean)
  Summary
}

DkSummariseSimulationFit <- function(Chains, Model, Scenario, Settings) {
  Names <- names(Scenario$Truth)
  DrawArray <- array(NA_real_, c(Settings$Sample, Settings$NChains, length(Names)),
                     dimnames = list(iteration = seq_len(Settings$Sample),
                                     chain = seq_len(Settings$NChains), variable = Names))
  for (Chain in seq_len(Settings$NChains)) {
    Draws <- Chains[[Chain]]$Draws
    Draws[, Model$Index$LogShape] <- exp(Draws[, Model$Index$LogShape, drop = FALSE])
    Draws[, Model$Index$LogitRho] <- plogis(Draws[, Model$Index$LogitRho])
    if (any(!is.finite(Draws))) stop("Nonfinite natural-scale parameter draws.")
    DrawArray[, Chain, ] <- Draws
  }
  Diagnostics <- posterior::summarise_draws(posterior::as_draws_array(DrawArray),
                                            "rhat", "ess_bulk", "ess_tail", "mcse_mean")
  Combined <- matrix(DrawArray, ncol = length(Names), dimnames = list(NULL, Names))
  Intervals <- apply(Combined, 2, quantile, probs = c(0.025, 0.975), type = 8, names = FALSE)
  Match <- match(Names, Diagnostics$variable)
  Summary <- data.frame(
    Parameter = Names, Truth = as.numeric(Scenario$Truth),
    Estimate = colMeans(Combined), PosteriorSd = apply(Combined, 2, sd),
    Lower95 = Intervals[1, ], Upper95 = Intervals[2, ],
    Rhat = Diagnostics$rhat[Match], EssBulk = Diagnostics$ess_bulk[Match],
    EssTail = Diagnostics$ess_tail[Match], McseMean = Diagnostics$mcse_mean[Match],
    NDraws = Settings$NChains*Settings$Sample, row.names = NULL
  )
  Summary$Covered <- Summary$Lower95 <= Summary$Truth & Summary$Upper95 >= Summary$Truth
  Summary$IntervalWidth <- Summary$Upper95 - Summary$Lower95
  if (any(!is.finite(as.matrix(Summary[, c("Estimate", "PosteriorSd", "Lower95", "Upper95")])))) {
    stop("Nonfinite posterior summary; the fit cannot contribute to parameter performance.")
  }
  DkApplyDiagnosticFlags(Summary, Settings)
}

DkMakeTasks <- function(Study, Settings, StudyDirectory, Signature) {
  Selected <- if (Settings$Stage == "Pilot" && !identical(Settings$PilotScenarios, "All")) {
    Settings$PilotScenarios
  } else names(Study$Scenarios)
  Replicates <- if (Settings$Stage == "Pilot") Settings$PilotReplicates else Settings$FullReplicates
  DkWithSeed(Settings$Seed, {
    ScenarioSeed <- get(".Random.seed", envir = .GlobalEnv)
    Tasks <- list()
    for (Name in names(Study$Scenarios)) {
      CurrentSeed <- ScenarioSeed
      if (Name %in% Selected) {
        for (Replicate in seq_len(Replicates)) {
          DataSeed <- CurrentSeed
          ChainSeeds <- vector("list", Settings$NChains)
          for (Chain in seq_len(Settings$NChains)) {
            CurrentSeed <- parallel::nextRNGSubStream(CurrentSeed)
            ChainSeeds[[Chain]] <- CurrentSeed
          }
          Directory <- file.path(StudyDirectory, "Replicates", Name, sprintf("Replicate%04d", Replicate))
          Tasks[[length(Tasks) + 1]] <- list(Scenario = Name, Replicate = Replicate,
                                             DataSeed = DataSeed, ChainSeeds = ChainSeeds,
                                             Directory = Directory, StudySignature = Signature)
          CurrentSeed <- parallel::nextRNGSubStream(CurrentSeed)
        }
      }
      ScenarioSeed <- parallel::nextRNGStream(ScenarioSeed)
    }
    Tasks
  })
}

DkReadTaskResult <- function(Task) {
  Result <- DkReadRds(file.path(Task$Directory, "Result.rds"))
  if (is.null(Result)) return(NULL)
  if (!identical(Result$StudySignature, Task$StudySignature) ||
      !identical(Result$Scenario, Task$Scenario) || !identical(Result$Replicate, Task$Replicate)) {
    stop("Mismatched replicate checkpoint: ", Task$Directory)
  }
  Result
}

RunDkSimulationReplicate <- function(Task, Context) {
  Settings <- Context$Settings
  Scenario <- Context$Study$Scenarios[[Task$Scenario]]
  Cached <- DkReadTaskResult(Task)
  if (!is.null(Cached) && (isTRUE(Cached$Success) || !Settings$RetryFailures)) return(Cached)
  dir.create(Task$Directory, recursive = TRUE, showWarnings = FALSE)
  Warnings <- character(0)
  Start <- proc.time()[[3]]
  Label <- paste(Task$Scenario, sprintf("replicate %04d", Task$Replicate))
  Result <- tryCatch(withCallingHandlers(DkWithSeed(Task$DataSeed, {
    DataFile <- file.path(Task$Directory, "Data.rds")
    Data <- DkReadRds(DataFile)
    if (is.null(Data)) {
      Y <- DkGenerateOutcomes(Scenario, Settings)
      Data <- list(Y = Y, Design = Scenario$Design, DataSeed = Task$DataSeed,
                   ChainSeeds = Task$ChainSeeds, StudySignature = Task$StudySignature)
      DkSaveRds(Data, DataFile)
    }
    if (!identical(Data$StudySignature, Task$StudySignature) ||
        !identical(Data$DataSeed, Task$DataSeed) || !identical(Data$ChainSeeds, Task$ChainSeeds) ||
        !identical(Data$Design, Scenario$Design)) stop("The saved dataset does not match its task.")
    Model <- DkMakeCompleteModel(Data$Y, Data$Design, Settings)
    Signature <- DkHashObject(list(Study = Task$StudySignature, Y = Data$Y,
                                   Design = Data$Design, Seeds = Task$ChainSeeds))
    # Initialization uses this simulated dataset and the specified priors only.
    Initial <- InitialDkParameters(Model, Settings)
    Chains <- lapply(seq_len(Settings$NChains), function(Chain) {
      RunDkSimulationChain(Chain, Model, Settings, Initial, Signature,
                           Task$Directory, Task$ChainSeeds[[Chain]], Label)
    })
    Summary <- DkSummariseSimulationFit(Chains, Model, Scenario, Settings)
    list(Success = TRUE, Summary = Summary,
         Acceptance = do.call(rbind, lapply(Chains, function(Chain) Chain$Acceptance)),
         SamplingSeconds = sum(vapply(Chains, function(Chain) Chain$RuntimeSeconds, numeric(1))),
         FitSignature = Signature, Error = "")
  }), warning = function(Warning) {
    Warnings <<- unique(c(Warnings, conditionMessage(Warning)))
    invokeRestart("muffleWarning")
  }), error = function(Error) {
    list(Success = FALSE, Summary = NULL, Acceptance = NULL,
         SamplingSeconds = NA_real_, FitSignature = NA_character_, Error = conditionMessage(Error))
  })
  Result$Scenario <- Task$Scenario
  Result$Replicate <- Task$Replicate
  Result$StudySignature <- Task$StudySignature
  Result$Warnings <- Warnings
  Result$AttemptSeconds <- proc.time()[[3]] - Start
  Result$FinishedUtc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  DkSaveRds(Result, file.path(Task$Directory, "Result.rds"))
  DkWriteLines(c(paste(Result$FinishedUtc, Label, if (Result$Success) "completed" else "FAILED"),
                 if (!Result$Success) Result$Error), file.path(Task$Directory, "Progress.txt"))
  if (Settings$PrintProgress) {
    cat(Label, if (Result$Success) "completed." else paste("FAILED:", Result$Error), "\n")
  }
  Result
}

# ---- Performance measures and manuscript tables ---------------------------

DkSafeMinimum <- function(Values) {
  Values <- Values[is.finite(Values)]
  if (length(Values) == 0) NA_real_ else min(Values)
}

DkSafeMaximum <- function(Values) {
  Values <- Values[is.finite(Values)]
  if (length(Values) == 0) NA_real_ else max(Values)
}

DkMcseAverage <- function(Values) {
  if (length(Values) < 2) NA_real_ else sd(Values)/sqrt(length(Values))
}

DkCoverageWilson <- function(Covered) {
  N <- length(Covered)
  if (N == 0) return(c(NA_real_, NA_real_))
  P <- mean(Covered)
  Z <- qnorm(0.975)
  Center <- (P + Z^2/(2*N))/(1 + Z^2/N)
  Radius <- Z*sqrt(P*(1 - P)/N + Z^2/(4*N^2))/(1 + Z^2/N)
  100*c(Center - Radius, Center + Radius)
}

DkCollectReports <- function(Tasks, Study, Settings, Results = NULL) {
  if (is.null(Results)) Results <- lapply(Tasks, DkReadTaskResult)
  if (length(Results) != length(Tasks)) stop("The result snapshot does not match the task list.")
  StatusRows <- SummaryRows <- AcceptanceRows <- WarningRows <- list()
  for (TaskIndex in seq_along(Tasks)) {
    Task <- Tasks[[TaskIndex]]
    Result <- Results[[TaskIndex]]
    Scenario <- Study$Scenarios[[Task$Scenario]]
    Success <- !is.null(Result) && isTRUE(Result$Success)
    Summary <- if (Success) DkApplyDiagnosticFlags(Result$Summary, Settings) else NULL
    State <- if (is.null(Result)) "Pending" else if (Success) "Completed" else "Failed"
    StatusRows[[length(StatusRows) + 1]] <- data.frame(
      Scenario = Task$Scenario, Replicate = Task$Replicate, N = Scenario$N,
      Setting = Scenario$Setting, Status = State,
      DiagnosticFlag = if (Success) any(Summary$DiagnosticFlag) else NA,
      FlaggedParameters = if (Success) sum(Summary$DiagnosticFlag) else NA_real_,
      MaximumRhat = if (Success) DkSafeMaximum(Summary$Rhat) else NA_real_,
      MinimumEssBulk = if (Success) DkSafeMinimum(Summary$EssBulk) else NA_real_,
      MinimumEssTail = if (Success) DkSafeMinimum(Summary$EssTail) else NA_real_,
      NonfiniteDiagnostics = if (Success) sum(!is.finite(as.matrix(Summary[, c("Rhat", "EssBulk", "EssTail", "McseMean")]))) else NA_real_,
      SamplingSeconds = if (!is.null(Result)) Result$SamplingSeconds else NA_real_,
      AttemptSeconds = if (!is.null(Result)) Result$AttemptSeconds else NA_real_,
      WarningCount = if (!is.null(Result)) length(Result$Warnings) else 0,
      Error = if (!is.null(Result)) Result$Error else "",
      FinishedUtc = if (!is.null(Result)) Result$FinishedUtc else "",
      stringsAsFactors = FALSE
    )
    if (Success) {
      SummaryRows[[length(SummaryRows) + 1]] <- cbind(
        data.frame(Scenario = Task$Scenario, Replicate = Task$Replicate, N = Scenario$N), Summary)
      AcceptanceRows[[length(AcceptanceRows) + 1]] <- cbind(
        data.frame(Scenario = Task$Scenario, Replicate = Task$Replicate), Result$Acceptance)
    }
    if (!is.null(Result) && length(Result$Warnings) > 0) {
      WarningRows[[length(WarningRows) + 1]] <- data.frame(
        Scenario = Task$Scenario, Replicate = Task$Replicate, Warning = Result$Warnings)
    }
  }
  Status <- do.call(rbind, StatusRows)
  Summaries <- if (length(SummaryRows) > 0) do.call(rbind, SummaryRows) else data.frame(
    Scenario = character(), Replicate = numeric(), Parameter = character())
  Acceptance <- if (length(AcceptanceRows) > 0) do.call(rbind, AcceptanceRows) else data.frame(
    Scenario = character(), Replicate = numeric(), Chain = numeric(), Block = character(),
    Acceptance = numeric(), ProposalScale = numeric())
  Warnings <- if (length(WarningRows) > 0) do.call(rbind, WarningRows) else data.frame(
    Scenario = character(), Replicate = numeric(), Warning = character())
  PerformanceRows <- DiagnosticRows <- list()
  for (Name in unique(Status$Scenario)) {
    Scenario <- Study$Scenarios[[Name]]
    Current <- Status[Status$Scenario == Name, , drop = FALSE]
    Completed <- sum(Current$Status == "Completed")
    Attempted <- sum(Current$Status != "Pending")
    Failed <- sum(Current$Status == "Failed")
    Flagged <- sum(Current$DiagnosticFlag, na.rm = TRUE)
    DiagnosticRows[[length(DiagnosticRows) + 1]] <- data.frame(
      Scenario = Name, Setting = Scenario$Setting, N = Scenario$N, Planned = nrow(Current),
      Attempted = Attempted, Completed = Completed, Failed = Failed,
      Pending = sum(Current$Status == "Pending"), Flagged = Flagged,
      FlaggedPercent = if (Completed > 0) 100*Flagged/Completed else NA_real_,
      MaximumRhat = DkSafeMaximum(Current$MaximumRhat),
      MinimumEssBulk = DkSafeMinimum(Current$MinimumEssBulk),
      MinimumEssTail = DkSafeMinimum(Current$MinimumEssTail),
      ReplicatesWithNonfiniteDiagnostics = sum(Current$NonfiniteDiagnostics > 0, na.rm = TRUE),
      MeanSamplingSeconds = if (Completed > 0) mean(Current$SamplingSeconds[Current$Status == "Completed"]) else NA_real_
    )
    for (Parameter in Study$Metadata$Parameter) {
      Values <- Summaries[Summaries$Scenario == Name & Summaries$Parameter == Parameter, , drop = FALSE]
      Count <- nrow(Values)
      Truth <- unname(Scenario$Truth[Parameter])
      Errors <- Values$Estimate - Truth
      SquaredErrors <- Errors^2
      Rmse <- if (Count > 0) sqrt(mean(SquaredErrors)) else NA_real_
      RmseMcse <- if (Count < 2) NA_real_ else if (Rmse == 0) 0 else DkMcseAverage(SquaredErrors)/(2*Rmse)
      Coverage <- if (Count > 0) mean(Values$Covered) else NA_real_
      Wilson <- DkCoverageWilson(Values$Covered)
      Meta <- Study$Metadata[Study$Metadata$Parameter == Parameter, , drop = FALSE]
      PerformanceRows[[length(PerformanceRows) + 1]] <- data.frame(
        Scenario = Name, Setting = Scenario$Setting, N = Scenario$N,
        Parameter = Parameter, Latex = Meta$Latex, Description = Meta$Description,
        Family = Meta$Family, Visit = Meta$Visit, Month = Meta$Month, Truth = Truth,
        Planned = nrow(Current), Attempted = Attempted, Contributing = Count,
        Failed = Failed, FlaggedReplicates = Flagged,
        ParameterFlagged = if (Count > 0) sum(Values$DiagnosticFlag) else 0,
        MeanEstimate = if (Count > 0) mean(Values$Estimate) else NA_real_,
        Bias = if (Count > 0) mean(Errors) else NA_real_, BiasMcse = DkMcseAverage(Errors),
        Rmse = Rmse, RmseMcse = RmseMcse,
        EmpiricalSd = if (Count > 1) sd(Values$Estimate) else NA_real_,
        MeanPosteriorSd = if (Count > 0) mean(Values$PosteriorSd) else NA_real_,
        MeanPosteriorSdMcse = DkMcseAverage(Values$PosteriorSd),
        CoveragePercent = 100*Coverage,
        CoverageMcsePercent = if (Count > 0) 100*sqrt(Coverage*(1 - Coverage)/Count) else NA_real_,
        CoverageWilsonLowerPercent = Wilson[1], CoverageWilsonUpperPercent = Wilson[2],
        MeanIntervalWidth = if (Count > 0) mean(Values$IntervalWidth) else NA_real_,
        MeanIntervalWidthMcse = DkMcseAverage(Values$IntervalWidth),
        MeanMcmcMcse = if (Count > 0 && all(is.finite(Values$McseMean))) mean(Values$McseMean) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  list(Status = Status, Summaries = Summaries, Acceptance = Acceptance, Warnings = Warnings,
       Performance = do.call(rbind, PerformanceRows), Diagnostics = do.call(rbind, DiagnosticRows))
}

DkFormatNumber <- function(Value, Digits = 3) {
  if (!is.finite(Value)) return("---")
  if (abs(Value) < 0.5*10^(-Digits)) Value <- 0
  formatC(Value, format = "f", digits = Digits)
}

DkLatexNumber <- function(Value, Digits = 3, Before = "", After = "") {
  if (!is.finite(Value)) return("{---}")
  Text <- DkFormatNumber(Value, Digits)
  if (nzchar(Before)) Text <- paste0("{", Before, "}", Text)
  if (nzchar(After)) Text <- paste0(Text, "{", After, "}")
  Text
}

DkLatexNumberFormat <- function(Values, Digits = 3) {
  Values <- Values[is.finite(Values)]
  if (length(Values) == 0) return(paste0("1.", Digits))
  Text <- vapply(Values, DkFormatNumber, character(1), Digits = Digits)
  IntegerText <- sub("\\..*$", "", sub("^-", "", Text))
  paste0(if (any(startsWith(Text, "-"))) "-" else "", max(nchar(IntegerText)), ".", Digits)
}

DkWritePerformanceLatex <- function(Reports, Study, Settings, Directory) {
  Lines <- c("% Requires booktabs, amsmath, pdflscape, and siunitx.",
             "% Sample sizes are shown side by side, with one landscape table per marginal setting.",
             "% Coverage and its parenthesized MCSE occupy separate decimal-aligned columns.",
             "% Simulation Monte Carlo errors for all measures are in ParameterPerformance.csv.")
  J <- length(Settings$VisitMonths)
  for (SettingId in unique(Reports$Performance$Setting)) {
    Rows <- Reports$Performance[Reports$Performance$Setting == SettingId, , drop = FALSE]
    Diagnostics <- Reports$Diagnostics[Reports$Diagnostics$Setting == SettingId, , drop = FALSE]
    Sizes <- sort(unique(Rows$N))
    Parameters <- Study$Metadata$Parameter[Study$Metadata$Parameter %in% Rows$Parameter]
    SettingLabel <- if (SettingId == "Application") "application-like margins" else "greater concentration at zero"
    Prefix <- if (Settings$Stage == "Pilot") "Pilot simulation" else "Simulation"
    GroupHeaders <- vapply(Sizes, function(N) {
      paste0("\\multicolumn{5}{c}{$\\boldsymbol{n = ", N, "}$}")
    }, character(1))
    Rules <- vapply(seq_along(Sizes), function(Index) {
      First <- 3 + 5*(Index - 1)
      paste0("\\cmidrule(lr){", First, "-", First + 4, "}")
    }, character(1))
    Columns <- c("l", paste0("S[table-format = ", DkLatexNumberFormat(Rows$Truth), "]"))
    for (N in Sizes) {
      Current <- Rows[Rows$N == N, , drop = FALSE]
      Columns <- c(Columns,
                   paste0("S[table-format = ", DkLatexNumberFormat(Current$Bias), "]"),
                   paste0("S[table-format = ", DkLatexNumberFormat(Current$Rmse), "]"),
                   paste0("S[table-format = ", DkLatexNumberFormat(Current$CoveragePercent, 1), "]"),
                   "@{\\extracolsep{0pt}\\;}",
                   paste0("S[table-format = ", DkLatexNumberFormat(Current$CoverageMcsePercent, 1),
                          ", table-space-text-pre = {(}, table-space-text-post = {)}]"),
                   "@{\\extracolsep{\\fill}\\hspace{2\\tabcolsep}}",
                   paste0("S[table-format = ", DkLatexNumberFormat(Current$MeanIntervalWidth), "]"))
    }
    Headings <- c("\\textbf{Parameter}", "{\\textbf{Truth}}",
                  rep(c("{\\textbf{Bias}}", "{\\textbf{RMSE}}",
                        "\\multicolumn{2}{c}{\\shortstack{\\textbf{Coverage (\\%)}\\\\\\textbf{(MCSE)}}}",
                        "{\\textbf{Mean width}}"), length(Sizes)))
    Lines <- c(Lines, "", "\\begin{landscape}", "\\begin{table}[p]", "\\centering\\small",
               paste0("\\caption{", Prefix, ": complete-data parameter recovery with ", SettingLabel, ".}"),
               paste0("\\label{tab:dk-recovery-", tolower(Settings$Stage), "-", tolower(SettingId), "}"),
               "\\setlength{\\tabcolsep}{5pt}", "\\renewcommand{\\arraystretch}{1.12}",
               "\\sisetup{table-align-text-pre = false, table-align-text-post = false, group-digits = false, output-decimal-marker = {.}, round-mode = none}",
               paste0("\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}", paste(Columns, collapse = " "), "@{}}"),
               "\\toprule", paste0(" & & ", paste(GroupHeaders, collapse = " & "), " \\\\"),
               paste(Rules, collapse = ""), paste0(paste(Headings, collapse = " & "), " \\\\"), "\\midrule")
    for (Index in seq_along(Parameters)) {
      ParameterRows <- Rows[Rows$Parameter == Parameters[Index], , drop = FALSE]
      if (nrow(ParameterRows) != length(Sizes) || anyDuplicated(ParameterRows$N) ||
          length(unique(ParameterRows$Truth)) != 1) stop("Cannot combine this parameter across sample sizes.")
      Cells <- c(paste0("$", ParameterRows$Latex[1], "$"),
                 DkLatexNumber(ParameterRows$Truth[1]))
      for (N in Sizes) {
        Row <- ParameterRows[ParameterRows$N == N, , drop = FALSE]
        Coverage <- if (is.finite(Row$CoveragePercent)) {
          c(DkLatexNumber(Row$CoveragePercent, 1),
            DkLatexNumber(Row$CoverageMcsePercent, 1, "(", ")"))
        } else c("{---}", "{---}")
        Cells <- c(Cells, DkLatexNumber(Row$Bias), DkLatexNumber(Row$Rmse), Coverage,
                   DkLatexNumber(Row$MeanIntervalWidth))
      }
      Lines <- c(Lines, paste0(paste(Cells, collapse = " & "), " \\\\"))
      if (Index %% J == 0 && Index < length(Parameters)) Lines <- c(Lines, "\\addlinespace[3pt]")
    }
    Counts <- vapply(Sizes, function(N) {
      D <- Diagnostics[Diagnostics$N == N, , drop = FALSE]
      if (nrow(D) != 1) stop("Expected one diagnostic summary per sample size and setting.")
      paste0("$n = ", N, "$: ", D$Completed, "/", D$Planned, " completed fits, ",
             D$Failed, " failed, ", D$Pending, " pending, and ", D$Flagged, " flagged")
    }, character(1))
    Note <- paste0(paste(Counts, collapse = "; "), ". Flagged fits remain included. ",
                   "Bias and root mean squared error (RMSE) refer to posterior means. Coverage and mean width refer to 95\\% equal-tailed credible intervals. ",
                   "Parentheses give the simulation Monte Carlo standard error (MCSE) of coverage in percentage points. ",
                   "MCSEs for the other measures and Wilson intervals for coverage are in the accompanying CSV. ",
                   "Each fit uses all ", format(Settings$NChains*Settings$Sample, big.mark = ",", scientific = FALSE, trim = TRUE),
                   " post-warmup draws. Shape and correlation are summarized on their natural scales. ",
                   "Visits $j = 1,\\ldots,", J, "$ correspond to months ", paste(Settings$VisitMonths, collapse = ", "), ".")
    if (Settings$Stage == "Pilot") Note <- paste(Note, "Pilot results assess the simulation and sampler; they do not establish repeated-sampling performance.")
    Lines <- c(Lines, "\\bottomrule", "\\end{tabular*}", "\\par\\smallskip",
               "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright", Note,
               "\\end{minipage}", "\\end{table}", "\\end{landscape}")
  }
  DkWriteLines(Lines, file.path(Directory, "ParameterPerformance.tex"))
  CountColumns <- c("N", "Planned", "Completed", "Failed", "Flagged")
  Columns <- vapply(CountColumns, function(Name) {
    paste0("S[table-format = ", DkLatexNumberFormat(Reports$Diagnostics[[Name]], 0), "]")
  }, character(1))
  Lines <- c("% Requires booktabs, amsmath, and siunitx.", "\\begin{table}[tbp]", "\\centering\\small",
             paste0("\\caption{", if (Settings$Stage == "Pilot") "Pilot simulation" else "Simulation", " completion and sampling diagnostics.}"),
             paste0("\\label{tab:dk-simulation-diagnostics-", tolower(Settings$Stage), "}"),
             "\\sisetup{table-align-text-pre = false, table-align-text-post = false, group-digits = false, output-decimal-marker = {.}, round-mode = none}",
             paste0("\\begin{tabular}{@{}l ", paste(Columns, collapse = " "), "@{}}"), "\\toprule",
             "\\textbf{Margins} & {$\\boldsymbol{n}$} & {\\textbf{Planned}} & {\\textbf{Completed}} & {\\textbf{Failed}} & {\\textbf{Flagged}} \\\\",
             "\\midrule")
  for (Index in seq_len(nrow(Reports$Diagnostics))) {
    D <- Reports$Diagnostics[Index, ]
    Lines <- c(Lines, paste0(if (D$Setting == "Application") "Application-like" else "Stronger floor", " & ",
                             D$N, " & ", D$Planned, " & ", D$Completed, " & ", D$Failed, " & ", D$Flagged, " \\\\"))
  }
  Lines <- c(Lines, "\\bottomrule", "\\end{tabular}", "\\par\\smallskip",
             "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright",
             paste0("A completed fit is flagged if any parameter has rank-normalized split $\\widehat R > ",
                    Settings$RhatThreshold, "$; bulk or tail effective sample size below ", Settings$EssThreshold,
                    "; or a nonfinite diagnostic. Flagged fits remain in all performance summaries. ",
                    "Failed fits are reported separately and are not replaced by newly simulated datasets. ",
                    "Planned minus completed minus failed gives the number pending."),
             "\\end{minipage}", "\\end{table}")
  DkWriteLines(Lines, file.path(Directory, "SimulationDiagnostics.tex"))
}

DkExportSimulationReports <- function(Tasks, Study, Settings, Directory, Results = NULL) {
  dir.create(Directory, recursive = TRUE, showWarnings = FALSE)
  Reports <- DkCollectReports(Tasks, Study, Settings, Results)
  DkWriteCsv(Reports$Performance, file.path(Directory, "ParameterPerformance.csv"))
  DkWriteCsv(Reports$Diagnostics, file.path(Directory, "SimulationDiagnostics.csv"))
  DkWriteCsv(Reports$Summaries, file.path(Directory, "ReplicateParameterSummaries.csv"))
  DkWriteCsv(Reports$Status, file.path(Directory, "ReplicateStatus.csv"))
  DkWriteCsv(Reports$Acceptance, file.path(Directory, "SamplerAcceptance.csv"))
  DkWriteCsv(Reports$Warnings, file.path(Directory, "Warnings.csv"))
  DkWriteCsv(Reports$Status[Reports$Status$Status == "Failed", , drop = FALSE], file.path(Directory, "FailureLog.csv"))
  DkWritePerformanceLatex(Reports, Study, Settings, Directory)
  DkSaveRds(Settings, file.path(Directory, "RunSettings.rds"))
  Reports
}

# ---- Reproducibility, dynamic cores, and the single entry point -------------

DkStatisticalFunctionNames <- function() {
  c("Softplus", "LogOneMinusExp", "LogAddExp", "LogDiffExp", "LogNegativeLogOneMinusExp",
    "DkLogShapes", "DkLogSurvival", "DkCdf", "DkInterval", "DkPmf", "DkQuantileFromLogSurvival",
    "DkQuantile", "PitToGaussian", "MakeContinuousAr1", "BuildCopulaGeometry", "CopulaLogBySubject",
    "ComputeDkPrior", "EvaluateDkState", "UpdateDkPit", "RegularizeProposalCovariance",
    "InitialDkParameters", "BuildDkBlocks", "UpdateBlockMoments", "AdaptDkBlock",
    "DkWithSeed", "DkPrepareStudy", "DkMakeTasks", "DkParameterMetadata",
    "DkMakeCompleteModel", "DkGenerateOutcomes", "RunDkSimulationChain", "DkSummariseSimulationFit")
}

DkWorkerFunctionNames <- function() {
  unique(c(DkStatisticalFunctionNames(), "DkHashObject", "DkCommitFile", "DkSaveRds", "DkReadRds",
           "DkWriteLines", "DkReadTaskResult", "DkApplyDiagnosticFlags", "RunDkSimulationReplicate"))
}

DkFunctionFingerprint <- function(Specification) {
  # Hash UTF-8 text in binary mode so this small compatibility guard is
  # independent of R's serialization header and Windows line endings.
  Path <- tempfile("DkFunction")
  on.exit(unlink(Path), add = TRUE)
  Text <- paste(c(Specification$Formals, "", Specification$Body), collapse = "\n")
  writeBin(charToRaw(enc2utf8(Text)), Path)
  unname(tools::md5sum(Path))
}

DkCheckpointStateSpecification <- function(Specification) {
  # The only sampler edit in 1.1.1 caches the unchanged PIT Jacobian during
  # parameter proposals. Same-seed comparisons against 1.1.0 give identical
  # draws and acceptance decisions. Preserve the previous checkpoint identity
  # ONLY for this exact audited function. Any further edit changes its hash
  # and therefore creates a new study as usual; no settings are disregarded.
  if (!identical(DkFunctionFingerprint(Specification), "0b06d94e97f850be3f4793c19fc310ff")) return(Specification)
  StateBody <- parse(text = Specification$Body, keep.source = FALSE)[[1]]
  CachedStatement <- parse(text = c("if (is.null(Current)) {",
                                    "State$LogJacobianBySubject <- rowSums(-Softplus(-PitEta) - Softplus(PitEta), na.rm = TRUE)",
                                    "}"), keep.source = FALSE)[[1]]
  Position <- which(vapply(as.list(StateBody), identical, logical(1), y = CachedStatement))
  if (length(Position) != 1) stop("Cannot verify the state evaluator's checkpoint compatibility.")
  StateBody[[Position]] <- CachedStatement[[3]][[2]]
  Specification$Body <- deparse(StateBody)
  if (!identical(DkFunctionFingerprint(Specification), "420b68091282ed8b412767e2c1d9684e")) {
    stop("The restored checkpoint specification is not the verified 1.1.0 evaluator.")
  }
  Specification
}

DkStudySpecification <- function(Study, Settings) {
  ExecutionSettings <- c("Stage", "OutputRoot", "PilotReplicates", "FullReplicates", "PilotScenarios",
                         "UseParallel", "CoreFraction", "RetryFailures", "InstallMissingPackages",
                         "PrintProgress", "ProgressEvery", "RhatThreshold", "EssThreshold")
  Functions <- setNames(lapply(DkStatisticalFunctionNames(), function(Name) {
    Function <- get(Name, envir = environment(DkStudySpecification), inherits = TRUE)
    list(Formals = deparse(formals(Function)), Body = deparse(body(Function)))
  }), DkStatisticalFunctionNames())
  Functions$EvaluateDkState <- DkCheckpointStateSpecification(Functions$EvaluateDkState)
  list(ProgramVersion = "1.1.0", Settings = Settings[setdiff(names(Settings), ExecutionSettings)],
       Designs = Study$Designs, BaselinePool = Study$BaselinePool,
       BaselineCenter = Study$BaselineCenter, Scenarios = Study$Scenarios,
       RVersion = R.version.string, Platform = R.version$platform,
       PackageVersions = vapply(c("HSAUR2", "posterior"), function(Name) as.character(utils::packageVersion(Name)), character(1)),
       StatisticalFunctions = Functions)
}

DkChooseWorkers <- function(Settings, Pending) {
  Available <- as.numeric(parallelly::availableCores())
  Limit <- max(1, floor(Settings$CoreFraction*Available))
  Connections <- as.numeric(parallelly::availableCores(constraints = "connections-16"))
  Workers <- if (Pending == 0) 0 else if (Settings$UseParallel) min(Pending, Limit, Connections) else 1
  list(Available = Available, Limit = Limit, Workers = Workers)
}

DkRunTaskQueue <- function(Cluster, Tasks, Pending, OnResult,
                           OnHeartbeat = NULL, HeartbeatSeconds = 30) {
  if (length(Pending) == 0) return(invisible(0))
  if (length(Cluster) == 0) stop("The task queue needs at least one worker.")
  # These are the same dispatch/receive primitives used by parallel's own
  # load-balanced scheduler. Socket polling allows master-side reports and
  # heartbeats without a barrier between groups of datasets, also on Windows.
  SendCall <- getFromNamespace("sendCall", "parallel")
  ReceiveResult <- getFromNamespace("recvOneResult", "parallel")
  WorkerFunction <- function(Task) RunDkSimulationReplicate(Task, DkWorkerContext)
  environment(WorkerFunction) <- .GlobalEnv
  Active <- rep(NA_integer_, length(Cluster))
  Next <- 1
  Finished <- 0
  Submit <- function(Worker) {
    if (Next > length(Pending)) return(invisible(NULL))
    TaskIndex <- Pending[Next]
    SendCall(Cluster[[Worker]], WorkerFunction, list(Tasks[[TaskIndex]]), tag = TaskIndex)
    Active[Worker] <<- TaskIndex
    Next <<- Next + 1
    invisible(NULL)
  }
  for (Worker in seq_len(min(length(Cluster), length(Pending)))) Submit(Worker)
  LastHeartbeat <- proc.time()[[3]]
  while (any(!is.na(Active))) {
    Workers <- which(!is.na(Active))
    Ready <- socketSelect(lapply(Cluster[Workers], function(Node) Node$con), timeout = 1)
    if (any(Ready)) {
      Received <- tryCatch(ReceiveResult(Cluster), error = function(Error) {
        stop("A simulation worker disconnected: ", conditionMessage(Error),
             ". Saved results and completed chains can be resumed with the same script.", call. = FALSE)
      })
      Worker <- Received$node
      TaskIndex <- Received$tag
      if (length(Worker) != 1 || !Worker %in% Workers ||
          !identical(TaskIndex, Active[Worker])) stop("Unexpected task identity returned by a worker.")
      if (inherits(Received$value, "try-error")) {
        stop("Worker failed for ", Tasks[[TaskIndex]]$Scenario, ", replicate ",
             Tasks[[TaskIndex]]$Replicate, ": ", as.character(Received$value),
             ". Completed checkpoints are retained.", call. = FALSE)
      }
      Active[Worker] <- NA_integer_
      # Dispatch before exporting reports so this worker can continue fitting.
      Submit(Worker)
      OnResult(TaskIndex, Received$value)
      Finished <- Finished + 1
    }
    Now <- proc.time()[[3]]
    if (!is.null(OnHeartbeat) && Now - LastHeartbeat >= HeartbeatSeconds) {
      OnHeartbeat(sum(!is.na(Active)), max(0, length(Pending) - Next + 1))
      LastHeartbeat <- proc.time()[[3]]
    }
  }
  invisible(Finished)
}

DkWriteStudyDesign <- function(Study, Settings, Directory) {
  dir.create(Directory, recursive = TRUE, showWarnings = FALSE)
  for (Name in names(Study$Designs)) DkWriteCsv(Study$Designs[[Name]], file.path(Directory, paste0("FixedDesign", Name, ".csv")))
  TruthRows <- lapply(Study$Scenarios, function(Scenario) {
    cbind(data.frame(Scenario = Scenario$Name, Setting = Scenario$Setting, N = Scenario$N,
                     InterceptShift = Scenario$InterceptShift), Study$Metadata,
          Truth = as.numeric(Scenario$Truth))
  })
  DkWriteCsv(do.call(rbind, TruthRows), file.path(Directory, "ParameterTruths.csv"))
  Probabilities <- DkScenarioProbabilities(Study, Settings)
  DkWriteCsv(Probabilities, file.path(Directory, "ScenarioProbabilities.csv"))
  DkWriteCsv(data.frame(BaselineScore = Study$BaselinePool), file.path(Directory, "BaselinePool.csv"))
  Notes <- c(
    "SIMULATION STUDY 1: COMPLETE-DATA PARAMETER RECOVERY",
    "",
    "Aim: assess estimation of the original DK Gaussian copula parameters under a correctly specified model.",
    paste("Sample sizes:", paste(Settings$SampleSizes, collapse = ", "), "; equal treatment allocation."),
    paste("Visits (months):", paste(Settings$VisitMonths, collapse = ", "), "; upper score bound:", Settings$ScoreMaximum),
    paste("Baseline pool:", length(Study$BaselinePool), "BtheB participants with at least one recorded follow-up."),
    paste("Baseline center:", format(Study$BaselineCenter, digits = 16), "; scale:", Settings$BaselineScale),
    "Within each sample size, baseline scores are sampled once with replacement from that pool.",
    "A balanced treatment vector is randomized once. Both are then held fixed across all replicates.",
    "The same fixed design is used in the two marginal settings at each sample size.",
    "Performance therefore describes repeated outcome sampling conditional on these fixed designs.",
    "Application-like truths are fixed rounded application estimates, not draws from a posterior or prior.",
    paste("The stronger-floor setting adds", Settings$FloorInterceptShift, "to every intercept; all other truths are unchanged."),
    "ScenarioProbabilities.csv reports the resulting probabilities of zero, the maximum, and scores at most five.",
    "The endpoint probabilities average over the actual fixed baseline values in each assigned arm.",
    "",
    "Independent participant Gaussian vectors are generated with R[j,k] = rho^abs(t[j] - t[k]).",
    "The DK inverse CDF maps these vectors to complete integer outcomes. There are no missing values.",
    "The existing augmented likelihood and Metropolis-within-Gibbs sampler are fitted to each dataset.",
    "Initialization uses only the generated observations and specified priors; it does not use parameter truths.",
    "Priors are mutually independent across parameters and visits. Normal prior scales below are standard deviations.",
    paste("Intercept prior: Normal(0, SD =", Settings$InterceptPriorSd, ")."),
    paste("Treatment prior: Normal(0, SD =", Settings$TreatmentPriorSd, ")."),
    paste("Baseline coefficient prior: Normal(0, SD =", Settings$BaselinePriorSd, ")."),
    paste("Log shape prior: Normal(mean =", Settings$LogShapePriorMean, ", SD =", Settings$LogShapePriorSd, ")."),
    "Correlation prior: rho ~ Uniform(0, 1), where rho is the one-month latent correlation.",
    "The logit(rho) target adds log(rho) + log(1-rho) for the change of variables, evaluated stably.",
    "The initial correlation is 0.5 before chain-specific jitter; it is not a prior mean setting.",
    paste("Each replicate has", Settings$NChains, "chains; each has", Settings$Warmup, "warmup and", Settings$Sample, "saved draws."),
    "Adaptation stops at the end of warmup. All saved parameter draws are used, without thinning.",
    "The saved chain matrices use alpha, beta, gamma, log(phi), logit(rho) coordinates.",
    "Every draw of log(phi) and logit(rho) is transformed before computing natural-scale means, SDs, and intervals.",
    "",
    "Performance measures: bias and RMSE of posterior means; empirical SD; mean posterior SD;",
    "coverage and mean width of 95% equal-tailed credible intervals. Quantiles use R type 8.",
    "Simulation MCSEs: SD(errors)/sqrt(R); SD(squared errors)/(2*RMSE*sqrt(R));",
    "sqrt(coverage*(1-coverage)/R); SD(interval widths)/sqrt(R), where R is the contributing count.",
    "The RMSE MCSE is a delta-method approximation. Wilson intervals for coverage are also exported.",
    "When observed coverage is 0% or 100%, its plug-in MCSE is zero; use the Wilson interval to assess uncertainty.",
    "MCMC MCSEs for individual posterior means are distinct from simulation MCSEs across replicates.",
    "Flags never remove a completed fit. Actual failures and pending fits have explicit counts.",
    "A retry uses the identical simulated dataset and RNG streams, retaining any completed chains.",
    "If fits fail, performance measures are conditional on the fits that produced finite summaries.",
    "",
    "RNG: L'Ecuyer-CMRG. Scenarios use disjoint streams; each replicate reserves one substream",
    "for data generation and one for each chain. Scheduling and worker count do not change the draws.",
    "The pilot and full study share checkpoints when model, sampler, seeds, and software match.",
    "A new study fingerprint creates a separate folder automatically; older results remain available.",
    "Files named Chain1.rds, etc. contain all saved parameter draws and acceptance information.",
    "Latent augmentation values are not saved because this study does not perform imputation.",
    "",
    "The pilot is for checking generation, fitting, mixing, and runtime. Review its diagnostics before a full run.",
    paste("With", Settings$FullReplicates, "successful replicates, coverage near 95% has a simulation MCSE of",
          round(100*sqrt(0.95*0.05/Settings$FullReplicates), 2), "percentage points."),
    "Increasing the replication count extends an existing study. Changing MCMC lengths creates a new study folder.",
    "A failed or flagged pilot is a reason to investigate the setup; it is not evidence of parameter recovery.",
    "",
    "Software references:",
    "https://parallelly.futureverse.org/reference/availableCores.html",
    "https://mc-stan.org/posterior/",
    "https://github.com/r-lib/filelock"
  )
  DkWriteLines(Notes, file.path(Directory, "StudyDesign.txt"))
}

RunDkParameterSimulation <- function(Settings = DefaultDkSimulationSettings()) {
  DkCheckSimulationSettings(Settings)
  DkEnsurePackages(Settings)
  CheckDkDistribution()
  Study <- DkPrepareStudy(Settings)
  Specification <- DkStudySpecification(Study, Settings)
  Signature <- DkHashObject(Specification)
  StudyDirectory <- file.path(Settings$OutputRoot, paste0("Study", substr(Signature, 1, 12)))
  dir.create(StudyDirectory, recursive = TRUE, showWarnings = FALSE)
  StudyDirectory <- normalizePath(StudyDirectory, winslash = "/", mustWork = TRUE)
  Lock <- filelock::lock(file.path(StudyDirectory, "Run.lock"), timeout = 0)
  if (is.null(Lock)) stop("Another R process is running this study. Stop it or let it finish before resuming.")
  Cluster <- NULL
  WorkerPids <- integer(0)
  ThreadVariables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                       "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS")
  OldThreads <- Sys.getenv(ThreadVariables, unset = NA_character_)
  on.exit({
    if (!is.null(Cluster)) {
      # Stop local workers before releasing the study lock after an error or
      # interrupt, preventing writes from overlapping with a restarted run.
      if (length(WorkerPids) > 0) try(tools::pskill(WorkerPids), silent = TRUE)
      try(parallel::stopCluster(Cluster), silent = TRUE)
    }
    filelock::unlock(Lock)
    Missing <- is.na(OldThreads)
    if (any(Missing)) Sys.unsetenv(ThreadVariables[Missing])
    if (any(!Missing)) do.call(Sys.setenv, as.list(OldThreads[!Missing]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(setNames(rep("1", length(ThreadVariables)), ThreadVariables)))
  SpecificationFile <- file.path(StudyDirectory, "StudySpecification.rds")
  Existing <- DkReadRds(SpecificationFile)
  if (!is.null(Existing) && !identical(Existing, Specification)) stop("Study fingerprint collision; choose a different OutputRoot.")
  if (is.null(Existing)) DkSaveRds(Specification, SpecificationFile)
  DkWriteStudyDesign(Study, Settings, file.path(StudyDirectory, "Design"))
  Tasks <- DkMakeTasks(Study, Settings, StudyDirectory, Signature)
  ReportDirectory <- file.path(StudyDirectory, Settings$Stage)
  # Read checkpoints once, before starting workers. Updating this small result
  # snapshot avoids rereading every RDS file and never races a worker's writes.
  Results <- lapply(Tasks, DkReadTaskResult)
  Reports <- DkExportSimulationReports(Tasks, Study, Settings, ReportDirectory, Results)
  DkWriteLines(capture.output(utils::sessionInfo()), file.path(ReportDirectory, "SessionInfo.txt"))
  Pending <- which(vapply(Results, function(Result) {
    is.null(Result) || (!isTRUE(Result$Success) && Settings$RetryFailures)
  }, logical(1)))
  Cores <- DkChooseWorkers(Settings, length(Pending))
  RuntimeInformation <- list(Stage = Settings$Stage, AvailableCores = Cores$Available,
                             CoreFraction = Settings$CoreFraction, WorkerLimit = Cores$Limit,
                             Workers = Cores$Workers, PendingAtStart = length(Pending),
                             ImplementationVersion = "1.1.1",
                             StateEvaluator = list(Formals = deparse(formals(EvaluateDkState)),
                                                   Body = deparse(body(EvaluateDkState))),
                             StartedUtc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  DkSaveRds(RuntimeInformation, file.path(ReportDirectory, "RuntimeInformation.rds"))
  cat("\nDK parameter simulation:", Settings$Stage, "\n",
      "Available cores:", Cores$Available, ";", paste0(100*Settings$CoreFraction, "%"), "limit:", Cores$Limit,
      "; workers:", Cores$Workers, "\n",
      "Replicates scheduled:", length(Tasks), "; pending or retrying:", length(Pending), "\n",
      "Per replicate:", Settings$NChains, "chains x", Settings$Sample,
      "saved draws, plus", Settings$Warmup, "warmup iterations per chain.\n",
      "Report folder:", ReportDirectory, "\n\n")
  utils::flush.console()
  if (length(Pending) > 0) {
    Context <- list(Study = Study, Settings = Settings)
    Start <- proc.time()[[3]]
    SinceReport <- 0
    Status <- Reports$Status$Status
    SaveReports <- function() {
      Reports <<- DkExportSimulationReports(Tasks, Study, Settings, ReportDirectory, Results)
      SinceReport <<- 0
      cat("Progress:", sum(Reports$Diagnostics$Completed), "completed,",
          sum(Reports$Diagnostics$Failed), "failed,", sum(Reports$Diagnostics$Pending),
          "pending; elapsed", round((proc.time()[[3]] - Start)/60, 1), "minutes.\n")
      utils::flush.console()
    }
    OnResult <- function(TaskIndex, Result) {
      Task <- Tasks[[TaskIndex]]
      if (!is.list(Result) || !identical(Result$StudySignature, Task$StudySignature) ||
          !identical(Result$Scenario, Task$Scenario) || !identical(Result$Replicate, Task$Replicate) ||
          !is.logical(Result$Success) || length(Result$Success) != 1 || is.na(Result$Success)) {
        stop("The worker result does not match its scheduled dataset.")
      }
      Results[TaskIndex] <<- list(Result)
      Status[TaskIndex] <<- if (Result$Success) "Completed" else "Failed"
      SinceReport <<- SinceReport + 1
      if (SinceReport >= max(1, Cores$Workers)) SaveReports()
      invisible(NULL)
    }
    OnHeartbeat <- function(Running, Waiting) {
      if (Settings$PrintProgress) {
        cat("Status:", sum(Status == "Completed"), "completed,", sum(Status == "Failed"),
            "failed;", Running, "running,", Waiting, "waiting; elapsed",
            round((proc.time()[[3]] - Start)/60, 1), "minutes.\n")
        utils::flush.console()
      }
    }
    if (Settings$UseParallel) {
      if (Settings$PrintProgress) {
        cat("Starting", Cores$Workers, "workers...\n")
        utils::flush.console()
      }
      Cluster <- parallel::makePSOCKcluster(Cores$Workers, outfile = "", rscript_args = "--vanilla")
      Names <- DkWorkerFunctionNames()
      Functions <- mget(Names, envir = environment(RunDkParameterSimulation), inherits = TRUE)
      WorkerPids <- unlist(parallel::clusterCall(Cluster, function(Paths, Functions, Context) {
        .libPaths(Paths)
        for (Name in names(Functions)) {
          Function <- Functions[[Name]]
          environment(Function) <- .GlobalEnv
          assign(Name, Function, envir = .GlobalEnv)
        }
        assign("DkWorkerContext", Context, envir = .GlobalEnv)
        Sys.getpid()
      }, .libPaths(), Functions, Context), use.names = FALSE)
      if (Settings$PrintProgress) {
        cat("Workers ready. Cumulative reports refresh after every", Cores$Workers,
            "finished attempts; status updates every 30 seconds.\n")
        utils::flush.console()
      }
      DkRunTaskQueue(Cluster, Tasks, Pending, OnResult, OnHeartbeat)
    } else {
      for (TaskIndex in Pending) OnResult(TaskIndex, RunDkSimulationReplicate(Tasks[[TaskIndex]], Context))
    }
    if (SinceReport > 0) SaveReports()
  }
  if (!is.null(Cluster)) {
    parallel::stopCluster(Cluster)
    Cluster <- NULL
  }
  Reports <- DkExportSimulationReports(Tasks, Study, Settings, ReportDirectory, Results)
  RuntimeInformation$FinishedUtc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  DkSaveRds(RuntimeInformation, file.path(ReportDirectory, "RuntimeInformation.rds"))
  Failed <- sum(Reports$Diagnostics$Failed)
  Flagged <- sum(Reports$Diagnostics$Flagged)
  cat("\nFinished:", sum(Reports$Diagnostics$Completed), "completed fits;", Failed,
      "failed fits;", Flagged, "completed fits with diagnostic flags.\n",
      "LaTeX table:", file.path(ReportDirectory, "ParameterPerformance.tex"), "\n",
      "Diagnostics:", file.path(ReportDirectory, "SimulationDiagnostics.csv"), "\n")
  if (Settings$Stage == "Pilot") {
    cat("Review the pilot before changing Stage to 'Full'. Compatible pilot fits will be reused.\n")
  }
  if (Failed > 0 || Flagged > 0) {
    cat("Check ReplicateStatus.csv and FailureLog.csv. Flagged completed fits remain in the summaries.\n")
  }
  invisible(list(StudyDirectory = StudyDirectory, ReportDirectory = ReportDirectory,
                 Performance = Reports$Performance, Diagnostics = Reports$Diagnostics,
                 Complete = all(Reports$Status$Status != "Pending"), AllFitsUsable = Failed == 0))
}

# Source this file to run. Change the settings in DefaultDkSimulationSettings().
if (isTRUE(getOption("DkSimulationRun", TRUE))) {
  DkSimulationSettings <- DefaultDkSimulationSettings()
  DkSimulationResult <- RunDkParameterSimulation(DkSimulationSettings)
}