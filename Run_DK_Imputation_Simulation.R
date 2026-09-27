# Simulation study 2: DK reference-based imputation followed by BB analysis.
# Simulation version 1.0.0, 25 September 2026.
#
# Run from the simulation directory:
# source("Run_DK_Imputation_Simulation.R")
# This script is self-contained. Required packages are installed if missing.
# Default: full study, 100 replicates, n = 100, 50 imputations per rule.
# Fifty imputations use 13, 13, 12, and 12 draws from the four chains.
# Compatible pilot and full checkpoints are reused on subsequent runs.
#
# Parallelism is across replicates, using 60% of detected available cores.
# Each replicate fits one DK copula and 153 BB models under the default settings.
# All post-warmup DK draws contribute to posterior summaries and direct DK means.
# Fifty joint draws of parameters and retained latent values produce the
# completed datasets for each rule under the default settings.
#
# Storage includes compact results and selected imputations. Full glmmTMB
# objects are not saved. By default, full chain files are removed after
# summaries using all draws and selected joint draws have been saved successfully.
# Set KeepChainFiles = TRUE before running to retain the full chain files.
# Saved chain files support resuming interrupted fits from completed chains.
# BB results are checkpointed every ten fits and at the end of each rule.
# Sourcing a completed run refreshes reports from compatible saved results.
#
# SourceStudyDirectory is optional. An empty value regenerates Study 1's
# ApplicationN100 datasets from the same design and data seeds. A path to the
# original study folder instead loads Data.rds files and checks them against
# regenerated data. Study 1 outputs are accessed only for reading.
# Imputation uses a new posterior fitted to retained outcomes in each replicate;
# posterior fits from the complete data in Study 1 are not used for imputation.
#
# Main outputs in the printed Pilot/Full report folder:
#   ImputationPerformance.tex (active mean, all generating/imputation rules)
#   ImputationContrasts.tex (mean difference, matching rules)
#   ImputationPerformanceStandalone.tex (both tables, ready to compile)
#   TreatmentPerformance.csv (all visits, matching generating/analysis rules)
#   CrossRulePerformance.csv (all generating versus MI rule combinations)
#   PairedDifferences.csv (MI minus complete-data BB, paired within replicate)
#   ReplicateEstimates.csv, ReplicateStatus.csv, SamplingDiagnostics.csv
#   BbFitDiagnostics.csv, MissingnessSummary.csv, FailureLog.csv
#   PilotTiming.txt (runtime extrapolation using measured component times)
# Design/GeneratingTruth.csv gives the true standardized means and MCSEs.
# Simulation coverage MCSEs differ from imputation/MCMC error.
#
# Priors: independent alpha, beta, gamma ~ N(0, 5^2);
# log(phi) ~ N(0, 2^2); rho ~ Uniform(0, 1). Settings use standard deviations.
#
# To load functions without running: options(DkImputationRun = FALSE).
# Custom settings can be passed to RunDkImputationSimulation(Settings).

DefaultDkImputationSettings <- function() {
  list(
    Stage = "Full",
    OutputRoot = "DkImputationSimulation",
    SourceStudyDirectory = "",
    PilotReplicates = 5,
    FullReplicates = 100,
    N = 100,
    VisitMonths = c(2, 3, 5, 8),
    ScoreMaximum = 63,
    Tau = 0.5,
    BaselineScale = 10,
    InterceptTruth = c(-0.95, -1.11, -1.48, -1.74),
    TreatmentTruth = c(-0.38, -0.43, -0.28, -0.16),
    BaselineTruth = c(0.44, 0.60, 0.63, 0.52),
    ShapeTruth = c(2.32, 1.99, 2.03, 2.15),
    RhoTruth = 0.78,
    Rules = c(MAR = 0, Half = 0.5, J2R = 1),
    DropoutTargets = c(0, 0.25, 0.40, 0.45),
    DropoutBaselineCoefficient = log(1.25),
    DropoutPreviousCoefficient = log(1.5),
    DropoutPreviousCenter = 20,
    DropoutPreviousScale = 10,
    CalibrationSize = 200000,
    DesignSeed = 20260923,
    DataSeed = 20260924,
    SimulationSeed = 20260925,
    CalibrationSeed = 20261001,
    TruthSeed = 20261002,
    TruthBatchSize = 50000,
    TruthMinBatches = 40,
    TruthMaxBatches = 200,
    TruthMcseTolerance = 0.005,
    NChains = 4,
    Warmup = 5000,
    Sample = 10000,
    NImputations = 50,
    AdaptEvery = 100,
    PitSweeps = 2,
    InterceptPriorSd = 5,
    TreatmentPriorSd = 5,
    BaselinePriorSd = 5,
    LogShapePriorMean = 0,
    LogShapePriorSd = 2,
    RhoPrior = "Uniform",
    RhatThreshold = 1.01,
    EssThreshold = 400,
    MaxIterations = 2000,
    MaxEvaluations = 3000,
    MaxGradient = 0.01,
    QuadratureNodes = c(30, 60, 120, 240),
    QuadratureTolerance = 1e-6,
    ConfidenceLevel = 0.95,
    BbCheckpointEvery = 10,
    KeepChainFiles = FALSE,
    UseParallel = TRUE,
    CoreFraction = 0.60,
    RetryFailures = TRUE,
    InstallMissingPackages = TRUE,
    PrintProgress = TRUE,
    ProgressEvery = 1000
  )
}

MiCheckSettings <- function(Settings) {
  if (!Settings$Stage %in% c("Pilot", "Full")) stop("Stage must be 'Pilot' or 'Full'.")
  Integers <- c("PilotReplicates", "FullReplicates", "N", "ScoreMaximum", "NChains",
                "Warmup", "Sample", "NImputations", "AdaptEvery", "PitSweeps",
                "CalibrationSize", "TruthBatchSize", "TruthMinBatches", "TruthMaxBatches",
                "BbCheckpointEvery", "MaxIterations", "MaxEvaluations", "ProgressEvery")
  for (Name in Integers) {
    Value <- Settings[[Name]]
    if (length(Value) != 1 || !is.finite(Value) || Value < 1 || Value != floor(Value)) {
      stop(Name, " must be a positive integer.")
    }
  }
  if (Settings$N < 20 || Settings$N %% 2 != 0) stop("N must be even and at least 20.")
  if (Settings$NChains < 2 || Settings$Sample < 100 || Settings$ScoreMaximum < 2) {
    stop("Use at least two chains, 100 saved draws, and a score bound of two.")
  }
  if (Settings$NImputations/Settings$NChains < 2 ||
      Settings$NImputations > Settings$NChains*Settings$Sample) {
    stop("Use at least two imputations per chain and no more than the saved draw count.")
  }
  if (Settings$FullReplicates < Settings$PilotReplicates) stop("The full target is smaller than the pilot.")
  J <- length(Settings$VisitMonths)
  if (J < 2 || any(!is.finite(Settings$VisitMonths)) || any(diff(Settings$VisitMonths) <= 0)) {
    stop("Visit times must be finite and increasing.")
  }
  for (Name in c("InterceptTruth", "TreatmentTruth", "BaselineTruth", "ShapeTruth", "DropoutTargets")) {
    if (length(Settings[[Name]]) != J || any(!is.finite(Settings[[Name]]))) {
      stop(Name, " must have one finite value per visit.")
    }
  }
  if (any(Settings$ShapeTruth <= 0) || Settings$RhoTruth <= 0 || Settings$RhoTruth >= 1) {
    stop("Invalid generating shape or correlation.")
  }
  P <- Settings$DropoutTargets
  if (P[1] != 0 || any(diff(P) <= 0) || any(P < 0 | P >= 1)) {
    stop("Dropout targets must start at zero and increase to less than one.")
  }
  if (!identical(Settings$Rules, c(MAR = 0, Half = 0.5, J2R = 1))) {
    stop("This study uses MAR, Half, and J2R weights 0, 0.5, and 1.")
  }
  for (Name in c("BaselineScale", "DropoutPreviousScale", "TruthMcseTolerance", "InterceptPriorSd",
                 "TreatmentPriorSd", "BaselinePriorSd", "LogShapePriorSd", "MaxGradient",
                 "QuadratureTolerance", "EssThreshold")) {
    if (length(Settings[[Name]]) != 1 || !is.finite(Settings[[Name]]) || Settings[[Name]] <= 0) {
      stop(Name, " must be positive and finite.")
    }
  }
  for (Name in c("DropoutBaselineCoefficient", "DropoutPreviousCoefficient", "DropoutPreviousCenter",
                 "LogShapePriorMean", "RhatThreshold")) {
    if (length(Settings[[Name]]) != 1 || !is.finite(Settings[[Name]])) stop(Name, " must be finite.")
  }
  if (!identical(Settings$RhoPrior, "Uniform") || Settings$Tau != 0.5) stop("Use the median model and uniform rho prior.")
  if (Settings$TruthMinBatches < 5 || Settings$TruthMaxBatches < Settings$TruthMinBatches ||
      Settings$CalibrationSize %% Settings$N != 0 || Settings$TruthBatchSize %% Settings$N != 0) {
    stop("Truth batches need at least five batches and batch sizes divisible by N.")
  }
  for (Name in c("DesignSeed", "DataSeed", "SimulationSeed", "CalibrationSeed", "TruthSeed")) {
    Value <- Settings[[Name]]
    if (length(Value) != 1 || !is.finite(Value) || Value < 1 || Value >= .Machine$integer.max || Value != floor(Value)) {
      stop(Name, " must be a valid positive integer seed.")
    }
  }
  SeedNames <- c("DesignSeed", "DataSeed", "SimulationSeed", "CalibrationSeed", "TruthSeed")
  if (anyDuplicated(unlist(Settings[SeedNames]))) stop("Use distinct seeds for the five simulation stages.")
  if (length(Settings$CoreFraction) != 1 || !is.finite(Settings$CoreFraction) ||
      Settings$CoreFraction <= 0 || Settings$CoreFraction > 1) stop("CoreFraction must be in (0, 1].")
  if (Settings$ConfidenceLevel <= 0 || Settings$ConfidenceLevel >= 1) stop("Invalid confidence level.")
  Nodes <- Settings$QuadratureNodes
  if (length(Nodes) < 2 || any(!is.finite(Nodes)) || any(Nodes < 5 | Nodes != floor(Nodes)) || any(diff(Nodes) <= 0)) {
    stop("Quadrature nodes must be increasing integers, at least five.")
  }
  for (Name in c("KeepChainFiles", "UseParallel", "RetryFailures", "InstallMissingPackages", "PrintProgress")) {
    if (!is.logical(Settings[[Name]]) || length(Settings[[Name]]) != 1 || is.na(Settings[[Name]])) {
      stop(Name, " must be TRUE or FALSE.")
    }
  }
  for (Name in c("OutputRoot", "SourceStudyDirectory")) {
    if (!is.character(Settings[[Name]]) || length(Settings[[Name]]) != 1 || is.na(Settings[[Name]])) {
      stop(Name, " must be a folder name.")
    }
  }
  if (!nzchar(Settings$OutputRoot)) stop("OutputRoot must not be empty.")
  invisible(TRUE)
}

MiEnsurePackages <- function(Settings) {
  Required <- c("HSAUR2", "posterior", "parallelly", "filelock", "glmmTMB", "TMB", "Matrix")
  Missing <- Required[!vapply(Required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(Missing) > 0 && Settings$InstallMissingPackages) {
    utils::install.packages(Missing, repos = "https://cloud.r-project.org")
  }
  Missing <- Required[!vapply(Required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(Missing) > 0) stop("Install these packages: ", paste(Missing, collapse = ", "))
  invisible(Required)
}

# Distribution, copula sampler, and BB analysis routines.
# Statistical function definitions contribute to the study signature.

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

DkChooseWorkers <- function(Settings, Pending) {
  Available <- as.numeric(parallelly::availableCores())
  Limit <- max(1, floor(Settings$CoreFraction*Available))
  Connections <- as.numeric(parallelly::availableCores(constraints = "connections-16"))
  Workers <- if (Pending == 0) 0 else if (Settings$UseParallel) min(Pending, Limit, Connections) else 1
  list(Available = Available, Limit = Limit, Workers = Workers)
}

DkSafeMinimum <- function(Values) {
  Values <- Values[is.finite(Values)]
  if (length(Values) == 0) NA_real_ else min(Values)
}

DkSafeMaximum <- function(Values) {
  Values <- Values[is.finite(Values)]
  if (length(Values) == 0) NA_real_ else max(Values)
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

BbNormalQuadrature <- function(Nodes) {
  Jacobi <- matrix(0, Nodes, Nodes)
  Indices <- seq_len(Nodes - 1)
  Jacobi[cbind(Indices, Indices + 1)] <- sqrt(Indices)
  Jacobi[cbind(Indices + 1, Indices)] <- sqrt(Indices)
  Eigen <- eigen(Jacobi, symmetric = TRUE)
  Order <- order(Eigen$values)
  list(Nodes = Eigen$values[Order], Weights = Eigen$vectors[1, Order]^2)
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

MiRunChain <- function(Chain, Model, Settings, Initial, Signature,
                       Directory, Seed, TaskLabel) {
  CacheFile <- file.path(Directory, paste0("Chain", Chain, ".rds"))
  Cached <- DkReadRds(CacheFile)
  if (!is.null(Cached)) {
    if (!identical(Cached$Signature, Signature) || !identical(Cached$Seed, Seed)) {
      stop("Chain checkpoint does not match this dataset and RNG stream: ", CacheFile)
    }
    if (!identical(dim(Cached$Draws), c(as.integer(Settings$Sample), length(Model$ParameterNames))) ||
        any(!is.finite(Cached$Draws))) stop("Invalid draws in checkpoint: ", CacheFile)
    if (!identical(Cached$SelectedIterations, MiSelectedIterations(Settings, Chain)) ||
        !identical(dim(Cached$ZObserved), as.integer(c(length(MiSelectedIterations(Settings, Chain)), sum(Model$Observed))))) {
      stop("Invalid selected latent values in chain checkpoint.")
    }
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
  SelectedIterations <- MiSelectedIterations(Settings, Chain)
  ZObserved <- matrix(NA_real_, length(SelectedIterations), sum(Model$Observed))
  SelectedRows <- match(seq_len(Settings$Sample), SelectedIterations)
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
      if (!is.na(SelectedRows[Row])) ZObserved[SelectedRows[Row], ] <- State$Z[Model$Observed]
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
                 SelectedIterations = SelectedIterations, ZObserved = ZObserved,
                 Acceptance = Acceptance,
                 RuntimeSeconds = proc.time()[[3]] - StartTime)
  DkSaveRds(Result, CacheFile)
  if (Settings$PrintProgress) cat(TaskLabel, "chain", Chain, "completed in",
                                  round(Result$RuntimeSeconds, 1), "seconds.\n")
  Result
}

MiRunTaskQueue <- function(Cluster, Tasks, Pending, OnResult,
                           OnHeartbeat = NULL, HeartbeatSeconds = 30) {
  if (length(Pending) == 0) return(invisible(0))
  if (length(Cluster) == 0) stop("The task queue needs at least one worker.")
  # These are the same dispatch/receive primitives used by parallel's own
  # load-balanced scheduler. Socket polling allows master-side reports and
  # heartbeats without a barrier between groups of datasets, also on Windows.
  SendCall <- getFromNamespace("sendCall", "parallel")
  ReceiveResult <- getFromNamespace("recvOneResult", "parallel")
  WorkerFunction <- function(Task) MiRunReplicate(Task, MiWorkerContext)
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



# ---- Generating model, dropout calibration, and independent truth ----------

MiStudyDesign <- function(Settings) {
  DataEnvironment <- new.env(parent = baseenv())
  utils::data("BtheB", package = "HSAUR2", envir = DataEnvironment)
  Data <- DataEnvironment$BtheB
  Keep <- rowSums(!is.na(Data[, c("bdi.2m", "bdi.3m", "bdi.5m", "bdi.8m")])) > 0
  Pool <- Data$bdi.pre[Keep]
  Center <- mean(Pool)
  Design <- DkWithSeed(Settings$DesignSeed, {
    Baseline <- sample(Pool, Settings$N, replace = TRUE)
    Treatment <- sample(rep(c(0, 1), each = Settings$N/2))
    data.frame(Participant = seq_len(Settings$N), Treatment = Treatment,
               BaselineScore = Baseline, Baseline = (Baseline - Center)/Settings$BaselineScale)
  })
  if (qr(cbind(1, Design$Treatment, Design$Baseline))$rank < 3) stop("The generating design is rank deficient.")
  Truth <- c(Settings$InterceptTruth, Settings$TreatmentTruth, Settings$BaselineTruth,
             Settings$ShapeTruth, Settings$RhoTruth)
  names(Truth) <- DkParameterMetadata(Settings)$Parameter
  list(Design = Design, BaselinePool = Pool, BaselineCenter = Center,
       Metadata = DkParameterMetadata(Settings),
       Scenario = list(Name = paste0("ApplicationN", Settings$N), N = Settings$N,
                       Truth = Truth, Design = Design))
}

MiGenerateBase <- function(Baseline, Treatment, Settings, Z = NULL) {
  N <- length(Baseline)
  J <- length(Settings$VisitMonths)
  if (is.null(Z)) Z <- matrix(rnorm(N*J), N, J)%*%chol(MakeContinuousAr1(Settings$VisitMonths, Settings$RhoTruth))
  Y <- matrix(0, N, J, dimnames = list(NULL, paste0("Month", Settings$VisitMonths)))
  for (Visit in seq_len(J)) {
    Eta <- Settings$InterceptTruth[Visit] + Settings$TreatmentTruth[Visit]*Treatment +
      Settings$BaselineTruth[Visit]*Baseline
    Y[, Visit] <- DkQuantileFromLogSurvival(pnorm(Z[, Visit], lower.tail = FALSE, log.p = TRUE),
                                            Eta, log(Settings$ShapeTruth[Visit]), Settings$ScoreMaximum, Settings$Tau)
  }
  list(Y = Y, Z = Z)
}

MiDropoutPredictor <- function(Baseline, Previous, Settings) {
  Settings$DropoutBaselineCoefficient*Baseline +
    Settings$DropoutPreviousCoefficient*(Previous - Settings$DropoutPreviousCenter)/Settings$DropoutPreviousScale
}

MiCalibrateDropout <- function(Study, Settings) {
  DkWithSeed(Settings$CalibrationSeed, {
    Rows <- rep(seq_len(Settings$N), length.out = Settings$CalibrationSize)
    B <- Study$Design$Baseline[Rows]
    A <- Study$Design$Treatment[Rows]
    Y <- MiGenerateBase(B, A, Settings)$Y
    Survival <- rep(1, length(Rows))
    Intercepts <- rep(NA_real_, ncol(Y))
    Achieved <- rep(0, ncol(Y))
    for (Visit in seq.int(2, ncol(Y))) {
      Predictor <- MiDropoutPredictor(B, Y[, Visit - 1], Settings)
      Target <- Settings$DropoutTargets[Visit]
      Objective <- function(Value) mean(1 - Survival*plogis(-Value - Predictor)) - Target
      Intercepts[Visit] <- uniroot(Objective, c(-40, 40), tol = 1e-12)$root
      Survival <- Survival*plogis(-Intercepts[Visit] - Predictor)
      Achieved[Visit] <- mean(1 - Survival)
    }
    list(Intercepts = Intercepts, Summary = data.frame(Visit = seq_len(ncol(Y)),
                                                       Month = Settings$VisitMonths, TargetMissing = Settings$DropoutTargets,
                                                       CalibrationMissing = Achieved, HazardIntercept = Intercepts))
  })
}

MiGenerateDropout <- function(Y, Baseline, Intercepts, Settings) {
  N <- nrow(Y)
  J <- ncol(Y)
  Trigger <- rep(J + 1, N)
  Uniforms <- matrix(runif(N*(J - 1)), N, J - 1)
  for (Visit in seq.int(2, J)) {
    Probability <- plogis(Intercepts[Visit] + MiDropoutPredictor(Baseline, Y[, Visit - 1], Settings))
    Drop <- Trigger == J + 1 & Uniforms[, Visit - 1] < Probability
    Trigger[Drop] <- Visit
  }
  Trigger
}

MiCompleteTruths <- function(Base, Design, Trigger, Settings) {
  J <- ncol(Base$Y)
  Complete <- setNames(lapply(Settings$Rules, function(Omega) Base$Y), names(Settings$Rules))
  for (Rule in names(Settings$Rules)[-1]) {
    for (Visit in seq_len(J)) {
      Rows <- which(Design$Treatment == 1 & Trigger <= Visit)
      if (length(Rows) == 0) next
      Eta <- Settings$InterceptTruth[Visit] + (1 - Settings$Rules[[Rule]])*Settings$TreatmentTruth[Visit] +
        Settings$BaselineTruth[Visit]*Design$Baseline[Rows]
      Complete[[Rule]][Rows, Visit] <- DkQuantileFromLogSurvival(
        pnorm(Base$Z[Rows, Visit], lower.tail = FALSE, log.p = TRUE), Eta,
        log(Settings$ShapeTruth[Visit]), Settings$ScoreMaximum, Settings$Tau)
    }
  }
  Complete
}

MiDkMean <- function(Eta, LogShape, K, Tau = 0.5) {
  N <- max(length(Eta), length(LogShape))
  Eta <- rep_len(Eta, N)
  LogShape <- rep_len(LogShape, N)
  Result <- numeric(N)
  LogFractions <- log(-log(seq_len(K)/(K + 1)))
  for (Start in seq.int(1, N, by = 2048)) {
    Rows <- seq.int(Start, min(N, Start + 2047))
    Shapes <- DkLogShapes(Eta[Rows], LogShape[Rows], Tau)
    LogPower <- -exp(outer(Shapes$LogA, LogFractions, "+"))
    Result[Rows] <- rowSums(exp(-exp(LogNegativeLogOneMinusExp(LogPower) + Shapes$LogB)))
  }
  Result
}

MiExactMarMeans <- function(Study, Settings) {
  J <- length(Settings$VisitMonths)
  Means <- matrix(0, J, 2)
  for (Visit in seq_len(J)) {
    for (Arm in 0:1) {
      Eta <- Settings$InterceptTruth[Visit] + Settings$TreatmentTruth[Visit]*Arm +
        Settings$BaselineTruth[Visit]*Study$Design$Baseline
      Means[Visit, Arm + 1] <- mean(MiDkMean(Eta, log(Settings$ShapeTruth[Visit]), Settings$ScoreMaximum, Settings$Tau))
    }
  }
  Means
}

MiTruthBatch <- function(Study, Settings) {
  # Integrate D conditional on the whole MAR trajectory. Hazards use only
  # preceding outcomes, so each product is the probability of reaching a visit.
  # Exact MAR means act as a control variate; only the borrowing shift is sampled.
  B <- rep(Study$Design$Baseline, length.out = Settings$TruthBatchSize)
  Base <- MiGenerateBase(B, rep(1, length(B)), Settings)
  J <- ncol(Base$Y)
  Survival <- rep(1, length(B))
  Shift <- matrix(0, J, length(Settings$Rules), dimnames = list(NULL, names(Settings$Rules)))
  for (Visit in seq_len(J)) {
    if (Visit > 1) {
      Survival <- Survival*plogis(-Study$Dropout$Intercepts[Visit] - MiDropoutPredictor(B, Base$Y[, Visit - 1], Settings))
    }
    for (Rule in names(Settings$Rules)[-1]) {
      Eta <- Settings$InterceptTruth[Visit] + (1 - Settings$Rules[[Rule]])*Settings$TreatmentTruth[Visit] +
        Settings$BaselineTruth[Visit]*B
      Borrowed <- DkQuantileFromLogSurvival(pnorm(Base$Z[, Visit], lower.tail = FALSE, log.p = TRUE),
                                            Eta, log(Settings$ShapeTruth[Visit]), Settings$ScoreMaximum, Settings$Tau)
      Shift[Visit, Rule] <- mean((1 - Survival)*(Borrowed - Base$Y[, Visit]))
    }
  }
  Shift
}

MiPrepareTruth <- function(Study, Settings, Directory, Signature) {
  TruthNames <- c("TruthSeed", "TruthBatchSize")
  Key <- DkHashObject(list(Study = Signature, Settings = Settings[TruthNames],
                           Function = deparse(body(MiTruthBatch))))
  File <- file.path(Directory, paste0("Truth", substr(Key, 1, 12), ".rds"))
  Saved <- DkReadRds(File)
  Batches <- if (is.null(Saved)) list() else Saved$Batches
  if (!is.null(Saved) && !identical(Saved$Key, Key)) stop("Truth checkpoint does not match this study.")
  Seeds <- DkWithSeed(Settings$TruthSeed, {
    Seed <- get(".Random.seed", envir = .GlobalEnv)
    lapply(seq_len(Settings$TruthMaxBatches), function(Index) {
      Current <- Seed
      Seed <<- parallel::nextRNGStream(Seed)
      Current
    })
  })
  if (Settings$PrintProgress) cat("Calculating independent generating truth;", length(Batches), "saved batches.\n")
  repeat {
    Count <- length(Batches)
    if (Count >= Settings$TruthMinBatches) {
      Values <- do.call(cbind, lapply(Batches, as.vector))
      Mcse <- apply(Values, 1, sd)/sqrt(Count)
      if (max(Mcse) <= Settings$TruthMcseTolerance) break
    }
    if (Count >= Settings$TruthMaxBatches) stop("Truth precision is insufficient. Increase TruthMaxBatches and rerun.")
    Batches[[Count + 1]] <- DkWithSeed(Seeds[[Count + 1]], MiTruthBatch(Study, Settings))
    DkSaveRds(list(Key = Key, Batches = Batches), File)
    if (Settings$PrintProgress && (Count + 1) %% 5 == 0) {
      cat("Truth:", Count + 1, "batches of", Settings$TruthBatchSize, "trajectories completed.\n")
      utils::flush.console()
    }
  }
  Shift <- matrix(rowMeans(Values), ncol = length(Settings$Rules))
  Error <- matrix(Mcse, ncol = length(Settings$Rules))
  Exact <- MiExactMarMeans(Study, Settings)
  Rows <- list()
  for (Rule in seq_along(Settings$Rules)) {
    for (Visit in seq_len(nrow(Exact))) {
      Means <- c(Exact[Visit, 1], Exact[Visit, 2] + Shift[Visit, Rule])
      Rows[[length(Rows) + 1]] <- data.frame(
        Rule = names(Settings$Rules)[Rule], Visit = Visit, Month = Settings$VisitMonths[Visit],
        Summary = c("Reference", "Active", "Difference"),
        Truth = c(Means, Means[2] - Means[1]),
        TruthMCSE = c(0, Error[Visit, Rule], Error[Visit, Rule]),
        TruthBatches = Count, TruthTrajectories = Count*Settings$TruthBatchSize)
    }
  }
  Table <- do.call(rbind, Rows)
  DkWriteCsv(Table, file.path(Directory, "GeneratingTruth.csv"))
  Table
}

MiSelectedIterations <- function(Settings, Chain) {
  Count <- floor(Settings$NImputations/Settings$NChains) +
    as.numeric(Chain <= Settings$NImputations %% Settings$NChains)
  as.numeric(round(seq(1, Settings$Sample, length.out = Count)))
}

MiMakeTasks <- function(Settings, Directory, Signature) {
  R <- if (Settings$Stage == "Pilot") Settings$PilotReplicates else Settings$FullReplicates
  DataSeed <- DkWithSeed(Settings$DataSeed, get(".Random.seed", envir = .GlobalEnv))
  ReplicateSeed <- DkWithSeed(Settings$SimulationSeed, get(".Random.seed", envir = .GlobalEnv))
  Tasks <- vector("list", R)
  for (Replicate in seq_len(R)) {
    Current <- ReplicateSeed
    ChainSeeds <- vector("list", Settings$NChains)
    for (Chain in seq_len(Settings$NChains)) {
      Current <- parallel::nextRNGSubStream(Current)
      ChainSeeds[[Chain]] <- Current
    }
    Tasks[[Replicate]] <- list(Scenario = paste0("ApplicationN", Settings$N), Replicate = Replicate,
                               StudySignature = Signature, DataSeed = DataSeed, DropoutSeed = ReplicateSeed,
                               ChainSeeds = ChainSeeds, ImputationSeed = parallel::nextRNGSubStream(Current),
                               Directory = file.path(Directory, "Replicates", sprintf("Replicate%04d", Replicate)))
    # Study 1 reserved five successive substreams per replicate (four chains).
    for (Index in seq_len(5)) DataSeed <- parallel::nextRNGSubStream(DataSeed)
    ReplicateSeed <- parallel::nextRNGStream(ReplicateSeed)
  }
  Tasks
}

MiPrepareData <- function(Task, Study, Settings) {
  File <- file.path(Task$Directory, "Data.rds")
  Saved <- DkReadRds(File)
  if (!is.null(Saved)) {
    if (!identical(Saved$StudySignature, Task$StudySignature) || !identical(Saved$DataSeed, Task$DataSeed) ||
        !identical(Saved$DropoutSeed, Task$DropoutSeed)) stop("The saved data do not match this task.")
    return(Saved)
  }
  Base <- DkWithSeed(Task$DataSeed, MiGenerateBase(Study$Design$Baseline, Study$Design$Treatment, Settings))
  if (nzchar(Settings$SourceStudyDirectory)) {
    Source <- file.path(Settings$SourceStudyDirectory, "Replicates", Task$Scenario,
                        sprintf("Replicate%04d", Task$Replicate), "Data.rds")
    if (!file.exists(Source)) stop("Study 1 dataset not found: ", Source)
    Original <- readRDS(Source)
    if (!identical(Original$DataSeed, Task$DataSeed) || !identical(Original$Design, Study$Design) ||
        !isTRUE(all.equal(Original$Y, Base$Y, tolerance = 0))) {
      stop("Study 1 data differ from the configured design or seeds.")
    }
    Base$Y <- Original$Y
  }
  Trigger <- DkWithSeed(Task$DropoutSeed,
                        MiGenerateDropout(Base$Y, Study$Design$Baseline, Study$Dropout$Intercepts, Settings))
  Complete <- MiCompleteTruths(Base, Study$Design, Trigger, Settings)
  Observed <- outer(Trigger, seq_len(ncol(Base$Y)), ">")
  Y <- Base$Y
  Y[!Observed] <- NA_real_
  for (Value in Complete) {
    if (!identical(Value[Observed], Y[Observed])) stop("Generating truths have different retained scores.")
  }
  Saved <- list(StudySignature = Task$StudySignature, DataSeed = Task$DataSeed,
                DropoutSeed = Task$DropoutSeed, Trigger = Trigger, Y = Y, Complete = Complete,
                Design = Study$Design)
  DkSaveRds(Saved, File)
  Saved
}

MiMakeModel <- function(Data, Study, Settings) {
  Y <- Data$Y
  Observed <- !is.na(Y)
  J <- ncol(Y)
  Keys <- apply(Observed, 1, function(Row) paste(as.integer(Row), collapse = ""))
  Patterns <- lapply(unique(Keys), function(Key) {
    Subjects <- which(Keys == Key)
    list(Subjects = Subjects, Visits = which(Observed[Subjects[1], ]))
  })
  Index <- list(Intercept = seq_len(J), Treatment = J + seq_len(J), Baseline = 2*J + seq_len(J),
                LogShape = 3*J + seq_len(J), LogitRho = 4*J + 1)
  Names <- c(paste0("Intercept[", seq_len(J), "]"), paste0("Treatment[", seq_len(J), "]"),
             paste0("Baseline[", seq_len(J), "]"), paste0("LogShape[", seq_len(J), "]"), "LogitRho")
  list(Y = Y, Observed = Observed, Patterns = Patterns,
       ObservedCells = which(Observed, arr.ind = TRUE), MissingCells = which(!Observed, arr.ind = TRUE),
       Treatment = Study$Design$Treatment, Baseline = Study$Design$Baseline,
       NSubjects = nrow(Y), NVisits = J, Index = Index, ParameterNames = Names,
       Trigger = Data$Trigger, PostTrigger = !Observed, VisitMonths = Settings$VisitMonths,
       ScoreMaximum = Settings$ScoreMaximum, Tau = Settings$Tau)
}


# ---- Observed-data fitting and paired completed-data analyses --------------

MiTargetLabels <- function(Settings) {
  data.frame(Visit = rep(seq_along(Settings$VisitMonths), each = 3),
             Month = rep(Settings$VisitMonths, each = 3),
             Summary = rep(c("Reference", "Active", "Difference"), length(Settings$VisitMonths)))
}

MiDirectMeans <- function(Draws, Model) {
  Baselines <- sort(unique(Model$Baseline))
  Weights <- tabulate(match(Model$Baseline, Baselines), nbins = length(Baselines))/length(Model$Baseline)
  Result <- matrix(0, nrow(Draws), 3*Model$NVisits)
  Index <- Model$Index
  for (Visit in seq_len(Model$NVisits)) {
    for (Arm in 0:1) {
      Estimate <- numeric(nrow(Draws))
      # Deduplicated baseline values give exactly the same empirical average.
      for (Covariate in seq_along(Baselines)) {
        Eta <- Draws[, Index$Intercept[Visit]] + Arm*Draws[, Index$Treatment[Visit]] +
          Draws[, Index$Baseline[Visit]]*Baselines[Covariate]
        Estimate <- Estimate + Weights[Covariate]*MiDkMean(Eta, Draws[, Index$LogShape[Visit]], Model$ScoreMaximum, Model$Tau)
      }
      Result[, 3*(Visit - 1) + Arm + 1] <- Estimate
    }
    Result[, 3*Visit] <- Result[, 3*Visit - 1] - Result[, 3*Visit - 2]
  }
  Result
}

MiSummariseDirect <- function(Chains, Model, Settings) {
  Targets <- MiTargetLabels(Settings)
  Array <- array(NA_real_, c(Settings$Sample, Settings$NChains, nrow(Targets)),
                 dimnames = list(iteration = seq_len(Settings$Sample), chain = seq_len(Settings$NChains),
                                 variable = paste0(Targets$Summary, Targets$Visit)))
  for (Chain in seq_along(Chains)) Array[, Chain, ] <- MiDirectMeans(Chains[[Chain]]$Draws, Model)
  Combined <- matrix(Array, ncol = nrow(Targets))
  Intervals <- apply(Combined, 2, quantile, probs = c((1 - Settings$ConfidenceLevel)/2,
                                                      (1 + Settings$ConfidenceLevel)/2), type = 8, names = FALSE)
  Diagnostics <- posterior::summarise_draws(posterior::as_draws_array(Array),
                                            "rhat", "ess_bulk", "ess_tail", "mcse_mean")
  Diagnostics <- Diagnostics[match(dimnames(Array)$variable, Diagnostics$variable), ]
  Estimates <- data.frame(Rule = "MAR", Method = "ObservedDK", Targets,
                          Estimate = colMeans(Combined), Variance = apply(Combined, 2, var),
                          Lower = Intervals[1, ], Upper = Intervals[2, ], M = NA_real_,
                          IntervalType = "Credible")
  DiagnosticTable <- data.frame(Parameter = dimnames(Array)$variable,
                                Rhat = Diagnostics$rhat, EssBulk = Diagnostics$ess_bulk,
                                EssTail = Diagnostics$ess_tail, McseMean = Diagnostics$mcse_mean)
  list(Estimates = Estimates, Diagnostics = DiagnosticTable)
}

MiFitDk <- function(Task, Model, Study, Settings) {
  File <- file.path(Task$Directory, "DkFit.rds")
  Signature <- DkHashObject(list(Study = Task$StudySignature, Y = Model$Y, Design = Study$Design,
                                 Seeds = Task$ChainSeeds))
  Saved <- DkReadRds(File)
  if (!is.null(Saved)) {
    if (!identical(Saved$Signature, Signature)) stop("The DK fit checkpoint does not match the observed data.")
    return(Saved)
  }
  Started <- proc.time()[[3]]
  Initial <- InitialDkParameters(Model, Settings)
  Chains <- lapply(seq_len(Settings$NChains), function(Chain) {
    MiRunChain(Chain, Model, Settings, Initial, Signature, Task$Directory,
               Task$ChainSeeds[[Chain]], paste("Replicate", Task$Replicate))
  })
  PostStarted <- proc.time()[[3]]
  Summary <- DkSummariseSimulationFit(Chains, Model, Study$Scenario, Settings)
  Direct <- MiSummariseDirect(Chains, Model, Settings)
  Selected <- do.call(rbind, lapply(Chains, function(Chain) Chain$Draws[Chain$SelectedIterations, , drop = FALSE]))
  Latents <- do.call(rbind, lapply(Chains, function(Chain) Chain$ZObserved))
  DrawIndex <- do.call(rbind, lapply(Chains, function(Chain) {
    data.frame(Chain = Chain$Chain, PostWarmupIteration = Chain$SelectedIterations)
  }))
  DrawIndex$Imputation <- seq_len(nrow(DrawIndex))
  if (nrow(Selected) != Settings$NImputations || any(!is.finite(c(Selected, Latents)))) {
    stop("Invalid joint draws selected for imputation.")
  }
  Saved <- list(Signature = Signature, Summary = Summary, Direct = Direct,
                Selected = Selected, Latents = Latents, DrawIndex = DrawIndex,
                Acceptance = do.call(rbind, lapply(Chains, function(Chain) Chain$Acceptance)),
                SamplingSeconds = sum(vapply(Chains, function(Chain) Chain$RuntimeSeconds, numeric(1))),
                FitSeconds = sum(vapply(Chains, function(Chain) Chain$RuntimeSeconds, numeric(1))) +
                  proc.time()[[3]] - PostStarted, NPosteriorDraws = Settings$NChains*Settings$Sample)
  DkSaveRds(Saved, File)
  # The successful compact fit is sufficient for every remaining calculation.
  # Full raw chains are optional; all their draws were used in Summary/Direct.
  if (!Settings$KeepChainFiles) {
    Files <- file.path(Task$Directory, paste0("Chain", seq_len(Settings$NChains), ".rds"))
    unlink(Files)
  }
  Saved
}

MiImpute <- function(Task, Fit, Model, Study, Settings) {
  File <- file.path(Task$Directory, "Imputations.rds")
  Saved <- DkReadRds(File)
  if (!is.null(Saved)) {
    if (!identical(Saved$FitSignature, Fit$Signature) || !identical(Saved$Seed, Task$ImputationSeed)) {
      stop("The imputation checkpoint does not match the fit.")
    }
    return(Saved$Bundle)
  }
  Values <- DkWithSeed(Task$ImputationSeed, {
    Imputed <- setNames(lapply(Settings$Rules, function(Omega) {
      matrix(NA_integer_, nrow(Model$MissingCells), Settings$NImputations)
    }), names(Settings$Rules))
    for (Draw in seq_len(Settings$NImputations)) {
      Parameter <- Fit$Selected[Draw, ]
      Correlation <- MakeContinuousAr1(Model$VisitMonths, plogis(Parameter[Model$Index$LogitRho]))
      Latents <- SampleDkMissingGaussians(Model, Correlation, Fit$Latents[Draw, ])
      for (Rule in names(Settings$Rules)) {
        Imputed[[Rule]][, Draw] <- MapDkMissing(Parameter, Latents, Model, Settings$Rules[[Rule]])
      }
    }
    Imputed
  })
  Common <- Model$Treatment[Model$MissingCells[, 1]] == 0
  for (Rule in names(Values)) {
    if (!identical(Values[[Rule]][Common, , drop = FALSE], Values$MAR[Common, , drop = FALSE])) {
      stop("Reference imputations differ across rules.")
    }
  }
  # Verify the theoretical ordering conditional on each posterior beta draw.
  Active <- which(!Common)
  if (length(Active) > 0) {
    Visits <- Model$MissingCells[Active, 2]
    Sign <- t(sign(Fit$Selected[, Model$Index$Treatment[Visits], drop = FALSE]))
    if (any((Values$Half[Active, , drop = FALSE] - Values$MAR[Active, , drop = FALSE])*Sign > 0) ||
        any((Values$J2R[Active, , drop = FALSE] - Values$Half[Active, , drop = FALSE])*Sign > 0)) {
      stop("Paired imputations violate predictor ordering.")
    }
  }
  Bundle <- list(Data = data.frame(bdi.pre = Study$Design$BaselineScore,
                                   treatment = factor(ifelse(Study$Design$Treatment == 0, "TAU", "BtheB"),
                                                      levels = c("TAU", "BtheB"))),
                 ObservedScores = Model$Y, MissingCells = Model$MissingCells,
                 ImputedValues = Values, DrawIndex = Fit$DrawIndex,
                 VisitMonths = Settings$VisitMonths, ScoreMaximum = Settings$ScoreMaximum,
                 BaselineCenter = Study$BaselineCenter, BaselineScale = Settings$BaselineScale,
                 FitSignature = Fit$Signature, Rules = Settings$Rules)
  DkSaveRds(list(FitSignature = Fit$Signature, Seed = Task$ImputationSeed, Bundle = Bundle), File)
  Bundle
}

MiBbDiagnostic <- function(Fit, Method, Rule, Imputation) {
  Value <- function(Name, Default = NA_real_) if (is.null(Fit[[Name]])) Default else Fit[[Name]]
  data.frame(Method = Method, Rule = Rule, Imputation = Imputation, Success = Fit$Ok,
             Attempts = Value("Attempts"), RuntimeSeconds = Value("RuntimeSeconds"),
             Optimizer = Value("Optimizer", ""), Convergence = Value("Convergence"),
             PdHess = Value("PdHess", FALSE), MaxGradient = Value("MaxGradient"),
             QuadratureNodes = Value("Nodes"), QuadratureMeanError = Value("QuadratureMeanError"),
             QuadratureVarianceError = Value("QuadratureVarianceError"),
             RandomSd = Value("RandomSd"), Precision = Value("Precision"),
             Warnings = Value("Warnings", ""), Error = Value("Error", ""))
}

MiAnalyseBb <- function(Task, Data, Bundle, Settings) {
  Context <- BbPrepareContext(Bundle, Settings)
  Targets <- MiTargetLabels(Settings)
  File <- file.path(Task$Directory, "BbFits.rds")
  Saved <- DkReadRds(File)
  if (is.null(Saved)) Saved <- list(Signature = Task$StudySignature, Fits = list())
  if (!identical(Saved$Signature, Task$StudySignature)) stop("The BB checkpoint does not match this study.")
  SinceSave <- 0
  for (Method in c("CompleteBB", "ImputedBB")) {
    for (Rule in names(Settings$Rules)) {
      Current <- Context
      Count <- if (Method == "CompleteBB") 1 else Settings$NImputations
      if (Method == "CompleteBB") {
        Current$Bundle$ObservedScores <- Data$Complete[[Rule]]
        Current$Bundle$MissingCells <- matrix(numeric(), 0, 2)
        Current$Bundle$ImputedValues <- setNames(lapply(Settings$Rules, function(Omega) matrix(numeric(), 0, 1)), names(Settings$Rules))
      }
      Start <- NULL
      for (Imputation in seq_len(Count)) {
        Key <- paste(Method, Rule, Imputation, sep = ":")
        Fit <- Saved$Fits[[Key]]
        if (is.null(Fit) || (!isTRUE(Fit$Ok) && Settings$RetryFailures)) {
          Fit <- BbAnalyseOne(Imputation, Rule, Current, Settings, Start)
          Saved$Fits[[Key]] <- Fit
          SinceSave <- SinceSave + 1
        }
        if (isTRUE(Fit$Ok)) Start <- Fit$Start
        if (SinceSave >= Settings$BbCheckpointEvery || Imputation == Count) {
          DkSaveRds(Saved, File)
          SinceSave <- 0
          DkWriteLines(paste("Replicate", Task$Replicate, Method, Rule, Imputation, "of", Count),
                       file.path(Task$Directory, "Progress.txt"))
        }
      }
    }
  }
  Estimates <- list()
  Diagnostics <- list()
  ChainEstimates <- list()
  for (Method in c("CompleteBB", "ImputedBB")) {
    for (Rule in names(Settings$Rules)) {
      Count <- if (Method == "CompleteBB") 1 else Settings$NImputations
      Fits <- lapply(seq_len(Count), function(Index) Saved$Fits[[paste(Method, Rule, Index, sep = ":")]])
      Diagnostics[[length(Diagnostics) + 1]] <- do.call(rbind, lapply(seq_len(Count), function(Index) {
        MiBbDiagnostic(Fits[[Index]], Method, Rule, Index)
      }))
      if (!all(vapply(Fits, function(Fit) isTRUE(Fit$Ok), logical(1)))) next
      if (Method == "CompleteBB") {
        Fit <- Fits[[1]]
        Critical <- qnorm((1 + Settings$ConfidenceLevel)/2)
        Estimates[[length(Estimates) + 1]] <- data.frame(Rule = Rule, Method = Method, Targets,
                                                         Estimate = Fit$Estimate, Variance = Fit$Variance,
                                                         Lower = Fit$Estimate - Critical*sqrt(Fit$Variance), Upper = Fit$Estimate + Critical*sqrt(Fit$Variance),
                                                         M = 1, IntervalType = "Confidence")
      } else {
        Values <- do.call(rbind, lapply(Fits, function(Fit) Fit$Estimate))
        Variances <- do.call(rbind, lapply(Fits, function(Fit) Fit$Variance))
        Pooled <- BbPool(Values, Variances, Bundle$DrawIndex, Targets, Rule, Settings$ConfidenceLevel)
        Estimates[[length(Estimates) + 1]] <- data.frame(Rule = Rule, Method = Method, Targets,
                                                         Estimate = Pooled$Estimate, Variance = Pooled$TotalVariance,
                                                         Lower = Pooled$Lower, Upper = Pooled$Upper, M = Count, IntervalType = "Confidence")
        for (Chain in seq_len(Settings$NChains)) {
          Rows <- which(Bundle$DrawIndex$Chain == Chain)
          ChainEstimates[[length(ChainEstimates) + 1]] <- data.frame(Rule = Rule, Chain = Chain, Targets,
                                                                     Estimate = colMeans(Values[Rows, , drop = FALSE]))
        }
      }
    }
  }
  list(Estimates = if (length(Estimates)) do.call(rbind, Estimates) else NULL,
       Diagnostics = do.call(rbind, Diagnostics),
       ChainEstimates = if (length(ChainEstimates)) do.call(rbind, ChainEstimates) else NULL)
}

MiMissingness <- function(Data, Settings) {
  Rows <- list()
  for (Arm in c("All", "Reference", "Active")) {
    Selected <- if (Arm == "All") seq_len(nrow(Data$Y)) else which(Data$Design$Treatment == as.numeric(Arm == "Active"))
    for (Visit in seq_len(ncol(Data$Y))) {
      Rows[[length(Rows) + 1]] <- data.frame(Arm = Arm, Visit = Visit, Month = Settings$VisitMonths[Visit],
                                             ProportionMissing = mean(is.na(Data$Y[Selected, Visit])))
    }
  }
  do.call(rbind, Rows)
}

MiReadResult <- function(Task) {
  Result <- DkReadRds(file.path(Task$Directory, "Result.rds"))
  if (!is.null(Result) && (!identical(Result$StudySignature, Task$StudySignature) ||
                           !identical(Result$Replicate, Task$Replicate))) stop("The replicate checkpoint does not match its task.")
  Result
}

MiRunReplicate <- function(Task, Context) {
  Settings <- Context$Settings
  Study <- Context$Study
  Cached <- MiReadResult(Task)
  if (!is.null(Cached) && (isTRUE(Cached$Success) || !Settings$RetryFailures)) return(Cached)
  dir.create(Task$Directory, recursive = TRUE, showWarnings = FALSE)
  Started <- proc.time()[[3]]
  Warnings <- character()
  Result <- tryCatch(withCallingHandlers({
    Data <- MiPrepareData(Task, Study, Settings)
    Model <- MiMakeModel(Data, Study, Settings)
    Fit <- MiFitDk(Task, Model, Study, Settings)
    Bundle <- MiImpute(Task, Fit, Model, Study, Settings)
    Bb <- MiAnalyseBb(Task, Data, Bundle, Settings)
    Estimates <- rbind(Fit$Direct$Estimates, Bb$Estimates)
    Success <- all(Bb$Diagnostics$Success) && nrow(Estimates) == 7*3*Model$NVisits
    list(Success = Success, Estimates = Estimates, Parameters = Fit$Summary,
         TargetDiagnostics = Fit$Direct$Diagnostics, Acceptance = Fit$Acceptance,
         BbDiagnostics = Bb$Diagnostics, ChainEstimates = Bb$ChainEstimates,
         Missingness = MiMissingness(Data, Settings),
         SamplingSeconds = Fit$SamplingSeconds, DkFitSeconds = Fit$FitSeconds,
         BbSeconds = sum(Bb$Diagnostics$RuntimeSeconds),
         Error = if (Success) "" else "Some BB fits are unavailable; incomplete imputation sets were not pooled.")
  }, warning = function(Warning) {
    Warnings <<- unique(c(Warnings, conditionMessage(Warning)))
    invokeRestart("muffleWarning")
  }), error = function(Error) list(Success = FALSE, Error = conditionMessage(Error)))
  Result$Replicate <- Task$Replicate
  Result$Scenario <- Task$Scenario
  Result$StudySignature <- Task$StudySignature
  Result$Warnings <- Warnings
  Result$AttemptSeconds <- proc.time()[[3]] - Started
  Result$FinishedUtc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  DkSaveRds(Result, file.path(Task$Directory, "Result.rds"))
  if (Settings$PrintProgress) {
    cat("Replicate", Task$Replicate, if (Result$Success) "completed." else paste("needs review:", Result$Error), "\n")
    utils::flush.console()
  }
  Result
}


# ---- Cumulative performance, diagnostics, and decimal-aligned tables -------

MiMcse <- function(Values) {
  if (length(Values) < 2) NA_real_ else sd(Values)/sqrt(length(Values))
}

MiPerformance <- function(Values, Truth) {
  N <- nrow(Values)
  if (N == 0) return(NULL)
  Error <- Values$Estimate - Truth
  Covered <- Values$Lower <= Truth & Values$Upper >= Truth
  Squared <- Error^2
  Rmse <- sqrt(mean(Squared))
  EmpiricalVariance <- if (N > 1) var(Values$Estimate) else NA_real_
  P <- mean(Covered)
  Z <- qnorm(0.975)
  WilsonCenter <- (P + Z^2/(2*N))/(1 + Z^2/N)
  WilsonRadius <- Z*sqrt(P*(1 - P)/N + Z^2/(4*N^2))/(1 + Z^2/N)
  Width <- Values$Upper - Values$Lower
  data.frame(N = N, Truth = Truth, MeanEstimate = mean(Values$Estimate),
             Bias = mean(Error), BiasMCSE = MiMcse(Error), RMSE = Rmse,
             RmseMCSE = if (Rmse > 0) MiMcse(Squared)/(2*Rmse) else 0,
             EmpiricalSD = sqrt(EmpiricalVariance), MeanSE = mean(sqrt(Values$Variance)),
             RootMeanVariance = sqrt(mean(Values$Variance)), MeanVariance = mean(Values$Variance),
             EmpiricalVariance = EmpiricalVariance,
             VarianceRatio = if (is.finite(EmpiricalVariance) && EmpiricalVariance > 0) mean(Values$Variance)/EmpiricalVariance else NA_real_,
             Coverage = 100*P, CoverageMCSE = 100*sqrt(P*(1 - P)/N),
             CoverageWilsonLower = 100*(WilsonCenter - WilsonRadius),
             CoverageWilsonUpper = 100*(WilsonCenter + WilsonRadius),
             MeanWidth = mean(Width), WidthMCSE = MiMcse(Width))
}

MiBindResults <- function(Results, Field) {
  Values <- lapply(Results, function(Result) {
    if (is.null(Result) || is.null(Result[[Field]]) || nrow(Result[[Field]]) == 0) return(NULL)
    cbind(data.frame(Replicate = Result$Replicate), Result[[Field]])
  })
  Values <- Filter(Negate(is.null), Values)
  if (length(Values)) do.call(rbind, Values) else data.frame()
}

MiCollectReports <- function(Tasks, Results, Truth, Settings) {
  Status <- do.call(rbind, lapply(seq_along(Tasks), function(Index) {
    Task <- Tasks[[Index]]
    Result <- Results[[Index]]
    Value <- function(Name, Default = NA_real_) if (is.null(Result[[Name]])) Default else Result[[Name]]
    data.frame(Replicate = Task$Replicate,
               Status = if (is.null(Result)) "Pending" else if (isTRUE(Result$Success)) "Completed" else "Failed",
               SamplingSeconds = Value("SamplingSeconds"), DkFitSeconds = Value("DkFitSeconds"),
               BbSeconds = Value("BbSeconds"), AttemptSeconds = Value("AttemptSeconds"),
               Error = Value("Error", ""))
  }))
  Estimates <- MiBindResults(Results, "Estimates")
  Parameters <- MiBindResults(Results, "Parameters")
  if (nrow(Parameters)) Parameters <- DkApplyDiagnosticFlags(Parameters, Settings)
  TargetDiagnostics <- MiBindResults(Results, "TargetDiagnostics")
  Bb <- MiBindResults(Results, "BbDiagnostics")
  Missingness <- MiBindResults(Results, "Missingness")
  Chains <- MiBindResults(Results, "ChainEstimates")
  Performance <- list()
  Cross <- list()
  Paired <- list()
  if (nrow(Estimates)) {
    Combinations <- unique(Estimates[, c("Rule", "Method", "Visit", "Month", "Summary")])
    for (Index in seq_len(nrow(Combinations))) {
      Group <- Combinations[Index, ]
      Current <- Estimates[Estimates$Rule == Group$Rule & Estimates$Method == Group$Method &
                             Estimates$Visit == Group$Visit & Estimates$Summary == Group$Summary, ]
      Target <- Truth[Truth$Rule == Group$Rule & Truth$Visit == Group$Visit & Truth$Summary == Group$Summary, ]
      if (nrow(Target) != 1) stop("A generating target is missing or duplicated.")
      Performance[[length(Performance) + 1]] <- cbind(Group, TruthMCSE = Target$TruthMCSE,
                                                      MiPerformance(Current, Target$Truth))
      if (Group$Method == "ImputedBB") {
        for (Generating in names(Settings$Rules)) {
          Target <- Truth[Truth$Rule == Generating & Truth$Visit == Group$Visit & Truth$Summary == Group$Summary, ]
          Cross[[length(Cross) + 1]] <- cbind(data.frame(GeneratingRule = Generating), Group,
                                              MiPerformance(Current, Target$Truth))
        }
        Complete <- Estimates[Estimates$Rule == Group$Rule & Estimates$Method == "CompleteBB" &
                                Estimates$Visit == Group$Visit & Estimates$Summary == Group$Summary, ]
        Match <- match(Current$Replicate, Complete$Replicate)
        Keep <- which(!is.na(Match))
        Difference <- Current$Estimate[Keep] - Complete$Estimate[Match[Keep]]
        if (length(Difference)) {
          Paired[[length(Paired) + 1]] <- cbind(Group, data.frame(NPaired = length(Difference),
                                                                  MeanDifference = mean(Difference), MCSE = MiMcse(Difference),
                                                                  RootMeanSquareDifference = sqrt(mean(Difference^2))))
        }
      }
    }
  }
  Bind <- function(Values) if (length(Values)) do.call(rbind, Values) else data.frame()
  Warnings <- do.call(rbind, lapply(Results, function(Result) {
    if (is.null(Result) || length(Result$Warnings) == 0) return(NULL)
    data.frame(Replicate = Result$Replicate, Warning = Result$Warnings)
  }))
  if (is.null(Warnings)) Warnings <- data.frame(Replicate = numeric(), Warning = character())
  Diagnostics <- data.frame(Planned = length(Tasks), Completed = sum(Status$Status == "Completed"),
                            Failed = sum(Status$Status == "Failed"), Pending = sum(Status$Status == "Pending"),
                            MaxRhat = if (nrow(Parameters)) DkSafeMaximum(Parameters$Rhat) else NA_real_,
                            MinBulkESS = if (nrow(Parameters)) DkSafeMinimum(Parameters$EssBulk) else NA_real_,
                            MinTailESS = if (nrow(Parameters)) DkSafeMinimum(Parameters$EssTail) else NA_real_,
                            BbFits = nrow(Bb), BbFailed = if (nrow(Bb)) sum(!Bb$Success) else 0)
  MissingSummary <- if (nrow(Missingness)) aggregate(ProportionMissing ~ Arm + Visit + Month, Missingness, mean) else data.frame()
  list(Status = Status, Estimates = Estimates, Parameters = Parameters, TargetDiagnostics = TargetDiagnostics,
       Bb = Bb, Missingness = Missingness, MissingSummary = MissingSummary, Chains = Chains,
       Performance = Bind(Performance), Cross = Bind(Cross), Paired = Bind(Paired),
       Diagnostics = Diagnostics, Warnings = Warnings)
}

MiLatexNumber <- function(Value, Digits = 2) {
  if (!is.finite(Value)) return("{--}")
  Value <- round(Value, Digits)
  if (Value == 0) Value <- 0
  formatC(Value, format = "f", digits = Digits, decimal.mark = ".")
}

MiNumberFormat <- function(Values, Digits = 2) {
  Values <- Values[is.finite(Values)]
  if (!length(Values)) return(paste0("1.", Digits))
  Integers <- max(nchar(formatC(floor(abs(round(Values, Digits))), format = "f", digits = 0)))
  paste0(if (any(round(Values, Digits) < 0)) "-" else "", Integers, ".", Digits)
}

MiManuscriptRows <- function(Reports, Settings, Truth, Summary, CrossRules) {
  Metrics <- c("N", "Truth", "MeanEstimate", "Bias", "RMSE", "Coverage", "CoverageMCSE", "MeanWidth")
  Visit <- length(Settings$VisitMonths)
  Rows <- list()
  AddRow <- function(Generating, Rule, Method) {
    UseCross <- CrossRules && Method == "ImputedBB"
    Data <- if (UseCross) Reports$Cross else Reports$Performance
    Current <- data.frame()
    if (!is.null(Data) && nrow(Data)) {
      Keep <- Data$Rule == Rule & Data$Method == Method &
        Data$Visit == Visit & Data$Summary == Summary
      if (UseCross) Keep <- Keep & Data$GeneratingRule == Generating
      Current <- Data[Keep, , drop = FALSE]
    }
    if (nrow(Current) > 1) stop("Duplicate performance rows for a table entry.")
    Values <- if (nrow(Current)) Current[1, Metrics, drop = FALSE] else
      as.data.frame(as.list(setNames(rep(NA_real_, length(Metrics)), Metrics)))
    if (!is.null(Truth) && nrow(Truth)) {
      Target <- Truth[Truth$Rule == Generating & Truth$Visit == Visit &
                        Truth$Summary == Summary, , drop = FALSE]
      if (nrow(Target) != 1) stop("A table's generating truth is missing or duplicated.")
      Values$Truth <- Target$Truth
    }
    cbind(data.frame(GeneratingRule = Generating, Rule = Rule, Method = Method), Values)
  }
  for (Generating in names(Settings$Rules)) {
    Rows[[length(Rows) + 1]] <- AddRow(Generating, Generating, "CompleteBB")
    if (Generating == "MAR") {
      Rows[[length(Rows) + 1]] <- AddRow(Generating, "MAR", "ObservedDK")
    }
    for (Rule in if (CrossRules) names(Settings$Rules) else Generating) {
      Rows[[length(Rows) + 1]] <- AddRow(Generating, Rule, "ImputedBB")
    }
  }
  do.call(rbind, Rows)
}

MiManuscriptTable <- function(Data, Settings, Active) {
  Columns <- c("MeanEstimate", "Bias", "RMSE", "Coverage", "CoverageMCSE", "MeanWidth")
  Digits <- c(2, 2, 2, 1, 1, 2)
  Formats <- vapply(seq_along(Columns), function(Index) {
    MiNumberFormat(Data[[Columns[Index]]], Digits[Index])
  }, character(1))
  Specs <- paste0("S[table-format = ", Formats, "]")
  Specs[5] <- paste0("S[table-format = ", Formats[5],
                     ", table-space-text-pre = {(}, table-space-text-post = {)}]")
  RowEnd <- paste0(" ", intToUtf8(c(92, 92)))
  Count <- if (Active) 8 else 7
  Month <- format(tail(Settings$VisitMonths, 1), trim = TRUE, scientific = FALSE)
  Caption <- if (Active) {
    paste0("Standardized mean in the active arm at month~", Month, " in the imputation simulation")
  } else {
    paste0("Standardized mean difference at month~", Month,
           " under matching generating and imputation assumptions")
  }
  Caption <- paste0(Caption, if (Settings$Stage == "Pilot") " (pilot)" else "",
                    " ($n=", Settings$N, "$).")
  Header <- c("\\textbf{Analysis}",
              if (Active) "\\shortstack{\\textbf{Imputation}\\\\\\textbf{rule}}" else character(),
              "{\\textbf{Mean}}", "{\\textbf{Bias}}", "{\\textbf{RMSE}}",
              "\\multicolumn{2}{c}{\\shortstack{\\textbf{Coverage (\\%)}\\\\\\textbf{(MCSE)}}}",
              "{\\textbf{Width}}")
  Lines <- c("% Generated from saved results; requires booktabs and siunitx.",
             "\\begin{table}[tbp]", paste0("\\caption{", Caption, "}"),
             paste0("\\label{", if (Active) "tab:imputation-simulation" else "tab:imputation-contrasts", "}"),
             "\\centering\\small", "\\begingroup",
             paste0("\\setlength{\\tabcolsep}{", if (Active) 3 else 4, "pt}"),
             "\\renewcommand{\\arraystretch}{1.08}",
             "\\sisetup{group-digits = false, round-mode = none, table-align-text-pre = false, table-align-text-post = false}",
             paste0("\\begin{tabular}{@{}", if (Active) "ll" else "l",
                    paste(Specs[1:4], collapse = ""), "@{\\hspace{\\fontdimen2\\font}}",
                    Specs[5], Specs[6], "@{}}"),
             "\\toprule", paste0(paste(Header, collapse = " & "), RowEnd), "\\midrule")
  Previous <- ""
  for (Generating in names(Settings$Rules)) {
    Group <- Data[Data$GeneratingRule == Generating, , drop = FALSE]
    if (nzchar(Previous)) Lines <- c(Lines, "\\addlinespace[5pt]")
    RuleLabel <- if (Generating == "Half") "Halfway" else Generating
    Targets <- unique(Group$Truth[is.finite(Group$Truth)])
    if (length(Targets) > 1) stop("Generating truths disagree within a table block.")
    TargetText <- if (length(Targets)) paste0("$", MiLatexNumber(Targets), "$") else "--"
    Title <- paste0(if (Active) "Generating assumption: " else "", RuleLabel,
                    if (Active) "; true mean " else "; true difference ", TargetText)
    Lines <- c(Lines, paste0("\\multicolumn{", Count, "}{@{}l}{\\textbf{", Title, "}}", RowEnd),
               "\\addlinespace[2pt]")
    PreviousMethod <- ""
    for (Index in seq_len(nrow(Group))) {
      Current <- Group[Index, , drop = FALSE]
      Label <- switch(Current$Method, CompleteBB = "Complete data, BB",
                      ObservedDK = "Observed data, DK", ImputedBB = "DK imputation, BB")
      if (Active && Current$Method == "ImputedBB" && PreviousMethod == "ImputedBB") Label <- ""
      Rule <- if (Current$Method %in% c("CompleteBB", "ObservedDK")) {
        "---"
      } else if (Current$Rule == "Half") {
        "Halfway"
      } else {
        Current$Rule
      }
      Numbers <- vapply(seq_along(Columns), function(Column) {
        MiLatexNumber(Current[[Columns[Column]]], Digits[Column])
      }, character(1))
      if (is.finite(Current$Coverage) && Current$Coverage %in% c(0, 100)) {
        Numbers[5] <- "{--}"
      } else if (is.finite(Current$CoverageMCSE)) {
        Numbers[5] <- paste0("{(}", Numbers[5], "{)}")
      }
      Lines <- c(Lines, paste0(paste(c(Label, if (Active) Rule else character(), Numbers),
                                     collapse = " & "), RowEnd))
      PreviousMethod <- Current$Method
    }
    Previous <- Generating
  }
  Counts <- Data$N[is.finite(Data$N) & Data$N > 0]
  if (!length(Counts)) {
    CountText <- "Results will appear as replicates finish."
  } else if (length(Counts) == nrow(Data) && length(unique(Counts)) == 1) {
    CountText <- paste0("Based on ", Counts[1], " replicates, with ", Settings$NImputations,
                        if (Active) " imputations per rule." else " imputations under the generating rule.")
  } else {
    CountText <- paste0("Available estimates per row range from ", min(Counts), " to ", max(Counts),
                        ", with ", Settings$NImputations, " imputations per rule. Row counts are in ",
                        "TreatmentPerformance.csv and CrossRulePerformance.csv; unavailable entries are shown as dashes.")
  }
  Level <- MiLatexNumber(100*Settings$ConfidenceLevel, 0)
  Notes <- if (Active) {
    paste0(CountText, " Means are standardized to the same ", Settings$N, " fixed baseline values. ",
           "Complete-data analyses use all generated scores. Retained outcomes are shared across generating assumptions, ",
           "so a given imputation rule yields the same estimates but is evaluated against a different truth in each block. ",
           "Coverage and mean width refer to nominal ", Level, "\\% intervals: equal-tailed credible intervals for DK ",
           "and confidence intervals for BB, with Rubin pooling after imputation. ")
  } else {
    paste0("The contrast is active minus reference, standardized to the same fixed baseline values. ",
           CountText, " Coverage and mean width refer to nominal ", Level,
           "\\% intervals: equal-tailed credible intervals for DK and confidence intervals for BB, ",
           "with Rubin pooling after imputation. ")
  }
  Notes <- paste0(Notes, "Coverage Monte Carlo standard errors (MCSEs) are in percentage points. ",
                  if (any(Data$Coverage %in% c(0, 100))) {
                    "Boundary coverage MCSEs are shown as dashes; Wilson intervals in the accompanying CSV quantify uncertainty. "
                  } else "",
                  "RMSE, root mean squared error.")
  c(Lines, "\\bottomrule", "\\end{tabular}", "\\endgroup", "\\par\\smallskip",
    "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright", Notes,
    if (Settings$Stage == "Pilot") "Pilot performance estimates are preliminary." else character(),
    "\\end{minipage}", "\\end{table}")
}

MiWriteTable <- function(Reports, Settings, Directory, Truth = NULL) {
  Active <- MiManuscriptTable(MiManuscriptRows(Reports, Settings, Truth, "Active", TRUE), Settings, TRUE)
  Contrast <- MiManuscriptTable(MiManuscriptRows(Reports, Settings, Truth, "Difference", FALSE), Settings, FALSE)
  DkWriteLines(Active, file.path(Directory, "ImputationPerformance.tex"))
  DkWriteLines(Contrast, file.path(Directory, "ImputationContrasts.tex"))
  Standalone <- c("\\documentclass[11pt]{article}", "\\usepackage[margin=0.75in]{geometry}",
                  "\\usepackage[T1]{fontenc}", "\\usepackage{booktabs,siunitx}",
                  "\\setlength{\\belowcaptionskip}{6pt}",
                  "\\begin{document}", Active, "\\clearpage", Contrast, "\\end{document}")
  DkWriteLines(Standalone, file.path(Directory, "ImputationPerformanceStandalone.tex"))
  invisible(c(Active = file.path(Directory, "ImputationPerformance.tex"),
              Contrast = file.path(Directory, "ImputationContrasts.tex"),
              Standalone = file.path(Directory, "ImputationPerformanceStandalone.tex")))
}


MiWriteReports <- function(Tasks, Results, Study, Settings, Truth, Directory, Cores) {
  dir.create(Directory, recursive = TRUE, showWarnings = FALSE)
  Reports <- MiCollectReports(Tasks, Results, Truth, Settings)
  Files <- c(Status = "ReplicateStatus", Estimates = "ReplicateEstimates", Parameters = "SamplingDiagnostics",
             TargetDiagnostics = "DkMeanDiagnostics", Bb = "BbFitDiagnostics", Missingness = "ReplicateMissingness",
             MissingSummary = "MissingnessSummary", Chains = "MiEstimatesByChain", Performance = "TreatmentPerformance",
             Cross = "CrossRulePerformance", Paired = "PairedDifferences", Diagnostics = "SimulationDiagnostics",
             Warnings = "Warnings")
  for (Name in names(Files)) DkWriteCsv(Reports[[Name]], file.path(Directory, paste0(Files[[Name]], ".csv")))
  DkWriteCsv(Reports$Status[Reports$Status$Status == "Failed", c("Replicate", "Error"), drop = FALSE],
             file.path(Directory, "FailureLog.csv"))
  MiWriteTable(Reports, Settings, Directory, Truth)
  DkSaveRds(Settings, file.path(Directory, "RunSettings.rds"))
  Timings <- Reports$Status[Reports$Status$Status == "Completed", ]
  Time <- Timings$DkFitSeconds + Timings$BbSeconds
  Time <- Time[is.finite(Time) & Time > 0]
  FullWorkers <- if (Settings$UseParallel) min(Cores$Limit, Settings$FullReplicates) else 1
  TimingText <- c(paste("Completed pilot/full replicates available:", length(Time)),
                  paste("Available cores:", Cores$Available, "; full-run worker limit:", FullWorkers))
  if (length(Time)) {
    Remaining <- max(0, Settings$FullReplicates - nrow(Timings))
    Hours <- Remaining*mean(Time)/FullWorkers/3600
    TimingText <- c(TimingText, paste("Mean DK plus BB processing time per replicate (minutes):", round(mean(Time)/60, 2)),
                    paste("Approximate remaining full-run hours:", round(Hours, 2)),
                    "This projection assumes similar processing times and full worker utilization. It excludes setup and interruptions.")
  }
  DkWriteLines(TimingText, file.path(Directory, "PilotTiming.txt"))
  Reports
}


# ---- Study identity, execution, and restart handling ----------------------

MiStatisticalNames <- function() {
  c("DkWithSeed", "Softplus", "LogOneMinusExp", "LogAddExp", "LogDiffExp",
    "LogNegativeLogOneMinusExp", "DkLogShapes", "DkLogSurvival", "DkCdf", "DkInterval",
    "DkPmf", "DkQuantileFromLogSurvival", "DkQuantile", "PitToGaussian", "MakeContinuousAr1",
    "BuildCopulaGeometry", "CopulaLogBySubject", "ComputeDkPrior", "EvaluateDkState",
    "UpdateDkPit", "RegularizeProposalCovariance", "InitialDkParameters", "BuildDkBlocks",
    "UpdateBlockMoments", "AdaptDkBlock", "MiRunChain", "DkSummariseSimulationFit",
    "SampleDkMissingGaussians", "MapDkMissing", "BbNormalQuadrature", "BbPrepareContext",
    "BbMarginalMoments", "BbStandardize", "BbAnalyseOne", "BbPool", "BbMcse",
    "MiStudyDesign", "MiGenerateBase", "MiDropoutPredictor", "MiCalibrateDropout",
    "MiGenerateDropout", "MiCompleteTruths", "MiDkMean", "MiSelectedIterations", "MiMakeTasks",
    "MiPrepareData", "MiMakeModel", "MiTargetLabels", "MiDirectMeans", "MiSummariseDirect",
    "MiFitDk", "MiImpute", "MiAnalyseBb", "MiRunReplicate")
}

MiWorkerNames <- function() {
  unique(c(MiStatisticalNames(), "DkHashObject", "DkCommitFile", "DkSaveRds", "DkReadRds",
           "DkWriteLines", "DkApplyDiagnosticFlags", "MiBbDiagnostic", "MiMissingness", "MiReadResult"))
}

MiSpecification <- function(Study, Settings) {
  Execution <- c("Stage", "OutputRoot", "SourceStudyDirectory", "PilotReplicates", "FullReplicates",
                 "KeepChainFiles", "UseParallel", "CoreFraction", "RetryFailures", "InstallMissingPackages",
                 "PrintProgress", "ProgressEvery", "RhatThreshold", "EssThreshold", "BbCheckpointEvery",
                 "TruthSeed", "TruthBatchSize", "TruthMinBatches", "TruthMaxBatches", "TruthMcseTolerance")
  Functions <- setNames(lapply(MiStatisticalNames(), function(Name) {
    Function <- get(Name, envir = environment(MiSpecification), inherits = TRUE)
    list(Formals = deparse(formals(Function)), Body = deparse(body(Function)))
  }), MiStatisticalNames())
  list(Version = "DkImputationSimulation1.0.0", Settings = Settings[setdiff(names(Settings), Execution)],
       Design = Study$Design, BaselineCenter = Study$BaselineCenter, Dropout = Study$Dropout$Intercepts,
       Functions = Functions, RVersion = R.version.string, Platform = R.version$platform,
       Packages = vapply(c("HSAUR2", "posterior", "glmmTMB", "TMB", "Matrix"),
                         function(Name) as.character(utils::packageVersion(Name)), character(1)))
}

MiWriteDesign <- function(Study, Settings, Directory, Truth) {
  DkWriteCsv(Study$Design, file.path(Directory, "FixedDesign.csv"))
  DkWriteCsv(cbind(Study$Metadata, Truth = unname(Study$Scenario$Truth)), file.path(Directory, "ParameterTruths.csv"))
  DkWriteCsv(Study$Dropout$Summary, file.path(Directory, "DropoutCalibration.csv"))
  Lines <- c("Simulation study 2: imputation and treatment comparisons.",
             paste("N =", Settings$N, "; visits:", paste(Settings$VisitMonths, collapse = ", "), "; maximum:", Settings$ScoreMaximum),
             "The balanced treatment assignment and baseline sample are generated once and then fixed.",
             "Both treatment means are standardized to all fixed baseline covariates, with the same weights.",
             paste("Baseline centering value:", format(Study$BaselineCenter, digits = 12), "; scale:", Settings$BaselineScale),
             "All first-visit scores are observed. Missingness thereafter is monotone.",
             paste("Cumulative missingness calibration targets:", paste(Settings$DropoutTargets, collapse = ", ")),
             "At visit j, among participants still observed, dropout has logit probability",
             "  HazardIntercept[j] + DropoutBaselineCoefficient*B +",
             "  DropoutPreviousCoefficient*(PreviousScore - DropoutPreviousCenter)/DropoutPreviousScale.",
             paste("Baseline coefficient:", Settings$DropoutBaselineCoefficient,
                   "; previous-score coefficient:", Settings$DropoutPreviousCoefficient),
             "Intercepts are calibrated once on independent generated trajectories and then fixed.",
             "Dropout uses only observed history, never the current or future outcome.",
             "MAR, halfway, and J2R complete trajectories have the same retained scores and triggers.",
             "Only post-trigger active outcomes use their rule's DK quantile; they share generating Gaussian draws.",
             "One observed-data DK posterior therefore supports all rules in each replicate.",
             "The complete-data posterior from Study 1 is never used for the incomplete data.",
             "True MAR means use the finite survival sum. The reference-arm mean is unchanged under every rule.",
             "Other truths use exact MAR means plus an independent Monte Carlo borrowing shift.",
             "Dropout times are integrated conditional on each simulated MAR history; all fixed covariate values",
             "occur equally often. Means are averaged over histories and dropout, not conditional on realized triggers.",
             paste("Maximum achieved truth MCSE:", format(max(Truth$TruthMCSE), digits = 5)),
             "Four types of uncertainty are distinct: truth integration, between-replicate simulation error,",
             "posterior sampling error, and finite-imputation error. The report's MCSE columns concern replicates.",
             "Normal(0, SD = 5) priors apply to regression coefficients by default; log(shape) has SD = 2.",
             "rho is uniform on (0, 1), with its Jacobian included on the sampled logit coordinate.",
             paste("Configured SDs:", Settings$InterceptPriorSd, Settings$TreatmentPriorSd,
                   Settings$BaselinePriorSd, Settings$LogShapePriorSd),
             paste("MCMC:", Settings$NChains, "chains;", Settings$Warmup, "warmup and", Settings$Sample, "saved iterations each."),
             "All post-warmup draws enter posterior summaries and direct DK means.",
             paste(Settings$NImputations, "joint draws are selected at evenly spaced iterations within chains, with counts differing by at most one."),
             "Parameter draws are paired with retained latent values from that same iteration.",
             "Missing Gaussian draws are reused across rules. Recorded scores are never changed.",
             "The BB model has visit-specific intercepts and treatment effects, a common baseline coefficient,",
             "a common precision, and a normal participant random intercept. It is a working analysis model.",
             "Standardized BB means integrate the normal random intercept with checked Gauss-Hermite quadrature.",
             "Delta-method variances include fixed coefficients and the log random intercept SD, with cross covariances.",
             "Complete BB intervals use the large-sample normal approximation; MI intervals use Rubin's t approximation.",
             "Observed DK means have posterior equal-tailed intervals. The observed DK comparator targets MAR only.",
             "The active-mean table compares every MI rule with all three generating truths and retains complete-data benchmarks.",
             "The companion table reports treatment differences under matching rules. CrossRulePerformance.csv",
             "retains all rule comparisons without treating deliberate sensitivity shifts as recovery failures.",
             "Rubin variance is also compared with repeated-sampling variance. Coverage is not assumed to equal 95%.",
             "Every selected completed dataset must have a successful BB fit before its rule is pooled.",
             "Diagnostic flags do not exclude draws or replicates. Failures are recorded; retries use the same seeds.",
             "Each performance row reports its available replicate count; paired comparisons require both estimates.",
             "Both manuscript tables are restricted to the final visit; CSV outputs retain all visits.",
             "Pilot and full stages share checkpoints when scientific, sampler, and software settings agree.",
             "Changing report formatting, target replication count, core fraction, or stage does not invalidate fits.",
             "Only selected joint draws/completions, summaries, and compact BB results are retained by default.")
  DkWriteLines(Lines, file.path(Directory, "StudyDesign.txt"))
}

RunDkImputationSimulation <- function(Settings = DefaultDkImputationSettings()) {
  MiCheckSettings(Settings)
  MiEnsurePackages(Settings)
  CheckDkDistribution()
  dir.create(Settings$OutputRoot, recursive = TRUE, showWarnings = FALSE)
  Root <- normalizePath(Settings$OutputRoot, winslash = "/", mustWork = TRUE)
  Lock <- filelock::lock(file.path(Root, "Run.lock"), timeout = 0)
  if (is.null(Lock)) stop("Another R process is running this simulation folder.")
  Cluster <- NULL
  WorkerPids <- integer()
  ThreadNames <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")
  OldThreads <- Sys.getenv(ThreadNames, unset = NA_character_)
  on.exit({
    if (!is.null(Cluster)) {
      if (length(WorkerPids)) try(tools::pskill(WorkerPids), silent = TRUE)
      try(parallel::stopCluster(Cluster), silent = TRUE)
    }
    filelock::unlock(Lock)
    Missing <- is.na(OldThreads)
    if (any(Missing)) Sys.unsetenv(ThreadNames[Missing])
    if (any(!Missing)) do.call(Sys.setenv, as.list(OldThreads[!Missing]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(setNames(rep("1", length(ThreadNames)), ThreadNames)))
  if (Settings$PrintProgress) {
    cat("\nDK imputation simulation:", Settings$Stage, "\nPreparing the fixed design and generating truth.\n")
    utils::flush.console()
  }
  Study <- MiStudyDesign(Settings)
  CalibrationNames <- c("N", "VisitMonths", "ScoreMaximum", "Tau", "InterceptTruth", "TreatmentTruth",
                        "BaselineTruth", "ShapeTruth", "RhoTruth", "DropoutTargets", "DropoutBaselineCoefficient",
                        "DropoutPreviousCoefficient", "DropoutPreviousCenter", "DropoutPreviousScale", "CalibrationSize", "CalibrationSeed")
  CalibrationKey <- DkHashObject(list(Design = Study$Design, Settings = Settings[CalibrationNames],
                                      Function = deparse(body(MiCalibrateDropout)),
                                      Generator = deparse(body(MiGenerateBase)),
                                      Predictor = deparse(body(MiDropoutPredictor)),
                                      Quantile = deparse(body(DkQuantileFromLogSurvival)),
                                      Helpers = lapply(c("DkWithSeed", "MakeContinuousAr1", "DkLogShapes", "Softplus",
                                                         "LogOneMinusExp", "LogNegativeLogOneMinusExp"), function(Name) {
                                                           Function <- get(Name, envir = environment(RunDkImputationSimulation), inherits = TRUE)
                                                           list(Formals = deparse(formals(Function)), Body = deparse(body(Function)))
                                                         })))
  CalibrationDirectory <- file.path(Root, "Calibration")
  dir.create(CalibrationDirectory, showWarnings = FALSE)
  CalibrationFile <- file.path(CalibrationDirectory, paste0(CalibrationKey, ".rds"))
  Study$Dropout <- DkReadRds(CalibrationFile)
  if (is.null(Study$Dropout)) {
    Study$Dropout <- MiCalibrateDropout(Study, Settings)
    DkSaveRds(Study$Dropout, CalibrationFile)
  }
  Specification <- MiSpecification(Study, Settings)
  Signature <- DkHashObject(Specification)
  Directory <- file.path(Root, paste0("Study", substr(Signature, 1, 12)))
  DesignDirectory <- file.path(Directory, "Design")
  dir.create(DesignDirectory, recursive = TRUE, showWarnings = FALSE)
  Existing <- DkReadRds(file.path(Directory, "StudySpecification.rds"))
  if (!is.null(Existing) && !identical(Existing, Specification)) stop("Study identity mismatch. Use another OutputRoot.")
  if (is.null(Existing)) DkSaveRds(Specification, file.path(Directory, "StudySpecification.rds"))
  Truth <- MiPrepareTruth(Study, Settings, DesignDirectory, Signature)
  MiWriteDesign(Study, Settings, DesignDirectory, Truth)
  Tasks <- MiMakeTasks(Settings, Directory, Signature)
  Results <- lapply(Tasks, MiReadResult)
  Pending <- which(vapply(Results, function(Result) is.null(Result) ||
                            (!isTRUE(Result$Success) && Settings$RetryFailures), logical(1)))
  if (!length(Pending)) cat("No replicates are pending; refreshing reports from saved results.\n")
  Cores <- DkChooseWorkers(Settings, length(Pending))
  ReportDirectory <- file.path(Directory, Settings$Stage)
  Reports <- MiWriteReports(Tasks, Results, Study, Settings, Truth, ReportDirectory, Cores)
  DkWriteLines(capture.output(utils::sessionInfo()), file.path(ReportDirectory, "SessionInfo.txt"))
  cat("Available cores:", Cores$Available, ";", paste0(round(100*Settings$CoreFraction), "%"),
      "limit:", Cores$Limit, "; workers:", Cores$Workers, "\n",
      "Replicates scheduled:", length(Tasks), "; already completed:", sum(Reports$Status$Status == "Completed"),
      "; pending or retrying:", length(Pending), "\n",
      "Per replicate:", Settings$NChains, "DK chains;", Settings$NImputations, "imputations per rule;",
      3*Settings$NImputations + 3, "BB fits.\n",
      "Report folder:", ReportDirectory, "\n", sep = " ")
  utils::flush.console()
  Started <- proc.time()[[3]]
  if (length(Pending)) {
    Context <- list(Study = Study, Settings = Settings)
    Status <- Reports$Status$Status
    SinceReport <- 0
    SaveReports <- function() {
      Reports <<- MiWriteReports(Tasks, Results, Study, Settings, Truth, ReportDirectory, Cores)
      SinceReport <<- 0
      cat("Progress:", sum(Status == "Completed"), "completed,", sum(Status == "Failed"), "failed,",
          sum(Status == "Pending"), "pending; elapsed", round((proc.time()[[3]] - Started)/60, 1), "minutes.\n")
      utils::flush.console()
    }
    OnResult <- function(Index, Result) {
      if (!is.list(Result) || !identical(Result$StudySignature, Signature) ||
          !identical(Result$Replicate, Tasks[[Index]]$Replicate)) stop("A worker returned an unexpected task.")
      Results[Index] <<- list(Result)
      Status[Index] <<- if (isTRUE(Result$Success)) "Completed" else "Failed"
      SinceReport <<- SinceReport + 1
      if (SinceReport >= max(1, Cores$Workers)) SaveReports()
    }
    OnHeartbeat <- function(Running, Waiting) {
      if (Settings$PrintProgress) {
        cat("Status:", sum(Status == "Completed"), "completed,", sum(Status == "Failed"), "failed;",
            Running, "running,", Waiting, "waiting; elapsed", round((proc.time()[[3]] - Started)/60, 1), "minutes.\n")
        utils::flush.console()
      }
    }
    if (Settings$UseParallel && Cores$Workers > 1) {
      Cluster <- parallel::makePSOCKcluster(Cores$Workers, outfile = "", rscript_args = "--vanilla")
      Functions <- mget(MiWorkerNames(), envir = environment(RunDkImputationSimulation), inherits = TRUE)
      WorkerPids <- unlist(parallel::clusterCall(Cluster, function(Paths, Functions, Context) {
        .libPaths(Paths)
        for (Name in names(Functions)) {
          Function <- Functions[[Name]]
          environment(Function) <- .GlobalEnv
          assign(Name, Function, envir = .GlobalEnv)
        }
        assign("MiWorkerContext", Context, envir = .GlobalEnv)
        Sys.getpid()
      }, .libPaths(), Functions, Context), use.names = FALSE)
      MiRunTaskQueue(Cluster, Tasks, Pending, OnResult, OnHeartbeat)
    } else {
      for (Index in Pending) OnResult(Index, MiRunReplicate(Tasks[[Index]], Context))
    }
    if (SinceReport > 0) SaveReports()
  }
  if (!is.null(Cluster)) {
    parallel::stopCluster(Cluster)
    Cluster <- NULL
  }
  Reports <- MiWriteReports(Tasks, Results, Study, Settings, Truth, ReportDirectory, Cores)
  cat("\nFinished:", Reports$Diagnostics$Completed, "completed;", Reports$Diagnostics$Failed,
      "failed.\nActive-mean table:", file.path(ReportDirectory, "ImputationPerformance.tex"),
      "\nTreatment-difference table:", file.path(ReportDirectory, "ImputationContrasts.tex"),
      "\nStandalone tables:", file.path(ReportDirectory, "ImputationPerformanceStandalone.tex"), "\n")
  if (Settings$Stage == "Pilot") cat("Review the pilot and PilotTiming.txt before setting Stage = 'Full'. Compatible results will be reused.\n")
  if (Reports$Diagnostics$Failed > 0) cat("See FailureLog.csv and BbFitDiagnostics.csv. Rerunning retries the same datasets.\n")
  invisible(list(StudyDirectory = Directory, ReportDirectory = ReportDirectory, Truth = Truth,
                 Performance = Reports$Performance, Diagnostics = Reports$Diagnostics,
                 AllSuccessful = Reports$Diagnostics$Completed == length(Tasks)))
}

if (isTRUE(getOption("DkImputationRun", TRUE))) {
  DkImputationSettings <- DefaultDkImputationSettings()
  DkImputationResult <- RunDkImputationSimulation(DkImputationSettings)
}