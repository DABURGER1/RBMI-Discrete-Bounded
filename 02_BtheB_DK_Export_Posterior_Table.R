# Export posterior summaries for the BtheB DK Gaussian copula imputation model.
# Uses base R and all saved post-warmup draws, without refitting or imputation.
# Sourcing this file exports the table from the default application output folder.
#
# Default export:
# source("02_BtheB_DK_Export_Posterior_Table.R")
#
# Export from an existing result object:
# DkPosteriorTable <- ExportDkPosteriorTable(BthebDkResult)
#
# Export from an explicit result path:
# DkPosteriorTable <- ExportDkPosteriorTable(
#   Result = "BtheB_DK_RBMI/ModelOutputs/BtheB_DK_ApplicationResult.rds"
# )
# ChainDirectory can specify the folder containing the saved chain files.
#
# Outputs: a LaTeX table, a standalone LaTeX document, and an unrounded CSV.
# Table entries use three decimal places by default.
# The exported table requires booktabs, tabularx, and siunitx and supports \input.
# Summaries describe the fitted model under the assigned treatment, shared across rules.

DkTableNumberFormat <- function(Values, Digits) {
  Precision <- max(Digits)
  Rounded <- round(Values, Digits)
  Rounded <- Rounded[is.finite(Rounded)]
  if (length(Rounded) == 0) return(paste0("1.", Precision))
  IntegerDigits <- max(nchar(format(floor(abs(Rounded)), trim = TRUE, scientific = FALSE)))
  paste0(if (any(Rounded < 0)) "-" else "", IntegerDigits, ".", Precision)
}


# Interval spacing retains at least one ordinary word space in every row.
# Padding accounts for rounded values, signs, and displayed precision;
# LaTeX measures the required space in the current font.
DkTableIntervalSpacing <- function(Lower, Upper) {
  Fraction <- function(Value) {
    if (grepl(".", Value, fixed = TRUE)) sub("^[^.]*[.]", "", Value) else ""
  }
  LowerPrecision <- nchar(vapply(Lower, Fraction, character(1)))
  UpperInteger <- sub("[.].*$", "", Upper)
  UpperDigits <- nchar(sub("^-", "", UpperInteger))
  UpperSigned <- any(startsWith(UpperInteger, "-"))
  Padding <- unique(vapply(seq_along(Lower), function(Index) {
    paste0(
      if (LowerPrecision[Index] == 0 && max(LowerPrecision) > 0) "." else "",
      strrep("0", max(LowerPrecision) - LowerPrecision[Index]),
      if (UpperSigned && !startsWith(UpperInteger[Index], "-")) "-" else "",
      strrep("0", max(UpperDigits) - UpperDigits[Index])
    )
  }, character(1)))
  Lines <- c(
    "% Space interval bounds from their displayed contents, retaining decimal alignment.",
    "\\ifdefined\\DkCrIGap\\else\\newlength{\\DkCrIGap}\\fi",
    "\\ifdefined\\DkCrIPadding\\else\\newlength{\\DkCrIPadding}\\fi",
    "\\ifdefined\\DkCrICandidate\\else\\newlength{\\DkCrICandidate}\\fi",
    "\\setlength{\\DkCrIPadding}{\\maxdimen}"
  )
  for (Value in Padding) {
    Lines <- c(Lines,
               if (nzchar(Value)) paste0("\\settowidth{\\DkCrICandidate}{$", Value, "$}") else
                 "\\setlength{\\DkCrICandidate}{0pt}",
               "\\ifdim\\DkCrICandidate<\\DkCrIPadding",
               "  \\setlength{\\DkCrIPadding}{\\DkCrICandidate}\\fi"
    )
  }
  c(Lines, "\\setlength{\\DkCrIGap}{\\dimexpr\\fontdimen2\\font - \\DkCrIPadding\\relax}")
}

ExportDkPosteriorTable <- function(
    Result = file.path("BtheB_DK_RBMI", "ModelOutputs",
                       "BtheB_DK_ApplicationResult.rds"),
    OutputDirectory = NULL, ChainDirectory = NULL,
    Digits = 3, AnchorDigits = 3, RhoDigits = 3) {
  
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
    # Avoid displaying a negative zero after rounding.
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
    "\\mbox{\\tablenum[table-format = ", CommonFormats["Lower95"],
    ", table-space-text-pre = {[}, table-space-text-post = {,}]{", Lower(Row), "}",
    "\\hspace{\\DkCrIGap}\\tablenum[table-format = ", CommonFormats["Upper95"],
    ", table-space-text-post = {]}]{", Upper(Row), "}}"
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
    DkTableIntervalSpacing(
      mapply(Number, CommonRows$Lower95, CommonRows$Digits, USE.NAMES = FALSE),
      mapply(Number, CommonRows$Upper95, CommonRows$Digits, USE.NAMES = FALSE)
    ),
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
    "\\begin{minipage}{\\linewidth}", "\\footnotesize",
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
    "\\captionsetup{justification=justified,singlelinecheck=false}",
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

DkPosteriorTable <- ExportDkPosteriorTable()