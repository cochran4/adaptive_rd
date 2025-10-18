# R/data_io.R
# Read, merge, and prepare variables needed for Pooled Cohort Equations (ASCVD).

# ---- helpers ----

.safe_read_xpt <- function(path) {
  stopifnot(file.exists(path))
  df <- foreign::read.xport(path)
  if (!"SEQN" %in% names(df)) stop("SEQN not found in: ", basename(path))
  tibble::as_tibble(df)
}

# ---- main API ----

# 1) Read all component files from a directory (your P_*.XPT names)
read_nhanes_2017_2018 <- function(dir) {
  files <- c(
    demo   = file.path(dir, "P_DEMO.XPT"),
    bpx    = file.path(dir, "P_BPXO.XPT"),
    diq    = file.path(dir, "P_DIQ.XPT"),
    bpq    = file.path(dir, "P_BPQ.XPT"),
    hdl    = file.path(dir, "P_HDL.XPT"),
    mcq    = file.path(dir, "P_MCQ.XPT"),    # not used by PCE, but ok to keep
    smq    = file.path(dir, "P_SMQ.XPT"),
    tchol  = file.path(dir, "P_TCHOL.XPT")
  )
  if (any(!file.exists(files))) {
    missing <- names(files)[!file.exists(files)]
    stop("Missing expected XPT files: ", paste(missing, collapse = ", "))
  }
  purrr::map(files, .safe_read_xpt)
}

# 2) Build cohort dataset with PCE covariates
build_cohort <- function(
    nhanes_list,
    seed      = 1L,
    age_min   = 40L,
    age_max   = 79L,
    impute    = TRUE,
    summarize = FALSE   
) {
  
  # Merge everything on SEQN
  nhanes <- Reduce(function(x, y) dplyr::full_join(x, y, by = "SEQN"), nhanes_list)
  
  # Keep only fields needed for PCE
  # DEMO: age, sex, race
  # BPXO: systolic readings
  # BPQ:  BP meds
  # SMQ:  smoking (use current smoking SMQ040 if available)
  # DIQ:  diabetes
  # HDL:  HDL-C
  # TCHOL: total cholesterol
  data <- nhanes %>%
    select(
      SEQN,
      RIDAGEYR, RIAGENDR, RIDRETH3,                # age, sex, race
      BPXOSY1, BPXOSY2, BPXOSY3,                   # systolic readings (some may be missing)
      BPQ050A,                                     # BP meds (1 yes, 2 no)
      SMQ020, SMQ040,                              # smoker items
      DIQ010,                                      # diabetes (1 yes, 2 no)
      LBDHDD,                                      # HDL cholesterol
      LBXTC                                        # total cholesterol
    )
  
  # Filter to PCE age range (focus on 40–79)
  data <- data %>%
    filter(!is.na(RIDAGEYR), RIDAGEYR >= age_min, RIDAGEYR <= age_max)
  
  # Type normalization
  num_vars <- c("RIDAGEYR","BPXOSY1","BPXOSY2","BPXOSY3","LBDHDD","LBXTC")
  data <- data %>%
    mutate(across(all_of(num_vars), ~ suppressWarnings(as.numeric(.x))))
  
  # Average SBP across available readings (same visit)
  data <- data %>%
    mutate(
      SBP_mean = rowMeans(select(., BPXOSY1, BPXOSY2, BPXOSY3), na.rm = TRUE)
    )
  
  # Derive PCE covariates
  # Sex: 1=Male, 2=Female
  # Race: RIDRETH3 — 3 White, 4 Black; others -> "Other"
  # Current smoker: SMQ040 (1 every day, 2 some days = current; 3 not at all = non-current)
  # Diabetes: DIQ010 (1 yes, 2 no)
  # BP treated: BPQ050A (1 yes, 2 no)
  data <- data %>%
    mutate(
      sex = case_when(
        RIAGENDR == 1 ~ "Male",
        RIAGENDR == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      race = case_when(
        RIDRETH3 == 3 ~ "White",
        RIDRETH3 == 4 ~ "Black",
        TRUE ~ "Other"
      ),
      smoker_current = case_when(
        SMQ040 %in% c(1, 2) ~ 1L,
        SMQ040 %in% c(3)    ~ 0L,
        TRUE                ~ NA_integer_
      ),
      diabetes = case_when(
        DIQ010 == 1 ~ 1L,
        DIQ010 == 2 ~ 0L,
        TRUE        ~ NA_integer_
      ),
      bp_treated = case_when(
        BPQ050A == 1 ~ 1L,
        BPQ050A == 2 ~ 0L,
        TRUE         ~ NA_integer_
      )
    ) %>%
    # rename to standardized labels
    rename(
      person_id = SEQN,
      age     = RIDAGEYR,
      sbp     = SBP_mean,  
      hdl_c   = LBDHDD,
      tot_chol= LBXTC
    )
  
  # Minimal column set for PCE (with new names)
  pce_vars <- c("person_id","age","sex","race","sbp","smoker_current",
                "diabetes","bp_treated","hdl_c","tot_chol")
  
  data_pce <- data %>% select(all_of(pce_vars))

  if (!isTRUE(impute)) {
    out <- tibble::as_tibble(data_pce)
  } else {
    # Impute missing values
    set.seed(seed)
    imp <- mice::mice(data_pce, printFlag = FALSE)
    out <- mice::complete(imp) %>%
      dplyr::mutate(
        smoker_current = as.integer(round(smoker_current)),
        diabetes       = as.integer(round(diabetes)),
        bp_treated     = as.integer(round(bp_treated)),
        age            = as.numeric(age),
        sbp            = as.numeric(sbp),
        hdl_c          = as.numeric(hdl_c),
        tot_chol       = as.numeric(tot_chol)
      ) %>%
      tibble::as_tibble()
  }
  
  # Build terms needed to run risk prediction model
  out <- .build_terms(out)
  
  
  # If no summary requested, return the modeling tibble (back-compat)
  if (!isTRUE(summarize)) {
    return(out)
  }
  
  # Concise summaries on the PCE-ready variables (not the derived terms)
  summ <- summarize_pce(data_pce)
  
  # Return both
  list(
    data    = out,     # with PCE terms (for risk/modeling)
    summary = summ     # compact baseline tables (for reporting)
  )
}

# --------------------------------------------------------------------
# Helper: Safe natural logarithm that avoids -Inf for non-positive values
# Clamps values to at least min_val before taking log.
# --------------------------------------------------------------------
.safe_ln <- function(x, min_val = 1e-6) log(pmax(x, min_val))

# --------------------------------------------------------------------
# Helper: Map (sex, race) to one of the four PCE groups
# Treat "Other" as White by convention (PCE supports White/Black only).
# Returns a character vector of group labels:
#   "White_F", "White_M", "Black_F", "Black_M"
# --------------------------------------------------------------------
.pce_group <- function(sex, race) {
  dplyr::case_when(
    sex == "Female" & race == "White" ~ "White_F",
    sex == "Male"   & race == "White" ~ "White_M",
    sex == "Female" & race == "Other" ~ "Other_F",
    sex == "Male"   & race == "Other" ~ "Other_M",
    sex == "Female" & race == "Black" ~ "Black_F",
    sex == "Male"   & race == "Black" ~ "Black_M",
    TRUE ~ NA_character_  # triggers a hard stop later if any NA present
  )
}

# --------------------------------------------------------------------
# Build the model terms required by the PCE from standardized column names.
# The names here must match the coefficient names in `betas` for each group.
#
# Input: df (tibble/data.frame) with columns listed at the top of this file.
# Output: tibble with one row per person, columns:
#   group, ln_age, ln_age_sqrd, ln_total_cholest, ln_hdlC,
#   ln_treated_BP, ln_untreated_BP,
#   smoker, diabetes,
#   ln_age_totcholest, ln_age_hdlC, ln_age_BP, ln_age_ln_untreated_BP, ln_age_smoker
# --------------------------------------------------------------------
.build_terms <- function(df) {
  df %>%
    transmute(
      # Group label used to pick the right coefficient set
      group = .pce_group(sex, race),
      
      # Core logs per PCE
      ln_age           = .safe_ln(age),
      ln_age_sqrd      = (.safe_ln(age))^2,
      ln_total_cholest = .safe_ln(tot_chol),
      ln_hdlC          = .safe_ln(hdl_c),
      
      # SBP enters via treated / untreated components in PCE
      ln_sbp           = .safe_ln(sbp),
      ln_treated_BP    = ifelse(bp_treated == 1, ln_sbp, 0),
      ln_untreated_BP  = ifelse(bp_treated == 0, ln_sbp, 0),
      
      # Binary covariates (coded 0/1)
      smoker   = as.numeric(smoker_current == 1),
      diabetes = as.numeric(diabetes == 1)
    ) %>%
    mutate(
      # Interactions used in PCE
      ln_age_totcholest      = ln_age * ln_total_cholest,
      ln_age_hdlC            = ln_age * ln_hdlC,
      ln_age_BP              = ln_age * ln_treated_BP,
      ln_age_ln_untreated_BP = ln_age * ln_untreated_BP,
      ln_age_smoker          = ln_age * smoker
    ) %>%
    # ln_sbp was a helper to build treated/untreated terms; drop it now
    select(-ln_sbp) %>%
    as_tibble()
}

# ---- Concise summaries for PCE-ready columns (with missingness) ----
summarize_pce <- function(
    df,
    cont_vars = c("age","sbp","hdl_c","tot_chol"),
    cat_vars  = c("sex","race","smoker_current","diabetes","bp_treated")
) {
  # Continuous: mean, sd, min, max, N, Missing
  continuous_tbl <- df |>
    dplyr::summarise(dplyr::across(
      dplyr::all_of(cont_vars),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x,   na.rm = TRUE),
        min  = ~min(.x,  na.rm = TRUE),
        max  = ~max(.x,  na.rm = TRUE),
        N    = ~length(.x),
        Missing = ~sum(is.na(.x))
      )
    )) |>
    # Use regex to split at the LAST underscore only
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to   = c("Variable", ".value"),
      names_pattern = "^(.*)_(mean|sd|min|max|N|Missing)$"
    ) |>
    dplyr::mutate(
      Mean = round(mean, 1),
      SD   = round(sd,   1),
      Min  = round(min,  1),
      Max  = round(max,  1)
    ) |>
    dplyr::select(Variable, N, Missing, Mean, SD, Min, Max)
  
  # Categorical: counts and percents; NA -> "Missing"
  cat_list <- lapply(cat_vars, function(v) {
    vals <- df[[v]]
    # Explicit "Missing" label for NA
    lvl  <- ifelse(is.na(vals), "Missing", as.character(vals))
    tibble::tibble(Level = lvl) |>
      dplyr::count(Level, name = "n") |>
      dplyr::mutate(
        Variable = v,
        Percent  = round(100 * n / sum(n), 1)
      ) |>
      dplyr::relocate(Variable)
  })
  categorical_tbl <- dplyr::bind_rows(cat_list)
  
  list(continuous = continuous_tbl, categorical = categorical_tbl)
}