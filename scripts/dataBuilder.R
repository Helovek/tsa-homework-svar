#!/usr/bin/env Rscript

# Настройки окружения и путей
options(scipen = 999)

processed_dir <- file.path("data", "processed")
diagnostics_path <- file.path(processed_dir, "build_diagnostics.txt")
local_library_dir <- file.path(processed_dir, "r_lib")

ensure_processed_dir <- function() {
  if (!dir.exists(processed_dir)) {
    dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

configure_library_paths <- function() {
  ensure_processed_dir()
  if (!dir.exists(local_library_dir)) {
    dir.create(local_library_dir, recursive = TRUE, showWarnings = FALSE)
  }
  user_library_dir <- Sys.getenv("R_LIBS_USER")
  if (nzchar(user_library_dir) && !dir.exists(user_library_dir)) {
    dir.create(user_library_dir, recursive = TRUE, showWarnings = FALSE)
  }
  candidate_paths <- c(user_library_dir, local_library_dir, .libPaths())
  candidate_paths <- candidate_paths[nzchar(candidate_paths)]
  existing_paths <- candidate_paths[dir.exists(candidate_paths)]
  .libPaths(unique(normalizePath(existing_paths, winslash = "/", mustWork = FALSE)))
}

configure_library_paths()

write_diagnostics <- function(messages) {
  ensure_processed_dir()
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  lines <- c(
    sprintf("build time: %s", timestamp),
    "messages:",
    paste0("- ", messages)
  )
  writeLines(lines, diagnostics_path, useBytes = TRUE)
}

stop_with_diagnostics <- function(...) {
  message_text <- paste(..., collapse = "")
  write_diagnostics(message_text)
  stop(message_text, call. = FALSE)
}

# Проверка зависимостей
required_packages <- c(
  "readxl",
  "openxlsx",
  "dplyr",
  "tidyr",
  "lubridate",
  "stringr",
  "zoo",
  "seasonal",
  "x13binary"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop_with_diagnostics(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them manually before running scripts/dataBuilder.R."
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(zoo)
  library(seasonal)
  library(x13binary)
})

series_names <- c(
  "gdp",
  "cpi",
  "wage",
  "m2",
  "key_rate",
  "p_im",
  "p_cm",
  "overdue_debt"
)

raw_dir <- file.path("data", "raw")

# Пути к исходным файлам
file_paths <- list(
  gdp = file.path(raw_dir, "gdp", "gdp_real_sa_iep.xlsx"),
  m2 = file.path(raw_dir, "m2", "M2-M2_SA.xlsx"),
  key_rate = file.path(raw_dir, "key_rate", "inflation_key_rate.xlsx"),
  overdue_post = file.path(raw_dir, "overdue_debt", "01_05_D_Debt_subj.xlsx"),
  overdue_pre_a = file.path(raw_dir, "overdue_debt", "302-09.xlsx"),
  overdue_pre_b = file.path(raw_dir, "overdue_debt", "302-10.xlsx"),
  p_cm_old = file.path(raw_dir, "p_cm", "43561_C.xls"),
  p_cm_new = file.path(raw_dir, "p_cm", "57609_B.xls"),
  p_im_old = file.path(raw_dir, "p_im", "43561_DF+DG+DJ.xls"),
  p_im_new = file.path(raw_dir, "p_im", "57609_19+20+24+25.xls"),
  p_im_weights = file.path(raw_dir, "p_im", "vesa_2006-2026.xlsx"),
  cpi = file.path(raw_dir, "cpi", "cpi_q.xlsx"),
  wage = file.path(raw_dir, "wage", "ind_04-2026.xlsx")
)

missing_files <- names(file_paths)[!file.exists(unlist(file_paths, use.names = FALSE))]
if (length(missing_files) > 0) {
  stop_with_diagnostics(
    "Missing source files: ",
    paste(missing_files, collapse = ", ")
  )
}

# Базовые преобразования типов
coerce_numeric <- function(x) {
  if (inherits(x, "Date")) {
    return(as.numeric(x))
  }
  if (inherits(x, "POSIXt")) {
    return(as.numeric(as.Date(x)))
  }
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00a0", "")
  x <- stringr::str_replace_all(x, "\\s+", "")
  x <- stringr::str_replace_all(x, ",", ".")
  x <- stringr::str_replace_all(x, "%", "")
  x[x %in% c("", "NA", "NaN", "-", "--")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

parse_dmy_date <- function(x) {
  if (inherits(x, "Date")) {
    return(as.Date(x))
  }
  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }
  if (is.numeric(x)) {
    return(as.Date(x, origin = "1899-12-30"))
  }
  suppressWarnings(as.Date(as.character(x), format = "%d.%m.%Y"))
}

parse_month_date <- function(x) {
  if (inherits(x, "Date")) {
    return(as.Date(format(x, "%Y-%m-01")))
  }
  if (inherits(x, "POSIXt")) {
    return(as.Date(format(as.Date(x), "%Y-%m-01")))
  }
  if (is.numeric(x)) {
    x <- as.Date(x, origin = "1899-12-30")
    return(as.Date(format(x, "%Y-%m-01")))
  }

  char_x <- as.character(x)
  out <- as.Date(rep(NA_character_, length(char_x)))
  month_year_mask <- stringr::str_detect(char_x, "^\\d{1,2}\\.\\d{4}$")

  if (any(month_year_mask, na.rm = TRUE)) {
    parsed <- suppressWarnings(as.Date(
      paste0("01.", char_x[month_year_mask]),
      format = "%d.%m.%Y"
    ))
    out[month_year_mask] <- as.Date(format(parsed, "%Y-%m-01"))
  }

  numeric_x <- suppressWarnings(as.numeric(char_x))
  numeric_mask <- !month_year_mask & !is.na(numeric_x)
  if (any(numeric_mask)) {
    out[numeric_mask] <- as.Date(
      format(as.Date(numeric_x[numeric_mask], origin = "1899-12-30"), "%Y-%m-01")
    )
  }

  text_mask <- !month_year_mask & !numeric_mask & !is.na(char_x)
  if (any(text_mask)) {
    parsed <- suppressWarnings(as.Date(
      paste0("01.", char_x[text_mask]),
      format = "%d.%m.%Y"
    ))
    out[text_mask] <- as.Date(format(parsed, "%Y-%m-01"))
  }

  out
}

# Работа с датами и кварталами
quarter_end_date <- function(year_value, quarter_value) {
  as.Date(sprintf(
    "%04d-%s",
    year_value,
    c("03-31", "06-30", "09-30", "12-31")[quarter_value]
  ))
}

excel_col_to_index <- function(col_name) {
  letters_vec <- strsplit(toupper(col_name), "")[[1]]
  Reduce(function(acc, letter) {
    acc * 26 + match(letter, LETTERS)
  }, letters_vec, init = 0)
}

# Чтение и валидация рядов
read_sheet_frame <- function(path, sheet = 1, range = NULL, col_types = NULL) {
  args <- list(
    path = path,
    sheet = sheet,
    col_names = FALSE,
    .name_repair = "minimal"
  )
  if (!is.null(range)) {
    args$range <- range
  }
  if (!is.null(col_types)) {
    args$col_types <- col_types
  }
  as.data.frame(do.call(readxl::read_excel, args), stringsAsFactors = FALSE)
}

sort_unique_dates <- function(df, value_col) {
  df <- df[!is.na(df$date) & !is.na(df[[value_col]]), , drop = FALSE]
  df <- df[order(df$date), , drop = FALSE]
  if (any(duplicated(df$date))) {
    stop_with_diagnostics("Duplicate dates found in series: ", value_col)
  }
  rownames(df) <- NULL
  df
}

assert_monthly_regular <- function(df, value_col) {
  df <- sort_unique_dates(df, value_col)
  expected <- seq.Date(min(df$date), max(df$date), by = "month")
  if (length(expected) != nrow(df) || any(expected != df$date)) {
    missing_dates <- expected[!expected %in% df$date]
    stop_with_diagnostics(
      "Monthly gaps detected in series ",
      value_col,
      ": ",
      paste(format(missing_dates, "%Y-%m-%d"), collapse = ", ")
    )
  }
  df
}

assert_quarterly_regular <- function(df, value_col) {
  df <- sort_unique_dates(df, value_col)
  quarter_id <- lubridate::year(df$date) * 4L + lubridate::quarter(df$date)
  expected_id <- seq.int(min(quarter_id), max(quarter_id), by = 1L)
  if (length(expected_id) != nrow(df) || any(expected_id != quarter_id)) {
    missing_id <- expected_id[!expected_id %in% quarter_id]
    missing_year <- (missing_id - 1L) %/% 4L
    missing_quarter <- missing_id - missing_year * 4L
    missing_dates <- quarter_end_date(missing_year, missing_quarter)
    stop_with_diagnostics(
      "Quarterly gaps detected in series ",
      value_col,
      ": ",
      paste(format(missing_dates, "%Y-%m-%d"), collapse = ", ")
    )
  }
  df
}

month_to_quarter <- function(df, value_col) {
  df <- assert_monthly_regular(df, value_col)
  quarter_df <- df %>%
    mutate(
      year_value = lubridate::year(date),
      quarter_value = lubridate::quarter(date),
      month_weight = lubridate::days_in_month(date)
    ) %>%
    group_by(year_value, quarter_value) %>%
    summarise(
      n_months = n(),
      value = stats::weighted.mean(.data[[value_col]], month_weight),
      .groups = "drop"
    )
  quarter_df <- quarter_df[quarter_df$n_months == 3L, c("year_value", "quarter_value", "value"), drop = FALSE]
  quarter_df <- data.frame(
    date = quarter_end_date(quarter_df$year_value, quarter_df$quarter_value),
    value = quarter_df$value,
    stringsAsFactors = FALSE
  )
  names(quarter_df)[2] <- value_col
  assert_quarterly_regular(quarter_df, value_col)
}

# Сезонная корректировка уровней
seasonal_adjust <- function(df, value_col, frequency) {
  df <- sort_unique_dates(df, value_col)
  ts_start <- if (frequency == 12L) {
    c(lubridate::year(df$date[1]), lubridate::month(df$date[1]))
  } else if (frequency == 4L) {
    c(lubridate::year(df$date[1]), lubridate::quarter(df$date[1]))
  } else {
    stop_with_diagnostics("Unsupported seasonal frequency: ", frequency)
  }

  ts_data <- stats::ts(df[[value_col]], start = ts_start, frequency = frequency)
  adjusted <- tryCatch(
    {
      fit <- seasonal::seas(ts_data, transform.function = "none")
      as.numeric(seasonal::final(fit))
    },
    error = function(e) {
      stop_with_diagnostics(
        "seasonal::seas failed for series ",
        value_col,
        ": ",
        e$message
      )
    }
  )
  df[[value_col]] <- adjusted
  df
}

# Вспомогательные функции для Excel-таблиц
extract_year_block <- function(sheet_df, row_id, start_col, year_value, step = 1L, months = 12L) {
  if (is.character(start_col)) {
    start_col <- excel_col_to_index(start_col)
  }
  col_ids <- start_col + step * (0:(months - 1L))
  if (max(col_ids) > ncol(sheet_df)) {
    stop_with_diagnostics(
      "Requested columns exceed worksheet width at row ",
      row_id,
      ", year ",
      year_value
    )
  }
  values <- coerce_numeric(unlist(sheet_df[row_id, col_ids], use.names = FALSE))
  data.frame(
    date = seq.Date(as.Date(sprintf("%04d-01-01", year_value)), by = "month", length.out = months),
    value = values,
    stringsAsFactors = FALSE
  )
}

extract_years_until_empty <- function(sheet_df, row_id, start_col, first_year, step = 3L, months = 12L) {
  if (is.character(start_col)) {
    start_col <- excel_col_to_index(start_col)
  }
  out <- list()
  year_value <- first_year
  current_col <- start_col
  repeat {
    col_ids <- current_col + step * (0:(months - 1L))
    if (max(col_ids) > ncol(sheet_df)) {
      break
    }
    values <- coerce_numeric(unlist(sheet_df[row_id, col_ids], use.names = FALSE))
    if (all(is.na(values))) {
      break
    }
    out[[length(out) + 1L]] <- data.frame(
      date = seq.Date(as.Date(sprintf("%04d-01-01", year_value)), by = "month", length.out = months),
      value = values,
      stringsAsFactors = FALSE
    )
    year_value <- year_value + 1L
    current_col <- current_col + step * months
  }
  dplyr::bind_rows(out)
}

build_level_from_indices <- function(index_df, base_date = as.Date("2013-12-01"), base_value = 100) {
  index_df <- sort_unique_dates(index_df, "index")
  index_df <- index_df[index_df$date > base_date, , drop = FALSE]
  if (nrow(index_df) == 0L) {
    stop_with_diagnostics("No monthly index observations found after technical base date.")
  }
  levels <- numeric(nrow(index_df))
  current_level <- base_value
  for (i in seq_len(nrow(index_df))) {
    current_level <- current_level * index_df$index[i] / 100
    levels[i] <- current_level
  }
  out <- data.frame(
    date = index_df$date,
    level = levels,
    stringsAsFactors = FALSE
  )
  assert_monthly_regular(out, "level")
}

write_series_csv <- function(df, series_name) {
  ensure_processed_dir()
  utils::write.csv(
    df,
    file = file.path(processed_dir, paste0(series_name, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

log_if_needed <- function(df, series_name) {
  if (series_name == "key_rate") {
    return(df)
  }
  values <- df[[series_name]]
  if (any(!is.finite(values) | values <= 0)) {
    stop_with_diagnostics(
      "Cannot take logarithm of non-positive values in series ",
      series_name
    )
  }
  df[[series_name]] <- log(values)
  df
}

row_texts <- function(sheet_df) {
  apply(sheet_df, 1, function(row_values) {
    row_values <- as.character(row_values)
    row_values <- row_values[!is.na(row_values) & nzchar(trimws(row_values))]
    paste(row_values, collapse = " | ")
  })
}

parse_month_token <- function(x) {
  lower <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
  month_map <- c(
    "январь" = 1L, "янв" = 1L,
    "февраль" = 2L, "фев" = 2L,
    "март" = 3L, "мар" = 3L,
    "апрель" = 4L, "апр" = 4L,
    "май" = 5L,
    "июнь" = 6L, "июн" = 6L,
    "июль" = 7L, "июл" = 7L,
    "август" = 8L, "авг" = 8L,
    "сентябрь" = 9L, "сен" = 9L,
    "октябрь" = 10L, "окт" = 10L,
    "ноябрь" = 11L, "ноя" = 11L,
    "декабрь" = 12L, "дек" = 12L
  )
  if (lower %in% names(month_map)) {
    return(unname(month_map[lower]))
  }
  numeric_token <- suppressWarnings(as.integer(stringr::str_extract(lower, "(?<!\\d)(1[0-2]|0?[1-9])(?!\\d)")))
  if (!is.na(numeric_token) && numeric_token >= 1L && numeric_token <= 12L) {
    return(numeric_token)
  }
  NA_integer_
}

parse_year_token <- function(x) {
  token <- suppressWarnings(as.integer(stringr::str_extract(as.character(x), "(19|20)\\d{2}")))
  if (is.na(token) || token < 2000L || token > 2100L) {
    return(NA_integer_)
  }
  token
}

detect_monthly_columns <- function(sheet_df, max_header_rows = 10L) {
  header_n <- min(nrow(sheet_df), max_header_rows)
  header_df <- as.data.frame(
    lapply(sheet_df[seq_len(header_n), , drop = FALSE], as.character),
    stringsAsFactors = FALSE
  )
  header_matrix <- as.matrix(header_df)
  month_counts <- vapply(
    seq_len(nrow(header_matrix)),
    function(i) {
      sum(!is.na(vapply(header_matrix[i, ], parse_month_token, integer(1))), na.rm = TRUE)
    },
    integer(1)
  )
  month_row <- which.max(month_counts)
  if (length(month_row) == 0L || month_counts[month_row] < 12L) {
    stop_with_diagnostics("Failed to detect monthly header row in cpi.xls.")
  }
  year_counts <- vapply(
    seq_len(month_row),
    function(i) {
      sum(!is.na(vapply(header_matrix[i, ], parse_year_token, integer(1))), na.rm = TRUE)
    },
    integer(1)
  )
  year_row <- which.max(year_counts)
  if (length(year_row) == 0L || year_counts[year_row] == 0L) {
    stop_with_diagnostics("Failed to detect yearly header row in cpi.xls.")
  }
  year_values <- vapply(header_matrix[year_row, ], parse_year_token, integer(1))
  year_values <- zoo::na.locf(year_values, na.rm = FALSE)
  month_values <- vapply(header_matrix[month_row, ], parse_month_token, integer(1))
  valid_cols <- which(!is.na(year_values) & !is.na(month_values))
  if (length(valid_cols) < 12L) {
    stop_with_diagnostics("Too few monthly columns detected in cpi.xls.")
  }
  data.frame(
    col_id = valid_cols,
    date = as.Date(sprintf("%04d-%02d-01", year_values[valid_cols], month_values[valid_cols])),
    stringsAsFactors = FALSE
  )
}

detect_cpi_monthly_columns <- function(sheet_df) {
  year_row <- as.character(unlist(sheet_df[3, , drop = TRUE], use.names = FALSE))
  month_row <- as.character(unlist(sheet_df[4, , drop = TRUE], use.names = FALSE))
  metric_row <- stringr::str_to_lower(as.character(unlist(sheet_df[5, , drop = TRUE], use.names = FALSE)))

  year_values <- vapply(year_row, parse_year_token, integer(1))
  year_values <- zoo::na.locf(year_values, na.rm = FALSE)
  month_values <- vapply(month_row, parse_month_token, integer(1))

  valid_cols <- which(
    !is.na(year_values) &
      !is.na(month_values) &
      stringr::str_detect(metric_row, stringr::fixed("к предыдущему месяцу"))
  )

  if (length(valid_cols) == 0L) {
    stop_with_diagnostics("Failed to detect monthly CPI columns in cpi.xls.")
  }

  data.frame(
    col_id = valid_cols,
    date = as.Date(sprintf("%04d-%02d-01", year_values[valid_cols], month_values[valid_cols])),
    stringsAsFactors = FALSE
  )
}

select_best_row <- function(sheet_df, required_patterns) {
  text_rows <- row_texts(sheet_df)
  text_rows_lower <- stringr::str_to_lower(text_rows)
  required_patterns <- stringr::str_to_lower(required_patterns)
  matches <- which(vapply(
    text_rows_lower,
    function(x) {
      all(vapply(
        required_patterns,
        function(pattern_text) stringr::str_detect(x, stringr::fixed(pattern_text)),
        logical(1)
      ))
    },
    logical(1)
  ))
  if (length(matches) == 0L) {
    return(integer(0))
  }
  numeric_counts <- vapply(
    matches,
    function(i) sum(!is.na(coerce_numeric(unlist(sheet_df[i, ], use.names = FALSE)))),
    integer(1)
  )
  matches[which.max(numeric_counts)]
}

# Построение отдельных рядов
build_gdp <- function() {
  sheet_df <- readxl::read_excel(file_paths$gdp, sheet = "data")
  out <- data.frame(
    date = parse_dmy_date(sheet_df[[1]]),
    gdp = coerce_numeric(sheet_df[[8]]),
    stringsAsFactors = FALSE
  )
  out <- sort_unique_dates(out, "gdp")
  out <- out[out$date >= as.Date("2013-12-31"), , drop = FALSE]
  assert_quarterly_regular(out, "gdp")
}

build_m2 <- function() {
  sheet_df <- readxl::read_excel(file_paths$m2, sheet = "M2 data")
  monthly <- data.frame(
    date = parse_month_date(sheet_df[[1]]),
    m2 = coerce_numeric(sheet_df[[3]]),
    stringsAsFactors = FALSE
  )
  monthly <- sort_unique_dates(monthly, "m2")
  month_to_quarter(monthly, "m2")
}

build_key_rate <- function() {
  sheet_df <- readxl::read_excel(file_paths$key_rate)
  monthly <- data.frame(
    date = parse_month_date(sheet_df[[1]]),
    key_rate = coerce_numeric(sheet_df[[2]]),
    stringsAsFactors = FALSE
  )
  monthly <- sort_unique_dates(monthly, "key_rate")
  month_to_quarter(monthly, "key_rate")
}

build_overdue_debt <- function() {
  pre_a <- readxl::read_excel(file_paths$overdue_pre_a)
  pre_b <- readxl::read_excel(file_paths$overdue_pre_b)

  pre_monthly <- bind_rows(
    data.frame(
      date = parse_month_date(pre_a[[1]]),
      component = "a",
      value = coerce_numeric(pre_a[[2]]),
      stringsAsFactors = FALSE
    ),
    data.frame(
      date = parse_month_date(pre_b[[1]]),
      component = "b",
      value = coerce_numeric(pre_b[[2]]),
      stringsAsFactors = FALSE
    )
  ) %>%
    filter(!is.na(date), !is.na(value)) %>%
    group_by(date) %>%
    summarise(overdue_debt = sum(value), .groups = "drop")

  post_sheet <- read_sheet_frame(
    file_paths$overdue_post,
    sheet = "в т.ч. просроч. итого"
  )
  post_dates <- parse_month_date(unlist(post_sheet[2, -1, drop = TRUE], use.names = FALSE))
  post_values <- coerce_numeric(unlist(post_sheet[3, -1, drop = TRUE], use.names = FALSE))
  post_monthly <- data.frame(
    date = post_dates,
    overdue_debt = post_values,
    stringsAsFactors = FALSE
  )
  post_monthly <- sort_unique_dates(post_monthly, "overdue_debt")

  cutoff <- min(post_monthly$date, na.rm = TRUE)
  monthly <- bind_rows(
    pre_monthly %>% filter(date < cutoff),
    post_monthly %>% filter(date >= cutoff)
  )
  monthly <- sort_unique_dates(monthly, "overdue_debt")
  monthly <- seasonal_adjust(monthly, "overdue_debt", 12L)
  month_to_quarter(monthly, "overdue_debt")
}

build_p_cm <- function() {
  old_sheet <- read_sheet_frame(file_paths$p_cm_old, sheet = "Данные")
  new_sheet <- read_sheet_frame(file_paths$p_cm_new, sheet = "Данные")

  old_index <- bind_rows(
    extract_year_block(old_sheet, 8, "O", 2014),
    extract_year_block(old_sheet, 8, "AA", 2015),
    extract_year_block(old_sheet, 8, "AM", 2016)
  ) %>%
    transmute(date, index = value)

  new_index_pre2023 <- bind_rows(
    extract_year_block(new_sheet, 6, "D", 2017, step = 1),
    extract_year_block(new_sheet, 6, "P", 2018, step = 1),
    extract_year_block(new_sheet, 6, "AB", 2019, step = 1),
    extract_year_block(new_sheet, 6, "AN", 2020, step = 1),
    extract_year_block(new_sheet, 6, "AZ", 2021, step = 1),
    extract_year_block(new_sheet, 6, "BL", 2022, step = 1)
  ) %>%
    transmute(date, index = value)

  new_index_post2023 <- extract_years_until_empty(
    new_sheet,
    row_id = 7,
    start_col = "BX",
    first_year = 2023,
    step = 1
  ) %>%
    transmute(date, index = value)

  monthly_index <- bind_rows(
    old_index,
    new_index_pre2023,
    new_index_post2023
  ) %>%
    filter(!is.na(date), !is.na(index)) %>%
    arrange(date)

  monthly_level <- build_level_from_indices(monthly_index)
  monthly_level <- seasonal_adjust(monthly_level, "level", 12L)
  quarter_level <- month_to_quarter(monthly_level, "level")
  names(quarter_level)[2] <- "p_cm"
  quarter_level
}

build_component_monthly_index <- function(sheet_df, row_id, year_col_map, step = 1L) {
  pieces <- lapply(names(year_col_map), function(year_name) {
    extract_year_block(
      sheet_df = sheet_df,
      row_id = row_id,
      start_col = unname(year_col_map[[year_name]]),
      year_value = as.integer(year_name),
      step = step
    )
  })
  bind_rows(pieces) %>% transmute(date, index = value)
}

build_component_monthly_index_until_empty <- function(sheet_df, row_id, start_col, first_year, step = 3L) {
  extract_years_until_empty(
    sheet_df = sheet_df,
    row_id = row_id,
    start_col = start_col,
    first_year = first_year,
    step = step
  ) %>%
    transmute(date, index = value)
}

find_code_row <- function(sheet_df, code_col, code_value) {
  code_values <- stringr::str_trim(as.character(sheet_df[[code_col]]))
  match(as.character(code_value), code_values)
}

build_p_im <- function() {
  old_sheet <- read_sheet_frame(file_paths$p_im_old, sheet = "Данные")
  new_sheet <- read_sheet_frame(file_paths$p_im_new, sheet = "Данные")
  weight_old <- read_sheet_frame(file_paths$p_im_weights, sheet = "2")
  weight_new <- read_sheet_frame(file_paths$p_im_weights, sheet = "4")

  old_year_map <- c("2014" = "O", "2015" = "AA", "2016" = "AM")
  new_year_map <- c(
  "2017" = "D",
  "2018" = "P",
  "2019" = "AB",
  "2020" = "AN",
  "2021" = "AZ",
  "2022" = "BL"
  )

  component_series <- bind_rows(
    build_component_monthly_index(old_sheet, 8, old_year_map) %>% mutate(component = "DF"),
    build_component_monthly_index(old_sheet, 9, old_year_map) %>% mutate(component = "DG"),
    build_component_monthly_index(old_sheet, 10, old_year_map) %>% mutate(component = "DJ"),

    build_component_monthly_index(new_sheet, 6, new_year_map, step = 1) %>% mutate(component = "19"),
    build_component_monthly_index(new_sheet, 7, new_year_map, step = 1) %>% mutate(component = "20"),
    build_component_monthly_index(new_sheet, 8, new_year_map, step = 1) %>% mutate(component = "24"),
    build_component_monthly_index(new_sheet, 9, new_year_map, step = 1) %>% mutate(component = "25"),

    build_component_monthly_index_until_empty(new_sheet, 10, "BX", 2023, step = 1) %>% mutate(component = "19"),
    build_component_monthly_index_until_empty(new_sheet, 11, "BX", 2023, step = 1) %>% mutate(component = "20"),
    build_component_monthly_index_until_empty(new_sheet, 12, "BX", 2023, step = 1) %>% mutate(component = "24"),
    build_component_monthly_index_until_empty(new_sheet, 13, "BX", 2023, step = 1) %>% mutate(component = "25")
  ) %>%
    filter(!is.na(date), !is.na(index))

  weight_old_rows <- c(
    DF = find_code_row(weight_old, 2, "DF"),
    DG = find_code_row(weight_old, 2, "DG"),
    DJ = find_code_row(weight_old, 2, "DJ")
  )
  weight_new_rows <- c(
    `19` = find_code_row(weight_new, 2, "19"),
    `20` = find_code_row(weight_new, 2, "20"),
    `24` = find_code_row(weight_new, 2, "24"),
    `25` = find_code_row(weight_new, 2, "25")
  )
  if (any(is.na(weight_old_rows)) || any(is.na(weight_new_rows))) {
    stop_with_diagnostics("Failed to locate one or more p_im weight rows by industry code.")
  }

  weight_old_df <- bind_rows(
    data.frame(year = 2014:2016, component = "DF", weight = coerce_numeric(unlist(weight_old[weight_old_rows["DF"], 5:7], use.names = FALSE))),
    data.frame(year = 2014:2016, component = "DG", weight = coerce_numeric(unlist(weight_old[weight_old_rows["DG"], 5:7], use.names = FALSE))),
    data.frame(year = 2014:2016, component = "DJ", weight = coerce_numeric(unlist(weight_old[weight_old_rows["DJ"], 5:7], use.names = FALSE)))
  )

  weight_new_years <- 2017:2026
  weight_new_df <- bind_rows(
    data.frame(year = weight_new_years, component = "19", weight = coerce_numeric(unlist(weight_new[weight_new_rows["19"], 3:12], use.names = FALSE))),
    data.frame(year = weight_new_years, component = "20", weight = coerce_numeric(unlist(weight_new[weight_new_rows["20"], 3:12], use.names = FALSE))),
    data.frame(year = weight_new_years, component = "24", weight = coerce_numeric(unlist(weight_new[weight_new_rows["24"], 3:12], use.names = FALSE))),
    data.frame(year = weight_new_years, component = "25", weight = coerce_numeric(unlist(weight_new[weight_new_rows["25"], 3:12], use.names = FALSE)))
  )

  monthly_index <- component_series %>%
    mutate(year = lubridate::year(date)) %>%
    left_join(bind_rows(weight_old_df, weight_new_df), by = c("year", "component")) %>%
    group_by(date) %>%
    filter(!is.na(weight), !is.na(index)) %>%
    mutate(weight = weight / sum(weight)) %>%
    summarise(index = sum(weight * index), .groups = "drop") %>%
    arrange(date)

  monthly_level <- build_level_from_indices(monthly_index)
  monthly_level <- seasonal_adjust(monthly_level, "level", 12L)
  quarter_level <- month_to_quarter(monthly_level, "level")
  names(quarter_level)[2] <- "p_im"
  quarter_level
}

build_cpi <- function() {
  sheet_df <- read_sheet_frame(file_paths$cpi, sheet = "Данные")
  if (ncol(sheet_df) < 5L || nrow(sheet_df) < 9L) {
    stop_with_diagnostics("Unexpected layout in cpi_q.xlsx: too few rows or columns.")
  }

  quarter_labels <- stringr::str_trim(as.character(sheet_df[[2]]))
  metric_labels <- stringr::str_trim(as.character(sheet_df[[3]]))
  year_values <- coerce_numeric(sheet_df[[1]])
  quarter_map <- c("I квартал" = 1L, "II квартал" = 2L, "III квартал" = 3L, "IV квартал" = 4L)
  quarter_values <- unname(quarter_map[quarter_labels])

  quarter_rows <- which(
    !is.na(year_values) &
      !is.na(quarter_values) &
      metric_labels == "К предыдущему кварталу"
  )
  if (length(quarter_rows) == 0L) {
    stop_with_diagnostics("Failed to detect quarterly CPI rows in cpi_q.xlsx.")
  }

  cpi_index <- data.frame(
    date = quarter_end_date(year_values[quarter_rows], quarter_values[quarter_rows]),
    index = ifelse(
      year_values[quarter_rows] <= 2022,
      coerce_numeric(sheet_df[[4]][quarter_rows]),
      coerce_numeric(sheet_df[[5]][quarter_rows])
    ),
    stringsAsFactors = FALSE
  )
  cpi_index <- sort_unique_dates(cpi_index, "index")
  assert_quarterly_regular(cpi_index, "index")

  base_date <- as.Date("2013-12-31")
  if (!base_date %in% cpi_index$date) {
    stop_with_diagnostics("CPI base quarter 2013Q4 is missing in cpi_q.xlsx.")
  }

  chained_rows <- cpi_index[cpi_index$date > base_date, , drop = FALSE]
  levels <- numeric(nrow(chained_rows))
  current_level <- 100
  for (i in seq_len(nrow(chained_rows))) {
    current_level <- current_level * chained_rows$index[i] / 100
    levels[i] <- current_level
  }

  quarter_level <- data.frame(
    date = chained_rows$date,
    cpi = levels,
    stringsAsFactors = FALSE
  )
  quarter_level <- seasonal_adjust(quarter_level, "cpi", 4L)
  assert_quarterly_regular(quarter_level, "cpi")
}

build_wage <- function(cpi_quarterly) {
  sheet_name <- readxl::excel_sheets(file_paths$wage)
  sheet_name <- sheet_name[stringr::str_trim(sheet_name) == "4.2"]
  if (length(sheet_name) != 1L) {
    stop_with_diagnostics("Failed to locate sheet 4.2 in ind_04-2026.xlsx.")
  }

  wage_block <- read_sheet_frame(
    file_paths$wage,
    sheet = sheet_name,
    range = "A7:F34"
  )

  year_values <- coerce_numeric(wage_block[[1]])
  quarterly_values <- wage_block[, 3:6, drop = FALSE]
  nominal_quarterly <- bind_rows(lapply(seq_along(year_values), function(i) {
    data.frame(
      date = quarter_end_date(year_values[i], 1:4),
      wage_nominal = coerce_numeric(unlist(quarterly_values[i, ], use.names = FALSE)),
      stringsAsFactors = FALSE
    )
  }))
  nominal_quarterly <- sort_unique_dates(nominal_quarterly, "wage_nominal")
  nominal_quarterly <- seasonal_adjust(nominal_quarterly, "wage_nominal", 4L)

  merged <- nominal_quarterly %>%
    inner_join(cpi_quarterly, by = "date")

  # CPI is chained to 2013M12 = 100, so division by CPI/100 gives real wages in base-period prices.
  real_quarterly <- merged %>%
    transmute(
      date = date,
      wage = wage_nominal / (cpi / 100)
    )
  assert_quarterly_regular(real_quarterly, "wage")
}

# Сборка итогового набора
build_all_series <- function() {
  cpi_quarterly <- build_cpi()
  series_list <- list(
    gdp = build_gdp(),
    cpi = cpi_quarterly,
    wage = build_wage(cpi_quarterly),
    m2 = build_m2(),
    key_rate = build_key_rate(),
    p_im = build_p_im(),
    p_cm = build_p_cm(),
    overdue_debt = build_overdue_debt()
  )
  names(series_list) <- series_names
  series_list
}

# Запись результатов
write_outputs <- function(series_list) {
  final_series_list <- lapply(names(series_list), function(name) {
    final_df <- log_if_needed(series_list[[name]], name)
    write_series_csv(final_df, name)
    final_df
  })
  names(final_series_list) <- names(series_list)

  merged <- Reduce(
    function(left_df, right_df) merge(left_df, right_df, by = "date", all = FALSE),
    final_series_list
  )
  merged <- merged[stats::complete.cases(merged), , drop = FALSE]
  if (nrow(merged) == 0L) {
    stop_with_diagnostics("No common quarterly sample across the final SVAR variables.")
  }
  assert_quarterly_regular(merged[, c("date", "gdp"), drop = FALSE], "gdp")
  utils::write.csv(
    merged,
    file = file.path(processed_dir, "svar_dataset.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

# Точка входа
main <- function() {
  ensure_processed_dir()
  series_list <- build_all_series()
  write_outputs(series_list)
  if (file.exists(diagnostics_path)) {
    file.remove(diagnostics_path)
  }
  message("dataBuilder.R completed successfully.")
}

if (!identical(tolower(Sys.getenv("DATA_BUILDER_SKIP_MAIN")), "true")) {
  main()
}
