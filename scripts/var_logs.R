.libPaths(c(normalizePath("data/processed/r_lib", winslash = "/", mustWork = FALSE), .libPaths()))

library(knitr)
library(vars)

raw_df <- read.csv("data/processed/svar_dataset.csv")

if ("date" %in% names(raw_df)) {
  raw_df$date <- as.Date(raw_df$date)
}

cols <- c(
  "p_cm",
  "m2",
  "key_rate",
  "p_im",
  "cpi",
  "wage",
  "gdp",
  "overdue_debt"
)

df <- raw_df[, cols]

criteria_names <- c("AIC(n)", "HQ(n)", "SC(n)", "FPE(n)")

pick_final_lag <- function(selection) {
  votes <- table(selection[criteria_names])
  winners <- as.integer(names(votes)[votes == max(votes)])

  # ponytail: if criteria tie, keep the smaller lag as the conservative default.
  min(winners)
}

selection <- VARselect(df, lag.max = 4, type = "const")

criteria_table <- as.data.frame(t(selection$criteria))
criteria_table$lag <- seq_len(nrow(criteria_table))
criteria_table <- criteria_table[, c("lag", criteria_names)]

selected_lags <- selection$selection[criteria_names]
final_lag <- pick_final_lag(selected_lags)

dir.create("tables/output", showWarnings = FALSE, recursive = TRUE)

criteria_table_md <- data.frame(
  lag = as.character(criteria_table$lag),
  AIC = sprintf("%.4f", criteria_table[["AIC(n)"]]),
  HQ = sprintf("%.4f", criteria_table[["HQ(n)"]]),
  SC = sprintf("%.4f", criteria_table[["SC(n)"]]),
  FPE = format(signif(criteria_table[["FPE(n)"]], 4), scientific = TRUE, trim = TRUE),
  stringsAsFactors = FALSE
)

criteria_table_tex <- criteria_table_md

criteria_table_md[criteria_table$lag == final_lag, ] <- lapply(
  criteria_table_md[criteria_table$lag == final_lag, , drop = FALSE],
  function(x) paste0("**", x, "**")
)

criteria_table_tex[criteria_table$lag == final_lag, ] <- lapply(
  criteria_table_tex[criteria_table$lag == final_lag, , drop = FALSE],
  function(x) paste0("\\textbf{", x, "}")
)

# Markdown-таблица строится через knitr::kable по аналогии с остальными таблицами проекта.
criteria_md <- knitr::kable(
  criteria_table_md,
  format = "pipe",
  escape = FALSE,
  align = c("r", "r", "r", "r", "r"),
  caption = "Информационные критерии VAR для спецификации const"
)

writeLines(
  criteria_md,
  "tables/output/var_ic_const.md",
  useBytes = TRUE
)

# LaTeX-таблица также строится через knitr::kable; escape = FALSE нужен для жирного выделения лучшего лага.
criteria_tex <- knitr::kable(
  criteria_table_tex,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = c("r", "r", "r", "r", "r"),
  caption = "Информационные критерии VAR для спецификации const",
  label = "tab:var_ic_const"
)

writeLines(
  criteria_tex,
  "tables/output/var_ic_const.tex",
  useBytes = TRUE
)

cat("\n==============================\n")
cat("Specification: const\n")
cat("==============================\n")
print(
  data.frame(
    lag = criteria_table$lag,
    AIC = round(criteria_table[["AIC(n)"]], 4),
    HQ = round(criteria_table[["HQ(n)"]], 4),
    SC = round(criteria_table[["SC(n)"]], 4),
    FPE = signif(criteria_table[["FPE(n)"]], 4)
  ),
  row.names = FALSE
)

cat("\nSelected lag by criterion:\n")
print(selected_lags)
cat("Final lag:", final_lag, "\n")
cat("Recommended lag for svar.R (type = \"const\"):", final_lag, "\n")

# =====================
# appendix: VAR summary tables
# captions are placed below tables
# =====================

dir.create("tables/output", showWarnings = FALSE, recursive = TRUE)

var_symbols <- c(
  p_cm = "P_{cm}",
  m2 = "M2",
  key_rate = "R",
  p_im = "P_{im}",
  cpi = "P_y",
  wage = "W",
  gdp = "y",
  overdue_debt = "OvDeb"
)

var_labels <- paste0("$", var_symbols, "$")
names(var_labels) <- names(var_symbols)

term_labels <- c(
  const = "Константа",
  setNames(
    paste0("${", var_symbols[cols], "}_{t-1}$"),
    paste0(cols, ".l1")
  )
)
model_var_appendix <- VAR(df, p = 1, type = "const")

write_table_bottom_caption <- function(table_df,
                                       path,
                                       caption,
                                       label,
                                       align,
                                       size = "\\small",
                                       resize = FALSE) {
  tabular <- knitr::kable(
    table_df,
    format = "latex",
    booktabs = FALSE,
    escape = FALSE,
    align = align
  )

  if (resize) {
    body <- c(
      "\\resizebox{\\textwidth}{!}{%",
      as.character(tabular),
      "}%"
    )
  } else {
    body <- as.character(tabular)
  }

  out <- c(
    "\\begin{table}[H]",
    "\\centering",
    size,
    body,
    paste0("\\caption{", caption, "}"),
    paste0("\\label{tab:", label, "}"),
    "\\end{table}"
  )

  writeLines(out, path, useBytes = TRUE)
}

# ---- information criteria ----

criteria_table_tex <- data.frame(
  lag = as.character(criteria_table$lag),
  AIC = sprintf("%.4f", criteria_table[["AIC(n)"]]),
  HQ = sprintf("%.4f", criteria_table[["HQ(n)"]]),
  SC = sprintf("%.4f", criteria_table[["SC(n)"]]),
  FPE = format(
    signif(criteria_table[["FPE(n)"]], 4),
    scientific = TRUE,
    trim = TRUE
  ),
  stringsAsFactors = FALSE
)

criteria_table_tex[criteria_table$lag == final_lag, ] <- lapply(
  criteria_table_tex[criteria_table$lag == final_lag, , drop = FALSE],
  function(x) paste0("\\textbf{", x, "}")
)

write_table_bottom_caption(
  table_df = criteria_table_tex,
  path = "tables/output/var_ic_const.tex",
  caption = "Информационные критерии VAR для спецификации с константой",
  label = "var_ic_const",
  align = c("r", "r", "r", "r", "r")
)

# ---- roots / stability ----

roots_table <- data.frame(
  "Номер корня" = seq_along(roots(model_var_appendix)),
  "Модуль" = sprintf("%.4f", roots(model_var_appendix)),
  check.names = FALSE
)

write_table_bottom_caption(
  table_df = roots_table,
  path = "tables/output/var_roots_const.tex",
  caption = "Модули корней приведённой VAR(1)",
  label = "var_roots_const",
  align = c("r", "r")
)

# ---- equation fit summary ----

fit_table <- do.call(
  rbind,
  lapply(names(model_var_appendix$varresult), function(eq_name) {
    eq_summary <- summary(model_var_appendix$varresult[[eq_name]])

    data.frame(
      "Уравнение" = unname(var_labels[eq_name]),
      "$R^2$" = sprintf("%.4f", eq_summary$r.squared),
      "Скорр. $R^2$" = sprintf("%.4f", eq_summary$adj.r.squared),
      "Ст. ошибка регрессии" = sprintf("%.4f", eq_summary$sigma),
      check.names = FALSE
    )
  })
)

write_table_bottom_caption(
  table_df = fit_table,
  path = "tables/output/var_equation_fit_const.tex",
  caption = "Качество подгонки уравнений приведённой VAR(1)",
  label = "var_equation_fit_const",
  align = c("l", "r", "r", "r")
)
# ---- residual correlation matrix ----

resid_cor <- round(cor(resid(model_var_appendix)), 4)

resid_cor_table <- as.data.frame(resid_cor, check.names = FALSE)

# Явно заменяем технические имена переменных на LaTeX-обозначения
rownames(resid_cor_table) <- unname(var_labels[rownames(resid_cor_table)])
colnames(resid_cor_table) <- unname(var_labels[colnames(resid_cor_table)])

resid_cor_table <- data.frame(
  "Уравнение" = rownames(resid_cor_table),
  resid_cor_table,
  check.names = FALSE,
  row.names = NULL
)

write_table_bottom_caption(
  table_df = resid_cor_table,
  path = "tables/output/var_residual_cor_const.tex",
  caption = "Корреляционная матрица остатков приведённой VAR(1)",
  label = "var_residual_cor_const",
  align = c("l", rep("r", length(cols))),
  size = "\\scriptsize",
  resize = TRUE
)

# ---- full coefficient summary by equations ----

coef_tables_tex <- c()

for (eq_name in names(model_var_appendix$varresult)) {
  coef_matrix <- coef(summary(model_var_appendix$varresult[[eq_name]]))

  terms_raw <- rownames(coef_matrix)
  terms_print <- ifelse(
    terms_raw %in% names(term_labels),
    unname(term_labels[terms_raw]),
    terms_raw
  )

  coef_table <- data.frame(
    "Регрессор" = terms_print,
    "Оценка" = sprintf("%.4f", coef_matrix[, "Estimate"]),
    "Ст. ошибка" = sprintf("%.4f", coef_matrix[, "Std. Error"]),
    "$t$-стат." = sprintf("%.3f", coef_matrix[, "t value"]),
    "$p$-знач." = sprintf("%.3f", coef_matrix[, "Pr(>|t|)"]),
    check.names = FALSE
  )

  tabular <- knitr::kable(
    coef_table,
    format = "latex",
    booktabs = FALSE,
    escape = FALSE,
    align = c("l", "r", "r", "r", "r")
  )

  coef_tables_tex <- c(
    coef_tables_tex,
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    as.character(tabular),
    paste0(
      "\\caption{Оценки коэффициентов уравнения ",
      unname(var_labels[eq_name]),
      " в приведённой VAR(1)}"
    ),
    paste0("\\label{tab:var_coef_", eq_name, "}"),
    "\\end{table}",
    ""
  )
}

writeLines(
  coef_tables_tex,
  "tables/output/var_coefficients_const.tex",
  useBytes = TRUE
)