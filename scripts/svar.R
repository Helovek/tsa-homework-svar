# install.packages(c('vars', 'svars', 'knitr', 'kableExtra'))

.libPaths(c(normalizePath("data/processed/r_lib", winslash = "/", mustWork = FALSE), .libPaths()))

library(ggplot2)
library(dplyr)
library(tibble)
library(vars)
library(knitr)

raw_df <- read.csv('data/processed/svar_dataset.csv')

if ("date" %in% names(raw_df)) {
  raw_df$date <- as.Date(raw_df$date)
}

cols <- c(
  'p_cm',
  'm2',
  'key_rate',
  'p_im',
  'cpi',
  'wage',
  'gdp', 
  'overdue_debt'
)

var_names  <- c(
  "p_cm",
  "m2",
  "key_rate",
  "p_im",
  "cpi",
  "wage",
  "gdp",
  "overdue_debt"
)

eq_names <- c(
  "Pcm",
  "MD",
  "MS",
  "Pim",
  "Py",
  "W",
  "y",
  "OvDeb"
)

df <- raw_df[,cols]

# =====================
# base models
# =====================
model_none  <- VAR(df, p = 1, type = "none")
model_const <- VAR(df, p = 1, type = "const")
model_both  <- VAR(df, p = 1, type = "both")

roots(model_none)
roots(model_const)
roots(model_both)

model <- model_const
u <- resid(model)

# =====================
# reduced VAR
# =====================
make_A <- function(b) {
  A <- matrix(0, 8, 8)

  colnames(A) <- var_names
  rownames(A) <- eq_names

  i <- 1

  # P_cm bounds
  A['Pcm',] <- b[i:(i+7)]
  i <- i + 8
  # MD bounds
  a_md <- b[i]
  b_md <- b[(i+1)]
  i <- i + 2
  A['MD', ] <- c(0, a_md, b_md, 0, -a_md, 0, -a_md, 0)
  # MS bounds
  A['MS', c('p_cm', 'm2', 'key_rate', 'overdue_debt')] <- b[i:(i+3)]
  i <- i + 4
  # slow sector
  # P_im
  A['Pim', c('p_cm', 'p_im', 'cpi', 'wage', 'gdp', 'overdue_debt')] <- b[i:(i+5)]
  i <- i + 6
  # P_y bounds
  A['Py', c('p_cm', 'cpi', 'wage', 'gdp', 'overdue_debt')] <- b[i:(i+4)]
  i <- i + 5

  # wage bounds
  A['W', c('p_cm', 'wage', 'gdp', 'overdue_debt')] <- b[i:(i+3)]
  i <- i + 4

  # y bounds
  A['y', c('p_cm', 'gdp', 'overdue_debt')] <- b[i:(i+2)]
  i <- i + 3

  # OvDeb bounds
  A['OvDeb', c('p_cm', 'overdue_debt')] <- b[i:(i+1)]
  
  return(A)
}

likelihood <- function(b, u) {
  A <- make_A(b)

  T <- nrow(u)
  k <- ncol(u)
  det_A <- determinant(A, logarithm = TRUE)
  eps <- t(A %*% t(u))
  eps_quad <- sum(eps^2)

  if (!is.finite(det_A$modulus) || det_A$sign == 0) {
    return(1e12)
  }
  if (any(!is.finite(eps))) {
    return(1e12)
  }

  loglik <-
    T * as.numeric(det_A$modulus) -
    0.5 * eps_quad -
    0.5 * T * k * log(2 * pi)

  return(-loglik/T)
} 


n_par <- 34
b0 <- rep(0.05, n_par)

n_par <- 34
b0 <- rep(0.05, n_par)

b0[1]  <- 1   
b0[9]  <- 1   
b0[13] <- 1   
b0[16] <- 1   
b0[22] <- 1   
b0[27] <- 1   
b0[31] <- 1   
b0[34] <- 1   

fit_nm <- optim(
  par = b0,
  fn = likelihood,
  u = u,
  method = "Nelder-Mead",
  control = list(
    maxit = 50000,
    reltol = 1e-10
  )
)

fit_bfgs <- optim(
  par = fit_nm$par,
  fn = likelihood,
  u = u,
  method = "BFGS",
  control = list(
    maxit = 10000,
    reltol = 1e-12
  ),
  hessian = TRUE
)

fit <- fit_bfgs

fit$convergence
fit$value

A_hat <- make_A(fit$par)

# A_hat matrix
round(A_hat, 4)

Sigma_u <- crossprod(u) / nrow(u)
eps_hat <- t(A_hat %*% t(u))

# =====================
# output directories
# =====================
dir.create("tables/output", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/output", showWarnings = FALSE, recursive = TRUE)

format_quarter <- function(date_value) {
  date_parts <- as.POSIXlt(date_value)
  sprintf("%dQ%d", date_parts$year + 1900, date_parts$mon %/% 3 + 1)
}

format_table_number <- function(x) {
  ifelse(is.finite(x), sprintf("%.3f", round(x, 3)), "Inf")
}

write_tabular_only <- function(df_input, file_path, align = NULL, escape = TRUE) {
  tabular_text <- knitr::kable(
    df_input,
    format = "latex",
    booktabs = TRUE,
    escape = escape,
    align = align
  )

  tabular_text <- sub("^\\s*\\n+", "", tabular_text)

  writeLines(tabular_text, file_path, useBytes = TRUE)
}

sample_period_label <- if ("date" %in% names(raw_df)) {
  paste0(
    format_quarter(min(raw_df$date, na.rm = TRUE)),
    "--",
    format_quarter(max(raw_df$date, na.rm = TRUE))
  )
} else {
  paste0("1--", nrow(raw_df))
}

var_specification_df <- data.frame(
  item = c(
    "frequency",
    "sample period",
    "observations after lags",
    "lag order",
    "variables",
    "transformations",
    "model type"
  ),
  value = c(
    "quarterly",
    sample_period_label,
    as.character(nrow(u)),
    "p = 1",
    "Pcm, M2, R, Pim, Py, W, y, OvDeb",
    "log-levels for level variables; R without log",
    "VAR(1) / SVAR(1)"
  ),
  stringsAsFactors = FALSE
)

write_tabular_only(
  var_specification_df,
  "tables/output/var_specification.tex",
  align = c("l", "l")
)

var_selection <- vars::VARselect(
  df,
  lag.max = 8,
  type = model$type
)

var_lag_selection_df <- data.frame(
  lag = seq_len(ncol(var_selection$criteria)),
  AIC = format_table_number(var_selection$criteria["AIC(n)", ]),
  HQ = format_table_number(var_selection$criteria["HQ(n)", ]),
  SC = format_table_number(var_selection$criteria["SC(n)", ]),
  FPE = format_table_number(var_selection$criteria["FPE(n)", ]),
  stringsAsFactors = FALSE
)

write_tabular_only(
  var_lag_selection_df,
  "tables/output/var_lag_selection.tex",
  align = c("r", "r", "r", "r", "r")
)

message(
  "Lag selection note: final specification remains p = 1 ",
  "due to compactness and SVAR stability on a short quarterly sample."
)

var_roots_modulus <- as.numeric(roots(model, modulus = TRUE))
var_stable <- all(var_roots_modulus < 1)

var_roots_df <- data.frame(
  root_id = c(as.character(seq_along(var_roots_modulus)), "stable"),
  modulus = c(
    format_table_number(var_roots_modulus),
    if (var_stable) "TRUE" else "FALSE"
  ),
  stringsAsFactors = FALSE
)

write_tabular_only(
  var_roots_df,
  "tables/output/var_roots.tex",
  align = c("l", "r")
)

serial_diag <- try(
  vars::serial.test(model, lags.pt = 8, type = "PT.asymptotic"),
  silent = TRUE
)
arch_diag <- try(
  vars::arch.test(model, lags.multi = 5),
  silent = TRUE
)
normality_diag <- try(
  vars::normality.test(model),
  silent = TRUE
)

if (
  !inherits(serial_diag, "try-error") &&
  !inherits(arch_diag, "try-error") &&
  !inherits(normality_diag, "try-error")
) {
  var_residual_diagnostics_df <- data.frame(
    test = c("Portmanteau", "ARCH", "Normality (JB)"),
    statistic = format_table_number(c(
      as.numeric(serial_diag$serial$statistic),
      as.numeric(arch_diag$arch.mul$statistic),
      as.numeric(normality_diag$jb.mul$JB$statistic)
    )),
    p_value = format_table_number(c(
      as.numeric(serial_diag$serial$p.value),
      as.numeric(arch_diag$arch.mul$p.value),
      as.numeric(normality_diag$jb.mul$JB$p.value)
    )),
    stringsAsFactors = FALSE
  )

  write_tabular_only(
    var_residual_diagnostics_df,
    "tables/output/var_residual_diagnostics.tex",
    align = c("l", "r", "r")
  )
}

# Подписи для русскоязычной работы.
# В таблицах оставляем математические обозначения, а на графиках — русские названия шоков.
var_labels <- c(
  p_cm = "$P_{cm}$",
  m2 = "$M2$",
  key_rate = "$R$",
  p_im = "$P_{im}$",
  cpi = "$P_y$",
  wage = "$W$",
  gdp = "$y$",
  overdue_debt = "$OvDeb$"
)

eq_labels <- c(
  Pcm = "$P_{cm}$",
  MD = "$MD$",
  MS = "$MS$",
  Pim = "$P_{im}$",
  Py = "$P_y$",
  W = "$W$",
  y = "$y$",
  OvDeb = "$OvDeb$"
)

shock_labels_ru <- c(
  Pcm = "Сырьевые цены",
  MD = "Денежный спрос",
  MS = "Денежная политика",
  Pim = "Промежуточные цены",
  Py = "Потребительские цены",
  W = "Заработная плата",
  y = "Выпуск",
  OvDeb = "Просроченная задолженность"
)

# =====================
# block 1: A_hat table
# =====================
A_hat_table <- data.frame(
  "Уравнение" = eq_labels[rownames(A_hat)],
  round(A_hat, 4),
  check.names = FALSE,
  row.names = NULL
)

A_hat_tabular <- knitr::kable(
  A_hat_table,
  format = "latex",
  booktabs = FALSE,
  escape = FALSE,
  align = c("l", rep("r", length(var_names)))
)

A_hat_md <- knitr::kable(
  A_hat_table,
  format = "pipe",
  escape = FALSE,
  align = c("l", rep("r", length(var_names))),
  caption = "Оценённая матрица текущих структурных связей $\\hat A$"
)

A_hat_tex <- paste(
  "\\begin{table}[H]",
  "\\centering",
  "\\scriptsize",
  "\\resizebox{\\textwidth}{!}{%",
  as.character(A_hat_tabular),
  "}%",
  "\\caption{Оценённая матрица текущих структурных связей $\\hat A$}",
  "\\label{tab:a_hat}",
  "\\end{table}",
  sep = "\n"
)

writeLines(
  A_hat_tex,
  "tables/output/A_hat.tex",
  useBytes = TRUE
)

writeLines(
  A_hat_md,
  "tables/output/A_hat.md",
  useBytes = TRUE
)

# =====================
# block 2: structural shocks
# =====================
colnames(eps_hat) <- eq_names

if ("date" %in% names(raw_df)) {
  shock_dates <- tail(raw_df$date, nrow(eps_hat))
} else {
  shock_dates <- seq_len(nrow(eps_hat))
}

eps_df <- as.data.frame(eps_hat) %>%
  dplyr::mutate(date = shock_dates) %>%
  dplyr::select(date, dplyr::everything())

write.csv(
  eps_df,
  "data/processed/structural_shocks.csv",
  row.names = FALSE
)

# Диагностические таблицы по структурным шокам.
write.csv(
  round(cov(eps_hat), 4),
  "tables/output/structural_shocks_cov.csv",
  row.names = TRUE
)

write.csv(
  round(cor(eps_hat), 4),
  "tables/output/structural_shocks_cor.csv",
  row.names = TRUE
)

write.csv(
  data.frame(
    shock = names(apply(eps_hat, 2, sd)),
    sd = as.numeric(apply(eps_hat, 2, sd))
  ),
  "tables/output/structural_shocks_sd.csv",
  row.names = FALSE
)

# =====================
# block 3: impulse response functions
# =====================

# Impact matrix
impact_mat <- solve(A_hat)

rownames(impact_mat) <- var_names
colnames(impact_mat) <- eq_names


# VAR lag matrices
get_var_lag_matrices <- function(model, var_names) {
  p_lags <- model$p
  k <- length(var_names)

  Phi <- array(
    NA_real_,
    dim = c(k, k, p_lags),
    dimnames = list(
      response = var_names,
      lagged_var = var_names,
      lag = paste0("lag_", 1:p_lags)
    )
  )

  model_coefs <- coef(model)

  for (lag_id in 1:p_lags) {
    for (response in var_names) {
      response_coefs <- model_coefs[[response]][, "Estimate"]

      for (lagged_var in var_names) {
        coef_name <- paste0(lagged_var, ".l", lag_id)

        if (!coef_name %in% names(response_coefs)) {
          stop(paste("Coefficient not found:", coef_name))
        }

        Phi[response, lagged_var, lag_id] <- response_coefs[coef_name]
      }
    }
  }

  return(Phi)
}

compute_reduced_multipliers <- function(Phi_input, horizon) {
  k_local <- dim(Phi_input)[1]
  p_lags <- dim(Phi_input)[3]

  Psi_local <- array(
    0,
    dim = c(horizon + 1, k_local, k_local),
    dimnames = list(
      horizon = 0:horizon,
      response = dimnames(Phi_input)$response,
      innovation = dimnames(Phi_input)$lagged_var
    )
  )

  Psi_local[1, , ] <- diag(k_local)

  for (h in 1:horizon) {
    Psi_h <- matrix(0, nrow = k_local, ncol = k_local)

    for (lag_id in 1:min(h, p_lags)) {
      Psi_h <- Psi_h + Phi_input[, , lag_id] %*% Psi_local[h - lag_id + 1, , ]
    }

    Psi_local[h + 1, , ] <- Psi_h
  }

  Psi_local
}

compute_structural_irf <- function(Phi_input, impact_input, horizon, response_names, shock_names) {
  Psi_local <- compute_reduced_multipliers(Phi_input, horizon)

  irf_local <- array(
    NA_real_,
    dim = c(horizon + 1, length(response_names), length(shock_names)),
    dimnames = list(
      horizon = 0:horizon,
      response = response_names,
      shock = shock_names
    )
  )

  for (h in 0:horizon) {
    irf_local[h + 1, , ] <- Psi_local[h + 1, , ] %*% impact_input
  }

  list(
    Psi = Psi_local,
    irf = irf_local
  )
}

Phi <- get_var_lag_matrices(model, var_names)

# Reduced-form multipliers
H <- 24
k <- length(var_names)

irf_main <- compute_structural_irf(
  Phi_input = Phi,
  impact_input = impact_mat,
  horizon = H,
  response_names = var_names,
  shock_names = eq_names
)

Psi <- irf_main$Psi
irf_array <- irf_main$irf

# Bootstrap confidence intervals
set.seed(123)

B <- 200
alpha <- 0.10

p_lags <- model$p
model_type <- model$type

u_centered <- scale(
  resid(model),
  center = TRUE,
  scale = FALSE
)

model_coefs <- coef(model)
df_mat <- as.matrix(df)

make_boot_sample <- function() {
  y_boot <- df_mat
  boot_resid <- u_centered[
    sample(seq_len(nrow(u_centered)), size = nrow(u_centered), replace = TRUE),
    ,
    drop = FALSE
  ]

  for (t in (p_lags + 1):nrow(y_boot)) {
    for (response in var_names) {
      response_coefs <- model_coefs[[response]][, "Estimate"]
      y_hat <- 0

      if ("const" %in% names(response_coefs)) {
        y_hat <- y_hat + response_coefs["const"]
      }

      if ("trend" %in% names(response_coefs)) {
        y_hat <- y_hat + response_coefs["trend"] * t
      }

      for (lag_id in 1:p_lags) {
        for (lagged_var in var_names) {
          coef_name <- paste0(lagged_var, ".l", lag_id)
          y_hat <- y_hat +
            response_coefs[coef_name] *
            y_boot[t - lag_id, lagged_var]
        }
      }

      y_boot[t, response] <- y_hat + boot_resid[t - p_lags, response]
    }
  }

  as.data.frame(y_boot)
}

compute_irf_from_model <- function(model_boot) {
  Phi_boot <- get_var_lag_matrices(model_boot, var_names)
  compute_structural_irf(
    Phi_input = Phi_boot,
    impact_input = impact_mat,
    horizon = H,
    response_names = var_names,
    shock_names = eq_names
  )$irf
}

irf_boot_array <- array(
  NA_real_,
  dim = c(H + 1, k, length(eq_names), B),
  dimnames = list(
    horizon = 0:H,
    response = var_names,
    shock = eq_names,
    bootstrap = 1:B
  )
)

for (b in 1:B) {
  boot_df <- make_boot_sample()

  model_boot <- try(
    VAR(
      boot_df,
      p = p_lags,
      type = model_type
    ),
    silent = TRUE
  )

  if (!inherits(model_boot, "try-error")) {
    irf_boot_array[, , , b] <- compute_irf_from_model(model_boot)
  }
}

irf_lower_array <- apply(
  irf_boot_array,
  c(1, 2, 3),
  quantile,
  probs = alpha / 2,
  na.rm = TRUE
)

irf_upper_array <- apply(
  irf_boot_array,
  c(1, 2, 3),
  quantile,
  probs = 1 - alpha / 2,
  na.rm = TRUE
)

# Long format
irf_df <- expand.grid(
  horizon = 0:H,
  response = var_names,
  shock = eq_names
) %>%
  dplyr::mutate(
    irf = as.vector(irf_array),
    irf_lower = as.vector(irf_lower_array),
    irf_upper = as.vector(irf_upper_array)
  )

write.csv(
  irf_df,
  "data/processed/irf_all_shocks.csv",
  row.names = FALSE
)

# Plot labels
var_plot_labels <- c(
  p_cm = "P[cm]",
  m2 = "M2",
  key_rate = "R",
  p_im = "P[im]",
  cpi = "P[y]",
  wage = "W",
  gdp = "y",
  overdue_debt = "OvDeb"
)

eq_plot_labels <- c(
  Pcm = "P[cm]",
  MD = "MD",
  MS = "MS",
  Pim = "P[im]",
  Py = "P[y]",
  W = "W",
  y = "y",
  OvDeb = "OvDeb"
)

irf_full_df <- irf_df %>%
  dplyr::mutate(
    response_plot = factor(
      response,
      levels = var_names,
      labels = var_plot_labels[var_names]
    ),
    shock_plot = factor(
      shock,
      levels = eq_names,
      labels = eq_plot_labels[eq_names]
    )
  )

# Full IRF grid
p_irf_full <- ggplot(
  irf_full_df,
  aes(x = horizon, y = irf)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25
  ) +
  geom_line(
    aes(y = irf_lower),
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_line(
    aes(y = irf_upper),
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_line(
    linewidth = 0.35
  ) +
  facet_grid(
    response_plot ~ shock_plot,
    scales = "free_y",
    switch = "y",
    labeller = label_parsed
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),

    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.35
    ),
    panel.grid = element_blank(),

    strip.background = element_blank(),
    strip.placement = "outside",

    strip.text.x = element_text(
      size = 8,
      margin = margin(b = 3)
    ),
    strip.text.y.left = element_text(
      size = 8,
      angle = 0,
      margin = margin(r = 4)
    ),

    panel.spacing.x = unit(0.12, "lines"),
    panel.spacing.y = unit(0.12, "lines"),

    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  "figures/output/irf_full_grid.pdf",
  p_irf_full,
  width = 8.5,
  height = 8.5
)

ggsave(
  "figures/output/irf_full_grid.png",
  p_irf_full,
  width = 8.5,
  height = 8.5,
  dpi = 400
)

p_irf_full_presentation <- ggplot(
  irf_full_df,
  aes(x = horizon, y = irf)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    colour = "grey45"
  ) +
  geom_ribbon(
    aes(ymin = irf_lower, ymax = irf_upper),
    inherit.aes = TRUE,
    fill = "grey70",
    alpha = 0.45
  ) +
  geom_line(
    linewidth = 0.55,
    colour = "black"
  ) +
  facet_grid(
    response_plot ~ shock_plot,
    scales = "free_y",
    switch = "y",
    labeller = label_parsed
  ) +
  scale_x_continuous(
    breaks = c(0, 10, 20),
    limits = c(0, H)
  ) +
  labs(
    x = "Horizon (quarters)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(
      size = 13,
      face = "bold",
      margin = margin(b = 5)
    ),
    strip.text.y.left = element_text(
      size = 13,
      face = "bold",
      angle = 0,
      margin = margin(r = 6)
    ),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12),
    axis.ticks = element_line(linewidth = 0.3),
    panel.spacing.x = unit(0.18, "lines"),
    panel.spacing.y = unit(0.18, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  "figures/output/irf_full_grid_presentation.png",
  p_irf_full_presentation,
  width = 12,
  height = 10,
  dpi = 300
)

irf_compact_response_order <- c("key_rate", "m2", "cpi", "gdp")
irf_compact_shock_order <- c("MS", "MD")

irf_ms_md_compact_df <- irf_df %>%
  dplyr::filter(
    response %in% irf_compact_response_order,
    shock %in% irf_compact_shock_order
  ) %>%
  dplyr::mutate(
    response_plot = factor(
      response,
      levels = irf_compact_response_order,
      labels = var_plot_labels[irf_compact_response_order]
    ),
    shock_plot = factor(
      shock,
      levels = irf_compact_shock_order,
      labels = eq_plot_labels[irf_compact_shock_order]
    )
  )

p_irf_ms_md_compact <- ggplot(
  irf_ms_md_compact_df,
  aes(x = horizon, y = irf)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    colour = "grey45"
  ) +
  geom_ribbon(
    aes(ymin = irf_lower, ymax = irf_upper),
    fill = "grey70",
    alpha = 0.45
  ) +
  geom_line(
    linewidth = 0.8,
    colour = "black"
  ) +
  facet_grid(
    response_plot ~ shock_plot,
    scales = "free_y",
    switch = "y",
    labeller = label_parsed
  ) +
  scale_x_continuous(
    breaks = c(0, 4, 8, 12, 16, 20, 24),
    limits = c(0, H)
  ) +
  labs(
    x = "Horizon (quarters)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(
      size = 14,
      face = "bold",
      margin = margin(b = 5)
    ),
    strip.text.y.left = element_text(
      size = 14,
      face = "bold",
      angle = 0,
      margin = margin(r = 7)
    ),
    axis.text = element_text(size = 11, colour = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.ticks = element_line(linewidth = 0.3),
    panel.spacing.x = unit(0.35, "lines"),
    panel.spacing.y = unit(0.35, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  "figures/output/irf_ms_md_compact.png",
  p_irf_ms_md_compact,
  width = 13,
  height = 6,
  dpi = 300
)

# Normalized MS shock
ms_key_rate_impact <- impact_mat["key_rate", "MS"]

irf_ms_scaled <- irf_df %>%
  dplyr::filter(shock == "MS") %>%
  dplyr::mutate(
    irf_scaled = irf / ms_key_rate_impact,
    irf_scaled_lower = irf_lower / ms_key_rate_impact,
    irf_scaled_upper = irf_upper / ms_key_rate_impact
  )

write.csv(
  irf_ms_scaled,
  "data/processed/irf_ms_scaled_key_rate_1pp.csv",
  row.names = FALSE
)

# =====================
# block 4: historical decomposition of forecast errors
# =====================

compute_historical_decomposition <- function(irf_input, shocks_input) {
  T_obs <- nrow(shocks_input)
  n_response <- dim(irf_input)[2]
  n_shocks <- dim(irf_input)[3]

  hd_array <- array(
    0,
    dim = c(T_obs, n_response, n_shocks),
    dimnames = list(
      time = seq_len(T_obs),
      response = dimnames(irf_input)$response,
      shock = dimnames(irf_input)$shock
    )
  )

  for (t in seq_len(T_obs)) {
    for (lag_h in 0:(t - 1)) {
      for (shock_id in seq_len(n_shocks)) {
        hd_array[t, , shock_id] <- hd_array[t, , shock_id] +
          irf_input[lag_h + 1, , shock_id] * shocks_input[t - lag_h, shock_id]
      }
    }
  }

  hd_array
}

simulate_zero_shock_path <- function(model, df_input, var_names) {
  y_path <- as.matrix(df_input[, var_names, drop = FALSE])
  model_coefs_local <- coef(model)
  p_lags_local <- model$p

  for (t in (p_lags_local + 1):nrow(y_path)) {
    y_next <- numeric(length(var_names))
    names(y_next) <- var_names

    for (response in var_names) {
      response_coefs <- model_coefs_local[[response]][, "Estimate"]
      y_hat <- 0

      if ("const" %in% names(response_coefs)) {
        y_hat <- y_hat + response_coefs["const"]
      }

      if ("trend" %in% names(response_coefs)) {
        y_hat <- y_hat + response_coefs["trend"] * t
      }

      for (lag_id in 1:p_lags_local) {
        for (lagged_var in var_names) {
          coef_name <- paste0(lagged_var, ".l", lag_id)
          y_hat <- y_hat + response_coefs[coef_name] * y_path[t - lag_id, lagged_var]
        }
      }

      y_next[response] <- y_hat
    }

    y_path[t, ] <- y_next
  }

  y_path
}

H_hd <- nrow(eps_hat) - 1

historical_irf_array <- compute_structural_irf(
  Phi_input = Phi,
  impact_input = impact_mat,
  horizon = H_hd,
  response_names = var_names,
  shock_names = eq_names
)$irf

historical_decomp_array <- compute_historical_decomposition(
  irf_input = historical_irf_array,
  shocks_input = eps_hat
)

historical_forecast_errors <- apply(
  historical_decomp_array,
  c(1, 2),
  sum
)

df_aligned <- tail(df, nrow(eps_hat))
zero_shock_path <- simulate_zero_shock_path(model, df, var_names)
zero_shock_aligned <- tail(zero_shock_path, nrow(eps_hat))

reconstruction_gap <- max(
  abs((as.matrix(df_aligned) - historical_forecast_errors) - zero_shock_aligned),
  na.rm = TRUE
)

if (!is.finite(reconstruction_gap) || reconstruction_gap > 1e-5) {
  warning(
    paste(
      "Historical decomposition check failed.",
      "Maximum reconstruction gap:",
      signif(reconstruction_gap, 4)
    )
  )
}

historical_decomp_df <- tibble::tibble(
  date = rep(shock_dates, times = length(var_names) * length(eq_names)),
  response = rep(rep(var_names, each = length(shock_dates)), times = length(eq_names)),
  shock = rep(eq_names, each = length(shock_dates) * length(var_names)),
  value = as.vector(historical_decomp_array)
)

historical_forecast_error_df <- tibble::tibble(
  date = rep(shock_dates, times = length(var_names)),
  response = rep(var_names, each = length(shock_dates)),
  shock = "ForecastErrors",
  value = as.vector(historical_forecast_errors)
)

historical_output_df <- dplyr::bind_rows(
  historical_forecast_error_df,
  historical_decomp_df
)

write.csv(
  historical_output_df,
  "data/processed/historical_decomposition.csv",
  row.names = FALSE
)

component_panel_labels <- c(
  ForecastErrors = "Forecast~errors",
  Pcm = "P[cm]",
  MD = "MD",
  MS = "MS",
  Pim = "P[im]",
  Py = "P[y]",
  W = "W",
  y = "y",
  OvDeb = "T[bk]"
)

shock_interval_df <- tibble::tibble(
  label = c(
    "2014K4-2015K4",
    "2020K2-2020K4",
    "2022K1-2022K4"
  ),
  start = as.Date(c("2014-10-01", "2020-04-01", "2022-01-01")),
  end = as.Date(c("2015-12-31", "2020-12-31", "2022-12-31"))
) %>%
  dplyr::mutate(
    midpoint = start + (end - start) / 2
  )

make_historical_decomposition_plot <- function(response_name) {
  plot_df <- historical_output_df %>%
    dplyr::filter(response == response_name) %>%
    dplyr::mutate(
      panel = factor(
        shock,
        levels = c("ForecastErrors", eq_names),
        labels = component_panel_labels[c("ForecastErrors", eq_names)]
      )
    )

  y_limit <- max(abs(plot_df$value), na.rm = TRUE)

  if (!is.finite(y_limit) || y_limit == 0) {
    y_limit <- 1
  }

  interval_label_df <- shock_interval_df %>%
    dplyr::mutate(
      panel = factor(
        unname(component_panel_labels["ForecastErrors"]),
        levels = levels(plot_df$panel)
      ),
      y = y_limit * 0.9
    )

  p_hd <- ggplot(
    plot_df,
    aes(x = date, y = value)
  ) +
    geom_rect(
      data = shock_interval_df,
      aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.25
    ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.25
    ) +
    geom_col(
      width = 75,
      fill = "black",
      colour = "black",
      linewidth = 0.1
    ) +
    facet_wrap(
      ~ panel,
      ncol = 1,
      scales = "fixed",
      labeller = label_parsed,
      axes = "all_x",
      axis.labels = "all_x"
    ) +
    coord_cartesian(
      ylim = c(-y_limit, y_limit)
    ) +
    geom_text(
      data = interval_label_df,
      aes(
        x = midpoint,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 2.2,
      vjust = 1,
      family = ""
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.35
      ),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 8,
        face = "bold"
      ),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      axis.ticks = element_line(linewidth = 0.25),
      panel.spacing.y = unit(0.12, "lines"),
      plot.margin = margin(8, 8, 8, 8)
    )

  if (inherits(plot_df$date, "Date")) {
    p_hd <- p_hd +
      scale_x_date(
        date_breaks = "2 years",
        date_labels = "%Y",
        expand = expansion(mult = c(0.01, 0.01))
      )
  }

  p_hd
}

compact_hd_line_colors <- c(
  total = "black",
  MS = "#B22222",
  MD = "#1F78B4",
  non_monetary = "#4D4D4D",
  monetary = "#7F0000",
  price_raw = "#E67E22",
  real_financial = "#2E8B57"
)

compact_hd_legend_labels <- c(
  total = "итого",
  MS = "MS",
  MD = "MD",
  non_monetary = "прочие",
  monetary = "ДКП",
  price_raw = "цены/сырьё",
  real_financial = "реал./фин."
)

make_compact_hd_plot <- function(plot_df, line_order, panel_order, facet_scales = "fixed") {
  y_limit <- max(abs(plot_df$value), na.rm = TRUE)

  if (!is.finite(y_limit) || y_limit == 0) {
    y_limit <- 1
  }

  p_compact_hd <- ggplot(
    plot_df,
    aes(x = date, y = value, colour = component)
  ) +
    geom_rect(
      data = shock_interval_df,
      aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.25
    ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      colour = "grey45"
    ) +
    geom_line(linewidth = 0.95) +
    facet_wrap(
      ~ panel,
      ncol = 2,
      scales = facet_scales
    ) +
    scale_color_manual(
      values = compact_hd_line_colors[line_order],
      breaks = line_order,
      labels = compact_hd_legend_labels[line_order]
    ) +
    labs(
      x = NULL,
      y = NULL,
      colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      strip.text = element_text(
        size = 15,
        face = "bold"
      ),
      axis.text = element_text(size = 11, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.4, "lines"),
      panel.spacing = unit(0.8, "lines"),
      plot.margin = margin(10, 10, 10, 10)
    )

  if (identical(facet_scales, "fixed")) {
    p_compact_hd <- p_compact_hd +
      coord_cartesian(
        ylim = c(-y_limit, y_limit)
      )
  }

  if (inherits(plot_df$date, "Date")) {
    p_compact_hd <- p_compact_hd +
      scale_x_date(
        date_breaks = "2 years",
        date_labels = "%Y",
        expand = expansion(mult = c(0.01, 0.01))
      )
  }

  p_compact_hd
}

money_hd_compact_df <- historical_output_df %>%
  dplyr::filter(
    response %in% c("m2", "key_rate"),
    shock %in% c("ForecastErrors", "MS", "MD", "Pcm", "Pim", "Py", "W", "y", "OvDeb")
  ) %>%
  dplyr::mutate(
    component = dplyr::case_when(
      shock == "ForecastErrors" ~ "total",
      shock == "MS" ~ "MS",
      shock == "MD" ~ "MD",
      TRUE ~ "non_monetary"
    ),
    value = dplyr::if_else(response == "m2", 100 * value, value),
    panel = dplyr::case_when(
      response == "m2" ~ "M2, log points × 100",
      response == "key_rate" ~ "Key rate, percentage points",
      TRUE ~ response
    )
  ) %>%
  dplyr::group_by(date, panel, component) %>%
  dplyr::summarise(value = sum(value), .groups = "drop") %>%
  dplyr::mutate(
    component = factor(
      component,
      levels = c("total", "MS", "MD", "non_monetary")
    ),
    panel = factor(
      panel,
      levels = c("M2, log points × 100", "Key rate, percentage points")
    )
  )

p_hd_money_block_compact <- make_compact_hd_plot(
  plot_df = money_hd_compact_df,
  line_order = c("total", "MS", "MD", "non_monetary"),
  panel_order = c("M2, log points × 100", "Key rate, percentage points"),
  facet_scales = "free_y"
)

ggsave(
  "figures/output/hd_money_block_compact.png",
  p_hd_money_block_compact,
  width = 13,
  height = 4.8,
  dpi = 300
)

prices_output_hd_compact_df <- historical_output_df %>%
  dplyr::filter(
    response %in% c("cpi", "gdp"),
    shock %in% c("ForecastErrors", "MS", "MD", "Pcm", "Pim", "Py", "W", "y", "OvDeb")
  ) %>%
  dplyr::mutate(
    component = dplyr::case_when(
      shock == "ForecastErrors" ~ "total",
      shock %in% c("MS", "MD") ~ "monetary",
      shock %in% c("Pcm", "Pim", "Py") ~ "price_raw",
      TRUE ~ "real_financial"
    ),
    panel = dplyr::case_when(
      response == "cpi" ~ "CPI",
      response == "gdp" ~ "GDP",
      TRUE ~ response
    )
  ) %>%
  dplyr::group_by(date, panel, component) %>%
  dplyr::summarise(value = sum(value), .groups = "drop") %>%
  dplyr::mutate(
    component = factor(
      component,
      levels = c("total", "monetary", "price_raw", "real_financial")
    ),
    panel = factor(
      panel,
      levels = c("CPI", "GDP")
    )
  )

prices_output_hd_check <- prices_output_hd_compact_df %>%
  dplyr::select(date, panel, component, value) %>%
  tidyr::pivot_wider(
    names_from = component,
    values_from = value
  ) %>%
  dplyr::mutate(
    gap = monetary + price_raw + real_financial - total
  )

prices_output_hd_max_gap <- max(abs(prices_output_hd_check$gap), na.rm = TRUE)

if (is.finite(prices_output_hd_max_gap) && prices_output_hd_max_gap > 1e-8) {
  warning(
    paste(
      "Compact HD aggregation check failed.",
      "Maximum gap:",
      signif(prices_output_hd_max_gap, 4)
    )
  )
}

p_hd_prices_output_compact <- make_compact_hd_plot(
  plot_df = prices_output_hd_compact_df,
  line_order = c("total", "monetary", "price_raw", "real_financial"),
  panel_order = c("CPI", "GDP")
)

ggsave(
  "figures/output/hd_prices_output_compact.png",
  p_hd_prices_output_compact,
  width = 13,
  height = 4.8,
  dpi = 300
)

shock_series_panel_labels <- c(
  Pcm = "P[cm]",
  MD = "MD",
  MS = "MS",
  Pim = "P[im]",
  W = "W",
  Py = "P[y]",
  y = "y",
  OvDeb = "T[bk]"
)

make_structural_shocks_plot <- function(shock_subset, file_stub) {
  plot_df <- eps_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(shock_subset),
      names_to = "shock",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      panel = factor(
        shock,
        levels = shock_subset,
        labels = shock_series_panel_labels[shock_subset]
      )
    )

  y_limit <- max(abs(plot_df$value), na.rm = TRUE)

  if (!is.finite(y_limit) || y_limit == 0) {
    y_limit <- 1
  }

  top_panel_label <- unname(shock_series_panel_labels[shock_subset[1]])

  interval_label_df <- shock_interval_df %>%
    dplyr::mutate(
      panel = factor(
        top_panel_label,
        levels = levels(plot_df$panel)
      ),
      y = y_limit * 0.9
    )

  p_shocks <- ggplot(
    plot_df,
    aes(x = date, y = value)
  ) +
    geom_rect(
      data = shock_interval_df,
      aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.25
    ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.25
    ) +
    geom_col(
      width = 75,
      fill = "black",
      colour = "black",
      linewidth = 0.1
    ) +
    facet_wrap(
      ~ panel,
      ncol = 1,
      scales = "fixed",
      labeller = label_parsed,
      axes = "all_x",
      axis.labels = "all_x"
    ) +
    coord_cartesian(
      ylim = c(-y_limit, y_limit)
    ) +
    geom_text(
      data = interval_label_df,
      aes(
        x = midpoint,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 2.2,
      vjust = 1,
      family = ""
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.35
      ),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 8,
        face = "bold"
      ),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      axis.ticks = element_line(linewidth = 0.25),
      panel.spacing.y = unit(0.12, "lines"),
      plot.margin = margin(8, 8, 8, 8)
    )

  if (inherits(plot_df$date, "Date")) {
    p_shocks <- p_shocks +
      scale_x_date(
        date_breaks = "1 year",
        date_labels = "%Y",
        expand = expansion(mult = c(0.01, 0.01))
      )
  }

  ggsave(
    file.path("figures/output", paste0(file_stub, ".pdf")),
    p_shocks,
    width = 8.5,
    height = 7.2
  )

  ggsave(
    file.path("figures/output", paste0(file_stub, ".png")),
    p_shocks,
    width = 8.5,
    height = 7.2,
    dpi = 400
  )
}

article_hd_responses <- c("m2", "key_rate", "cpi", "gdp")

for (response_name in article_hd_responses) {
  plot_object <- make_historical_decomposition_plot(response_name)
  file_stub <- paste0("historical_decomposition_", response_name)

  ggsave(
    file.path("figures/output", paste0(file_stub, ".pdf")),
    plot_object,
    width = 8.5,
    height = 12
  )

  ggsave(
    file.path("figures/output", paste0(file_stub, ".png")),
    plot_object,
    width = 8.5,
    height = 12,
    dpi = 400
  )
}

make_structural_shocks_plot(
  shock_subset = c("Pcm", "MD", "MS", "Pim"),
  file_stub = "structural_shocks_8a"
)

make_structural_shocks_plot(
  shock_subset = c("W", "Py", "y", "OvDeb"),
  file_stub = "structural_shocks_8b"
)

# =====================
# block 5: policy implications
# =====================

# =====================
# structural lag matrices
# =====================

Gamma <- array(
  NA_real_,
  dim = c(length(eq_names), length(var_names), model$p),
  dimnames = list(
    equation = eq_names,
    lagged_var = var_names,
    lag = paste0("lag_", 1:model$p)
  )
)

for (lag_id in 1:model$p) {
  Gamma[, , lag_id] <- A_hat %*% Phi[, , lag_id]
}

# =====================
# counterfactual policy systems
# =====================

make_policy_system <- function(policy_type) {
  A0_policy <- A_hat
  Gamma_policy <- Gamma
  
  if (policy_type == "fixed_key_rate") {
    # Фиксированная R:
    # уравнение MS заменяется условием на текущую key_rate.
    # Систематическая реакция MS на лаги выключается.
    A0_policy["MS", ] <- 0
    A0_policy["MS", "key_rate"] <- A_hat["MS", "key_rate"]
    
    Gamma_policy["MS", , ] <- 0
    
    return(
      list(
        A0 = A0_policy,
        Gamma = Gamma_policy
      )
    )
  }
  
  if (policy_type == "fixed_m2") {
    # Фиксированный M2:
    # уравнение MS заменяется условием на текущий M2.
    # Систематическая реакция MS на лаги выключается.
    A0_policy["MS", ] <- 0
    A0_policy["MS", "m2"] <- A_hat["MS", "m2"]
    
    Gamma_policy["MS", , ] <- 0
    
    return(
      list(
        A0 = A0_policy,
        Gamma = Gamma_policy
      )
    )
  }
  
  stop(paste("Unknown policy type:", policy_type))
}

# =====================
# counterfactual IRF function
# =====================

compute_counterfactual_irf <- function(policy_system) {
  A0_policy <- policy_system$A0
  Gamma_policy <- policy_system$Gamma
  
  Phi_policy <- array(
    NA_real_,
    dim = c(length(var_names), length(var_names), model$p),
    dimnames = list(
      response = var_names,
      lagged_var = var_names,
      lag = paste0("lag_", 1:model$p)
    )
  )
  
  for (lag_id in 1:model$p) {
    Phi_policy[, , lag_id] <- solve(A0_policy) %*% Gamma_policy[, , lag_id]
  }
  
  impact_policy <- solve(A0_policy)
  rownames(impact_policy) <- var_names
  colnames(impact_policy) <- eq_names
  
  Psi_policy <- array(
    0,
    dim = c(H + 1, k, k),
    dimnames = list(
      horizon = 0:H,
      response = var_names,
      innovation = var_names
    )
  )
  
  Psi_policy[1, , ] <- diag(k)
  
  for (h in 1:H) {
    Psi_h <- matrix(0, nrow = k, ncol = k)
    
    for (lag_id in 1:min(h, model$p)) {
      Psi_h <- Psi_h + Phi_policy[, , lag_id] %*% Psi_policy[h - lag_id + 1, , ]
    }
    
    Psi_policy[h + 1, , ] <- Psi_h
  }
  
  irf_policy <- array(
    NA_real_,
    dim = c(H + 1, k, length(eq_names)),
    dimnames = list(
      horizon = 0:H,
      response = var_names,
      shock = eq_names
    )
  )
  
  for (h in 0:H) {
    irf_policy[h + 1, , ] <- Psi_policy[h + 1, , ] %*% impact_policy
  }
  
  return(irf_policy)
}

# =====================
# compute counterfactual IRFs
# =====================

irf_fixed_key_rate <- compute_counterfactual_irf(
  make_policy_system("fixed_key_rate")
)

irf_fixed_m2 <- compute_counterfactual_irf(
  make_policy_system("fixed_m2")
)

# =====================
# helper: IRF array to long format
# =====================

irf_to_long <- function(irf_input, policy_label) {
  expand.grid(
    horizon = 0:H,
    response = var_names,
    shock = eq_names
  ) %>%
    dplyr::mutate(
      policy = policy_label,
      irf = as.vector(irf_input)
    )
}

# =====================
# common plotting function
# =====================

make_policy_plot <- function(plot_df) {
  ggplot(
    plot_df,
    aes(
      x = horizon,
      y = irf,
      linetype = policy
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.25
    ) +
    geom_line(
      linewidth = 0.35
    ) +
    scale_linetype_manual(
      values = c(
        "Исходная политика" = "solid",
        "Контрфактическая политика" = "dashed"
      )
    ) +
    facet_grid(
      response_plot ~ shock_plot,
      scales = "free_y",
      switch = "y",
      labeller = label_parsed
    ) +
    labs(
      x = NULL,
      y = NULL,
      linetype = NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.35
      ),
      panel.grid = element_blank(),
      
      strip.background = element_blank(),
      strip.placement = "outside",
      
      strip.text.x = element_text(
        size = 8,
        margin = margin(b = 3)
      ),
      strip.text.y.left = element_text(
        size = 8,
        angle = 0,
        margin = margin(r = 4)
      ),
      
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      
      panel.spacing.x = unit(0.12, "lines"),
      panel.spacing.y = unit(0.12, "lines"),
      
      plot.margin = margin(8, 8, 8, 8)
    )
}

make_policy_plot_presentation <- function(plot_df) {
  policy_levels <- levels(plot_df$policy)

  ggplot(
    plot_df,
    aes(
      x = horizon,
      y = irf,
      colour = policy,
      linetype = policy
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      colour = "grey45"
    ) +
    geom_line(
      linewidth = 0.6
    ) +
    scale_colour_manual(
      values = stats::setNames(c("black", "#B22222"), policy_levels)
    ) +
    scale_linetype_manual(
      values = stats::setNames(c("solid", "dashed"), policy_levels)
    ) +
    facet_grid(
      response_plot ~ shock_plot,
      scales = "free_y",
      switch = "y",
      labeller = label_parsed
    ) +
    scale_x_continuous(
      breaks = c(0, 10, 20),
      limits = c(0, H)
    ) +
    labs(
      x = "Horizon (quarters)",
      y = NULL,
      colour = NULL,
      linetype = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(
        size = 13,
        face = "bold",
        margin = margin(b = 5)
      ),
      strip.text.y.left = element_text(
        size = 13,
        face = "bold",
        angle = 0,
        margin = margin(r = 6)
      ),
      axis.text = element_text(size = 10, colour = "black"),
      axis.title.x = element_text(size = 12, margin = margin(t = 8)),
      axis.ticks = element_line(linewidth = 0.3),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      panel.spacing.x = unit(0.18, "lines"),
      panel.spacing.y = unit(0.18, "lines"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

make_policy_plot_compact <- function(plot_df) {
  policy_levels <- levels(plot_df$policy)

  ggplot(
    plot_df,
    aes(
      x = horizon,
      y = irf,
      colour = policy,
      linetype = policy
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      colour = "grey45"
    ) +
    geom_line(
      linewidth = 0.8
    ) +
    scale_colour_manual(
      values = stats::setNames(c("black", "#B22222"), policy_levels)
    ) +
    scale_linetype_manual(
      values = stats::setNames(c("solid", "dashed"), policy_levels)
    ) +
    facet_grid(
      response_plot ~ shock_plot,
      scales = "free_y",
      switch = "y",
      labeller = label_parsed
    ) +
    scale_x_continuous(
      breaks = c(0, 4, 8, 12, 16, 20, 24),
      limits = c(0, H)
    ) +
    labs(
      x = "Horizon (quarters)",
      y = NULL,
      colour = NULL,
      linetype = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(
        size = 14,
        face = "bold",
        margin = margin(b = 5)
      ),
      strip.text.y.left = element_text(
        size = 14,
        face = "bold",
        angle = 0,
        margin = margin(r = 7)
      ),
      axis.text = element_text(size = 11, colour = "black"),
      axis.title.x = element_text(size = 12, margin = margin(t = 8)),
      axis.ticks = element_line(linewidth = 0.3),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      panel.spacing.x = unit(0.35, "lines"),
      panel.spacing.y = unit(0.35, "lines"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# =====================
# Figure: fixed R policy
# =====================

policy_fixed_key_rate_df <- dplyr::bind_rows(
  irf_to_long(irf_array, "Исходная политика"),
  irf_to_long(irf_fixed_key_rate, "Контрфактическая политика")
) %>%
  dplyr::mutate(
    response_plot = factor(
      response,
      levels = var_names,
      labels = var_plot_labels[var_names]
    ),
    shock_plot = factor(
      shock,
      levels = eq_names,
      labels = eq_plot_labels[eq_names]
    ),
    policy = factor(
      policy,
      levels = c(
        "Исходная политика",
        "Контрфактическая политика"
      )
    )
  )

p_policy_fixed_key_rate <- make_policy_plot(policy_fixed_key_rate_df)

p_policy_fixed_key_rate_presentation <- make_policy_plot_presentation(
  policy_fixed_key_rate_df
)

policy_fixed_key_rate_compact_df <- policy_fixed_key_rate_df %>%
  dplyr::filter(
    response %in% irf_compact_response_order,
    shock %in% irf_compact_shock_order
  ) %>%
  dplyr::mutate(
    response_plot = factor(
      response,
      levels = irf_compact_response_order,
      labels = var_plot_labels[irf_compact_response_order]
    ),
    shock_plot = factor(
      shock,
      levels = irf_compact_shock_order,
      labels = eq_plot_labels[irf_compact_shock_order]
    )
  )

p_policy_fixed_key_rate_compact <- make_policy_plot_compact(
  policy_fixed_key_rate_compact_df
)

ggsave(
  "figures/output/policy_fixed_key_rate.pdf",
  p_policy_fixed_key_rate,
  width = 8.5,
  height = 8.5
)

ggsave(
  "figures/output/policy_fixed_key_rate.png",
  p_policy_fixed_key_rate,
  width = 8.5,
  height = 8.5,
  dpi = 400
)

ggsave(
  "figures/output/policy_fixed_key_rate_presentation.png",
  p_policy_fixed_key_rate_presentation,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  "figures/output/policy_fixed_key_rate_compact.png",
  p_policy_fixed_key_rate_compact,
  width = 13,
  height = 6,
  dpi = 300
)

# =====================
# Figure: fixed M2 policy
# =====================

policy_fixed_m2_df <- dplyr::bind_rows(
  irf_to_long(irf_array, "Исходная политика"),
  irf_to_long(irf_fixed_m2, "Контрфактическая политика")
) %>%
  dplyr::mutate(
    response_plot = factor(
      response,
      levels = var_names,
      labels = var_plot_labels[var_names]
    ),
    shock_plot = factor(
      shock,
      levels = eq_names,
      labels = eq_plot_labels[eq_names]
    ),
    policy = factor(
      policy,
      levels = c(
        "Исходная политика",
        "Контрфактическая политика"
      )
    )
  )

p_policy_fixed_m2 <- make_policy_plot(policy_fixed_m2_df)

p_policy_fixed_m2_presentation <- make_policy_plot_presentation(
  policy_fixed_m2_df
)

policy_fixed_m2_compact_df <- policy_fixed_m2_df %>%
  dplyr::filter(
    response %in% irf_compact_response_order,
    shock %in% irf_compact_shock_order
  ) %>%
  dplyr::mutate(
    response_plot = factor(
      response,
      levels = irf_compact_response_order,
      labels = var_plot_labels[irf_compact_response_order]
    ),
    shock_plot = factor(
      shock,
      levels = irf_compact_shock_order,
      labels = eq_plot_labels[irf_compact_shock_order]
    )
  )

p_policy_fixed_m2_compact <- make_policy_plot_compact(
  policy_fixed_m2_compact_df
)

ggsave(
  "figures/output/policy_fixed_m2.pdf",
  p_policy_fixed_m2,
  width = 8.5,
  height = 8.5
)

ggsave(
  "figures/output/policy_fixed_m2.png",
  p_policy_fixed_m2,
  width = 8.5,
  height = 8.5,
  dpi = 400
)

ggsave(
  "figures/output/policy_fixed_m2_presentation.png",
  p_policy_fixed_m2_presentation,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  "figures/output/policy_fixed_m2_compact.png",
  p_policy_fixed_m2_compact,
  width = 13,
  height = 6,
  dpi = 300
)
