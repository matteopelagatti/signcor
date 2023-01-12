library(tidyverse)
library(ggplot2)
library(pcaPP)


n <- 5000
nsim <- 10000
with_kendall <- TRUE


# ---- Sim Stochastic Volatility model
rho <- 0.7

smpl_cor <- matrix(0, n, nsim)
my_cor   <- matrix(0, n, nsim)
if (with_kendall) {
  kend_cor <- matrix(0, n, nsim)
}

cat("Example 1\n")

set.seed(2020)
for (i in 1:nsim){
  if (i %% 100 == 0) {
    cat("simulation", i , " ")
    print(Sys.time())
  }
  log_var <- apply(matrix(rnorm(2 * n, sd = 0.04), n, 2), 2, cumsum)
  xi1 <- rnorm(n)
  xi2 <- rho * xi1 + sqrt(1 - rho^2) * rnorm(n)
  x1 <- xi1 * exp(log_var[, 1]/2)
  x2 <- xi2 * exp(log_var[, 2]/2)
  
  x1_ <- x1 - cumsum(x1)/(1:n)
  x2_ <- x2 - cumsum(x2)/(1:n)

  x1_med <- x1 - sapply(1:n, function(i) median(x1[1:i]))
  x2_med <- x2 - sapply(1:n, function(i) median(x2[1:i]))

  smpl_cov  <- cumsum(x1_*x2_)/(1:n)
  smpl_var1 <- cumsum(x1_^2)/(1:n)
  smpl_var2 <- cumsum(x2_^2)/(1:n)
  smpl_cor[, i]  <- smpl_cov / sqrt(smpl_var1 * smpl_var2)
  
  smpl_sign_cov <- cumsum(sign(x1_med*x2_med))/(1:n)
  my_cor[, i]   <- sin(pi / 2 * smpl_sign_cov)
  if (with_kendall) {
    kend_cor[-(1:10), i] <- sin(
      pi / 2 * c(sapply(11:n, function(k) cor.fk(x1[1:k], x2[1:k])))
    )
  }
}


smpl_cor_q <- t(apply(smpl_cor[-(1:10),], 1, quantile, prob = c(0.05, 0.5, 0.95)))
my_cor_q   <- t(apply(my_cor[-(1:10),], 1, quantile, prob = c(0.05, 0.5, 0.95)))
my_asy_sd <- sqrt((1 - rho^2) * (pi^2/4 - asin(rho)^2))
my_cor_th_q <- cbind(lwr = rho + my_asy_sd * qnorm(0.05) / sqrt(11:n),
                     med = rho,
                     upr = rho + my_asy_sd * qnorm(0.95) / sqrt(11:n))
if (with_kendall) {
  kend_cor_q <- t(apply(kend_cor[-(1:10), ], 1, quantile, prob = c(0.05, 0.5, 0.95)))
#  kend_cor_q <- rbind(matrix(NA, 10, 3), kend_cor_q)
}

df <- data.frame(n = 11:n, estimator = "SMPL COR (MC)", quantile = "q(0.05)", value = smpl_cor_q[, 1])
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SMPL COR (MC)", quantile = "q(0.50)", value = smpl_cor_q[, 2]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SMPL COR (MC)", quantile = "q(0.95)", value = smpl_cor_q[, 3]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (MC)", quantile = "q(0.05)", value = my_cor_q[, 1]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (MC)", quantile = "q(0.50)", value = my_cor_q[, 2]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (MC)", quantile = "q(0.95)", value = my_cor_q[, 3]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (ASY)", quantile = "q(0.05)", value = my_cor_th_q[, 1]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (ASY)", quantile = "q(0.50)", value = my_cor_th_q[, 2]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (ASY)", quantile = "q(0.95)", value = my_cor_th_q[, 3]))
if (with_kendall) {
  df <- df %>% rbind(data.frame(n = 11:n, estimator = "KEND COR (MC)", quantile = "q(0.05)", value = kend_cor_q[, 1]))
  df <- df %>% rbind(data.frame(n = 11:n, estimator = "KEND COR (MC)", quantile = "q(0.50)", value = kend_cor_q[, 2]))
  df <- df %>% rbind(data.frame(n = 11:n, estimator = "KEND COR (MC)", quantile = "q(0.95)", value = kend_cor_q[, 3]))
}

# df %>%
#   filter(n >= 100) %>%
#   ggplot(aes(x = n, y = value, linetype = quantile, color = estimator)) +
#   geom_line() + geom_hline(yintercept = rho, linetype = 2) +
#   theme_bw(base_family = "serif")


df %>%
  filter(n >= 100) %>%
  mutate(`n (in thousands)` = n/1000) %>%
  ggplot(aes(x = `n (in thousands)`, y = value, linetype = quantile, color = quantile)) +
  geom_line() + geom_hline(yintercept = rho, linetype = 2) +
  facet_grid(cols = vars(fct_relevel(estimator, c("SIGN COR (ASY)",
                                                  "SIGN COR (MC)",
                                                  "KEND COR (MC)",
                                                  "SMPL COR (MC)")))) +
  theme_bw(base_family = "serif")

ggsave("rho_estim_sv_facet_demean.pdf", width = 6, height = 3, family = "Times")


# Stationary CCC-GARCH
rho <- 0.8

smpl_cor <- matrix(0, n, nsim)
my_cor   <- matrix(0, n, nsim)
if (with_kendall) {
  kend_cor <- matrix(0, n, nsim)
}


cccsim <- function(nsim, n, omega, alpha, beta, R) {
  m <- nrow(R)
  eps <- array(0, c(n, m, nsim))
  cva <- array(0, c(n, m, nsim))
  v <- omega / (1 - alpha - beta)
  s <- sqrt(v)
  for (i in seq(nsim)) {
    xi <- mvnfast::rmvt(n, numeric(m), R, df = 5)
    eps[1, , i] <- xi[1, ] * s
    cva[1, , i] <- v
    for (t in 2:n) {
      cva[t, , i] <- omega + alpha * eps[t-1, , i]^2 + beta * cva[t-1, , i]
      eps[t, , i] <- xi[t, ] * sqrt(cva[t, , i])
    }
  }
  list(eps = eps, cva = cva)
}

# garchfit <- function(x, init, var1) {
#   eps2 <- x^2
#   n <- length(x)
#   obj <- function(par) {
#     vart <- stats::filter(x = par[1]^2 + par[2]^2*eps2,
#                           filter = par[3]^2,
#                           method = "recursive",
#                           init = var1)
#     mean(log(vart[-n]) + eps2[-1]/vart[-n])
#   }
#   opt <- optim(par = sqrt(init), fn = obj, method = "BFGS")
#   pars <- opt$par^2
#   x / sqrt(c(var1,
#          stats::filter(pars[1] + pars[2]*eps2, pars[3], "recursive", init = var1))[-n])
# }

cat("Example 2")

set.seed(2020)
R <- matrix(c(1, rho, rho, 1), 2, 2)
out <- cccsim(nsim, n, c(0.002, 0.001), c(0.05, 0.06), c(0.89, 0.90), R)

# zeta1 <- apply(out$eps[, 1, ], 2, garchfit, init = c(0.002, 0.05, 0.89), var1 = 0.03333)
# zeta2 <- apply(out$eps[, 2, ], 2, garchfit, init = c(0.001, 0.06, 0.90), var1 = 0.025)
# zeta  <- array(cbind(zeta1, zeta2), c(n, nsim, 2))

sgn_cor <- apply(out$eps, 3,
                 function(x) sin(pi/2 * 
                                   cumsum(sign((x[, 1] - sapply(1:n, function(i) median(x[1:i, 1]))) * 
                                               (x[, 2] - sapply(1:n, function(i) median(x[1:i, 2])))
                                               )
                                          ) / (1:nrow(x))
                                 )
                 )
smp_cor <- apply(out$eps, 3,
                 function(x) {
                   x1 <- x[, 1] - cumsum(x[, 1])/(1:nrow(x))
                   x2 <- x[, 2] - cumsum(x[, 2])/(1:nrow(x))
                   cumsum(x1*x2)/sqrt(cumsum(x1*x1)*cumsum(x2*x2))
                 }
                 )
if (with_kendall) {
  for (i in 1:nsim) {
    if (i %% 100 == 0) cat("nsim =", i, "\n")
    x1 <- out$eps[, 1, i]
    x2 <- out$eps[, 2, i]
    kend_cor[-(1:10), i] <- sin(
      pi / 2 * c(sapply(11:n, function(k) cor.fk(x1[1:k], x2[1:k])))
    )
  }
}  
# ccc_cor <- apply(zeta, 2, function(x) cumsum(Rfast::rowprods(x))/sqrt(cumsum(x[, 1]^2)*cumsum(x[, 2]^2)))

smpl_cor_q <- t(apply(smp_cor[-(1:10),], 1, quantile, prob = c(0.05, 0.5, 0.95)))
my_cor_q   <- t(apply(sgn_cor[-(1:10),], 1, quantile, prob = c(0.05, 0.5, 0.95)))
my_asy_sd <- sqrt((1 - rho^2) * (pi^2/4 - asin(rho)^2))
my_cor_th_q <- cbind(lwr = rho + my_asy_sd * qnorm(0.05) / sqrt(11:n),
                     med = rho,
                     upr = rho + my_asy_sd * qnorm(0.95) / sqrt(11:n))
if (with_kendall) {
  kend_cor_q <- t(apply(kend_cor[-(1:10),], 1, quantile, prob = c(0.05, 0.5, 0.95)))
  #  kend_cor_q <- rbind(matrix(NA, 10, 3), kend_cor_q)
}

df <- data.frame(n = 11:n, estimator = "SMPL COR (MC)", quantile = "q(0.05)", value = smpl_cor_q[, 1])
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SMPL COR (MC)", quantile = "q(0.50)", value = smpl_cor_q[, 2]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SMPL COR (MC)", quantile = "q(0.95)", value = smpl_cor_q[, 3]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (MC)", quantile = "q(0.05)", value = my_cor_q[, 1]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (MC)", quantile = "q(0.50)", value = my_cor_q[, 2]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (MC)", quantile = "q(0.95)", value = my_cor_q[, 3]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (ASY)", quantile = "q(0.05)", value = my_cor_th_q[, 1]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (ASY)", quantile = "q(0.50)", value = my_cor_th_q[, 2]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (ASY)", quantile = "q(0.95)", value = my_cor_th_q[, 3]))
if (with_kendall) {
  df <- df %>% rbind(data.frame(n = 11:n, estimator = "KEND COR (MC)", quantile = "q(0.05)", value = kend_cor_q[, 1]))
  df <- df %>% rbind(data.frame(n = 11:n, estimator = "KEND COR (MC)", quantile = "q(0.50)", value = kend_cor_q[, 2]))
  df <- df %>% rbind(data.frame(n = 11:n, estimator = "KEND COR (MC)", quantile = "q(0.95)", value = kend_cor_q[, 3]))
}

# df %>%
#   filter(n >= 100) %>%
#   ggplot(aes(x = n, y = value, linetype = quantile, color = estimator)) +
#   geom_line() + geom_hline(yintercept = rho, linetype = 2) +
#   theme_bw(base_family = "serif")


df %>%
  filter(n >= 100) %>%
  mutate(`n (in thousands)` = n/1000) %>%
  ggplot(aes(x = `n (in thousands)`, y = value, linetype = quantile, color = quantile)) +
  geom_line() + geom_hline(yintercept = rho, linetype = 2) +
  facet_grid(cols = vars(fct_relevel(estimator, c("SIGN COR (ASY)",
                                                  "SIGN COR (MC)",
                                                  "KEND COR (MC)",
                                                  "SMPL COR (MC)")))) +
  theme_bw(base_family = "serif")

ggsave("rho_estim_garch_facet_demean.pdf", width = 6, height = 3, family = "Times")


# ---- Sim two-variance model
rho <- -0.9

cat("Example 3")

set.seed(2020)
smpl_cor <- matrix(0, n, nsim)
my_cor   <- matrix(0, n, nsim)
if (with_kendall) {
  kend_cor <- matrix(0, n, nsim)
}

for (i in 1:nsim){
  if (i %% 100 == 0) {
    cat("simulation", i , " ")
    print(Sys.time())
  }
  x1 <- rnorm(n)
  x2 <- rho * x1 + sqrt(1 - rho^2) * rnorm(n)
  x2[seq(1, n, 2)] <- x2[seq(1, n, 2)]*5
  
  x1_ <- x1 - cumsum(x1)/(1:n)
  x2_ <- x2 - cumsum(x2)/(1:n)

  x1_med <- x1 - sapply(1:n, function(i) median(x1[1:i]))
  x2_med <- x2 - sapply(1:n, function(i) median(x2[1:i]))
  
  smpl_cov  <- cumsum(x1_*x2_)/(1:n)
  smpl_var1 <- cumsum(x1_^2)/(1:n)
  smpl_var2 <- cumsum(x2_^2)/(1:n)
  smpl_cor[, i]  <- smpl_cov / sqrt(smpl_var1 * smpl_var2)
  
  smpl_sign_cov <- cumsum(sign(x1_med*x2_med))/(1:n)
  my_cor[, i]   <- sin(pi / 2 * smpl_sign_cov)
  if (with_kendall) {
    kend_cor[-(1:10), i] <- sin(
      pi / 2 * c(sapply(11:n, function(k) cor.fk(x1[1:k], x2[1:k])))
    )
  }
}


smpl_cor_q <- t(apply(smpl_cor[-(1:10),], 1, quantile, prob = c(0.05, 0.5, 0.95)))
my_cor_q   <- t(apply(my_cor[-(1:10),], 1, quantile, prob = c(0.05, 0.5, 0.95)))
my_asy_sd <- sqrt((1 - rho^2) * (pi^2/4 - asin(rho)^2))
my_cor_th_q <- cbind(lwr = rho + my_asy_sd * qnorm(0.05) / sqrt(11:n),
                     med = rho,
                     upr = rho + my_asy_sd * qnorm(0.95) / sqrt(11:n))
if (with_kendall) {
  kend_cor_q <- t(apply(kend_cor[-(1:10),], 1, quantile, prob = c(0.05, 0.5, 0.95)))
  #  kend_cor_q <- rbind(matrix(NA, 10, 3), kend_cor_q)
}

df <- data.frame(n = 11:n, estimator = "SMPL COR (MC)", quantile = "q(0.05)", value = smpl_cor_q[, 1])
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SMPL COR (MC)", quantile = "q(0.50)", value = smpl_cor_q[, 2]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SMPL COR (MC)", quantile = "q(0.95)", value = smpl_cor_q[, 3]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (MC)", quantile = "q(0.05)", value = my_cor_q[, 1]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (MC)", quantile = "q(0.50)", value = my_cor_q[, 2]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (MC)", quantile = "q(0.95)", value = my_cor_q[, 3]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (ASY)", quantile = "q(0.05)", value = my_cor_th_q[, 1]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (ASY)", quantile = "q(0.50)", value = my_cor_th_q[, 2]))
df <- df %>% rbind(data.frame(n = 11:n, estimator = "SIGN COR (ASY)", quantile = "q(0.95)", value = my_cor_th_q[, 3]))
if (with_kendall) {
  df <- df %>% rbind(data.frame(n = 11:n, estimator = "KEND COR (MC)", quantile = "q(0.05)", value = kend_cor_q[, 1]))
  df <- df %>% rbind(data.frame(n = 11:n, estimator = "KEND COR (MC)", quantile = "q(0.50)", value = kend_cor_q[, 2]))
  df <- df %>% rbind(data.frame(n = 11:n, estimator = "KEND COR (MC)", quantile = "q(0.95)", value = kend_cor_q[, 3]))
}

# df %>%
#   filter(n >= 100) %>%
#   ggplot(aes(x = n, y = value, linetype = quantile, color = estimator)) +
#   geom_line() + geom_hline(yintercept = rho, linetype = 2) +
#   theme_bw(base_family = "serif")


df %>%
  filter(n >= 100) %>%
  mutate(`n (in thousands)` = n/1000) %>%
  ggplot(aes(x = `n (in thousands)`, y = value, linetype = quantile, color = quantile)) +
  geom_line() + geom_hline(yintercept = rho, linetype = 2) +
  facet_grid(cols = vars(fct_relevel(estimator, c("SIGN COR (ASY)",
                                                  "SIGN COR (MC)",
                                                  "KEND COR (MC)",
                                                  "SMPL COR (MC)")))) +
  theme_bw(base_family = "serif")

ggsave("rho_estim_mix_facet_demean.pdf", width = 6, height = 3, family = "Times")


# Kendall's tau
microbenchmark::microbenchmark(
  tau1 <- Kendall::Kendall(x1, x2),
  tau2 <- pcaPP::cor.fk(x1, x2),
  nu1  <- mean(sign(x1*x2)),
  times = 10
)

vrho <- seq(0, 0.95, 0.05)
n <- 100
tau <- nu <- matrix(0, nsim, length(vrho))
set.seed(2021)
for (i in seq_along(vrho)) {
  for (s in 1:nsim) {
    z1 <- rnorm(n)
    z2 <- vrho[i] * z1 + sqrt(1 - vrho[i]^2) * rnorm(n)
    tau[s, i] <- pcaPP::cor.fk(z1, z2)
    nu[s, i] <- mean(sign(z1*z2))
  }
}

colnames(tau) <- vrho
colnames(nu)  <- vrho

df_values <- rbind(
  cbind(as.data.frame(tau), statistic = "tau"),
  cbind(as.data.frame(nu), statistic = "nu")) %>%
  pivot_longer(colnames(tau), names_to = "rho")
df_values$estimate <- sin(df_values$value * pi / 2)

mse <- function(x, m) mean((x-m)^2)

df_values %>% group_by(statistic, rho) %>%
  summarise(mse_phi = mse(estimate, as.numeric(first(rho)))) ->
  df_mse

df_mse %>%
  pivot_wider(names_from = statistic, values_from = mse_phi) %>%
  mutate(eff = tau/nu, rho = as.numeric(rho)) %>%
  ggplot(aes(x = rho, y = eff)) +
  geom_line()

ggplot(df_values, aes(x = rho, y = estimate, fill = statistic)) +
  geom_boxplot(outlier.shape = NA)

# ---- speed
library(microbenchmark)
vn <- 10^(1:8)
spd <- data.frame(expr = character(0), time = numeric(0), n = integer(0))
for (n in vn) {
  w1 <- rnorm(n)
  w2 <- rnorm(n)
  spd <- rbind(spd,
               cbind(microbenchmark(
                       nu  <- mean(sign(w1*w2)),
                       tau <- pcaPP::cor.fk(w1, w2),
                     times = 10),
                     n = n))
}

spd$expr <- trimws(substr(spd$expr, 1, 3))

spd %>% group_by(n, expr) %>% summarise(mean_time = mean(time)) %>%
  ggplot(aes(x = n, y = mean_time, col = expr)) + geom_line() + geom_point() +
  scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')

spd %>% group_by(n, expr) %>% summarise(mean_time = mean(time)) %>%
  pivot_wider(names_from = expr, values_from = mean_time) %>%
  mutate(relative_speed = nu/tau) %>%
  ggplot(aes(x = n, y = relative_speed)) +
  geom_line() + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')
