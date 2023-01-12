library(tidyverse)
library(rmgarch)
library(xts)

#functions
dynscor <- function(X, nobs_init = 100, transf = function(x) sin(pi*x/2),
                    hold_out = 0) {
  Y <- sign(as.matrix(X))
  k <- ncol(Y)
  n <- nrow(Y)
  P <- matrix(as.integer(0), n, k*(k-1)/2,
              dimnames = list(NULL, paste("V", 1:(k*(k-1)/2))))
  cnt <- 1
  for (i in 1:(k-1)) for (j in (i+1):k) {
    P[, cnt] <- as.integer((Y[, i]*Y[, j] + 1)/2)
    colnames(P)[cnt] <- paste(colnames(Y)[i], colnames(Y)[j], sep = "_")
    cnt <- cnt + 1
  }
  Pmean <- colMeans(P)
  Pmean0 <- if (nobs_init > 1 && nobs_init <= n) colMeans(P[1:nobs_init, , drop = FALSE]) else Pmean
  wrange <- seq(n - hold_out)
  wP <- P[wrange, ]
  obj <- function(alpha) {
    if (alpha == 0) {
      p <- rep(Pmean, each = n)
    } else {
      p <- stats::filter(rbind(Pmean0, wP*alpha),
                         1 - alpha, "recursive")[wrange, ]
    }
    -mean(wP*log(p) + (1 - wP)*log(1-p))
  }
  opt <- optimize(obj, interval = c(0, 1))
  scor <- transf(stats::filter(rbind(Pmean0, P*opt$minimum),
                               1 - opt$minimum, "recursive")*2 - 1)
  colnames(scor) <- colnames(P)
  list(alpha = opt$minimum, logLik = -n*opt$minimum, holdout = hold_out,
       scor = if (class(X)[1] == "xts") xts(scor[-nrow(scor), ], time(X)) else scor)
}

dcc <- function(X, udist = c("norm", "std"), mdist = c("mvnorm", "mvt"),
                include_mean = FALSE, model = "iGARCH", hold_out = 0) {
  k <- ncol(X)
  garch11spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = include_mean),
                            variance.model = list(garchOrder = c(1, 1), model = model),
                            distribution.model = udist[1])
  
  dccspec <- dccspec(uspec = multispec(replicate(k, garch11spec)),
                     dccOrder = c(1, 1),
                     distribution = mdist[1])
  
  dccfit(dccspec, data = X, out.sample = hold_out)
}

matrify <- function(x, on_diag = 1) {
  m <- (1 + sqrt(1 + 8*length(x)))/2
  if (m != floor(m)) stop("wrong dimensions")
  M <- matrix(NA, m, m)
  M[lower.tri(M)] <- x
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  diag(M) <- on_diag
  M
}

gmv_portfolio <- function(S, long = FALSE) {
  if (long) {
    m <- dim(S)[1]
    A <- cbind(1, diag(m))
    b <- c(1, rep(0, m))
    qpout <- quadprog::solve.QP(S, rep(0, m), A, b, 1)
    return(qpout$solution * (qpout$solution >= 0))
  }
  iS <- chol2inv(chol(S))
  cs <- colSums(iS)
  cs / sum(cs)
}

sim2series <- function(n, df = 5, cor_type = c("sin", "step", "rw", "sigmoid")) {
  v1 <- exp(cumsum(rnorm(n, sd = 0.01)))/1000
  v2 <- exp(cumsum(rnorm(n, sd = 0.01)))/1000
  crs <- switch(cor_type[1],
                sin     = sin((1:n)*pi*2/n),
                step    = -1.5 * ((1:n) < 1000) + 0.9,
                rw      = 1/(1 + exp(-cumsum(rnorm(n, sd = 0.05)))),
                sigmoid = 1/(1 + exp(-((1:n)/n - 0.5)*5))
                )
  z1 <- rt(n, df)
  z2 <- z1 * crs + sqrt(1-crs^2) * rt(n, df)
  X <- cbind(z1*v1, z2*v2)
  list(ret = X, cor = crs)
}

plot_cor <- function(n = 2000, df = 5, cor_types = c("sin", "step", "rw", "sigmoid")) {
  dt <- tibble()
  for (cor_type in cor_types){
    tmp <- sim2series(n, df, cor_type)
    dcc <- dcc(tmp$ret)
    dyn <- dynscor(tmp$ret)
    dt <- rbind(dt,
                tibble(shape = cor_type,
                       time = 1:n,
                       actual = tmp$cor,
                       dcc_cor = rcor(dcc)[1, 2, ],
                       sign_cor = dyn$scor[-(n+1)])
          )
  }
  dt <- pivot_longer(dt, actual:sign_cor, names_to = "type", values_to = "correlation")
  ggplot(dt, aes(x = time, y = correlation, color = type, linetype = type)) +
    geom_line() +
    facet_wrap(vars(shape)) +
    theme_bw(base_family = "serif")
}

make_sims <- function(nsim, n, df, cor_type,
                      udist = c("norm", "std"), mdist = c("mvnorm", "mvt")) {
  mae <- matrix(0, nsim, 2, dimnames = list(NULL, c("DCC", "DSCOR")))
  mse <- matrix(0, nsim, 2, dimnames = list(NULL, c("DCC", "DSCOR")))
  systime <- matrix(0, nsim, 2, dimnames = list(NULL, c("DCC", "DSCOR")))
  for (i in seq(nsim)) {
    sms <- sim2series(n, df, cor_type)
    systime[i, 1] <- system.time(dcc <- rcor(dcc(sms$ret, udist, mdist))[1, 2, ])[2]
    systime[i, 2] <- system.time(sco <- dynscor(sms$ret)$scor[-(n+1)])[2]
    dcc_dif <- dcc - sms$cor
    sco_dif <- sco - sms$cor
    mae[i, ] <- c(mean(abs(dcc_dif)), mean(abs(sco_dif)))
    mse[i, ] <- c(mean((dcc_dif)^2), mean((sco_dif)^2))
  }
  list(mae = mae, mse = mse, systime = systime)
}

# script

set.seed(400)
dt <- tibble()
for (cor_type in c("sin", "step", "rw", "sigmoid")) {
  for (df in c(8, 4, 2)) {
    cat(cor_type, df, "\n")
    ex1 <- make_sims(nsim = 10, n = 2000, df = df, cor_type = cor_type)
    lex1 <- lapply(ex1, colMeans)
    dt <- rbind(dt,
                tibble(cor_type = cor_type,
                       df = df,
                       mae_dcc = lex1$mae[1],
                       mae_sign_cor = lex1$mae[2],
                       mse_dcc = lex1$mse[1],
                       mse_sign_cor = lex1$mse[2],
                ))
  }
}

dt %>%
  mutate(mae_ratio = mae_sign_cor / mae_dcc,
         mse_ratio = mse_sign_cor / mse_dcc) %>%
  relocate(mae_ratio, .after = mae_sign_cor) %>%
  relocate(mse_ratio, .after = mse_sign_cor) %>%
  arrange(desc(df), cor_type) %>%
  write.csv("cor_estimate.csv")
  

set.seed(444)
plot_cor(); ggsave("cor_estimate.pdf", width = 7, height = 6)

# ---- Application to 10 SP100 stocks

# -- data preparation
pri <- readxl::read_xlsx("SP100_30y.xlsx", sheet = "Foglio1")
names(pri) <- stringr::str_remove(names(pri), "-US")
pri$DATE <- as.Date(pri$DATE, format = "%m/%d/%Y")

ret <- pri %>%
  transmute_at(.vars = vars(RTX:DD), .funs = function(x) 100*(x - lag(x))/lag(x)) %>%
  mutate(DATE = pri$DATE) %>%
  relocate(DATE) %>%
  slice(-1)
  
# -- sample selection
vstocks <- 1:10
vobs <- 5044:7812 # from 2010-01-04 to 2020-12-31
hout <- 505 # last two years

any(is.na(ret[vobs, 1 + vstocks]))

# -- univariate GARCH part
uspec <- ugarchspec(variance.model = list(model = "sGARCH"),
                    mean.model = list(armaOrder = c(0, 0),
                                      include.mean = TRUE))
gout <- list()
gfcs <- matrix(0, hout, length(vstocks),
               dimnames = list(NULL, names(ret[, -1])[vstocks]))
for (nm in names(ret[, -1])[vstocks]) {
  gout[[nm]] <- ugarchfit(uspec, xts(pull(ret[vobs, ], nm), ret$DATE[vobs]), out.sample = hout)
  gfcs[, nm] <- ugarchforecast(gout[[nm]], n.ahead = 1, n.roll = hout) %>%
    sigma() %>%
    as.numeric() %>%
    `[`(1:hout)
}

# -- dynamic sign correlation part
dscor <- dynscor(xts(ret[vobs, 1 + vstocks], ret$DATE[vobs]), hold_out = hout)
dscor_fcs <- dscor$scor["2019-01-01/"]

# -- DCC
fitted_dcc <- dcc(xts(ret[vobs, 1 + vstocks], ret$DATE[vobs]),
                  include_mean = TRUE, model = "sGARCH",
                  hold_out = hout)

fcs_dcc <- dccforecast(fitted_dcc, n.roll = hout)

# ---- Portfolios
our_p <- matrix(0, hout, length(vstocks))
dcc_p <- matrix(0, hout, length(vstocks))
our_p0 <- matrix(0, hout, length(vstocks))
dcc_p0 <- matrix(0, hout, length(vstocks))

for (t in seq(hout)) {
  S <- Matrix::nearPD(matrify(dscor_fcs[t, ]), corr = TRUE)$mat *
    outer(sqrt(gfcs[t, ]), sqrt(gfcs[t, ]))
  our_p[t, ] <- gmv_portfolio(S)
  our_p0[t, ] <- gmv_portfolio(S, TRUE)
  
  S <- fcs_dcc@mforecast$H[[t]][,, 1]
  dcc_p[t, ] <- gmv_portfolio(S)
  dcc_p0[t, ] <- gmv_portfolio(S, TRUE)
}

our_p_ret <- rowSums(our_p * as.matrix(tail(ret[vobs, 1 + vstocks], hout)))
our_p0_ret <- rowSums(our_p0 * as.matrix(tail(ret[vobs, 1 + vstocks], hout)))
dcc_p_ret <- rowSums(dcc_p * as.matrix(tail(ret[vobs, 1 + vstocks], hout)))
dcc_p0_ret <- rowSums(dcc_p0 * as.matrix(tail(ret[vobs, 1 + vstocks], hout)))

plot(cumprod(1 + our_p_ret/100), type = "l")
plot(cumprod(1 + our_p0_ret/100), type = "l")
plot(cumprod(1 + dcc_p_ret/100), type = "l")
plot(cumprod(1 + dcc_p0_ret/100), type = "l")

var(our_p_ret)
var(dcc_p_ret)
var(our_p0_ret)
var(dcc_p0_ret)

# PLOTS
# variances
vol <- data.frame(
  DATE = tail(ret$DATE[vobs], hout),
  PORTFOLIO = "unconstrained",
  SCOR = sqrt(cumsum(our_p_ret^2)/(1:hout)*252),
  DCC  = sqrt(cumsum(dcc_p_ret^2)/(1:hout)*252)
)

vol <- vol %>% rbind(data.frame(
  DATE = tail(ret$DATE[vobs], hout),
  PORTFOLIO = "no short",
  SCOR = sqrt(cumsum(our_p0_ret^2)/(1:hout)*252),
  DCC  = sqrt(cumsum(dcc_p0_ret^2)/(1:hout)*252)
))

vol <- vol %>% pivot_longer(cols = SCOR:DCC,
                            names_to = "MODEL",
                            values_to = "VOLATILITY")

vol %>%
  filter(DATE >= "2019-02-01") %>%
  ggplot(aes(x = DATE, y = VOLATILITY, color = MODEL, linetype = PORTFOLIO)) +
  geom_line() +
  scale_x_date(minor_breaks = as.Date(
      paste0(rep(c(2019, 2020), each = 12), "-", 1:12, "-01")
    )
  ) +
  theme_bw(base_family = "serif")

ggsave("portfolio_vol.pdf", width = 7, height = 4)

# sharpe ratios
sharpe <- data.frame(
  DATE = tail(ret$DATE[vobs], hout),
  PORTFOLIO = "unconstrained",
  SCOR = (cumsum(our_p_ret)/(1:hout)*252/sqrt(cumsum(our_p_ret^2)/(1:hout)*252)),
  DCC  = (cumsum(dcc_p_ret)/(1:hout)*252/sqrt(cumsum(dcc_p_ret^2)/(1:hout)*252))
)

sharpe <- sharpe %>% rbind(
  data.frame(
    DATE = tail(ret$DATE[vobs], hout),
    PORTFOLIO = "no short",
    SCOR = (cumsum(our_p0_ret)/(1:hout)*252/sqrt(cumsum(our_p0_ret^2)/(1:hout)*252)),
    DCC  = (cumsum(dcc_p0_ret)/(1:hout)*252/sqrt(cumsum(dcc_p0_ret^2)/(1:hout)*252))
  )
)

sharpe <- sharpe %>% pivot_longer(cols = SCOR:DCC,
                            names_to = "MODEL",
                            values_to = "SHARPE RATIO")

sharpe %>%
  filter(DATE >= "2019-03-01") %>%
  ggplot(aes(x = DATE, y = `SHARPE RATIO`, color = MODEL, linetype = PORTFOLIO)) +
  geom_line() +
  scale_x_date(minor_breaks = as.Date(
      paste0(rep(c(2019, 2020), each = 12), "-", 1:12, "-01")
    )
  ) +
  theme_bw(base_family = "serif")

ggsave("portfolio_sharpe.pdf", width = 7, height = 4)
