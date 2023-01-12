Sys.setlocale("LC_ALL","English")

library(tidyverse)
library(ggplot2)
library(ggfortify)

load("SP100_constituents_5min.Rdata")
names(fulldt) %>% stringr::str_remove(".US") -> names(fulldt)
ret <- fulldt[2:nrow(fulldt), ] / as.matrix(fulldt[1:(nrow(fulldt)-1)]) - 1


# stock pair
stck1 <- "NFLX"
stck2 <- "MSFT"

two <- cbind(ret[, stck1], ret[, stck2])
two <- two[two[, 1] != 0 & two[, 2] != 0, ] # keep only if both not zero

#days <- as.Date(time(two))
library(Rfast)

nsim <- 10000
smpl_cor <- matrix(0, nsim, 2)
sign_cor <- matrix(0, nsim, 2)
trim_cor10 <- matrix(0, nsim, 2)
trim_cor05 <- matrix(0, nsim, 2)
set.seed(202102)
for (i in 1:nsim) {
  days <- logical(nrow(two))
  days[sample.int(nrow(two), nrow(two)/2)] <- TRUE
  smpl_cor[i, 1] <- cor(two[days, 1], two[days, 2])
  smpl_cor[i, 2] <- cor(two[!days, 1], two[!days, 2])
  sign_cor[i, 1] <- sin(mean(sign(two[days, 1]*two[days, 2]))*pi/2)
  sign_cor[i, 2] <- sin(mean(sign(two[!days, 1]*two[!days, 2]))*pi/2)
  trim_cor05[i, 1] <- cor(two[days & rowMaxs(abs(two), T) < 0.05, 1], two[days & rowMaxs(abs(two), T) < 0.05, 2])
  trim_cor05[i, 2] <- cor(two[!days & rowMaxs(abs(two), T) < 0.05, 1], two[!days & rowMaxs(abs(two), T) < 0.05, 2])
  trim_cor10[i, 1] <- cor(two[days & rowMaxs(abs(two), T) < 0.1, 1], two[days & rowMaxs(abs(two), T) < 0.1, 2])
  trim_cor10[i, 2] <- cor(two[!days & rowMaxs(abs(two), T) < 0.1, 1], two[!days & rowMaxs(abs(two), T) < 0.1, 2])
}

# df_diff <- tibble("smpl cor" = smpl_cor[, 1] - smpl_cor[, 2],
#                   "sign cor" = sign_cor[, 1] - sign_cor[, 2],
#                   "trim cor (5%)" = trim_cor05[, 1] - trim_cor05[, 2],
#                   "trim cor (10%)" = trim_cor10[, 1] - trim_cor10[, 2]) %>%
#   pivot_longer(cols = 1:4, names_to = "estimate", values_to = "difference")
# 
# 
# df_diff %>% ggplot(aes(x = difference, fill = estimate)) +
#   geom_density(alpha = 0.5) +
#   theme_bw(base_family = "serif") + ggtitle(paste(stck1, stck2, sep = "-"))

df_estim <- tibble("smpl cor" = as.vector(smpl_cor),
                   "sign cor" = as.vector(sign_cor),
                   "trim cor (5%)" = as.vector(trim_cor05),
                   "trim cor (10%)" = as.vector(trim_cor10)) %>%
  pivot_longer(cols = 1:4, names_to = "estimate", values_to = "value")

df_estim %>% ggplot(aes(x = value, fill = estimate)) +
  geom_density(alpha = 0.5) +
  theme_bw(base_family = "serif") + ggtitle(paste(stck1, stck2, sep = "-"))

ggsave(paste0(stck1, "_", stck2, ".pdf"), width = 4, height = 4)

# Plot one time series
# autoplot(ret[, 100:95], scales = "free_y") + theme_bw()
# ggsave("first_5_stocks.pdf", width = 6, height = 8)
