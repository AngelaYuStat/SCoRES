# Simulation for Multiplier Bootstrap
library(fda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(refund)

library(SCoRES)
data(pupil)
pupil_fpca <- SCoRES::prepare_pupil_fpca(pupil)
fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
                        s(seconds, by = use, k=30, bs = "cr") +
                        s(seconds, by = age, k = 30, bs = "cr") +
                        s(seconds, by = gender, k = 30, bs = "cr") +
                        s(id, by = Phi1, bs="re") +
                        s(id, by = Phi2, bs="re") +
                        s(id, by = Phi3, bs="re") +
                        s(id, by = Phi4, bs="re"),
                      method = "fREML", data = pupil_fpca, discrete = TRUE)

results_pupil_cma <- SCoRES::SCB_functional_outcome(
  data_df = pupil,
  object = fosr_mod,
  method = "cma",
  fitted = TRUE,
  alpha = 0.05,
  outcome = "percent_change",
  domain = "seconds",
  subset = c("use = 1"),
  id = "id")

results_pupil_multiplier <- SCoRES::SCB_functional_outcome(
  data_df = pupil,
  object = fosr_mod,
  method = "multiplier",
  fitted = TRUE,
  alpha = 0.05,
  outcome = "percent_change",
  domain = "seconds",
  subset = c("use = 1"),
  id = "id")

results <- list(mu_hat = results_pupil_cma$mu_hat,
                scb_low1 = results_pupil_cma$scb_low,
                scb_up1 = results_pupil_cma$scb_up,
                scb_low2 = results_pupil_multiplier$scb_low,
                scb_up2 = results_pupil_multiplier$scb_up,
                domain = results_pupil_cma$domain)

plot_data <- data.frame(
  x = results$domain,
  mu = results$mu_hat,
  scb_low1 = results$scb_low1,
  scb_up1 = results$scb_up1,
  scb_low2 = results$scb_low2,
  scb_up2 = results$scb_up2
)

ggplot(plot_data, aes(x = x)) +
  # 绘制均值函数曲线
  geom_line(aes(y = mu, color = "Estimated Mean"), linewidth = 1.2) +

  # 绘制第一条置信带的上边界（CMA）
  geom_line(aes(y = scb_up1, color = "CMA Upper Bound"), linewidth = 0.8, linetype = "solid") +

  # 绘制第一条置信带的下边界（CMA）
  geom_line(aes(y = scb_low1, color = "CMA Lower Bound"), linewidth = 0.8, linetype = "solid") +

  # 绘制第二条置信带的上边界（Multiplier）
  geom_line(aes(y = scb_up2, color = "Multiplier Upper Bound"), linewidth = 0.8, linetype = "dashed") +

  # 绘制第二条置信带的下边界（Multiplier）
  geom_line(aes(y = scb_low2, color = "Multiplier Lower Bound"), linewidth = 0.8, linetype = "dashed") +

  # 自定义颜色和线型
  scale_color_manual(
    name = NULL,
    values = c(
      "Estimated Mean" = "black",
      "CMA Upper Bound" = "blue",
      "CMA Lower Bound" = "blue",
      "Multiplier Upper Bound" = "red",
      "Multiplier Lower Bound" = "red"
    ),
    breaks = c("Estimated Mean", "CMA Upper Bound", "Multiplier Upper Bound"),
    labels = c("Estimated Mean", "CMA Bootstrap", "Multiplier Bootstrap")
  ) +

  # 添加标题和轴标签
  labs(
    title = "Comparison of Simultaneous Confidence Bands",
    x = "Domain",
    y = "Value"
  ) +

  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

library(tidyverse)
library(splines)

set.seed(235324523)

I = 100 # number of subjects
D = 50 # number of timepoints per subject


t = seq(0, 1, length.out = D)

## simulating from this model
# Y(t) = beta_0(t) + beta_1(t) X_i + b_i(t) + epsilon_i(t)


## define covariate (I'm using a binary one)

n_simulations = 30 # 模拟次数

cover_matrix <- matrix(NA, nrow = n_simulations, ncol = D)

for(sim in 1:n_simulations){
X = rbinom(I, 1, 0.6)

## set fixed effects
## I decided to let beta0 = 0, you can look at coverage of beta1
beta0 = 0
beta1 <- sin(2 * pi * 3 * t)


#plot(t, beta1)

## set random elements
t_basis = bs(t, df = 5)
bi_coef = matrix(rnorm(5 * I), nrow = I, ncol = 5)

bi = bi_coef %*% t(t_basis) # random effects


epsilon = matrix(rnorm(I * D, sd = 0.5), nrow = I, ncol = D)


## define Y
xb = matrix(NA, nrow = I, ncol = D)
for(i in 1:I){
  xb[i,] = beta0 + X[i] * beta1
}


Y = xb + bi + epsilon



# store outcome and covariates as dataframe
simdat = tibble(id = rep(1:I, each = D),
                t = rep(1:D, times = I),
                Y = as.vector(t(Y)),
                X = rep(X, each = D))

# 假设您的数据框为 df，包含以下变量：
# Y: 响应变量
# t: 时间变量
# X: 预测变量
# subject_id: 个体ID

mean_mod <- mgcv::gam(
  Y ~ s(t, k = 7, bs = "cr") +
    s(t, by = X, k = 7, bs = "cr"),
  data = simdat, method = "REML"
)

# Prepare residuals
resid_df <- simdat %>%
  filter(!is.na(Y)) %>%
  select(id, t) %>%
  mutate(resid = mean_mod$residuals) %>%
  pivot_wider(
    names_from = t,
    values_from = resid,
    names_prefix = "resid."
  )

resid_mat <- as.matrix(resid_df[, -1])
rownames(resid_mat) <- resid_df$id

# FPCA estimation
fpca_results <- fpca.face(
  Y = resid_mat,
  argvals = unique(simdat$t),
  knots = 15
)

# Create output dataset
eigenfunctions <- as.data.frame(fpca_results$efunctions)
colnames(eigenfunctions) <- paste0("Phi", seq(1, fpca_results$npc))
eigenfunctions$t <- unique(simdat$t)

simdat <- simdat %>%
  left_join(eigenfunctions, by = "t") %>%
  as_tibble() %>%
  arrange(id, t) %>%
  mutate(id = factor(id))

model <- mgcv::bam(Y ~ s(t, k=7, bs="cr") +
                       s(t, by = X, k=7, bs = "cr") +
                       s(id, by = Phi1, bs="re") +
                       s(id, by = Phi2, bs="re")+
                       s(id, by = Phi3, bs="re") +
                       s(id, by = Phi4, bs="re"),
                       method = "fREML", data = simdat, discrete = TRUE)
scb <- SCB_functional_outcome(
  simdat,
  object = model,
  method = "multiplier",
  outcome = "Y",
  domain = "t",
  fitted = FALSE,
  subset = c("X = 1"),
  id = "id"
)
covered <- (scb$scb_low <= beta1) & (beta1 <= scb$scb_up)
print(mean(covered, na.rm = TRUE))

cover_matrix[sim, ] <- covered
}
coverage_by_time <- colMeans(cover_matrix, na.rm = TRUE)
coverage_by_time
final_coverage <- mean(coverage_by_time)
final_coverage

# 快速诊断k=5和k=10模型的区别
model_k5 <- bam(Y ~ s(t, k = 5, bs = "cr") +
                  s(t, by = X, k = 5, bs = "cr"),
                method = "fREML", data = simdat)

model_k10 <- bam(Y ~ s(t, k = 10, bs = "cr") +
                   s(t, by = X, k = 10, "cr"),
                 method = "fREML", data = simdat)

# 比较拟合优度
cat("k=5 AIC:", AIC(model_k5), "\n")
cat("k=10 AIC:", AIC(model_k10), "\n")

# 比较参数估计
cat("k=5 edf:", sum(model_k5$edf), "\n")
cat("k=10 edf:", sum(model_k10$edf), "\n")

coverage_by_time <- colMeans(cover_results, na.rm = TRUE)
# 计算总体覆盖率
overall_coverage <- mean(cover_results, na.rm = TRUE)

  scb <- SCB_functional_outcome(
    data_long,
    object = gam_mod,
    method = "multiplier",
    outcome = "y",
    domain = "time",
    subset = c("group = 1"),
    id = "id"
  )

  # 计算当前模拟的覆盖率
  for (t in 1:n_time) {
    coverage_results[sim, t] <- mean(
      Y_matrix[, t] >= scb$scb_low[t] &
        Y_matrix[, t] <= scb$scb_up[t]
    )
  }
}

# 4. 计算最终覆盖率
final_coverage <- colMeans(coverage_results)

# 5. 输出结果
cat("各时间点平均SCB覆盖率:", round(mean(final_coverage)*100, 2), "%\n")

# 可视化各时间点覆盖率
plot(time_points, final_coverage, type = "l", ylim = c(0.8, 1),
     xlab = "Time", ylab = "Coverage Rate",
     main = "Group=1的SCB覆盖率（500次模拟）")
abline(h = 0.95, col = "red", lty = 2)
grid()

# 2. 初始化覆盖率存储
coverage_results <- matrix(NA, nrow = 15, ncol = n_time)
# 将数据转换为矩阵
Y_matrix <- matrix(group1_data$beta_hat, nrow = n_subjects, byrow = TRUE)
# 3. 进行500次模拟
for (sim in 1:15) {
  # 计算当前模拟的SCB
  scb <- SCB_functional_outcome(
    data_long,
    object = gam_mod,
    method = "multiplier",
    fitted = FALSE,
    outcome = "y",
    domain = "time",
    subset = c("group = 1"),
    id = "id"
  )

  # 计算当前模拟的覆盖率
  for (t in 1:n_time) {
    coverage_results[sim, t] <- mean(
      Y_matrix[, t] >= scb$scb_low[t] &
        Y_matrix[, t] <= scb$scb_up[t]
    )
  }
}

# 4. 计算最终覆盖率
final_coverage <- colMeans(coverage_results)

# 5. 输出结果
cat("各时间点平均SCB覆盖率:", round(mean(final_coverage)*100, 2), "%\n")

