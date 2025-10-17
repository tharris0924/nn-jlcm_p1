#! k=4
i <- 10; n <- 1000
  wkdir <- paste0("NNEM/Example_k4/k=4_Results_n=", n, "_c0.05")
  file <- paste0(wkdir, "/jlcmm_fit", i, n, ".RData")
  load(file)
  # wkdir <- paste0("NNEM/Example_k-3/k=3_Results_n=", n, "_c0.05")
  file <- paste0(wkdir, "/nnem_size", i, n, ".RData")
  load(file)

obj <- jlcmm_fit
nclasses <- ncol(obj$pprob)-2
coefs <- obj$best
coefend <- min(which(grepl('Weibull',names(coefs))))-1
coefs_multilogit <- matrix(coefs[1:coefend],nrow=nclasses-1)
tmpX <- cbind(1,result$x)
linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
exps=exp(cbind(linearval,0));
pi_logit <- exps/rowSums(exps)
          # True mixing proportions
          true_mixing_prop <- function(x) {
              pi1 <- (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4)) * as.numeric(x < 3)
              pi2 <- (1 - (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4))) * as.numeric(x < 0)
              pi3 <- (2 * exp(-0.1 * (x - 6)^4)) / (1 + exp(-0.1 * (x - 6)^4)) * as.numeric(x >= 3)
              pi4 <- 1 - (pi1 + pi2+ pi3)
              return(cbind(pi1, pi2, pi3, pi4))
            }
          
          pi_true <- true_mixing_prop(result$x)
      # library(dplyr)
      pi_true <- true_mixing_prop(result$x)
result$pj <- as.data.frame(result$pj)
#! true pis
df <- data.frame(pi_true, result$x)

i <- 1
tmpdf <- cbind(df[,i], df$result.x,i)
df_long <- data.frame(tmpdf)
for(i in 2:ncol(df)-1){
  tmpdf <- cbind(df[,i], df$result.x,i)
  df_long <- rbind(df_long,tmpdf)
}
# df_long <- df_long|>rename(`Props`=`V1`, `x`=`V2`)

colnames(df_long) <- c("Props","x","i")
df_long$i <- factor(df_long$i)
library(forcats)

# my_factor <- factor(c("A", "B", "C", "A", "B"))

# Recode with forcats (part of tidyverse)
df_long$i <- fct_recode(df_long$i,
                       "Component 1" = "1",
                       "Component 2" = "2",
                       "Component 3" = "3",
                       "Component 4" = "4",
                      # "JLCMM p1"="5",
                      # "JLCMM p2"="6",
                      # "JLCMM p3"="7",
                      # "JLCMM p4"="8",
                      # "NNEM p1"="9",
                      # "NNEM p2"="10",
                      # "NNEM p3"="11",
                      # "NNEM p4"="12"
                    )
library(ggthemes)
library(ggplot2)
library(scales)
cols <- hue_pal()(4)
p1 <- ggplot(data = df_long,aes(x=x,y=Props, col=factor(i))) + geom_line() + theme_minimal() + labs(title = "True Mixing Proportions",
       x = "x", 
       y = latex2exp::TeX("\\pi(x)")) +theme(legend.title=element_blank())+
      scale_color_manual(values = c(cols[1], cols[2], cols[3],cols[4]))
p1
df <- data.frame(pi_logit, result$x)

i <- 1
tmpdf <- cbind(df[,i], df$result.x,(i+4))
df_long <- data.frame(tmpdf)
for(i in 2:ncol(df)-1){
  tmpdf <- cbind(df[,i], df$result.x,(i+4))
  df_long <- rbind(df_long,data.frame(tmpdf))
}
# df_long <- df_long|>rename(`Props`=`V1`, `x`=`V2`)

colnames(df_long) <- c("Props","x","i")
df_long$i <- factor(df_long$i)
library(forcats)

# my_factor <- factor(c("A", "B", "C", "A", "B"))

# Recode with forcats (part of tidyverse)
df_long$i <- fct_recode(df_long$i,
                      #  "True p1" = "1",
                      #  "True p2" = "2",
                      #  "True p3" = "3",
                      #  "True p4" = "4",
                      "Component 1"="5",
                      "Component 2"="6",
                      "Component 3"="7",
                      "Component 4"="8",
                      # "NNEM p1"="9",
                      # "NNEM p2"="10",
                      # "NNEM p3"="11",
                      # "NNEM p4"="12"
                    )
library(ggthemes)
library(ggplot2)
library(scales)
cols <- hue_pal()(4)
p2 <- ggplot(data = df_long,aes(x=x,y=Props, col=factor(i))) + geom_line() + theme_minimal() + labs(title = "JLCMM",
       x = "x", 
       y = latex2exp::TeX("\\pi(x)")) +theme(legend.title=element_blank())+
      scale_color_manual(values = c(cols[1], cols[2], cols[3],cols[4]))

p2



df <- data.frame(result$pj, result$x)

i <- 1
tmpdf <- cbind(df[,i], df$result.x,(i+8))
df_long <- data.frame(tmpdf)
for(i in 2:ncol(df)-1){
  tmpdf <- cbind(df[,i], df$result.x,(i+8))
  df_long <- rbind(df_long,data.frame(tmpdf))
}
# df_long <- df_long|>rename(`Props`=`V1`, `x`=`V2`)

colnames(df_long) <- c("Props","x","i")
df_long$i <- factor(df_long$i)
library(forcats)

# my_factor <- factor(c("A", "B", "C", "A", "B"))

# Recode with forcats (part of tidyverse)
df_long$i <- fct_recode(df_long$i,
                      #  "True p1" = "1",
                      #  "True p2" = "2",
                      #  "True p3" = "3",
                      #  "True p4" = "4",
                      # "JLCMM p1"="5",
                      # "JLCMM p2"="6",
                      # "JLCMM p3"="7",
                      # "JLCMM p4"="8",
                      "Component 1"="9",
                      "Component 2"="10",
                      "Component 3"="11",
                      "Component 4"="12"
                    )
library(ggthemes)
library(ggplot2)
library(scales)
cols <- hue_pal()(4)
p3 <- ggplot(data = df_long,aes(x=x,y=Props, col=factor(i))) + geom_line() + theme_minimal() + labs(title = "NN-JLCM",
       x = "x", 
       y = latex2exp::TeX("\\pi(x)")) +theme(legend.title=element_blank())+
      scale_color_manual(values = c(cols[1], cols[2], cols[3],cols[4]))

p3


library(patchwork)
combined_plot_patchwork <- p1/ p2 / p3 +
  plot_annotation(
    title = "Mixing Proportions for Four Component Model",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )
combined_plot_patchwork

#!k=3
# df <- data.frame(true_mixing_prop(x), pi_logit, result$theta$pj,result$x)
# x <- as.matrix(x)

  i <- 1
  wkdir <- paste0("NNEM/Example_k-3/k=3_Results_n=", n, "_c0.05")
  file <- paste0(wkdir, "/jlcmm_fit", i, n, ".RData")
  load(file)
  wkdir <- paste0("NNEM/Example_k-3/k=3_Results_n=", n, "_c0.05")
  file <- paste0(wkdir, "/nnem_size", i, n, ".RData")
  load(file)

obj <- jlcmm_fit
nclasses <- ncol(obj$pprob)-2
coefs <- obj$best
coefend <- min(which(grepl('Weibull',names(coefs))))-1
coefs_multilogit <- matrix(coefs[1:coefend],nrow=nclasses-1)
tmpX <- cbind(1,result$x)
linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
exps=exp(cbind(linearval,0));
pi_logit <- exps/rowSums(exps)
      true_mixing_prop <- function(x) {
        pi1 <- (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4)) * as.numeric(x < 3)
        pi2 <- (1 - (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4))) *
          as.numeric(x < 0)
        pi3 <- (2 * exp(-0.1 * (x - 6)^4)) /
          (1 + exp(-0.1 * (x - 6)^4)) *
          as.numeric(x >= 3)
        pi4 <- 1 - (pi1 + pi3)
        return(cbind(pi1, pi3, pi4))
      }
      # library(dplyr)
      pi_true <- true_mixing_prop(result$x)
result$pj <- as.data.frame(result$pj)
df <- data.frame(pi_true, pi_logit, result$pj,result$x)

i <- 1
tmpdf <- cbind(df[,i], df$result.x,i)
df_long <- data.frame(tmpdf)
for(i in 2:ncol(df)-1){
  tmpdf <- cbind(df[,i], df$result.x,i)
  df_long <- rbind(df_long,tmpdf)
}
# df_long <- df_long|>rename(`Props`=`V1`, `x`=`V2`)

colnames(df_long) <- c("Props","x","i")
df_long$i <- factor(df_long$i)
library(forcats)

# my_factor <- factor(c("A", "B", "C", "A", "B"))

# Recode with forcats (part of tidyverse)
df_long$i <- fct_recode(df_long$i,
                       "True p1" = "1",
                       "True p2" = "2",
                       "True p3" = "3",
                      "JLCMM p1"="4",
                      "JLCMM p2"="5",
                      "JLCMM p3"="6",
                      "NN-JLCM p1"="7",
                      "NN-JLCM p2"="8",
                      "NN-JLCM p3"="9")
library(ggthemes)
library(ggplot2)
library(scales)
cols <- hue_pal()(3)
ggplot(data = df_long,aes(x=x,y=Props, col=factor(i))) + geom_line() + theme_minimal() + labs(title = "Mixing Proporitons Comparisons for Three Component Model",
       x = "x", 
       y = latex2exp::TeX("\\pi(x)")) +theme(legend.title=element_blank())+
      scale_color_manual(values = c(cols[1], cols[1], cols[1], cols[2], cols[2], cols[2], cols[3], cols[3], cols[3]))

df <- data.frame(true_mixing_prop(x), pi_logit, result$theta$pj,result$x)
# x <- as.matrix(x)

# pis <- 
# cols <- list(unlist(colnames(df)))
# # cols <- rep(list('true one', "true two", "jlcmm one", "jlcmm two", "nnem one", "nnem two"))
# tmpdf <- cbind(df[,i], df$x,i)

#! k=2
  i <- 1
  wkdir <- paste0("NNEM/Example2/k=2_Results_n=", n, "_c0.05")
  file <- paste0(wkdir, "/jlcmm_fit", i, n, ".RData")
  load(file)
  # wkdir <- paste0("NNEM/Example_k-3/k=3_Results_n=", n, "_c0.05")
  file <- paste0(wkdir, "/nnem_size", i, n, ".RData")
  load(file)
obj <- jlcmm_fit
nclasses <- ncol(obj$pprob)-2
coefs <- obj$best
coefend <- min(which(grepl('Weibull',names(coefs))))-1
coefs_multilogit <- matrix(coefs[1:coefend],nrow=nclasses-1)
tmpX <- cbind(1,result$x)
linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
exps=exp(cbind(linearval,0));
pi_logit <- exps/rowSums(exps)
# lines(x, pi_logit[,1])
# lines(x, pi_logit[,2])
# lines(x,sdsds$theta_orig$pj[,1])
# lines(x,sdsds$theta_orig$pj[,2])
# test_class <- t(apply(linearval, 1,function(x){ }))
# pis <- NULL
# true_p <- MSE_p <- numeric(3)
# pis_logit <- pis_nnem <- NULL
# for(k in seq_len(3)){
# MAB_pi <- true_pi <-nnems_pis <- matrix(ncol=B, nrow = 1000)
# for(i in seq_len(B)){
#   # print(k)
# if(!is.null(pij[i][[1]])){
# MAB_pi[,i] <- pij[i][[1]][,k]
# true_pi[,i] <- pij_true[i][[1]][,k]
# nnems_pis[,i] <- pij_nnem[i][[1]][,k]
# }}
# pis <- cbind(pis,rowMeans(true_pi, na.rm = T))
# pis_logit <- cbind(pis_logit, rowMeans(MAB_pi, na.rm = T))
# pis_nnem <- cbind(pis_nnem,  rowMeans(nnems_pis, na.rm = T))
# pi2 <- mean(abs(rowMeans(MAB_pi, na.rm = T) - rowMeans(true_pi, na.rm = T)));
# pi3 <- mean(rowMeans(((MAB_pi) - (true_pi))^2, na.rm = T))
# true_p[k] <- c(pi2)
# MSE_p[k] <- pi3

# }

      # library(dplyr)
      pi_true <- result$theta_orig$pj
result$pj <- as.data.frame(result$pj)
df <- data.frame(pi_true, pi_logit, result$pj,result$x)


# cols <- list(unlist(colnames(df)))
# # cols <- rep(list('true one', "true two", "jlcmm one", "jlcmm two", "nnem one", "nnem two"))
# tmpdf <- cbind(df[,i], df$x,i)

i <- 1
tmpdf <- cbind(df[,i], df$result.x,i)
df_long <- data.frame(tmpdf)
for(i in 2:ncol(df)-1){
  tmpdf <- cbind(df[,i], df$result.x,i)
  df_long <- rbind(df_long,tmpdf)
}
# df_long <- df_long|>rename(`Props`=`V1`, `x`=`V2`)

colnames(df_long) <- c("Props","x","i")
df_long$i <- factor(df_long$i)
library(forcats)

# my_factor <- factor(c("A", "B", "C", "A", "B"))

# Recode with forcats (part of tidyverse)
df_long$i <- fct_recode(df_long$i,
                       "True p1" = "1",
                       "True p2" = "2",
                      "JLCMM p1"="3",
                      "JLCMM p2"="4",
                      "NN-JLCM p1"="5",
                      "NN-JLCM p2"="6")
library(ggthemes)
library(ggplot2)
library(scales)
cols <- hue_pal()(3)
ggplot(data = df_long,aes(x=x,y=Props, col=factor(i))) + geom_line() + theme_minimal() + labs(title = "Mixing Proporitons Comparisons for Two Component Model",
       x = "x", 
       y = latex2exp::TeX("\\pi(x)")) +theme(legend.title=element_blank())+
      scale_color_manual(values = c(cols[1], cols[1], cols[2], cols[2], cols[3], cols[3]))

df <- data.frame(true_mixing_prop(x), pi_logit, result$theta$pj,result$x)

# Overlay true mixture density on the histogram of generated survival times
library(ggplot2)
library(dplyr)

# Use your existing parameters
shapes <- round(shapes_est)
scales <- round(scales_est)

# Calculate mixture proportions from your data
mixture_props <- table(component_assignments) / n
print("True mixture proportions:")
print(mixture_props)

# Function to calculate Weibull density
weibull_density <- function(t, shape, scale) {
  (shape/scale) * (t/scale)^(shape-1) * exp(-(t/scale)^shape)
}

# Function to calculate mixture density
mixture_density <- function(t, shapes, scales, props) {
  density_sum <- 0
  for(k in 1:length(shapes)) {
    density_sum <- density_sum + props[k] * weibull_density(t, shapes[k], scales[k])
  }
  return(density_sum)
}

# Create data frame for plotting
hist_data <- data.frame(times = times, 
                       component = factor(component_assignments),
                       censored = factor(delta, labels = c("Censored", "Observed")))

# Calculate density curve points
time_seq <- seq(0, max(times) * 1.2, length.out = 200)
mixture_curve <- sapply(time_seq, function(t) mixture_density(t, shapes, scales, mixture_props))

# Individual component densities
component_densities <- data.frame()
for(k in 1:4) {
  comp_density <- data.frame(
    time = time_seq,
    density = mixture_props[k] * weibull_density(time_seq, shapes[k], scales[k]),
    component = paste("Component", k)
  )
  component_densities <- rbind(component_densities, comp_density)
}

mixture_curve_df <- data.frame(time = time_seq, density = mixture_curve)
cols <- hue_pal()(4)

# Create the plot with ggplot2
p1 <- ggplot() +
  # Histogram of observed times
  geom_histogram(data = hist_data, aes(x = times, y = after_stat(density), fill = component), 
                 alpha = 0.7, bins = 30, color = "white", size = 0.3) +
  
  # Individual component densities (dashed lines)
  geom_line(data = component_densities, aes(x = time, y = density, color = component), 
            linetype = "dashed", size = 0.8, alpha = 0.8) +
  scale_x_continuous(limits=c(65,110))+
  # True mixture density (solid black line)
  geom_line(data = mixture_curve_df, aes(x = time, y = density), 
            color = "black", size = 1.5, alpha = 0.9) +
  
  labs(
    title = "Histogram of Generated Survival Times with True Mixture Density",
    subtitle = paste("4-component Weibull mixture with proportions:", 
                     paste(round(mixture_props, 3), collapse = ", ")),
    x = "Survival Time",
    y = "Density",
    fill = "True Component",
    color = "Component Density"
  ) +
  
  scale_fill_manual(values = c("1" = cols[1], "2" = cols[2], 
                              "3" = cols[3], "4" = cols[4])) +
  scale_color_manual(values = c("Component 1" = cols[1], "Component 2" = cols[2], 
                               "Component 3" = cols[3], "Component 4" = cols[4])) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(title = "True Component"),
    color = guide_legend(title = "Component Density", override.aes = list(linetype = "dashed"))
  )

print(p1)
# Plot true class-specific survival probabilities
library(ggplot2)
library(dplyr)

# Use your existing parameters
shapes <- round(shapes_est)
scales <- round(scales_est)

print("True Weibull parameters:")
print(data.frame(Class = 1:4, Shape = shapes, Scale = scales))

# Create time sequence for plotting
time_seq <- seq(0, 200, by = 1)

# Function to calculate Weibull survival probability
weibull_survival <- function(t, shape, scale) {
  exp(-(t/scale)^shape)
}

# Calculate survival probabilities for each class
survival_data <- data.frame()

for(class in 1:4) {
  class_survival <- data.frame(
    time = time_seq,
    survival_prob = weibull_survival(time_seq, shapes[class], scales[class]),
    class = paste("Class", class)
  )
  survival_data <- rbind(survival_data, class_survival)
}

# Create the plot
p1 <- ggplot(survival_data, aes(x = time, y = survival_prob, color = class)) +
  geom_line(size = 1.2) +
  labs(
    title = "True Class-Specific Survival Probabilities",
    x = "Time",
    y = "Survival Probability",
    color = "Latent Class"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(50, 110)) +
  scale_color_manual(values = c("Class 1" = cols[1], "Class 2" = cols[2], 
                               "Class 3" = cols[3], "Class 4" = cols[4])) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p1)

# Alternative: Base R plot
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
plot(time_seq, weibull_survival(time_seq, shapes[1], scales[1]), 
     type = "l", col = "blue", lwd = 2,
     ylim = c(0, 1), xlab = "Time", ylab = "Survival Probability",
     main = "True Class-Specific Survival Probabilities")

for(class in 2:4) {
  lines(time_seq, weibull_survival(time_seq, shapes[class], scales[class]), 
        col = c("orange", "green", "red")[class-1], lwd = 2)
}

legend("topright", 
       legend = paste("Class", 1:4, "(λ=", scales, ", κ=", shapes, ")"),
       col = c("blue", "orange", "green", "red"),
       lwd = 2, cex = 0.8)

# Additional analysis: median survival times
median_survival <- numeric(4)
for(class in 1:4) {
  # Median survival time for Weibull: scale * (ln(2))^(1/shape)
  median_survival[class] <- scales[class] * (log(2))^(1/shapes[class])
}

print("Median survival times by class:")
print(data.frame(Class = 1:4, 
                 Shape = shapes, 
                 Scale = scales, 
                 Median_Survival = round(median_survival, 2)))

# Hazard ratios (relative to Class 1)
cat("\nHazard characteristics:\n")
for(class in 1:4) {
  cat("Class", class, "- Shape:", shapes[class], 
      "(", ifelse(shapes[class] > 1, "increasing", 
                  ifelse(shapes[class] == 1, "constant", "decreasing")), 
      "hazard), Scale:", scales[class], "\n")
}


# Combined subplot: Mixture density, survival probabilities, and longitudinal profiles
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# Use your existing parameters
shapes <- round(shapes_est)
scales <- round(scales_est)
cols <- hue_pal()(4)

# Calculate mixture proportions from your data
mixture_props <- table(component_assignments) / n

# Function to calculate Weibull density
weibull_density <- function(t, shape, scale) {
  (shape/scale) * (t/scale)^(shape-1) * exp(-(t/scale)^shape)
}

# Function to calculate mixture density
mixture_density <- function(t, shapes, scales, props) {
  density_sum <- 0
  for(k in 1:length(shapes)) {
    density_sum <- density_sum + props[k] * weibull_density(t, shapes[k], scales[k])
  }
  return(density_sum)
}

# Function to calculate Weibull survival probability
weibull_survival <- function(t, shape, scale) {
  exp(-(t/scale)^shape)
}

# =============================================================================
# PLOT 1: Mixture Density Histogram
# =============================================================================

# Create data frame for plotting
hist_data <- data.frame(times = times, 
                       component = factor(component_assignments),
                       censored = factor(delta, labels = c("Censored", "Observed")))

# Calculate density curve points
time_seq <- seq(0, max(times) * 1.2, length.out = 200)
mixture_curve <- sapply(time_seq, function(t) mixture_density(t, shapes, scales, mixture_props))

# Individual component densities
component_densities <- data.frame()
for(k in 1:4) {
  comp_density <- data.frame(
    time = time_seq,
    density = mixture_props[k] * weibull_density(time_seq, shapes[k], scales[k]),
    component = paste("Component", k)
  )
  component_densities <- rbind(component_densities, comp_density)
}

mixture_curve_df <- data.frame(time = time_seq, density = mixture_curve)

# Create mixture density plot
p1 <- ggplot() +
  geom_histogram(data = hist_data, aes(x = times, y = after_stat(density), fill = component), 
                 alpha = 0.7, bins = 30, color = "white", size = 0.3) +
  geom_line(data = component_densities, aes(x = time, y = density, color = component), 
            linetype = "dashed", size = 0.8, alpha = 0.8) +
  geom_line(data = mixture_curve_df, aes(x = time, y = density), 
            color = "black", size = 1.5, alpha = 0.9) +
  scale_x_continuous(limits = c(65, 110)) +
  labs(
    title = "A) Mixture Density of Survival Times",
    x = "Survival Time",
    y = "Density",
    fill = "Component",
    color = "Component"
  ) +
  scale_fill_manual(values = c("1" = cols[1], "2" = cols[2], 
                              "3" = cols[3], "4" = cols[4])) +
  scale_color_manual(values = c("Component 1" = cols[1], "Component 2" = cols[2], 
                               "Component 3" = cols[3], "Component 4" = cols[4])) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# =============================================================================
# PLOT 2: Survival Probabilities
# =============================================================================

# Calculate survival probabilities for each class
time_seq_surv <- seq(50, 110, by = 1)
survival_data <- data.frame()

for(class in 1:4) {
  class_survival <- data.frame(
    time = time_seq_surv,
    survival_prob = weibull_survival(time_seq_surv, shapes[class], scales[class]),
    class = paste("Class", class)
  )
  survival_data <- rbind(survival_data, class_survival)
}

# Create survival plot
p2 <- ggplot(survival_data, aes(x = time, y = survival_prob, color = class)) +
  geom_line(size = 1.2) +
  labs(
    title = " Class-Specific Survival Probabilities",
    x = "Time",
    y = "Survival Probability",
    color = "Class"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(50, 110)) +
  scale_color_manual(values = c("Class 1" = cols[1], "Class 2" = cols[2], 
                               "Class 3" = cols[3], "Class 4" = cols[4])) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# =============================================================================
# PLOT 3: Longitudinal Profiles (from original code)
# =============================================================================

# Sample subjects for trajectory plot
sample_subjects <- sample(1:n_subjects, min(50, n_subjects))
data_sample <- datalong %>% filter(ID %in% sample_subjects)

p3 <- ggplot(data_sample, aes(x = measure_time, y = y, group = ID, color = factor(g))) +
  geom_line(alpha = 0.7, size = 0.5) +
  geom_point(alpha = 0.5, size = 0.8) +
  labs(
    title = " Sample Longitudinal Trajectories",
    x = "Time",
    y = "Response (Y)",
    color = "Class"
  ) +
  scale_color_manual(values = c("1" = cols[1], "2" = cols[2], 
                               "3" = cols[3], "4" = cols[4]),
                     labels = c("Class 1", "Class 2", "Class 3", "Class 4")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(title = "True Latent Class", ncol = 4))

# =============================================================================
# COMBINE PLOTS
# =============================================================================
library(tidyverse)
# # Create the combined plot
# combined_plot <- grid.arrange(
#   p1, p2, p3,
#   ncol = 2,
#   nrow = 2,
#   layout_matrix = rbind(c(1, 2),
#                        c(3, 3)),
#   top = textGrob("Joint Latent Class Model: True Parameter Visualization", 
#                  gp = gpar(fontsize = 16, fontface = "bold")),
#   bottom = textGrob(paste("4-component mixture | Proportions:", 
#                          paste(round(mixture_props, 3), collapse = ", "),
#                    gp = gpar(fontsize = 10))
# )

# Alternative: Using patchwork (if available)
library(patchwork)
combined_plot_patchwork <- (p2) / p3 +
  plot_annotation(
    title = "Joint Latent Class Survival Probabilities and Trajectories",
    subtitle = paste("4-component model | Proportions of True Classes:", 
                    paste(round(mixture_props, 3), collapse = ", ")),
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )
combined_plot_patchwork
# Print summary information
cat("\n=== TRUE MODEL PARAMETERS ===\n")
cat("Mixture proportions:\n")
print(round(mixture_props, 3))
cat("\nWeibull parameters:\n")
print(data.frame(Class = 1:4, Shape = shapes, Scale = scales))
cat("\nFixed effects (betas):\n")
print(betas_matrix)
cat("\nRandom effects covariance matrix:\n")
print(round(D_g_og, 3))
cat("\nResidual error SD:", residual_err, "\n")

# Save the plot (optional)
ggsave("joint_model_visualization.png", combined_plot, 
        width = 12, height = 10, dpi = 300)
