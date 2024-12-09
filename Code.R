#Q1

#1.1

# Read the CSV file
water_temp_data <- read.csv("/Users/ravishankar/Library/CloudStorage/OneDrive-DeakinUniversity/Deakin/SIT743/Assignment_2/AIMSWaterTempData.csv", header = FALSE)

# Plot the histogram
hist(water_temp_data[,1], main = "Histogram of Underwater Temperature", xlab = "Temperature", ylab = "Frequency", breaks = 30)

##########
#1.2
##########

# Fit a single Gaussian model
mu <- mean(water_temp_data[,1])
sigma <- sd(water_temp_data[,1])

# Plot the histogram and the single Gaussian density distribution
hist(water_temp_data[,1], main = "Histogram with Single Gaussian Density", xlab = "Temperature", ylab = "Frequency", breaks = 30, prob = TRUE)
curve(dnorm(x, mean = mu, sd = sigma), col = "blue", lwd = 2, add = TRUE)
legend("topright", legend = c("Single Gaussian"), col = c("blue"), lwd = 2)

#Printing MLE estimates

cat("Mean: ",mu)
cat("\n")
cat("Standard Deviation: ",sigma)
cat("\n")


#1.3

# Load the mixtools package
library(mixtools)

# Fit a mixture of Gaussians model with 3 components
mixture_model <- normalmixEM(water_temp_data[,1], k = 3)

# Extract the mixing coefficients, means, and standard deviations
mixing_coefficients <- mixture_model$lambda
means <- mixture_model$mu
standard_deviations <- sqrt(mixture_model$sigma)

# Create a sequence of x-values for plotting the density curves

x <- seq(min(water_temp_data[,1]), max(water_temp_data[,1]), length.out = 100)

# Calculate the density values for each Gaussian component

density_component1 <- dnorm(x, mean = means[1], sd = standard_deviations[1]) * mixing_coefficients[1]

density_component2 <- dnorm(x, mean = means[2], sd = standard_deviations[2]) * mixing_coefficients[2]

density_component3 <- dnorm(x, mean = means[3], sd = standard_deviations[3]) * mixing_coefficients[3]

# Calculate the combined density values

combined_density <- density_component1 + density_component2 + density_component3

# Plot the histogram and the mixture of Gaussians density distributions
hist(water_temp_data[,1], main = "Histogram and Mixture of Gaussians", xlab = "Temperature", col = "lightblue", prob = TRUE, ylim = c(0, max(combined_density)))

lines(x, density_component1, col = "red", lwd = 2)
lines(x, density_component2, col = "green", lwd = 2)
lines(x, density_component3, col = "blue", lwd = 2)
lines(x, combined_density, col = "black", lwd = 2)

# Print the mixing coefficients, means, and standard deviations

cat("Mixing Coefficients:\n")
print(mixing_coefficients)

cat("\nMeans:\n")
print(means)

cat("\nStandard Deviations:\n")
print(standard_deviations)


#1.4

# Plot the log-likelihood values over iterations
plot(mixture_model$all.loglik, main = "Log-Likelihood Plot", xlab = "Iteration", ylab = "Log-Likelihood")


#1.5 - Included in the PDF


#2

#2.1 - 2.6  - Included in the PDF

#2.7

#BiocManager::install("graph")
#BiocManager::install("ggm")
library(graph)
library(ggm)

# Define the Bayesian network structure
dag <- DAG(K ~ T + L, A ~ L + H, V ~ W, S ~ K + A + V)

# Plot the Bayesian network
plotGraph(dag, nodesize = 20, tcltk = FALSE, vc = "green")

# Q2.7 (a)
# W ⊥ K | S
dSep(dag, first = "W", second = "K", cond = "S") # False

# Q2.7 (b)
# H ⊥ V | {K,A}
dSep(dag, first = "H", second = "V", cond = c("K", "A")) # True

# Q2.8 (a)
# Markov blanket of A (Accident type)
library(bnlearn)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("Rgraphviz")

model_net = model2network("[T][L][H][W][K|T:L][A|L:H][V|W][S|K:A:V]")
plot(model_net)

# Markov blanket
mb(model_net, "A")

# Plot the Markov blanket of A in the network
graphviz.plot(model_net, highlight = list(nodes = mb(model_net, "A"), col = "green", fill = "red"))

#2.8 (b) - Included in the PDF

# Q2.8 (c)
# Find the Markov blanket of K (Kilowatt Power) and plot it
mb(model_net, "K")
graphviz.plot(model_net, highlight = list(nodes = mb(model_net, "K"), col = "blue", fill = "yellow"))

#2.9 - Included in the PDF


#Q3

#3.1 - 3.4 - Included in the PDF

#3.5 (a)

# BiocManager::install(c("gRain", "RBGL", "gRbase"))
# BiocManager::install(c("Rgraphviz"))

library("Rgraphviz")
library(RBGL)
#install.packages("gRain")
library(gRbase)
library(gRain)

#Assigning alpha=0.3, beta=0.6, lambda=0.5, gamma=0.2, sigma=0.3

alpha <- 0.3
beta <- 0.6
lambda <- 0.5
gamma <- 0.2
sigma <- 0.3

# CPT for A
A_cpt <- cptable(~ A, values = c(alpha, 1 - alpha), levels = c("0", "1"))

# CPT for B
B_cpt <- cptable(~ B, values = c(gamma, 1 - gamma), levels = c("0", "1"))

# CPT for C
C_cpt <- cptable(~ C | A:B, values = c(0.1, 0.4, 0.5, 0.5, 0.2, 0.3, 0.2, lambda, (0.8 - lambda), 0.1, 0.1, 0.8), levels = c("0", "1", "2"))

# CPT for D
D_cpt <- cptable(~ D | C, values = c(0.7, 0.3, sigma, (1 - sigma), 0.4, 0.6), levels = c("0", "1"))

# CPT for E
E_cpt <- cptable(~ E | C, values = c(0.4, 0.6, beta, (1 - beta), 0.8, 0.2), levels = c("0", "1"))

# Define the network structure
network_structure <- grain(compileCPT(list(A_cpt, B_cpt, C_cpt, D_cpt, E_cpt)))

# Plot the belief network
plot(network_structure$dag)

#Q3.5 (B)

compiled_nodes <- compileCPT(list(A_cpt, B_cpt, C_cpt, D_cpt, E_cpt))

# Probability table for A
print(compiled_nodes$A)

# Probability table for B
print(compiled_nodes$B)

# Probability table for C
print(compiled_nodes$C)

# Probability table for D
print(compiled_nodes$D)

# Probability table for E
print(compiled_nodes$E)

#Q3.5 (C)

#i)
# Marginal distribution of C
querygrain(network_structure, nodes = "C")

#ii)
# Joint distribution of A, B, and C
querygrain(network_structure, nodes = c("A", "B", "C"), type="joint")

#iii)
# Set the evidence for A=0, D=1, and E=0
network_with_evidence <- setEvidence(network_structure, nodes = c("A", "D", "E"), states = c("0", "1", "0"))

# Query the network for P(C=2 | A=0, D=1, E=0)
querygrain(network_with_evidence, nodes = "C")


library(bnlearn)

# Load the alarm dataset
data(alarm)
summary(alarm)

# Create the true network structure
true_network_structure <- paste0("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF][LVF]",
                                 "[STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA][HRSA|ERCA:HR][ANES]",
                                 "[APL][TPR|APL][ECO2|ACO2:VLNG][KINK][MINV|INT:VLNG][FIO2][PVS|FIO2:VALV]",
                                 "[SAO2|PVS:SHNT][PAP|PMB][PMB][SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC]",
                                 "[MVS][VMCH|MVS][VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG]",
                                 "[ACO2|VALV][CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]")

true_dag <- model2network(true_network_structure)
par(mfrow = c(1, 1))
graphviz.plot(true_dag, shape = "circle")

# Q4.1 Learn Bayesian network structures using hill-climbing algorithm with BIC and BDe scores
data_sizes <- c(500, 5000, 10000)

cat("Question 4.1\n")

for (size in data_sizes) {
  alarm_data <- alarm[1:size, ]
  
  # Learn network structure using hill-climbing with BIC score
  hc_bic_model <- hc(alarm_data, score = "bic")
  bic_score_value <- score(hc_bic_model, alarm_data, type = "bic")
  cat(sprintf("BIC score for sample size %d: %f\n", size, bic_score_value))
  graphviz.plot(hc_bic_model, main = paste("Hill-climbing with BIC score (Sample size:", size, ")"), shape = "circle")
  
  # Learn network structure using hill-climbing with BDe score
  hc_bde_model <- hc(alarm_data, score = "bde")
  bde_score_value <- score(hc_bde_model, alarm_data, type = "bde")
  cat(sprintf("BDe score for sample size %d: %f\n", size, bde_score_value))
  graphviz.plot(hc_bde_model, main = paste("Hill-climbing with BDe score (Sample size:", size, ")"), shape = "circle")
}

cat("\n")

# Q4.3 (a) Learn Bayesian network structures using the full dataset with BIC and BDe scores
cat("Question 4.3 (a)\n")

hc_bic_full_model <- hc(alarm, score = "bic")
bic_full_score <- score(hc_bic_full_model, alarm, type = "bic")
cat(sprintf("BIC score for full dataset: %f\n", bic_full_score))
graphviz.plot(hc_bic_full_model, main = "Hill-climbing with BIC score (Full dataset)", shape = "circle")

hc_bde_full_model <- hc(alarm, score = "bde")
bde_full_score <- score(hc_bde_full_model, alarm, type = "bde")
cat(sprintf("BDe score for full dataset: %f\n", bde_full_score))
graphviz.plot(hc_bde_full_model, main = "Hill-climbing with BDe score (Full dataset)", shape = "circle")

cat("\n")

# Q4.3 (b) Compare the learned networks with the true network structure
cat("Question 4.3 (b)\n")

bic_comparison <- compare(true_dag, hc_bic_full_model)
bde_comparison <- compare(true_dag, hc_bde_full_model)

cat("Comparison of BIC network with true network:\n")
print(bic_comparison)
graphviz.compare(true_dag, hc_bic_full_model, main = c("True Network", "BIC Network"))

cat("Comparison of BDe network with true network:\n")
print(bde_comparison)
graphviz.compare(true_dag, hc_bde_full_model, main = c("True Network", "BDe Network"))

cat("\n")

# Q4.3 (c) Fit the data to the BIC network and show CPD table for "HR"
cat("Question 4.3 (c)\n")

bic_fitted_model <- bn.fit(hc_bic_full_model, alarm)
hr_cpd_table <- bic_fitted_model$HR
cat("CPD table for 'HR' variable:\n")
print(hr_cpd_table)

cat("\n")

# Q4.3 (d) Find the probability of P(HR = "HIGH" | BP = "LOW", PRSS = "NORMAL")
cat("Question 4.3 (d)\n")

query_probability <- cpquery(bic_fitted_model, event = (HR == "HIGH"), evidence = (BP == "LOW" & PRSS == "NORMAL"))
cat(sprintf("P(HR = 'HIGH' | BP = 'LOW', PRSS = 'NORMAL'): %f\n", query_probability))


#Q5 - Included in the PDF


#Q6#
# Load necessary libraries
library(dplyr)
library(lubridate)
library(ggplot2)
library(corrplot)
library(readr)
library(bnlearn)  # Add this library for Bayesian network analysis
#install.packages("lavaan")
#library(lavaan)
#install.packages("semStructuralEquations")
#library(semStructuralEquations)

#Reading the crash dataset
crash_data <- read.csv("/Users/ravishankar/Library/CloudStorage/OneDrive-DeakinUniversity/Deakin/SIT743/Assignment_2/Crash_Reporting_-_Drivers_Data.csv")
colnames(crash_data)

columns_for_csv <- c("Crash.Date.Time", "Collision.Type", "Weather", "Surface.Condition", "Light", "Injury.Severity", 
                     "Driver.Distracted.By", "Speed.Limit", "Vehicle.Movement")

# Filter crash data
crash_data <- crash_data %>% 
  filter(Agency.Name == "Montgomery County Police") %>%
  select(columns_for_csv)

# Convert Crash.Date.Time to POSIXct format  
crash_data$Crash.Date.Time <- parse_date_time(crash_data$`Crash.Date.Time`, orders = "%m/%d/%Y %I:%M:%S %p", tz = "UTC")

# Filter dates and convert to Date format
crash_data <- crash_data %>%
  filter(Crash.Date.Time >= as.Date("2023-01-01") & Crash.Date.Time <= as.Date("2023-12-31")) %>%
  mutate(Crash.Date.Time = as.Date(Crash.Date.Time))

# Sort by Crash.Date.Time
crash_data <- crash_data %>% arrange(Crash.Date.Time)

#Reading the Weather CSV file
weather_data <- read.csv("/Users/ravishankar/Library/CloudStorage/OneDrive-DeakinUniversity/Deakin/SIT743/Assignment_2/montgomery_data.csv", header = TRUE)

#Dropping first column (Country Name)
weather_data <- weather_data[, -1]

# Converting imperial units to metric units
weather_data[, c("tempmax", "tempmin", "temp")] <- round((weather_data[, c("tempmax", "tempmin", "temp")] - 32) * 5/9, 2)
weather_data[, c("windspeed", "windgust")] <- round(weather_data[, c("windspeed", "windgust")] * 1.60934, 2)

# Converting datetime column to Date format
weather_data$datetime <- as.Date(weather_data$datetime)

# Merge crash data and weather data
crash_data <- left_join(crash_data, weather_data, by = c("Crash.Date.Time" = "datetime"))

# Inspect the structure of the dataset
str(crash_data)

# Convert categorical variables to factors
factor_variables <- c("Collision.Type", "Weather", "Surface.Condition", "Light", "Injury.Severity", "Driver.Distracted.By", 
                      "Speed.Limit", "Vehicle.Movement")
crash_data[factor_variables] <- lapply(crash_data[factor_variables], as.factor)

# Select 15-25 variables for the analysis
selected_variables <- c("Collision.Type", "Surface.Condition", "Light", "Injury.Severity", "Driver.Distracted.By", 
                        "Speed.Limit", "Vehicle.Movement", "tempmax", "tempmin", "temp", "humidity", "windgust", 
                        "windspeed", "cloudcover", "visibility", "uvindex")

crash_data_subset <- crash_data[, selected_variables]

# Perform exploratory analysis
par(mfrow = c(4, 2), mar = c(2, 2, 2, 2))
for (var in selected_variables) {
  if (is.numeric(crash_data_subset[[var]])) {
    hist(crash_data_subset[[var]], main = paste("Histogram of", var), xlab = var)
  } else {
    barplot(table(crash_data_subset[[var]]), main = paste("Barplot of", var), xlab = var, col = "lightblue")
  }
}

# Summary statistics for numeric variables
summary(crash_data_subset[, sapply(crash_data_subset, is.numeric)])

# Summary statistics for categorical variables
summary(crash_data_subset[, sapply(crash_data_subset, is.factor)])

# Perform Bayesian structure learning and parameter learning
# Example using bnlearn package (you can experiment with different algorithms)

# Continuous variables
continuous_vars <- c("tempmax", "tempmin", "temp", "humidity", "windgust", "windspeed", "cloudcover", "visibility", "uvindex")

# Categorical variables
categorical_vars <- c("Collision.Type", "Surface.Condition", "Light", "Injury.Severity", "Driver.Distracted.By", 
                      "Speed.Limit", "Vehicle.Movement")

# Convert continuous variables to discrete variables (example: using quantile-based discretization)
crash_data_discrete <- crash_data_subset
for (var in continuous_vars) {
  crash_data_discrete[[var]] <- crash_data_discrete[[var]] + runif(length(crash_data_discrete[[var]]), -1e-6, 1e-6)
  crash_data_discrete[[var]] <- cut(crash_data_discrete[[var]], breaks = quantile(crash_data_discrete[[var]], probs = seq(0, 1, 0.25), na.rm = TRUE), include.lowest = TRUE)
}


# Perform structure learning using different algorithms

# 1. Hill-Climbing (HC) algorithm
hc_bn <- hc(crash_data_discrete)

# 2. Tabu search algorithm
tabu_bn <- tabu(crash_data_discrete)

# 3. Max-Min Hill-Climbing (MMHC) algorithm
mmhc_bn <- mmhc(crash_data_discrete)

# 4. Hybrid Parents and Children (HPC) algorithm
hpc_bn <- hpc(crash_data_discrete)


# Plot the learned Bayesian network structures
par(mfrow = c(1, 1))
graphviz.plot(hc_bn, main = "Hill-Climbing")
graphviz.plot(tabu_bn, main = "Tabu")
graphviz.plot(mmhc_bn, main = "MMHC")
graphviz.plot(hpc_bn, main = "HPC")


# Perform parameter learning for each learned structure
hc_bn_fitted <- bn.fit(hc_bn, crash_data_discrete)
tabu_bn_fitted <- bn.fit(tabu_bn, crash_data_discrete)
mmhc_bn_fitted <- bn.fit(mmhc_bn, crash_data_discrete)
hpc_bn_fitted <- bn.fit(hpc_bn, crash_data_discrete)
bdeu_bn_fitted <- bn.fit(bdeu_bn, crash_data_discrete)

# Query the Bayesian networks

# Compute the probability: P(Injury Severity = "POSSIBLE INJURY" | Surface Condition = "WET")
query1 <- cpquery(hc_bn_fitted, event = (Injury.Severity == "POSSIBLE INJURY"), evidence = (Surface.Condition == "WET"))
print(query1)

# Check if "Injury Severity" is conditionally independent of "Speed Limit" given the learned network structure
query2 <- dsep(hc_bn, "Injury.Severity", "Speed.Limit")
print(query2)

# Find the Markov blanket of "Driver Distracted By"
markov_blanket <- mb(hc_bn, "Driver.Distracted.By")
print(markov_blanket)