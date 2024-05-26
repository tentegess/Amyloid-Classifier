# install.packages("e1071")
# install.packages("tidyverse")
# install.packages("caret")
# install.packages("here")
# install.packages("pROC")
library(readxl)
library(dplyr)
library(caret)
library(e1071)
library(here)
library(pROC)

# Load the data
data <- read_excel(here("waltzdb_export.xlsx"))
data <- select(data, Sequence, Classification)

# Separate numeric and categorical columns
numeric_columns <- sapply(data, is.numeric)
categorical_columns <- !numeric_columns

# Impute missing numeric data with median
numeric_data <- data[, numeric_columns]
preProcess_numeric <- preProcess(numeric_data, method = 'medianImpute')
data[, numeric_columns] <- predict(preProcess_numeric, newdata = numeric_data)

# Handling missing categorical data by imputing with the most common value (mode)
categorical_data <- data[, categorical_columns]
for(col in names(categorical_data)) {
    mode_value <- names(which.max(table(categorical_data[[col]])))
    categorical_data[[col]][is.na(categorical_data[[col]])] <- mode_value
}
data[, categorical_columns] <- categorical_data

data$Classification <- factor(data$Classification)
levels(data$Classification) <- make.names(levels(data$Classification))



pssm <- matrix(
  c(
    -0.26, -0.32, -0.27, -0.14, -0.43, -0.22,  # A
    -0.45, -0.41, -0.46, -0.33, -0.52, -0.35,  # R
    -0.40, -0.34, -0.49, -0.27, -0.46, -0.30,  # N
    -0.49, -0.43, -0.56, -0.41, -0.56, -0.36,  # D
    -0.09, -0.21, 0.03, -0.05, -0.17, -0.05,   # C
    -0.37, -0.30, -0.36, -0.34, -0.48, -0.32,  # Q
    -0.51, -0.41, -0.43, -0.30, -0.61, -0.39,  # E
    -0.23, -0.37, -0.46, -0.37, -0.30, -0.33,  # G
    -0.32, -0.26, -0.26, -0.30, -0.35, -0.25,  # H
    -0.06, -0.08, 0.26, 0.09, -0.06, -0.07,    # I
    -0.10, -0.18, 0.02, 0.04, -0.22, -0.13,    # L
    -0.39, -0.45, -0.51, -0.35, -0.59, -0.32,  # K
    -0.17, -0.25, -0.02, -0.10, -0.19, -0.18,  # M
    -0.13, -0.11, 0.05, -0.03, -0.13, -0.11,  # F
    -0.56, -0.38, -0.56, -0.51, -0.42, -0.45,  # P
    -0.37, -0.35, -0.41, -0.30, -0.48, -0.23,  # S
    -0.34, -0.33, -0.28, -0.23, -0.40, -0.23,  # T
    -0.17, -0.17, -0.09, -0.06, -0.12, -0.16,  # W
    -0.23, -0.11, -0.13, -0.06, -0.18, -0.15,  # Y
    -0.05, -0.14, 0.19, 0.14, -0.19, 0.01      # V
  ),
  nrow = 20, ncol = 6, byrow = TRUE
)

rownames(pssm) <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
colnames(pssm) <- 1:6

encode_sequence <- function(sequence, pssm) {
    score <- 0
    for (i in 1:nchar(sequence)) {
        aa <- substring(sequence, i, i)
        if (aa %in% rownames(pssm)) {
            score <- score + pssm[aa, i]
        } else {
            warning(paste("Amino acid", aa, "at position", i, "not in PSSM rownames. Score set to 0."))
            score <- score + 0
        }
    }
    return(score)
}

encoded_data <- do.call(rbind, lapply(data$Sequence, function(seq) encode_sequence(seq, pssm)))
encoded_data <- as.data.frame(encoded_data)
data <- cbind(data, encoded_data)
data <- data[, !(names(data) %in% 'Sequence')]
# Split data into training and testing sets
set.seed(123)
training_samples <- createDataPartition(data$Classification, p = 0.8, list = FALSE)
train_data <- data[training_samples, ]
test_data <- data[-training_samples, ]


# kroswalidacja

train_control <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = "final",
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

nu_svm_model <- list(
  type = "Classification",
  library = "e1071",
  loop = NULL,
  method = "svmNu",
  parameters = data.frame(parameter = c("nu", "sigma", "gamma"), class = rep("numeric", 3), label = c("nu", "sigma", "gamma")),
  grid = function(x, y, len = NULL, search = "grid") {
    expand.grid(sigma = seq(0.05, 0.1, length = 6), nu = seq(0.3, 0.7, length = 5), gamma = seq(0.06, 0.1, length = 5))
  },
  fit = function(x, y, wts, param, lev, last, weights, classProbs) {
    svm(x = x, y = y, probability = classProbs, type = "nu-classification", nu = param$nu, kernel = "radial", gamma = param$gamma)
  },
  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, newdata, probability = TRUE)
  },
  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    preds <- predict(modelFit, newdata, probability = TRUE)
    attr(preds, "probabilities")
  },
  levels = function(x) x$levels,
  tags = c("SVM", "Support Vector Machines", "nu-SVM", "Radial Basis Function", "RBF"),
  sort = function(x) x
)

# Train the nu-SVM model
svm_cv <- train(
  Classification ~ .,
  data = train_data,
  method = nu_svm_model,
  trControl = train_control,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  tuneLength = 10
)

print(svm_cv)

predictions <- predict(svm_cv, newdata=test_data, type = "prob")
prob_scores <- predictions[, "amyloid"]
class_labels <- factor(ifelse(prob_scores > 0.5, "amyloid", "non.amyloid"), levels = c("amyloid", "non.amyloid"))
roc_curve <- roc(response = as.factor(test_data$Classification), predictor = prob_scores)
conf_matrix <- confusionMatrix(class_labels, as.factor(test_data$Classification))
print(conf_matrix)
print(paste("AUC:", auc(roc_curve)))
plot(roc_curve)