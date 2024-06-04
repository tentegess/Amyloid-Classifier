# install.packages("e1071")
# install.packages("tidyverse")
# install.packages("caret")
# install.packages("here")
library(readxl)
library(dplyr)
library(caret)
library(e1071)
library(here)

# Load the data
data <- read_excel(here("waltzdb_export.xlsx"))
data <- select(data, Sequence, Classification, `TEM Staining`, `Th-T Binding`, WALTZ, TANGO)

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

# Encode Classification as binary
data$Classification <- ifelse(data$Classification == "amyloid", 1, 0)
data$Classification <- as.factor(data$Classification)

# One-hot encoding for "TEM Staining" using model.matrix
data <- cbind(data, model.matrix(~ `TEM Staining` - 1, data = data))

# Normalize numerical features
numeric_columns_to_normalize <- c("WALTZ", "TANGO")  # Add other numeric columns as needed
data[numeric_columns_to_normalize] <- scale(data[numeric_columns_to_normalize])

# Example of a simple numerical encoding for the 'Sequence'
# encode_sequence <- function(sequence) {
#   sapply(sequence, function(x) {
#     sum(utf8ToInt(x) - utf8ToInt("A") + 1)
#   })
# }
#data$Encoded_Sequence <- encode_sequence(data$Sequence)

amino_acids <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

encode_sequence <- function(sequence, amino_acids) {
  # Initialize a matrix of zeroes
  encoding <- matrix(0, nrow = 1, ncol = length(amino_acids), dimnames = list(NULL, amino_acids))
  
  # Fill the matrix: set 1 in the column corresponding to each amino acid in the sequence
  amino_acids_in_seq <- strsplit(sequence, "")[[1]]
  valid_aas <- amino_acids_in_seq[amino_acids_in_seq %in% amino_acids]
  encoding[1, valid_aas] <- 1
  
  return(encoding)
}

encoded_data <- do.call(rbind, lapply(data$Sequence, function(seq) encode_sequence(seq, amino_acids)))
encoded_data <- as.data.frame(encoded_data)
data <- cbind(data, encoded_data)
data <- data[, !(names(data) %in% 'Sequence')]
# Split data into training and testing sets
set.seed(123)
training_samples <- createDataPartition(data$Classification, p = 0.8, list = FALSE)
train_data <- data[training_samples, ]
test_data <- data[-training_samples, ]


svm_model <- svm(Classification ~ ., data = train_data, type = 'C-classification', kernel = 'radial')

summary(svm_model)

predictions <- predict(svm_model, test_data)

conf_matrix <- confusionMatrix(as.factor(predictions), as.factor(test_data$Classification))

# kroswalidacja
tune_grid <- expand.grid(
  sigma = seq(0.01, 0.1, length = 10),  
  C = 2^seq(-1, 1, length = 10))

train_control <- trainControl(
  method = "cv",          
  number = 10,     
  search = "grid")

svm_cv <- train(
  Classification ~ .,
  data = data,
  method = "svmRadial",
  metric = "Accuracy",
  trControl = train_control,
  tuneGrid = tune_grid,
  preProcess = c("center", "scale"),
  tuneLength = 10)


#wypisanie wyników
print(conf_matrix)

print(svm_cv)