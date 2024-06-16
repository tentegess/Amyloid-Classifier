library(readxl)
library(dplyr)
library(caret)
library(e1071)
library(here)
library(doParallel)
library(randomForest)

# Rejestrujemy klaster równoległy
numCores <- detectCores()
registerDoParallel(cores = numCores - 1)

# Load the data
data <- read_excel(here("waltzdb_export.xlsx"))
data <- select(data, Sequence, Classification)
data <- data[nchar(data$Sequence) == 6, ]

data$Classification <- factor(data$Classification)

amino_acids <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

#(format A1,R1 etc)
encode_sequence <- function(sequence, amino_acids) {
  if (nchar(sequence) != 6) {
    stop("Each sequence must be exactly 6 amino acids long.")
  }
  
  encoding <- matrix(0, nrow = 1, ncol = length(amino_acids) * 6,
                     dimnames = list(NULL, paste0(rep(1:6, each = length(amino_acids)), rep(amino_acids, times = 6))))
  
  amino_acids_in_seq <- strsplit(sequence, "")[[1]]
  for (i in seq_along(amino_acids_in_seq)) {
    aa <- amino_acids_in_seq[i]
    if (aa %in% amino_acids) {
      colname <- paste0(i, aa) 
      encoding[1, colname] <- 1
    }
  }
  
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

# Utworzenie kontroli treningowej dla RFE
ctrl <- rfeControl(functions=rfFuncs, method="cv", number=10)

# Wykonanie RFE
results <- rfe(train_data[, -1], train_data$Classification, sizes=c(1:ncol(train_data)-1), rfeControl=ctrl)

# Wyświetlenie wyników
print(results)

# Wybór optymalnej liczby zmiennych
optimal_vars <- predictors(results)
print(optimal_vars)


stopImplicitCluster()