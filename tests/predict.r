if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
if (!"caret" %in% installed.packages()) install.packages("caret")
if (!"e1071" %in% installed.packages()) install.packages("e1071")
if (!"here" %in% installed.packages()) install.packages("here")
if (!"doParallel" %in% installed.packages()) install.packages("doParallel")
if (!"mltools" %in% installed.packages()) install.packages("mltools")

library(readxl)
library(dplyr)
library(caret)
library(e1071)
library(here)
library(randomForest)
library(doParallel)
library(mltools)
library(pROC)

numCores <- detectCores()
registerDoParallel(cores = numCores - 1)

# Wczytujemy dane
data <- read_excel(here("waltzdb_export.xlsx"))
data <- select(data, Sequence, Classification)
data <- data[nchar(data$Sequence) == 6, ]

data$Classification <- factor(data$Classification)
levels(data$Classification) <- make.names(levels(data$Classification))
data_basic <- data

amino_acids <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

# Funkcja do kodowania sekwencji
encode_sequence <- function(sequence, amino_acids) {
  if (nchar(sequence) != 6) {
    stop("Sekwencja musi zawierać 6 aminokwasów")
  }
  
  #tworzenie dataframe'a
  encoding <- matrix(0, nrow = 1, ncol = length(amino_acids) * 6,
                     dimnames = list(NULL, paste0(rep(amino_acids, each = 6), rep(1:6, times = length(amino_acids)))))
  
  #Uzupełnianie dataframe'a na podstawie sekwencji
  amino_acids_in_seq <- strsplit(sequence, "")[[1]]
  for (i in seq_along(amino_acids_in_seq)) {
    aa <- amino_acids_in_seq[i]
    if (aa %in% amino_acids) {
      colname <- paste0(aa, i)
      encoding[1, colname] <- 1
    }
  }
  
  return(encoding)
}

# Kodowanie danych
encoded_data <- do.call(rbind, lapply(data$Sequence, function(seq) encode_sequence(seq, amino_acids)))
encoded_data <- as.data.frame(encoded_data)
data <- cbind(data, encoded_data)
data <- data[, !(names(data) %in% 'Sequence')]

# Dołaczanie danych fizyczno-chemicznych
sum_index_values <- function(sequence, amino_acids, aa_index) {
  # Rozdzielenie sekwencji na pojedyncze aminokwasy
  amino_acids_in_seq <- strsplit(sequence, "")[[1]]
  
  # Sumowanie wartości indeksu dla każdego aminokwasu w sekwencji
  sum_values <- sum(aa_index[amino_acids_in_seq])
  
  return(sum_values)
}


wartosci_do_wyciagniecia <- c("BAEK050101", "GEIM800107", "QIAN880121", "MANP780101", "PONP930101")
wyniki <- list()

for (wartosc in wartosci_do_wyciagniecia) {
  if (wartosc %in% names(aa_index1_dto)) {
    wyniki[[wartosc]] <- aa_index1_dto[[wartosc]]
  }
}

for (i in seq_along(wyniki)) {
  column_physic <- sapply(data_basic$Sequence, function(seq)
    sum_index_values(seq, names(wyniki[[i]]$I), wyniki[[i]]$I))

  column_physic_df <- as.data.frame(column_physic)
  names(column_physic_df) <- as.character(wyniki[[i]]$H)
  data <- cbind(data, column_physic_df)
}

#Podział danych na zestawy treningowe i testowe
set.seed(123)
training_samples <- createDataPartition(data$Classification, p = 0.8, list = FALSE)
train_data <- data[training_samples, ]
test_data <- data[-training_samples, ]


train_control <- trainControl(
  method = "cv",
  number = 10,
  search = "grid",
  allowParallel = TRUE,
  #classProbs = TRUE, 
  #summaryFunction = twoClassSummary 
)

tune_grid <- expand.grid(sigma = seq(0.01, 0.1, length = 10), C = 2^seq(-1, 1, length = 10))

svm_cv <- train(Classification ~ ., data = train_data, method = "svmRadial", metric = "Accuracy", trControl = train_control, tuneGrid = tune_grid, preProcess = c("center", "scale"), tuneLength = 10)

# Podsumowanie modelu SVM
summary(svm_cv)

# Predykcje i macierz pomyłek
predictions <- predict(svm_cv, test_data)
conf_matrix <- confusionMatrix(as.factor(predictions), as.factor(test_data$Classification), mode = "everything")
print(conf_matrix)
mcc_value <- mcc(preds = predictions, actuals = as.factor(test_data$Classification))
print(paste("MCC:", mcc_value))

# prob_predictions <- predict(svm_cv, test_data, type = "prob")
# roc_curve <- roc(response = test_data$Classification, predictor = prob_predictions[,"amyloid"])
# plot(roc_curve, main = "ROC Curve")
# auc_value <- auc(roc_curve)
# print(paste("AUC:", auc_value))



# set.seed(123) 
# index <- createDataPartition(data$Classification, p = 0.8, list = FALSE)
# train_data <- data[index, ]
# test_data <- data[-index, ]

# train_control <- trainControl(method = "cv", number = 10, search = "grid", allowParallel = TRUE,)

# tune_grid <- expand.grid(.mtry = c(1, sqrt(ncol(data)-1), (ncol(data)-1)/2))

# tuned_model_train <- train(x = train_data[, !names(train_data) %in% "Classification"], y = train_data$Classification, method = "rf",
#                            trControl = train_control, tuneGrid = tune_grid, ntree = 500)

# test_predictions <- predict(tuned_model_train, test_data[, !names(test_data) %in% "Classification"])
# test_confusion_matrix <- confusionMatrix(test_predictions, as.factor(test_data$Classification), mode = "everything")
# print(test_confusion_matrix)
# mcc_value <- mcc(preds = test_predictions, actuals = as.factor(test_data$Classification))
# print(paste("MCC:", mcc_value))

stopImplicitCluster()