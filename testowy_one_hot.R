library(readxl)
library(dplyr)
library(caret)
library(e1071)
library(here)
library(randomForest)
library(doParallel)

# Rejestrujemy klaster równoległy
numCores <- detectCores()
registerDoParallel(cores = numCores - 1)

# Wczytujemy dane
data <- read_excel(here("waltzdb_export.xlsx"))
data <- select(data, Sequence, Classification)
data <- data[nchar(data$Sequence) == 6, ]

data$Classification <- factor(data$Classification)
data_basic <- data

amino_acids <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

# Funkcja do kodowania sekwencji
encode_sequence <- function(sequence, amino_acids) {
  encoding <- matrix(0, nrow = 1, ncol = length(amino_acids), dimnames = list(NULL, amino_acids))
  amino_acids_in_seq <- strsplit(sequence, "")[[1]]
  valid_aas <- amino_acids_in_seq[amino_acids_in_seq %in% amino_acids]
  encoding[1, valid_aas] <- 1
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

for (i in 1:200) {
  column_physic <- sapply(data_basic$Sequence, function(seq)
    sum_index_values(seq, names(aa_index1_dto[[i]]$I), aa_index1_dto[[i]]$I))

  column_physic_df <- as.data.frame(column_physic)
  names(column_physic_df) <- as.character(aa_index1_dto[[i]]$H)
  data <- cbind(data, column_physic_df)
}

# for (i in seq_along(aa_index1_dto)) {
#   column_physic <- sapply(data_basic$Sequence, function(seq) 
#     sum_index_values(seq, names(aa_index1_dto[[i]]$I), aa_index1_dto[[i]]$I))
#   
#   column_physic_df <- as.data.frame(column_physic)
#   names(column_physic_df) <- as.character(aa_index1_dto[[i]]$H)
#   data <- cbind(data, column_physic_df)
# }


# Podział danych na zestawy treningowe i testowe
set.seed(123)
training_samples <- createDataPartition(data$Classification, p = 0.8, list = FALSE)
train_data <- data[training_samples, ]
test_data <- data[-training_samples, ]

# Trenowanie modelu SVM
train_control <- trainControl(method = "cv", number = 10, search = "grid", allowParallel = TRUE)
tune_grid <- expand.grid(sigma = seq(0.01, 0.1, length = 10), C = 2^seq(-1, 1, length = 10))

svm_cv <- train(Classification ~ ., data = train_data, method = "svmRadial", metric = "Accuracy", trControl = train_control, tuneGrid = tune_grid, preProcess = c("center", "scale"), tuneLength = 10)

# Podsumowanie modelu SVM
summary(svm_cv)

# Predykcje i macierz pomyłek
predictions <- predict(svm_cv, test_data)
conf_matrix <- confusionMatrix(as.factor(predictions), as.factor(test_data$Classification))

# Kroswalidacja
svm_model <- svm(Classification ~ ., data = train_data, type = 'C-classification', kernel = 'radial')
print(conf_matrix)

# Utworzenie kontroli treningowej dla RFE
ctrl <- rfeControl(functions=rfFuncs, method="cv", number=10, allowParallel = TRUE)

# Wybierz kolumny od 21 do końca
selected_data <- train_data[, 22:ncol(train_data)]

# Wykonanie RFE
results <- rfe(selected_data, train_data$Classification, sizes=c(1:ncol(selected_data)), rfeControl=ctrl)

# Wyświetlenie wyników RFE
print(results)


# Wybór optymalnej liczby zmiennych
optimal_vars <- predictors(results)
print(optimal_vars)

# Zakończenie równoległego przetwarzania
stopImplicitCluster()
