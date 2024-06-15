library(readxl)
library(dplyr)
library(caret)
library(e1071)
library(here)
library(randomForest)
library(doParallel)
library(stringr)


# Rejestrujemy klaster równoległy
num_cores <- detectCores()
registerDoParallel(cores = num_cores - 1)

# Wczytujemy dane
data <- read_excel(here("waltzdb_export.xlsx"))
data <- select(data, Sequence, Classification)
data <- data[nchar(data$Sequence) == 6, ]

data$Classification <- factor(data$Classification)
data_basic <- data

amino_acids_normal <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H",
                        "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

amino_acids <- c()
for (aa1 in amino_acids_normal) {
  for (aa2 in amino_acids_normal) {
    combination <- paste(aa1, aa2, sep = "")
    amino_acids <- c(amino_acids, combination)
  }
}


# Funkcja do kodowania sekwencji
encode_sequence <- function(sequence) {
  encoding <- matrix(0, nrow = 1, ncol = length(amino_acids) * 3,
                     dimnames = list(NULL, paste0(rep(1:3, each = length(amino_acids)), rep(amino_acids, times = 3))))
  # Rozdzielenie sekwencji na pary aminokwasów
  amino_acids_in_seq <- strsplit(sequence, "(?<=.{2})", perl = TRUE)[[1]]
  for (i in seq_along(amino_acids_in_seq)) {
    aa <- amino_acids_in_seq[i]
    if (aa %in% amino_acids) {
      colname <- paste0(i, aa)
      encoding[1, colname] <- 1
    }
  }
  
  return(encoding)
}

sequence <- data$Sequence[1]
# Kodowanie danych
encoded_data <- do.call(rbind, lapply(data$Sequence, function(seq) encode_sequence(seq)))
encoded_data <- as.data.frame(encoded_data)
data <- cbind(data, encoded_data)
data <- data[, !(names(data) %in% 'Sequence')]

# Dołaczanie danych fizyczno-chemicznych
sum_index_values <- function(sequence, aa_index) {
  # Rozdzielenie sekwencji na pary aminokwasów
  amino_acids_in_seq <- strsplit(sequence, "(?<=.{2})", perl = TRUE)[[1]]
  sum_values <- 0

  for (i in seq_along(amino_acids_in_seq)) {
    # Sumowanie wartości indeksu dla każdego aminokwasu w sekwencji
    sum_values <- sum_values + sum(aa_index[str_sub(amino_acids_in_seq[i], 2, 2), str_sub(amino_acids_in_seq[i], 1, 1)])
  }
  return(sum_values)
}

for (i in seq_along(aa_index_2_3_dto)) {
  column_physic <- sapply(data_basic$Sequence, function(seq)
                            sum_index_values(seq, aa_index_2_3_dto[[i]]$I))

  column_physic_df <- as.data.frame(column_physic)
  names(column_physic_df) <- as.character(aa_index_2_3_dto[[i]]$H)
  data <- cbind(data, column_physic_df)
}

# Podział danych na zestawy treningowe i testowe
set.seed(123)
training_samples <- createDataPartition(data$Classification, p = 0.8, list = FALSE)
train_data <- data[training_samples, ]
test_data <- data[-training_samples, ]

# Trenowanie modelu SVM
train_control <- trainControl(method = "cv", number = 10, search = "grid", allowParallel = TRUE)
tune_grid <- expand.grid(sigma = seq(0.01, 0.1, length = 10), C = 2^seq(-1, 1, length = 10))

# Dodanie dodatkowych cech
train_data_encoded <- cbind(train_data, encoded_data)
test_data_encoded <- cbind(test_data, encoded_data)

svm_cv <- train(Classification ~ ., data = train_data_encoded, method = "svmRadial", metric = "Accuracy", trControl = train_control, tuneGrid = tune_grid, preProcess = c("center", "scale"), tuneLength = 10)

# Podsumowanie modelu SVM
summary(svm_cv)

# Predykcje i macierz pomyłek
predictions <- predict(svm_cv, test_data_encoded)
conf_matrix <- confusionMatrix(as.factor(predictions), as.factor(test_data$Classification))

# Kroswalidacja
svm_model <- svm(Classification ~ ., data = train_data, type = 'C-classification', kernel = 'radial')
print(conf_matrix)

# Utworzenie kontroli treningowej dla RFE
ctrl <- rfeControl(functions=rfFuncs, method="cv", number=10, allowParallel = TRUE)

# Wybierz kolumny od 21 do końca
selected_data <- train_data[, 22:ncol(train_data)]

# Wykonanie RFE
results_list <- list()
cat("Rozpoczęto RFE...\n")
for (i in 1:ncol(selected_data)) {
  cat("Przetwarzanie: ", i, " na ", ncol(selected_data), "\n")
  results <- rfe(selected_data, train_data$Classification, sizes=i, rfeControl=ctrl)
  results_list[[i]] <- results
}
cat("Zakończono RFE.\n")


# Wyświetlenie wyników RFE
print(results_list)

# Wybór optymalnej liczby zmiennych
optimal_vars <- predictors(results_list)
print(optimal_vars)

# Zakończenie równoległego przetwarzania
stopImplicitCluster()
