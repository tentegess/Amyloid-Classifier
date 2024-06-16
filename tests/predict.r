if (!"readxl" %in% installed.packages()) install.packages("readxl")
if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
if (!"caret" %in% installed.packages()) install.packages("caret")
if (!"here" %in% installed.packages()) install.packages("here")
if (!"doParallel" %in% installed.packages()) install.packages("doParallel")
if (!"mltools" %in% installed.packages()) install.packages("mltools")
if (!"stringr" %in% installed.packages()) install.packages("stringr")
if (!"kernlab" %in% installed.packages()) install.packages("kernlab")

library(readxl)
library(dplyr)
library(caret)
library(here)
library(doParallel)
library(mltools)
library(stringr)

numCores <- detectCores()
registerDoParallel(cores = numCores - 1)

# Wczytanie danych
data <- read_excel(here("waltzdb_export.xlsx"))
data <- select(data, Sequence, Classification)
#wybranie tylko sekwencji o długości 6 znaków
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

#index 1
selected_features <- c("BAEK050101", "GEIM800107", "QIAN880121", "MANP780101", "PONP930101")
results <- list()

for (val in selected_features) {
  if (val %in% names(aa_index1_dto)) {
    results[[val]] <- aa_index1_dto[[val]]
  }
}

for (i in seq_along(results)) {
  column_physic <- sapply(data_basic$Sequence, function(seq)
    sum_index_values(seq, names(results[[i]]$I), results[[i]]$I))

  column_physic_df <- as.data.frame(column_physic)
  names(column_physic_df) <- as.character(results[[i]]$H)
  data <- cbind(data, column_physic_df)
}

#index 2 i 3
selected_features <- c("DOSZ010101", "MEHP950103", "MEHP950102", "QU_C930101", "DOSZ010103", "ZHAC000106", "THOP960101", "ZHAC000103", "LIWA970101", "ZHAC000101")
results <- list()

for (index in seq_along(aa_index_2_3_dto)){
  if (aa_index_2_3_dto[[index]]$H %in% selected_features){
    results[[aa_index_2_3_dto[[index]]$H]] <- aa_index_2_3_dto[[index]]$I
  }
}

sum_index_values_index_2_3 <- function(sequence, aa_index) {
  # Rozdzielenie sekwencji na pary aminokwasów
  amino_acids_in_seq <- strsplit(sequence, "")[[1]]
  sum_values <- 0
  # Sumowanie wartości indeksu dla każdego aminokwasu w sekwencji
  for (i in seq(1, length(amino_acids_in_seq) - 1)) {
    sum_values <- sum_values + sum(aa_index[str_sub(amino_acids_in_seq[i]), amino_acids_in_seq[i+1]])
  }
  return(sum_values)
}

for (i in seq_along(results)) {
  column_physic <- sapply(data_basic$Sequence, function(seq)
                            sum_index_values_index_2_3(seq, results[[i]]))

  column_physic_df <- as.data.frame(column_physic)
  names(column_physic_df) <- as.character(names(results)[i])
  data <- cbind(data, column_physic_df)
}

#Podział danych na zestawy treningowe i testowe
set.seed(123)
training_samples <- createDataPartition(data$Classification, p = 0.8, list = FALSE)
train_data <- data[training_samples, ]
test_data <- data[-training_samples, ]

# Dostrajanie modelu za pomocą wyszukiwania po siatce
train_control <- trainControl(
  method = "cv",
  number = 10,
  search = "grid",
  allowParallel = TRUE,
)

tune_grid <- expand.grid(sigma = seq(0.01, 0.1, length = 10), C = 2^seq(-1, 1, length = 10))

# trening modelu
svm_cv <- train(Classification ~ ., data = train_data, method = "svmRadial",
metric = "Accuracy", trControl = train_control,
tuneGrid = tune_grid, preProcess = c("center", "scale"), tuneLength = 10)

# Podsumowanie modelu SVM
summary(svm_cv)

# Predykcje i macierz pomyłek
predictions <- predict(svm_cv, test_data)
conf_matrix <- confusionMatrix(as.factor(predictions), as.factor(test_data$Classification), mode = "everything")
print(conf_matrix)

#obliczenie mcc
mcc_value <- mcc(preds = predictions, actuals = as.factor(test_data$Classification))
print(paste("MCC:", mcc_value))

stopImplicitCluster()