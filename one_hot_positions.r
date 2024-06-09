library(readxl)
library(dplyr)
library(caret)
library(e1071)
library(here)
library(randomForest)
library(doParallel)

data <- read_excel(here("waltzdb_export.xlsx"))
data <- select(data, Sequence, Classification)
data <- data[nchar(data$Sequence) == 6, ]

amino_acids <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

#one hot encoder z uwzględnieniem pozycji (format A1-6, R1-6 etc)
# encode_sequence <- function(sequence, amino_acids) {
#   if (nchar(sequence) != 6) {
#     stop("Sekwencja musi zawierać 6 aminokwasów")
#   }
  
#   #tworzenie dataframe'a
#   encoding <- matrix(0, nrow = 1, ncol = length(amino_acids) * 6,
#                      dimnames = list(NULL, paste0(rep(amino_acids, each = 6), rep(1:6, times = length(amino_acids)))))
  
#   #Uzupełnianie dataframe'a na podstawie sekwencji
#   amino_acids_in_seq <- strsplit(sequence, "")[[1]]
#   for (i in seq_along(amino_acids_in_seq)) {
#     aa <- amino_acids_in_seq[i]
#     if (aa %in% amino_acids) {
#       colname <- paste0(aa, i)
#       encoding[1, colname] <- 1
#     }
#   }
  
#   return(encoding)
# }

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

# Kodowanie danych
encoded_data <- do.call(rbind, lapply(data$Sequence, function(seq) encode_sequence(seq, amino_acids)))
encoded_data <- as.data.frame(encoded_data)
data <- cbind(data, encoded_data)