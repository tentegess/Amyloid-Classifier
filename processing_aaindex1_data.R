if (!"seqinr" %in% installed.packages()) install.packages("seqinr")

library('seqinr')
source("aa_index1_dto.R")

process_aaindex <- function(aaindex) {
  full_names <- c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His',
                  'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val')
  one_letter_codes <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L',
                        'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  
  # Tworzenie pustego wektora do przechowywania wyników
  amino_acids_dict <- vector(mode="list", length=length(full_names))
  
  # Przypisywanie pełnych nazw do jednoliterowych kodów
  for (i in 1:length(full_names)) {
    amino_acids_dict[[i]] <- one_letter_codes[i]
  }
  
  # Przypisanie nazw do listy
  names(amino_acids_dict) <- full_names
  
  for(i in 1:length(aaindex)) {
    aaindex[[i]]["A"] <- NULL
    aaindex[[i]]["R"] <- NULL
    aaindex[[i]]["T"] <- NULL
    aaindex[[i]]["J"] <- NULL
    aaindex[[i]]["C"] <- NULL
  }
  
  for (i in seq_along(aaindex)) {
    names(aaindex[[i]]$I) <- amino_acids_dict[names(aaindex[[i]]$I)]
  }
  
  for (i in seq_along(aaindex)) {
    aaindex[[i]]$I[is.na(aaindex[[i]]$I)] <- 0
  }
  
  return(aaindex)
}

data(aaindex)
aa_index1_dto <- process_aaindex(aaindex)
