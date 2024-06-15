source("aa_index_2_3_dto.R")

aa_names <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H",
              "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

process_data <- function(data_text, aa_index) {
  lines <- unlist(strsplit(data_text, "\n"))
  matrix_data <- NULL
  is_matrix <- FALSE
  is_matrix_square <- FALSE
  h_accession <- ""
  i <- 1

  for (line in lines) {
    line <- trimws(line)
    if (startsWith(line, "H ")) {
      h_accession <- sub("H ", "", line)
    } else if (startsWith(line, "M ")) {
      is_matrix <- TRUE
      y <- 1
      matrix_data <- matrix(NaN, nrow = length(aa_names),
                            ncol = length(aa_names))
      next
    } else if (is_matrix) {
      values <- as.numeric(unlist(strsplit(line, " +")))
      if (length(values) == length(aa_names) + 1) {
        if (!is_matrix_square) {
          is_matrix_square <- TRUE
          next
        }
        values <- values[-1]
      } else if (length(values) == length(aa_names) && y == 1) {
        is_matrix_square <- TRUE
      }
      if (!is.null(matrix_data)) {
        for(x in seq_along(values)) {
          if (is.na(values[x])) {
            values[x] <- 0
          }
          matrix_data[y, x] <- values[x]
        }
        y <- y + 1
        if (y > length(aa_names)) {
          if (!is_matrix_square){
            matrix_data[upper.tri(matrix_data)] <-
              t(matrix_data)[upper.tri(matrix_data)]
          }

          is_matrix <- FALSE
          is_matrix_square <- FALSE

          aa_index[[length(aa_index) + 1]] <- list(
            H = h_accession,
            I = matrix(matrix_data, nrow = length(aa_names), 
                       ncol = length(aa_names),
                       dimnames = list(aa_names, aa_names))
          )
          i <- i + 1
        }
      }
    }
  }
  return(aa_index)
}

data_text2 <- readLines("Amyloid-Classifier/aaindex2")
data_text3 <- readLines("Amyloid-Classifier/aaindex3")

aa_index_2_3_dto = list()
# index 2
aa_index_2_3_dto <- process_data(data_text2, aa_index_2_3_dto)

# index 3
aa_index_2_3_dto <- process_data(data_text3, aa_index_2_3_dto)