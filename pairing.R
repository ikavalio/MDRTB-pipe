if (!require("hash")) {
  install.packages("hash")
  library(hash)
}

start_index <- 0
finish_index <- 0

or <- function(x,y) {return(x|y)}

is.resistant <- function(status) {
  return (!is.na(status) && status == 2) # 2 indicates resistance
}


find.resistance.pairs <- function(snps_data, phenotypes, window, misses_border) {
  pairs <- hash() # key is a SNP with lower position
  pairs.df <- data.frame()  
  visited <- rep(FALSE, ncol(snps_data))
  
  snp_index <- 1
  while (snp_index < ncol(snps_data)-window) {
    current_index <- snp_index
    if (visited[current_index]) {
      snp_index <- snp_index + 1
      next
    } else {
      visited[current_index] <- TRUE
    }
    vec1 <- snps_data[,current_index]
        
    pair_index <- current_index + window
    vec2 <- snps_data[,pair_index]
    
    paired <- test.snps(vec1, vec2, phenotypes, misses_border)
  
    if (paired) {
      snp1 <- as.numeric(colnames(snps_data)[current_index])
      snp2 <- as.numeric(colnames(snps_data)[pair_index])
      pairs[snp1] <- snp2
      pairs.df <- rbind(pairs.df, c(snp1, snp2))      
      
      visited[pair_index] <- TRUE
      break
    }
    
    snp_index <- snp_index + 1      
  }
  if (nrow(pairs.df) > 0) {
    colnames(pairs.df) <- c("key1", "key2")
  }
  return(pairs.df)
}

test.snps <- function(snp_vec1, snp_vec2, phenotypes, misses_border) {
  paired <- FALSE
  misses <- 0
  for (i in range(1, length(snp_vec1))) {
    pos1 <- snp_vec1[i]
    pos2 <- snp_vec2[i]
    if ((pos1 | pos2) & xor(pos1, pos2)
        & is.resistant(phenotypes[i])) {
      paired <- TRUE
    } else if (pos1 & pos2
               & is.resistant(phenotypes[i])) {
      if (misses >= misses_border) {
        paired <- FALSE
      }
      misses <- misses + 1
    }
  }
  return(paired)
}

pair.vectors <- function(snp_vec1, snp_vec2, phenotypes) {
  result <- rep(0, length(snp_vec1))
  for (i in range(1, length(snp_vec1))) {
    pos1 <- snp_vec1[i]
    pos2 <- snp_vec2[i]
    if ((pos1 | pos2) & xor(pos1, pos2)
        & is.resistant(phenotypes[i])) {
      result[i] <- 1
    } else {
      result[i] <- pos1 
    }
  }
  return(result)
}
