
#' Update classification accuracy matrix with experiment results
#'
#' @param where existing prediction accuracy matrix or NULL (will be created automatically)
#' @param name row name
#' @param predicted vector of predicted values (should contain no more than 2 unique values)
#' @param actual vector of actual values (should contain no more than 2 unique values)
#' @param neg.value value that should be interpreted as negative outcome
#' @param pos.value value that should be interpreted as positive outcome
#' @return where matrix with appended experiment results

update.classification.stats <- function(where, name, predicted, actual, neg.value = "1", pos.value = "2") {
  stopifnot(length(predicted) == length(actual))
  stopifnot(length(actual) > 0)
  stopifnot(length(unique(predicted)) < 3)
  stopifnot(length(unique(actual)) < 3)
  
  levels <- c(neg.value, pos.value)
  ac.f <- factor(actual, levels = levels)
  pr.f <- factor(predicted, levels = levels)
  
  stat <- table(ac.f, pr.f) # stat[i, j], i - actual, j - predicted
  update.classification.table(where, name, stat, neg.value, pos.value)
}

#' Update classification accuracy matrix with experiment results (using tabulated frequency table)
#'
#' @param where existing prediction accuracy matrix or NULL (will be created automatically)
#' @param name row name
#' @param rtable 2 x 2 experiment contingency table 
#' @param neg.value value that should be interpreted as negative outcome
#' @param pos.value value that should be interpreted as positive outcome
#' @return where matrix with appended experiment results

update.classification.table <- function(where, name, rtable, neg.value = "1", pos.value = "2") {
  stopifnot(all(dim(rtable) == c(2, 2)))
  
  stat <- rtable # stat[i, j], i - actual, j - predicted
  tp <- stat[pos.value, pos.value]
  fp <- stat[neg.value, pos.value]
  fn <- stat[pos.value, neg.value]
  tn <- stat[neg.value, neg.value]
  
  prec <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  total <- sum(rtable)
  acc <- (tp + tn) / total
  
  row <- list()
  row[["drug"]] <- name
  row[["1"]] <- sum(stat[pos.value,])
  row[["0"]] <- sum(stat[neg.value,])
  row[["1/1"]] <- tp / total
  row[["1/0"]] <- fn / total
  row[["0/1"]] <- fp / total
  row[["0/0"]] <- tn / total
  row[["precission"]] <- prec
  row[["recall"]] <- recall
  row[["f1"]] <- 2 * prec * recall / (prec + recall)
  row[["accuracy"]] <- acc
  
  if(is.null(where)) {
    res <- data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    res <- rbind(where, row)
  }
  res
}

roc.create.data <- function(observed, predicted, pos.value = "2", neg.value = "1", thresholds = sort(predicted)) {
  observed <- factor(observed, levels = c(neg.value, pos.value))
  res <- sapply(thresholds, function(th) {
    p <- ifelse(predicted < th, neg.value, pos.value)
    p <- factor(p, levels = c(neg.value, pos.value))
    stat <- table(observed, p)
    tp <- stat[pos.value, pos.value]
    fp <- stat[neg.value, pos.value]
    fn <- stat[pos.value, neg.value]
    tn <- stat[neg.value, neg.value]
    tpr <- tp / (tp + fn)
    fpr <- fp / (fp + tn)
    
    c(x = fpr, y = tpr)
  })
  d <- data.frame(t(res))
  d <- d[order(d$x, d$y),]
  d
}