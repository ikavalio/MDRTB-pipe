
plot.chart.JMeterRes <- function(results, type, ...) {
  .assert(is.function(type), "type must be a valid function")
  options <- list(...)
  Reduce(function(l, r) do.call(type, c(list(result = results, init.plot = l), r)), options)
}
