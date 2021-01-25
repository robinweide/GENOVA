.main_data <- function(discovery) {
  UseMethod(".main_data")
}

.main_data.genomescore_discovery <- function(discovery) {
  as.data.frame(discovery[[1]])
}

.main_data.domainogram_discovery <- function(discovery) {
  as.data.frame(discovery)
}

.metric_name <- function(discovery) {
  UseMethod(".metric_name")
}

.metric_name.IS_discovery <- function(discovery) {
  "Insulation Score"
}

.metric_name.CS_discovery <- function(discovery) {
  "Compartment Score"
}

.metric_name.DI_discovery <- function(discovery) {
  "Directionality Index"
}

.metric_name.virtual4C_discovery <- function(discovery) {
  "Signal"
}