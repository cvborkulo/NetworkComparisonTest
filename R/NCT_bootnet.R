# Silly wrapper for backward competability:
NCT_bootnet <- function(...){
  warning("'NCT_bootnet' is deprecated. Please use 'NCT' instead.")
  NCT(...)
}