
#' The SmartMatrix Class
#'
#' @slot matrix
#' @slot meta.data Contains meta-information 
setClass(
  Class = 'SmartMatrix',
  slots = c(
    matrix = 'matrix',
    row.data = 'data.frame',
    col.data = 'data.frame'
  )
)

SmartMatrix = function(matrix, row.data = NULL, col.data = NULL){
  
  rownames(row.data) = rownames(matrix)
  rownames(col.data) = colnames(matrix)
  
  new("SmartMatrix", 
      matrix = matrix, 
      row.data = row.data,
      col.data = col.data)
}

# "[.Smartmatrix" <- function(x, i, j, ...) {
#   if (missing(x = i) && missing(x = j)) {
#     return(x)
#   }
#   if (missing(x = i)) {
#     i <- NULL
#   } else if (missing(x = j)) {
#     j <- colnames(x = x)
#   }
#   if (is.logical(x = i)) {
#     if (length(i) != nrow(x = x)) {
#       stop("Incorrect number of logical values provided to subset features")
#     }
#     i <- rownames(x = x)[i]
#   }
#   if (is.logical(x = j)) {
#     if (length(j) != ncol(x = x)) {
#       stop("Incorrect number of logical values provided to subset cells")
#     }
#     j <- colnames(x = x)[j]
#   }
#   if (is.numeric(x = i)) {
#     i <- rownames(x = x)[i]
#   }
#   if (is.numeric(x = j)) {
#     j <- colnames(x = x)[j]
#   }
#   return(subset.SmartMatrix(x = x, rows = i, columns = j, ...))
# }

setMethod(
  f = "[",
  signature = "SmartMatrix",
  definition = function(x, i, j, ..., drop = FALSE) {
      # if (missing(x = i) && missing(x = j)) {
      #   return(x)
      # }
      # if (missing(x = i)) {
      #   i <- NULL
      # } else if (missing(x = j)) {
      #   j <- colnames(x = x)
      # }
      # if (is.logical(x = i)) {
      #   if (length(i) != nrow(x = x)) {
      #     stop("Incorrect number of logical values provided to subset features")
      #   }
      #   i <- rownames(x = x)[i]
      # }
      # if (is.logical(x = j)) {
      #   if (length(j) != ncol(x = x)) {
      #     stop("Incorrect number of logical values provided to subset cells")
      #   }
      #   j <- colnames(x = x)[j]
      # }
      # if (is.numeric(x = i)) {
      #   i <- rownames(x = x)[i]
      # }
      # if (is.numeric(x = j)) {
      #   j <- colnames(x = x)[j]
      # }
    
      x@matrix = x@matrix[i,j]
      x@row.data = x@row.data[i,]
      x@col.data = x@col.data[j,]
    
      return(x)
    }
)

# subset.SmartMatrix = function(x, rows, columns){
#   x@matrix = x@matrix[rows,columns]
#   x@row.data = x@row.data[rows]
#   x@col.data = x@col.data[columns]
#   return(x)
# }


