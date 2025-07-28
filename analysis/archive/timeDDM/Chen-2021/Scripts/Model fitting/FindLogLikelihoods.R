# Function to recursively print the size of all fields in an object
printSizeOfFields <- function(obj, prefix="") {
  if (is.list(obj) || is.data.frame(obj)) {
    for (name in names(obj)) {
      next_prefix <- paste0(prefix, if (nzchar(prefix)) "$", name)
      printSizeOfFields(obj[[name]], next_prefix)
    }
  } else {
    size <- object.size(obj)
    cat(prefix, ": ", size, " bytes\n", sep="")
  }
}

# Example usage:
# my_object <- list(a = 1:10, b = data.frame(x = 1:100, y = rnorm(100)), c = "Hello")
# printSizeOfFields(my_object)
