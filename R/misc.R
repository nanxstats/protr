# Convert string to character vectors

s2c = function (string) {
  if(is.character(string) && length(string) == 1) {
    return(.Call('s2c', string))
    } else {
      warning('Wrong argument type in s2c(), NA returned')
      return(NA)
    }
}

