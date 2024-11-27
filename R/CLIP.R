#' Input Data from Clipboard Independent of Platform
#'
#' This function reads data from the clipboard, independent of the operating system platform.
#'
#' @return A data frame containing the data from the clipboard.
#'
#' @details
#' This function checks the operating system platform and reads data from the clipboard accordingly.
#' On macOS, it uses the `pbpaste` command to read data from the clipboard.
#' On Windows and Linux, it reads data from the clipboard using the `read.delim` function with 'clipboard' as the filename.
#'
#' @examples
#' CLIP()
#'
#' @export
CLIP <- function() {
  pasted <- character()
  if (Sys.info()[1] == 'Darwin') #macOS
    pasted <- read.table(pipe('pbpaste'), head = FALSE)
  else 
    pasted <- read.delim('clipboard', head = FALSE) #Windows and Linux. Other platforms untested.
  pasted
}
