#' @export
summary_module <- function(dat,at)
{

  dat <- summary_popsumm(dat,at)
  dat <- summary_misc(dat,at)

  return(dat)
  
}


