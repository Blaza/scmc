mpoly2function <- function(x, varorder = vars(x), compile_c = FALSE){
  ## argument checking
  stopifnot(is.character(varorder))
  stopifnot(is.mpoly(x))

  if(!setequal(varorder, vars(x))){
    stop("varorder must contain all of the variables of x.",
         call. = FALSE)
  }

  ## deal with constant polynomials
  if(is.constant(x)) return( function(.) unlist(x)[["coef"]] )

  ## general polynomials as a bunch of arguments
  mpoly_string <- print(x, stars = TRUE, silent = TRUE)

  #### Extracting used degrees of vars which should be precomputed
  varnames <- vars(x)
  degrees <- apply(do.call(rbind, exponents(x)), 2, max)

  init_block <- ""

  # if there are several degrees of a variable, we'll create its matrix of
  # precalculated pows, so we need to update the mpoly_string to reflect that
  for (var in varnames) {
    if (degrees[var] > 1) {

      if (!compile_c) {
        var_pows_init <- glue("{var}_pows <- matrix(numeric(length({var}) * {degrees[var]}),
                              ncol = {degrees[var]})\n{var}_pows[,1] <- {var}\n\n")
      } else {
        var_pows_init <- glue("double {var}_pows[{degrees[var]}];{var}_pows[0] = {var}[i];")
      }

      for (deg in 2:degrees[var]) {
        if (!compile_c) {
          var_pows_init <- glue(var_pows_init, "{var}_pows[,{deg}] <- {var}*{var}_pows[,{deg-1}]\n\n")
          mpoly_string <- gsub(paste0(var,"\\*\\*",deg),
                               sprintf("%s_pows[, %s]", var, deg),
                               mpoly_string)

        } else {
          var_pows_init <- glue(var_pows_init, "{var}_pows[{deg-1}] = {var}[i]*{var}_pows[{deg-2}];")
          mpoly_string <- gsub(paste0(var,"\\*\\*",deg),
                               sprintf("%s_pows[%s]", var, deg-1),
                               mpoly_string)
        }
      }

      init_block <- glue(init_block, var_pows_init, "\n\n")
    }
  }

  #### Generating function from mpoly_string and adding precomputing code for new vars
  ##### R version
  if (!compile_c) {
    fun_args <- paste(varorder, collapse = ", ")

    fun_string <- glue("
      function({fun_args}) {{
        {init_block}
        {mpoly_string}
      }}
    ")
    return(eval(parse(text = fun_string)))
  }

  ##### Rcpp version
  # we need the indexed version of mpoly_string to use in for loop in C
  # we add spaces to beginning and end of mpoly_string to ease regex
  mpoly_string <- paste0(" ", mpoly_string, " ")
  for (var in varnames) {
    mpoly_string <- gsub(paste0(" ", var, " "), paste0(" ", var, "[i] "), mpoly_string)
  }

  fun_args <- paste(sapply(varorder, function(var){
    paste0("NumericVector ", var)
  }), collapse = ", ")

  fun_string <- glue("
    NumericVector fun({fun_args}) {{
      int n = {varorder[1]}.size();
      NumericVector result(n);
      for (int i = 0; i < n; i++) {{
        {init_block}
        result[i] = {mpoly_string};
      }}
      return result;
    }}
  ")

  Rcpp::cppFunction(fun_string) # -> fun
  return(fun)
}
