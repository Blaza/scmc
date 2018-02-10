# function to create the poly string recursively
gen_poly_string <- function(node) {
  if (length(node$children) == 0)
    return(
      glue::glue(switch(as.character(node$coef),
                        "1" = "{node$name}", # we don't want 1 * x, but just x
                        "-1" = "- {node$name}",
                        "{node$coef} * {node$name}"))
    )

  # if we don't have the current node in the poly, just take the children
  if (is.null(node$coef) || node$coef == 0) {
    return(
      glue::glue(ifelse(length(node$children) == 1,
                      "{node$name} * {childsum}", # if only one child don't add ()
                      "{node$name} * ({childsum} )"),
               childsum = paste(sapply(node$children, gen_poly_string), collapse=" + "))
    )
  }

  # if we are at the root node, just add the coef and generate the rest of the poly string, so we
  # don't have 1*(...), as that is more multiplication and thus costly.
  if(node$name == "1") {
    return(
      glue::glue("{node$coef} + {childsum}",
                 childsum = paste(sapply(node$children, gen_poly_string), collapse=" + "))
    )
  }

  glue::glue("{node$name} * ( {node$coef} + {childsum} )",
             childsum = paste(sapply(node$children, gen_poly_string), collapse=" + "))
}

mpoly2function <- function(x, compile_c = FALSE){
  if (!is.matrix(x))
    x <- as.matrix(x)

  # remove zero terms
  x <- x[x[, "coef"] != 0, ]

  # get names of vars
  varnames <- head(colnames(x), -1)

  # Initialise the dataframe which will be converted to a data.tree object
  poly_tree_df <- data.frame(pathString = character(nrow(x)), coef = numeric(nrow(x)),
                             stringsAsFactors = FALSE)

  # populate the data.frame from the mpoly PLEASE EXPLAIN SOMEWHERE
  for (i in seq_along(poly_tree_df[,1])) {
    monom <- x[i, ]
    monom <- monom[c(head(monom, -1) > 0, TRUE)] # get non zero degrees and coef
    poly_tree_df$coef[i] <- monom["coef"]
    if (length(monom) == 1) {
      poly_tree_df$pathString[i] <- "1"
    } else {
      monom_vars <- names(monom)[-which(names(monom)=="coef")]
      varspath <- paste(collapse = "/", sapply(monom_vars, function(v){
        paste(rep(v, monom[v]), collapse = "/")
      }))

      poly_tree_df$pathString[i] <- paste0("1/", varspath)
    }
  }

  # create the data.tree object from the dataframe
  poly_tree <- as.Node(poly_tree_df)

  # generate poly string from tree recursively
  poly_string <- gen_poly_string(poly_tree)
  poly_string <- gsub("\\+\\ \\-", "-", poly_string)

  #### Generating function from mpoly_string and adding precomputing code for new vars
  ##### R version
  if (!compile_c) {
    fun_args <- paste(varnames, collapse = ", ")

    fun_string <- glue("
      function({fun_args}) {{
        {poly_string}
      }}
    ")
    return(eval(parse(text = fun_string)))
  }

  ##### Rcpp version
  # we need the indexed version of mpoly_string to use in for loop in C
  # we add spaces to beginning and end of mpoly_string to ease regex
  poly_string <- paste0(" ", poly_string, " ")
  for (var in varnames) {
    poly_string <- gsub(paste0(" ", var, " "), paste0(" ", var, "[i] "), poly_string)
  }

  fun_args <- paste(sapply(varnames, function(var){
    paste0("NumericVector ", var)
  }), collapse = ", ")

  fun_string <- glue("
    NumericVector fun({fun_args}) {{
      int n = {varnames[1]}.size();
      NumericVector result(n);
      for (int i = 0; i < n; i++) {{
        result[i] = {poly_string};
      }}
      return result;
    }}
  ")

  Rcpp::cppFunction(fun_string) # -> fun
  return(fun)
}
