#' Find or fit new gam if residuals are not already stored in temporary
#' environment
#'
#' @param design Design matrix
#' @param response Column index of response variable
#' @param model_on Vector of column index/indices of predictor variables, or
#' `"null"` for a space-only model
#' @param distribution String describing distribution/family of response
#' @param modeltype String describing model structure: one of "separated",
#' "tensor", or "both". Not relevant when `model_on = "null"`.
#' @param Env temporary environment storing residuals while running
#' [conditional_independencies]
#' @param spatial_str named list with either coords, or sites & nb
#'
#' @keywords internal
#' @returns residuals of the designated model, either via accessing the temporary
#' environment or by calling [fit_custom_gam]
#' @export
find_or_fit <- function(design,
                        response = a, #variable to fit
                        model_on = b, #variable(s) to predict with or "null" for null model
                        distribution = distributions[a], #distribution of a
                        Env = res_store,
                        modeltype = modeltype,
                        spatial_str = space #named list with either coords, or sites & nb
){

  UseMethod("find_or_fit")

}

#' @export
find_or_fit.discrete_dsgn <- function(design,
                                      response = a, #variable to fit
                                      model_on = b, #variable(s) to predict with or "null" for null model
                                      distribution = distributions[a], #distribution of a
                                      Env = res_store,
                                      modeltype = modeltype,
                                      spatial_str = list(sites = sites, nb = neighborhood)
){

  cond_key <- paste(model_on, collapse = ";")
  key <- paste(response, cond_key, sep = "_")

  if(exists(key, envir = Env)){

    #extract from storage
    return(Env[[key]])

  }else{ #otherwise fit the and store residuals for later

    #force all arguments to pass to fit_custom_gam()
    design <- design
    model_on <- model_on
    sites <- spatial_str[["sites"]]
    distribution <- force(distribution)
    nb <- spatial_str[["nb"]]
    modeltype <- modeltype

    #fit and extract residuals
    Env[[key]] <- fit_custom_gam(design = design,
                                 response = response,
                                 model_on = model_on,
                                 distribution = distribution,
                                 modeltype = modeltype,
                                 sites = sites,
                                 nb = nb)

    return(Env[[key]])

  }

}

#' @export
find_or_fit.coord_dsgn <- function(design,
                                   response = a, #variable to fit
                                   model_on = b, #variable(s) to predict with or "null" for null model
                                   distribution = distributions[a], #distribution of a
                                   Env = res_store,
                                   modeltype = modeltype,
                                   spatial_str = list(coords = coords)
){

  cond_key <- paste(model_on, collapse = ";")
  key <- paste(response, cond_key, sep = "_")

  if(exists(key, envir = Env)){

    #extract from storage
    return(Env[[key]])

  }else{ #otherwise fit the and store residuals for later

    #force all arguments to pass to fit_custom_gam()
    design <- design
    model_on <- model_on
    distribution <- force(distribution)
    modeltype <- modeltype
    coords <- spatial_str[["coords"]]

    #fit and extract residuals
    Env[[key]] <- fit_custom_gam(design = design,
                                 response = response,
                                 model_on = model_on,
                                 distribution = distribution,
                                 modeltype = modeltype,
                                 coords = coords)

    return(Env[[key]])

  }

}

