#' Fit spatial generalized additive models for conditional independence testing
#'
#' @param design Design matrix
#' @param response Column index of response variable
#' @param model_on Vector of column index/indices of predictor variables, or
#' `"null"` for a space-only model
#' @param distribution String describing distribution/family of response
#' @param modeltype String describing model structure: one of "separated",
#' "tensor", or "both". Not relevant when `model_on = "null"`.
#' @param ... other arguments specific to spatial structure. See details.
#'
#' @section Fitting discrete space models:
#' When a `discrete_dsgn` object is passed to `fit_custom_gam`, it fits a
#' \link[mgcv:smooth.construct.mrf.smooth.spec]{Gaussian Markov random field}
#' which is either:
#' * added to thin plate splines of variables indicated by the `model_on`
#' argument (`modeltype = "separated"`)
#'
#' * in a tensor product with variables indicated by the `model_on`
#' argument (`modeltype = "tensor"`)
#'
#' * both of the above (`modeltype = "both"`)
#'
#' The following arguments must also be provided for discrete space models:
#'
#' * `sites` String naming the column of `design` containing site names
#'
#' * `nb` List defining neighborhood structure. See the \link[mgcv:smooth.construct.mrf.smooth.spec]{smooth
#' specifications} for more details. See [build_neighborhood] to construct.
#'
#' @section Fitting coordinate space models:
#' When a `coord_dsgn` object is passed to `fit_custom_gam`, it fits a
#' 2-dimensional \link[mgcv:smooth.terms]{thin plate spline} (`bs = "tp"`), either additively,
#' in a tensor product, or both (see description of model types in above section).
#'
#' The following argument must also be provided for coordinate space models:
#'
#' * `coords` A vector of strings naming the two columns containing coordinate
#' values in `design`
#'
#'
#' @returns Vector of raw residuals
#' @export
#'
#' @examples
#'
#' #Provide a coord_dsgn object
#' head(example_coord_data)
#'
#' #Fit spatial model of Agr ~ Agr_tminus1 + Agr_tminus2 and return residuals
#' resid_543 <- fit_custom_gam(design = example_coord_data,
#'                 response = 1,
#'                 model_on = c(2, 3),
#'                 modeltype = "separated",
#'                 distribution = "betar",
#'                 coords = c("x", "y"))
#'
#' resid_543
#'
fit_custom_gam <- function(design,
                           response = a, #variable to fit
                           model_on = b, #variable(s) to predict with or "null" for null model
                           distribution = distribution, #distribution of a
                           modeltype = modeltype,
                           ...
){

  UseMethod("fit_custom_gam")

}

#' @export
fit_custom_gam.discrete_dsgn <- function(design,
                                         response = a, #variable to fit
                                         model_on = b, #variable(s) to predict with or "null" for null model
                                         distribution = distribution, #distribution of a
                                         modeltype = modeltype,
                                         sites = sites, #name of sites variable
                                         nb = nb #and neighborhood
){

  nb <- force(nb)
  if (response %in% model_on) {
    stop("Self-model detected: response appears in predictors")
  }

  if("null" %in% model_on){

    #Create MRF smooth
    formula <- paste0("s(", sites, ", bs = 'mrf', xt = list(nb = nb))")
    #Connect to response
    formula <- paste(colnames(design)[response], formula, sep = " ~ ")
    #Format
    formula <- as.formula(formula)

    fam_fun <- get(distribution)

    null_model <- suppressWarnings(gam(formula, family = fam_fun(),
                                       data = design,
                                       method = "REML"))
    res <- residuals(null_model, type = "response")

    return(res)

  }else{

    if(modeltype == "separated"){
      #Initialize formula
      formula <- paste0(colnames(design)[response], " ~ ")


      #Find the right column in design for each predictor
      vars_c_design <- colnames(design)[model_on]

      for(c_ind in seq(1, length(model_on))){

        #Combine into a single formula

        c_formula <- paste("s(", vars_c_design[c_ind], ", bs = 'tp')", sep = "")

        if(c_ind == 1){
          formula <- paste(formula, c_formula)
        }else{
          formula <- paste(formula, c_formula, sep = " + ")
        }

      }

      #Add spatial smooth
      formula <- paste(formula,
                       paste0("s(", sites, ", bs = 'mrf', xt = list(nb = nb))"),
                       sep = " + ")

      formula <- as.formula(formula)

      fam_fun <- get(distribution)

      #Fit GAM
      fit_bam <- suppressWarnings(gam(formula = formula, family = fam_fun(),
                                      data = design,
                                      method = "REML"))

      #Return residuals
      res <- residuals(fit_bam, type = "response")

      return(res)
    }else{
      if(modeltype == "tensor"){

        #Initialize formula
        formula <- paste0(colnames(design)[response], " ~ ")


        #Find the right column in design for each predictor
        vars_c_design <- colnames(design)[model_on]

        for(c_ind in seq(1, length(model_on))){

          #Combine into a single formula

          c_formula <- paste0("te(", vars_c_design[c_ind],
                              ",",
                              sites,
                              ", bs = c('tp', 'mrf'), xt = list(nb = nb))")

          if(c_ind == 1){
            formula <- paste(formula, c_formula)
          }else{
            formula <- paste(formula, c_formula, sep = " + ")
          }

        }

        formula <- as.formula(formula)

        fam_fun <- get(distribution)

        #Fit GAM
        fit_bam <- suppressWarnings(gam(formula = formula, family = fam_fun(),
                                        data = design,
                                        method = "REML"))

        #Return residuals
        res <- residuals(fit_bam, type = "response")

      }else{
        if(modeltype == "both"){

          #Initialize formula
          formula <- paste0(colnames(design)[response], " ~ ")

          #Find the right column in design for each predictor
          vars_c_design <- colnames(design)[model_on]

          for(c_ind in seq(1, length(model_on))){

            #Combine into a single formula

            c_formula <- paste0("s(", vars_c_design[c_ind], ", bs = 'tp') + ",
                                "te(", vars_c_design[c_ind],
                                ",",
                                sites,
                                ", bs = c('tp', 'mrf'), xt = list(nb = nb))")

            if(c_ind == 1){
              formula <- paste(formula, c_formula)
            }else{
              formula <- paste(formula, c_formula, sep = " + ")
            }

          }

          #Add spatial smooth
          formula <- paste(formula,
                           paste0("s(", sites, ", bs = 'mrf', xt = list(nb = nb))"),
                           sep = " + ")

          formula <- as.formula(formula)

          fam_fun <- get(distribution)

          #Fit GAM
          fit_bam <- suppressWarnings(gam(formula = formula, family = fam_fun(),
                                          data = design,
                                          method = "REML"))


          #Return residuals
          res <- residuals(fit_bam, type = "response")

        }else{
          stop(paste0("\nUnknown model type ", modeltype,
                      ". Options are 'separated', 'tensor', or 'both'."))}
      }
    }
  }
}

#' @export
fit_custom_gam.coord_dsgn <- function(design,
                                      response = a, #variable to fit
                                      model_on = b, #variable(s) to predict with or "null" for null model
                                      distribution = distribution, #distribution of a
                                      modeltype = modeltype,
                                      coords = c("x", "y") #vector of names of coordinate axes
){

  if (response %in% model_on) {
    stop("Self-model detected: response appears in predictors")
  }

  if("null" %in% model_on){

    #Create spatial smooth
    formula <- paste0("s(", coords[1], ", ", coords[2], ")")
    #Connect to response
    formula <- paste(colnames(design)[response], formula, sep = " ~ ")
    #Format
    formula <- as.formula(formula)

    fam_fun <- get(distribution)

    null_model <- suppressWarnings(gam(formula, family = fam_fun(),
                                       data = design,
                                       method = "REML"))
    res <- residuals(null_model, type = "response")

    return(res)

  }else{

    if(modeltype == "separated"){
      #Initialize formula
      formula <- paste0(colnames(design)[response], " ~ ")


      #Find the right column in design for each predictor
      vars_c_design <- colnames(design)[model_on]

      for(c_ind in seq(1, length(model_on))){

        #Combine into a single formula

        c_formula <- paste("s(", vars_c_design[c_ind], ", bs = 'tp')", sep = "")

        if(c_ind == 1){
          formula <- paste(formula, c_formula)
        }else{
          formula <- paste(formula, c_formula, sep = " + ")
        }

      }

      #Add spatial smooth
      formula <- paste(formula,
                       paste0("s(", coords[1], ", ", coords[2], ")"),
                       sep = " + ")

      formula <- as.formula(formula)

      fam_fun <- get(distribution)

      #Fit GAM
      fit_gam <- suppressWarnings(gam(formula = formula, family = fam_fun(),
                                      data = design,
                                      method = "REML"))

      #Return residuals
      res <- residuals(fit_gam, type = "response")

      return(res)
    }else{
      if(modeltype == "tensor"){

        #Initialize formula
        formula <- paste0(colnames(design)[response], " ~ ")


        #Find the right column in design for each predictor
        vars_c_design <- colnames(design)[model_on]

        for(c_ind in seq(1, length(model_on))){

          #Combine into a single formula

          c_formula <- paste0("te(", vars_c_design[c_ind],
                              ",",
                              coords[1], ", ", coords[2],
                              ", bs = c('tp'))")

          if(c_ind == 1){
            formula <- paste(formula, c_formula)
          }else{
            formula <- paste(formula, c_formula, sep = " + ")
          }

        }

        formula <- as.formula(formula)

        fam_fun <- get(distribution)

        #Fit GAM
        fit_gam <- suppressWarnings(gam(formula = formula, family = fam_fun(),
                                        data = design,
                                        method = "REML"))

        #Return residuals
        res <- residuals(fit_gam, type = "response")

      }else{
        if(modeltype == "both"){

          #Initialize formula
          formula <- paste0(colnames(design)[response], " ~ ")

          #Find the right column in design for each predictor
          vars_c_design <- colnames(design)[model_on]

          for(c_ind in seq(1, length(model_on))){

            #Combine into a single formula

            c_formula <- paste0("s(", vars_c_design[c_ind], ", bs = 'tp') + ",
                                "te(", vars_c_design[c_ind],
                                ",",
                                coords[1], ", ", coords[2],
                                ", bs = c('tp'))")

            if(c_ind == 1){
              formula <- paste(formula, c_formula)
            }else{
              formula <- paste(formula, c_formula, sep = " + ")
            }

          }

          #Add spatial smooth
          formula <- paste(formula,
                           paste0("s(", coords[1], ", ", coords[2], ")"),
                           sep = " + ")

          formula <- as.formula(formula)

          fam_fun <- get(distribution)

          #Fit GAM
          fit_gam <- suppressWarnings(gam(formula = formula, family = fam_fun(),
                                          data = design,
                                          method = "REML"))


          #Return residuals
          res <- residuals(fit_gam, type = "response")

        }else{
          stop(paste0("\nUnknown model type ", modeltype,
                      ". Options are 'separated', 'tensor', or 'both'."))}
      }
    }
  }
}
