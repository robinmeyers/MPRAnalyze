
#' estimate random effect variance
#' @param obj the MpraObject
#' @param rand.factor the factor on which to estimate random effect
#' @param lib.factor the factor associating each sample to a library. Can be a
#' factor or the name of a column in the object's colAnnot. If not provided, the
#' data is assumed to have been generated from a single library, and constant
#' library depth is set.
#' @return the MpraObject with estimated value for random effect variance
#'
#' @export
#'
#' @importFrom magrittr set_rownames `%>%`
#' @importFrom dplyr as_tibble mutate select group_by group_map left_join
#' @importFrom tidyselect all_of
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom DESeq2 vst
#' @importFrom lme4 lmer VarCorr
#' @importFrom fitdistrplus fitdist
#' @import invgamma
#'
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=NULL, nbatch=2, nbc=20)
#' obj <- MpraObject(dnaCounts = data$obs.dna,
#'                   rnaCounts = data$obs.rna,
#'                   colAnnot = data$annot)
#' obj <- estimateRandomEffectVariance(obj, lib.factor = "sample", rand.factor = "barcode")

estimateRandomEffectVariance <- function(obj, rand.factor, lib.factor) {

    dna_long <- obj@dnaCounts %>% as_tibble(rownames = "row") %>%
        pivot_longer(!row, names_to = "col", values_to = "count") %>%
        left_join(obj@dnaAnnot %>% as_tibble(rownames = "col"), by = "col") %>%
        mutate(element = paste(.[["row"]], .[[rand.factor]], sep = '.'))

    dna_df <- dna_long %>%
        pivot_wider(id_cols = element, names_from = all_of(lib.factor), values_from = "count")

    dna_mat <- dna_df %>% as.data.frame() %>% set_rownames(.$element) %>%
        select(all_of(unique(dna_long[[lib.factor]]))) %>%
        as.matrix()

    rna_long <- obj@rnaCounts %>% as_tibble(rownames = "row") %>%
        pivot_longer(!row, names_to = "col", values_to = "count") %>%
        left_join(obj@rnaAnnot %>% as_tibble(rownames = "col"), by = "col") %>%
        mutate(element = paste(.[["row"]], .[[rand.factor]], sep = '.'))

    rna_df <- rna_long %>%
        pivot_wider(id_cols = element, names_from = all_of(lib.factor), values_from = "count")

    rna_mat <- rna_df %>% as.data.frame() %>% set_rownames(.$element) %>%
        select(all_of(unique(rna_long[[lib.factor]]))) %>%
        as.matrix()

    rna_cols <- seq_along(colnames(rna_mat))
    dna_cols <- seq_along(colnames(dna_mat)) + ncol(rna_mat)
    full_mat <- cbind(rna_mat, dna_mat[rownames(rna_mat),])

    full_vst <- vst(full_mat, nsub = min(ceiling(0.5*nrow(full_mat)), 1000))

    activity_mat <- full_vst[,rna_cols] - apply(full_vst[,dna_cols], 1, mean)

    activity_long <- activity_mat %>% as_tibble(rownames = "element") %>%
        pivot_longer(!element, names_to = lib.factor, values_to = "activity") %>%
        left_join(rna_long, by = c(lib.factor, "element"))


    lmers <- activity_long %>%
        group_by(row) %>%
        group_map(~ lmer(activity ~ (1 | element), data = .x))

    randvars <- map_df(lmers, ~ as.data.frame(VarCorr(.)) %>% dplyr::filter(grp == "element"))

    invgamma_fit <- fitdistrplus::fitdist(randvars$sdcor, "invgamma")
    invgamma_params <- invgamma_fit$estimate

    # random_effects_estimate <- activity_long %>%
    #     group_by(row) %>%
    #     mutate(mean_row_activity = mean(activity, na.rm = T)) %>%
    #     group_by(across(all_of(c("row", rand.factor, "mean_row_activity")))) %>%
    #     summarise(mean_rand_factor_activity = mean(activity, na.rm = T)) %>%
    #     mutate(mu = mean_rand_factor_activity - mean_row_activity)
    #
    # obj@randEffVar <- var(random_effects_estimate$mu)

    obj@randEffVar <- invgamma_params

    return(obj)

}

