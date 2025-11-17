

#A small change to the orig. SCTransform function
#Seurat 4.1.0
#
#compared to SCTransform
#
#added 
#Misc(object = assay.out, slot = "vst.out") <- vst.out
#in line 227
#fixes the issue of missing misc data in the SCT assay slot (vst.out)

SCTransform_v2 <- function (object, assay = "RNA", new.assay.name = "SCT", reference.SCT.model = NULL, 
    do.correct.umi = TRUE, ncells = 5000, residual.features = NULL, 
    variable.features.n = 3000, variable.features.rv.th = 1.3, 
    vars.to.regress = NULL, do.scale = FALSE, do.center = TRUE, 
    clip.range = c(-sqrt(x = ncol(x = object[[assay]])/30), sqrt(x = ncol(x = object[[assay]])/30)), 
    conserve.memory = FALSE, return.only.var.genes = TRUE, seed.use = 1448145, 
    verbose = TRUE, ...) 
{
    if (!is.null(x = seed.use)) {
        set.seed(seed = seed.use)
    }
    assay <- assay %||% DefaultAssay(object = object)
    assay.obj <- GetAssay(object = object, assay = assay)
    umi <- GetAssayData(object = assay.obj, slot = "counts")
    cell.attr <- slot(object = object, name = "meta.data")
    vst.args <- list(...)
    if ("batch_var" %in% names(x = vst.args)) {
        if (!(vst.args[["batch_var"]] %in% colnames(x = cell.attr))) {
            stop("batch_var not found in seurat object meta data")
        }
    }
    if (!is.null(x = reference.SCT.model)) {
        if (inherits(x = reference.SCT.model, what = "SCTModel")) {
            reference.SCT.model <- SCTModel_to_vst(SCTModel = reference.SCT.model)
        }
        if (is.list(x = reference.SCT.model) & inherits(x = reference.SCT.model[[1]], 
            what = "SCTModel")) {
            stop("reference.SCT.model must be one SCTModel rather than a list of SCTModel")
        }
        if ("latent_var" %in% names(x = vst.args)) {
            stop("custom latent variables are not supported when reference.SCT.model is given")
        }
        if (reference.SCT.model$model_str != "y ~ log_umi") {
            stop("reference.SCT.model must be derived using default SCT regression formula, `y ~ log_umi`")
        }
    }
    if ("latent_var" %in% names(x = vst.args)) {
        known.attr <- c("umi", "gene", "log_umi", "log_gene", 
            "umi_per_gene", "log_umi_per_gene")
        if (!all(vst.args[["latent_var"]] %in% c(colnames(x = cell.attr), 
            known.attr))) {
            stop("latent_var values are not from the set of cell attributes sctransform calculates by default and cannot be found in seurat object meta data")
        }
    }
    if (any(!vars.to.regress %in% colnames(x = cell.attr))) {
        stop("problem with second non-regularized linear regression; not all variables found in seurat object meta data; check vars.to.regress parameter")
    }
    if (any(c("cell_attr", "verbosity", "return_cell_attr", "return_gene_attr", 
        "return_corrected_umi") %in% names(x = vst.args))) {
        warning("the following arguments will be ignored because they are set within this function:", 
            paste(c("cell_attr", "verbosity", "return_cell_attr", 
                "return_gene_attr", "return_corrected_umi"), 
                collapse = ", "), call. = FALSE, immediate. = TRUE)
    }
    vst.args[["umi"]] <- umi
    vst.args[["cell_attr"]] <- cell.attr
    vst.args[["verbosity"]] <- as.numeric(x = verbose) * 2
    vst.args[["return_cell_attr"]] <- TRUE
    vst.args[["return_gene_attr"]] <- TRUE
    vst.args[["return_corrected_umi"]] <- do.correct.umi
    vst.args[["n_cells"]] <- min(ncells, ncol(x = umi))
    residual.type <- vst.args[["residual_type"]] %||% "pearson"
    res.clip.range <- vst.args[["res_clip_range"]] %||% c(-sqrt(x = ncol(x = umi)), 
        sqrt(x = ncol(x = umi)))
    if (!is.null(reference.SCT.model)) {
        sct.method <- "reference.model"
    }
    else if (!is.null(x = residual.features)) {
        sct.method <- "residual.features"
    }
    else if (conserve.memory) {
        sct.method <- "conserve.memory"
    }
    else {
        sct.method <- "default"
    }
    vst.out <- switch(EXPR = sct.method, reference.model = {
        if (verbose) {
            message("Using reference SCTModel to calculate pearson residuals")
        }
        do.center <- FALSE
        do.correct.umi <- FALSE
        vst.out <- reference.SCT.model
        clip.range <- vst.out$arguments$sct.clip.range
        umi.field <- paste0("nCount_", assay)
        vst.out$cell_attr <- if (umi.field %in% colnames(x = object[[]])) {
            data.frame(log_umi = log10(x = object[[umi.field, 
                drop = T]]))
        } else {
            data.frame(log_umi = log10(x = CalcN(object = object[[assay]])$nCount))
        }
        all.features <- intersect(x = rownames(x = vst.out$gene_attr), 
            y = rownames(x = umi))
        vst.out$gene_attr <- vst.out$gene_attr[all.features, 
            ]
        vst.out$model_pars_fit <- vst.out$model_pars_fit[all.features, 
            ]
        vst.out
    }, residual.features = {
        if (verbose) {
            message("Computing residuals for the ", length(x = residual.features), 
                " specified features")
        }
        return.only.var.genes <- TRUE
        do.correct.umi <- FALSE
        vst.args[["return_corrected_umi"]] <- FALSE
        vst.args[["residual_type"]] <- "none"
        vst.out <- do.call(what = "vst", args = vst.args)
        vst.out$gene_attr$residual_variance <- NA_real_
        vst.out
    }, conserve.memory = {
        return.only.var.genes <- TRUE
        vst.args[["residual_type"]] <- "none"
        vst.out <- do.call(what = "vst", args = vst.args)
        feature.variance <- get_residual_var(vst_out = vst.out, 
            umi = umi, residual_type = residual.type, res_clip_range = res.clip.range)
        vst.out$gene_attr$residual_variance <- NA_real_
        vst.out$gene_attr[names(x = feature.variance), "residual_variance"] <- feature.variance
        vst.out
    }, default = {
        vst.out <- do.call(what = "vst", args = vst.args)
        vst.out
    })
    feature.variance <- vst.out$gene_attr[, "residual_variance"]
    names(x = feature.variance) <- rownames(x = vst.out$gene_attr)
    if (verbose) {
        message("Determine variable features")
    }
    feature.variance <- sort(x = feature.variance, decreasing = TRUE)
    if (!is.null(x = variable.features.n)) {
        top.features <- names(x = feature.variance)[1:min(variable.features.n, 
            length(x = feature.variance))]
    }
    else {
        top.features <- names(x = feature.variance)[feature.variance >= 
            variable.features.rv.th]
    }
    vst.out <- switch(EXPR = sct.method, reference.model = {
        if (is.null(x = residual.features)) {
            residual.features <- top.features
        }
        residual.features <- Reduce(f = intersect, x = list(residual.features, 
            rownames(x = umi), rownames(x = vst.out$model_pars_fit)))
        residual.feature.mat <- get_residuals(vst_out = vst.out, 
            umi = umi[residual.features, , drop = FALSE], verbosity = as.numeric(x = verbose) * 
                2)
        vst.out$gene_attr <- vst.out$gene_attr[residual.features, 
            ]
        ref.residuals.mean <- vst.out$gene_attr[, "residual_mean"]
        vst.out$y <- sweep(x = residual.feature.mat, MARGIN = 1, 
            STATS = ref.residuals.mean, FUN = "-")
        vst.out
    }, residual.features = {
        residual.features <- intersect(x = residual.features, 
            y = rownames(x = vst.out$gene_attr))
        residual.feature.mat <- get_residuals(vst_out = vst.out, 
            umi = umi[residual.features, , drop = FALSE], verbosity = as.numeric(x = verbose) * 
                2)
        vst.out$y <- residual.feature.mat
        vst.out$gene_attr$residual_mean <- NA_real_
        vst.out$gene_attr$residual_variance <- NA_real_
        vst.out$gene_attr[residual.features, "residual_mean"] <- rowMeans2(x = vst.out$y)
        vst.out$gene_attr[residual.features, "residual_variance"] <- RowVar(x = vst.out$y)
        vst.out
    }, conserve.memory = {
        vst.out$y <- get_residuals(vst_out = vst.out, umi = umi[top.features, 
            ], residual_type = residual.type, res_clip_range = res.clip.range, 
            verbosity = as.numeric(x = verbose) * 2)
        vst.out$gene_attr$residual_mean <- NA_real_
        vst.out$gene_attr[top.features, "residual_mean"] = rowMeans2(x = vst.out$y)
        if (do.correct.umi & residual.type == "pearson") {
            vst.out$umi_corrected <- correct_counts(x = vst.out, 
                umi = umi, verbosity = as.numeric(x = verbose) * 
                  2)
        }
        vst.out
    }, default = {
        if (return.only.var.genes) {
            vst.out$y <- vst.out$y[top.features, ]
        }
        vst.out
    })
    if (do.correct.umi & residual.type == "pearson") {
        if (verbose) {
            message("Place corrected count matrix in counts slot")
        }
        assay.out <- CreateAssayObject(counts = vst.out$umi_corrected, 
            )
        vst.out$umi_corrected <- NULL
    }
    else {
        assay.out <- CreateAssayObject(counts = umi)
    }
    VariableFeatures(object = assay.out) <- residual.features %||% 
        top.features
    assay.out <- SetAssayData(object = assay.out, slot = "data", 
        new.data = log1p(x = GetAssayData(object = assay.out, 
            slot = "counts")))
    scale.data <- vst.out$y
    scale.data[scale.data < clip.range[1]] <- clip.range[1]
    scale.data[scale.data > clip.range[2]] <- clip.range[2]
    scale.data <- ScaleData(scale.data, features = NULL, vars.to.regress = vars.to.regress, 
        latent.data = cell.attr[, vars.to.regress, drop = FALSE], 
        model.use = "linear", use.umi = FALSE, do.scale = do.scale, 
        do.center = do.center, scale.max = Inf, block.size = 750, 
        min.cells.to.block = 3000, verbose = verbose)
    assay.out <- SetAssayData(object = assay.out, slot = "scale.data", 
        new.data = scale.data)
    vst.out$y <- NULL
    vst.out$arguments$sct.clip.range <- clip.range
    vst.out$arguments$sct.method <- sct.method
    Misc(object = assay.out, slot = "vst.out") <- vst.out
    assay.out <- as(object = assay.out, Class = "SCTAssay") # here, the misc slot from the previous line gets deleted
    assay.out <- SCTAssay(assay.out, assay.orig = assay)
    Misc(object = assay.out, slot = "vst.out") <- vst.out # added this to fix the issue with saving misc data
    slot(object = slot(object = assay.out, name = "SCTModel.list")[[1]], 
        name = "umi.assay") <- assay
    object[[new.assay.name]] <- assay.out
    if (verbose) {
        message(paste("Set default assay to", new.assay.name))
    }
    DefaultAssay(object = object) <- new.assay.name
    object <- LogSeuratCommand(object = object)
    return(object)
}


#When Seurat is NOT loaded, these two lines make the magic
#Now SCTransform_v2 can be invoked properly
e <- loadNamespace("Seurat")
environment(SCTransform_v2) <- e


