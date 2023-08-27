##' @title sc_dim
##' @rdname sc-dim
##' @param object Seurat or SingleCellExperiment object
##' @param dims selected dimensions (must be a two-length vector) that are used in visualization
##' @param reduction reduction method, default is NULL and will use the default setting store in the object
##' @param cells selected cells to plot (default is all cells)
##' @param slot slot to pull expression data from (e.g., 'count' or 'data')
##' @param mapping aesthetic mapping
##' @param ... additional parameters pass to 'scattermore::geom_scattermore()'
##' @return dimension reduction plot
##' @seealso
##'  [geom_scattermore][scattermore::geom_scattermore]; 
##' @export
sc_dim <- function(object, 
                    dims=c(1,2), reduction=NULL, colour_by="label",
                    cells=NULL, slot = "data", assay=NULL, 
                    mapping = NULL, ...) {

    if (inherits(object, "SingleCellExperiment")) {
        d <- get_dim_data_sc(object = object, features = NULL,
                    dims = dims, reduction = reduction, 
                    cells = cells, slot = slot, assay=assay,
                    colour_by=colour_by)
    } else {
        d <- get_dim_data(object = object, features = NULL,
                    dims = dims, reduction = reduction, 
                    cells = cells, slot = slot)
    }

    default_mapping <- aes(color=!!sym(colour_by))
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }               
    p <- sc_dim_internal(d, mapping, ...)
    return(p)
}

##' @importFrom tidydr theme_dr
sc_dim_internal <- function(data, mapping, ...) {
    dims <- names(data)[1:2]
    ggplot(data, aes(.data[[dims[1]]], .data[[dims[2]]])) + 
        sc_geom_point(mapping, ...) + 
        theme_dr()
} 

##' @importFrom SeuratObject FetchData
get_dim_data <- function(object, features = NULL, 
                    dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data") {
    if (is.null(cells)) {
        cells <- colnames(object)
    }
    if (!is.null(dims)) {
        if (is.null(reduction)) {
            reduction <- SeuratObject::DefaultDimReduc(object)
        }
        dims <- paste0(SeuratObject::Key(object = object[[reduction]]), dims)
    }
    SeuratObject::FetchData(object, vars = c(dims, "ident", 
        features), cells = cells, slot = slot)
}

get_dim_data_sc <- function(object, features = NULL, 
                    dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data", assay=NULL,
                    colour_by="label") {
    if (is.null(cells)) {
        cells <- colnames(object)
    }
    if (is.null(assay)) {
        assay <- assayNames(object)[1]
    }
    if (!is.null(dims)) {
        if (is.null(reduction)) {
            reduction <- SingleCellExperiment::reducedDims(pbmc_small) |>
                names() |> _[1]
        }
        reducedMat <- SingleCellExperiment::reducedDims(object)[[reduction]] |>
            as.data.frame()
        dims <- reducedMat |> colnames() |> _[dims]
    }
    if (is.null(features)) {
        feats <- NULL
        bind_df <- cbind(reducedMat[, dims],
            object |> colData() |> as.data.frame())
    } else {
        feats <- assay(object, assay) |> t() |>
            as.data.frame() |> select(features)
        bind_df <- cbind(reducedMat[, dims],
            feats,
            object |> colData() |> as.data.frame())
    }
    bind_df[cells, ] |> mutate(ident=!!sym(colour_by))
}

