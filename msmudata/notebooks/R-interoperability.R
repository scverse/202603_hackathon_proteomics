library(QFeatures)
library(anndataR)
library(reticulate)

## reticulate::py_install("mudata")

md <- import("mudata")

mdata <- md$read_h5mu("minimal-zarr/minimal1.h5mu")

qf <- QFeatures(
    list(precursors = mdata$mod[["precursors"]]$as_SingleCellExperiment(),
         proteins = mdata$mod[["proteins"]]$as_SingleCellExperiment(),
         genes = mdata$mod[["genes"]]$as_SingleCellExperiment()))


## Needed to add the second AssayLink
rowData(qf[["genes"]]) <- DataFrame(genes = rownames(qf[["genes"]]))


qf <- qf |>
    addAssayLink(from = "precursors",
                 to = "proteins",
                 varFrom = "proteins",
                 varTo = "proteins") |>
    addAssayLink(from = "proteins",
                 to = "genes",
                 varFrom = "genes",
                 varTo = "genes")


plot(qf)
rownames(qf["gene0", , ])
rownames(qf["gene1", , ])
rownames(qf["protein1", , ])
rownames(qf["protein2", , ])
