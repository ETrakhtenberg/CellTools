#Read reference and query seurat object files into the R 

reference <- readRDS(“/path-to-referencedata”) 
query <- readRDS(“/path-to-querydata”) 

#Denoise a reference Seurat Object

reference<- DeNoise(ref, lab = "type", coef = 1.5)

#' @param ref: Reference Seurat object with pre-calculated UMAP
#' @param lab: The reference label to treat as clusters. Default is "seurat_clusters"
#' @param coef: The coefficient to multiply by the IQR. Default is 1.5
#' @return: Returns a Seurat object with a new slot called "denoise"


# Find the optimal k-weight parameter for KNN
# This function iterates through a select range of k-weight parameters and stores the optimal parameter in the Seurat object
reference <- OptiK(reference, lab = "type", range = c(5,50), dims = 50, perc = 0.2)
#' @param data Reference Seurat object 
#' @param lab The reference label to treat as clusters. Default is "seurat_clusters"
#' @param range The range of k-weight values to loop through. Default is 10-100.
#' @param dims The number of PCs to consider when discovering anchors
#' @param perc The percent of cells to subsample from the refeerence for the test query. Default is 20%.
#' @param num If perc = F, OptiK will look for a specific number of cells specified here.
#' @param seed Set the seed for reproducible results. Default = 1984.
#' @return Returns a Seurat object with a new slot where the k-weight is stored: misc$CellTools$opti_k


#' Map query Seurat object to reference Seurat object
#' This function is a wrapper for Seurat's TransferData, IntegrateEmbeddings, and ProjectUMAP functions, meant to be run after
Result<- MapTo(ref, query, lab = "type", dims = 50)
#' @param ref Reference Seurat object 
#' @param query Wuery Seurat object
#' @param lab The reference label to treat as clusters. Default is "seurat_clusters"
#' @param k The k-weight to use when weighting anchors. Default is opti
#' @param get.anchors A logical that when true will first run Seurat's FindTransferAnchors. Default is TRUE.
#' @param dims Number of PCs to use when get.anchors = T. Default is 30
#' @return Returns a Seurat object with a new slot where the k-weight is stored: misc$CellTools$opti_k


   
#' Optimized Tracing  of Query to Reference Data
#' This function examines the global transcriptomic properties of single cell data.
correlation<-GlobalProbe(dat,lab="type")
#' @param dat A Seurat object
#' @param lab The reference @meta.data slot on which to compare global parameters
#' @return Global correlation analyses of the clusters provided in lab parameter























