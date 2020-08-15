
args <- commandArgs(trailingOnly = TRUE)
run_id <- args[2]
data_name <- args[3]
data_type <- args[4]
fam <- args[5]
type <- args[6]
num_workers <- args[7]

.id <- paste(fam, type, "residuals", sep="_")

# Number of observations
# ns <- c(100, 1000, 5000, 10000, 25000, 50000, 75000, 100000)
ns <- c(1000, 5000)

suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(DelayedArray)
  library(BiocSingular)
  library(BiocParallel)
  library(TENxBrainData)
  library(scry)
})

if(!file.exists(here("data_timing"))){
  dir.create(here("data_timing"))
}

# Benchmark for in-memory data (includes time to load data from disk)
if(data_type == "inmem"){
  time_realized_null <- matrix(nrow=length(ns), ncol=5)
  time_realized_pca  <- matrix(nrow=length(ns), ncol=5)
  for(i in seq_along(ns)){
    tenx <- TENxBrainData()[,seq_len(ns[i])]
    object <- tenx
  
    time_start <- proc.time()
    assay(object, "counts") <- as.matrix(assay(object, "counts"))
    object <- nullResiduals(object, assay=assay,
                                   fam = fam, type=type)
    time_end <- proc.time()
    time_realized_null[i,] <- c(time_end - time_start)
    
    time_start <- proc.time()
    pcs_out <- runPCA(assay(object, .id), rank = 10, 
                      BSPARAM = FastAutoParam(), get.rotation = FALSE,
                      BPPARAM = MulticoreParam(workers = workers))
    time_end <- proc.time()
    time_realized_pca[i,] <- c(time_end - time_start)
  }

  dat <- data.frame("null" = time_realized_null[,3], 
                                  "pca"  = time_realized_pca[,3])
  colnames(dat) <- c("user", "system", "elapsed")
  readr::write_csv(dat, path = here::here("data_timing", 
                          paste0(paste("time", run_id, data_name, 
                                       fam, type, data_type, sep="_"), 
                                 ".csv")))
}



# Benchmark for on-disk data (`DelayedArray`) (does not include time to load data from disk)
if(data_type == "ondisk"){
  time_da_null <- matrix(nrow=length(ns), ncol=5)
  time_da_pca <- matrix(nrow=length(ns), ncol=5)
  for(i in seq_along(ns)){
    object <- TENxBrainData()[,seq_len(ns[i])]
  
    time_start <- proc.time()
    object <- nullResiduals(object, assay=assay,
                            fam = fam, type=type)
    time_end <- proc.time()
    time_da_null[i,] <- c(time_end - time_start)
    
    time_start <- proc.time()
    pcs_out <- runPCA(assay(object, .id), rank = 10, 
                      BSPARAM = FastAutoParam(), get.rotation = FALSE,
                      BPPARAM = MulticoreParam(workers = workers))
    time_end <- proc.time()
    time_da_pca[i,] <- c(time_end - time_start)
  }
  
  dat <- data.frame("null" = time_da_null[,3], 
                    "pca"  = time_da_pca[,3])
  colnames(dat) <- c("user", "system", "elapsed")
  readr::write_csv(dat, path = here::here("data_timing", 
                                          paste0(paste("time", run_id, data_name, 
                                                       fam, type, data_type, sep="_"), 
                                                 ".csv")))
}

