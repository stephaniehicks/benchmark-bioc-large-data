
args <- commandArgs(trailingOnly = TRUE)
run_id <- args[2]
data_name <- args[3]
data_type <- args[4]
fam <- args[5]
model_type <- args[6]
num_workers <- args[7]
mode <- args[8]

.id <- paste(fam, model_type, "residuals", sep="_")

# Number of observations
n_obs <- c(100, 1000, 5000, 10000, 25000, 50000, 75000, 100000)
# n_obs <- c(1000, 5000)

suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(DelayedArray)
  library(BiocSingular)
  library(BiocParallel)
  library(TENxBrainData)
  library(scry)
  library(readr)
})

if(!file.exists(here("output_time"))){
  dir.create(here("output_time"))
}

if(!file.exists(here("output_mem"))){
  dir.create(here("output_mem"))
}

if(mode == "time"){ 
  # Benchmark for in-memory data (includes time to load data from disk)
  if(data_type == "inmem"){
    time_realized_null <- matrix(nrow=length(n_obs), ncol=5)
    time_realized_pca  <- matrix(nrow=length(n_obs), ncol=5)
    for(i in seq_along(n_obs)){
      tenx <- TENxBrainData()[,seq_len(n_obs[i])]
      object <- tenx
      
      time_start <- proc.time()
      assay(object, "counts") <- as(assay(object, "counts"), "matrix")
      object <- nullResiduals(object, assay="counts",
                              fam = fam, type = model_type)
      time_end <- proc.time()
      time_realized_null[i,] <- c(time_end - time_start)
      
      time_start <- proc.time()
      pcs_out <- runPCA(assay(object, .id), rank = 10, 
                        BSPARAM = FastAutoParam(), get.rotation = FALSE,
                        BPPARAM = MulticoreParam(workers = num_workers))
      time_end <- proc.time()
      time_realized_pca[i,] <- c(time_end - time_start)
      rm(object)
    }

    dat <- data.frame("run_id" = run_id,
                      "data_name" = data_name, 
                      "fam" = fam, 
                      "model_type" = model_type, 
                      "num_workers" = num_workers, 
                      "data_type" = data_type,
                      "n_obs" = n_obs,
                      "null" = time_realized_null[,3], 
                      "pca"  = time_realized_pca[,3])
    write_csv(dat, 
        path = here("output_time", 
                  paste0(paste("time", run_id, data_name, fam, model_type, 
                               paste0(num_workers, "workers"), data_type, sep="_"), 
                         ".csv")))
  }

  # Benchmark for on-disk data (`DelayedArray`) (does not include time to load data from disk)
  if(data_type == "ondisk"){
    time_da_null <- matrix(nrow=length(n_obs), ncol=5)
    time_da_pca  <- matrix(nrow=length(n_obs), ncol=5)
    for(i in seq_along(n_obs)){
      object <- TENxBrainData()[,seq_len(n_obs[i])]
   
      time_start <- proc.time()
        object <- nullResiduals(object, assay="counts",
                                fam = fam, type = model_type)
      time_end <- proc.time()
      time_da_null[i,] <- c(time_end - time_start)
    
      time_start <- proc.time()
      pcs_out <- runPCA(assay(object, .id), rank = 10, 
                        BSPARAM = FastAutoParam(), get.rotation = FALSE,
                        BPPARAM = MulticoreParam(workers = num_workers))
      time_end <- proc.time()
      time_da_pca[i,] <- c(time_end - time_start)
    }
  
    dat <- data.frame("run_id" = run_id,
                      "data_name" = data_name, 
                      "fam" = fam, 
                      "model_type" = model_type, 
                      "num_workers" = num_workers, 
                      "data_type" = data_type,
                      "n_obs" = n_obs,
                      "null" = time_da_null[,3], 
                      "pca"  = time_da_pca[,3])
    write_csv(dat, 
        path = here("output_time", 
              paste0(paste("time", run_id, data_name, fam, model_type, 
                           paste0(num_workers, "workers"), data_type, sep="_"), 
                     ".csv")))
  }
}



if(mode == "mem"){ 
  if(!file.exists(here("output_mem","tmp_files"))){
    dir.create(here("output_mem", "tmp_files"))
  }
  
  if(data_type == "inmem"){
    max_mem_realized <- vector(length = length(n_obs))

    now <- format(Sys.time(), "%b%d%H%M%OS3")
     for(i in seq_along(n_obs)){
      tenx <- TENxBrainData()[,seq_len(n_obs[i])]
      object <- tenx
      out_name <- here("output_mem", "tmp_files",
                       paste0(paste("mem", run_id, data_name, fam, model_type, 
                                    paste0(num_workers, "workers"), data_type, 
                                    n_obs[i], now, sep="_"), 
                              ".out"))
      
      Rprof(filename = out_name, append = FALSE, memory.profiling = TRUE)
      assay(object, "counts") <- as(assay(object, "counts"), "matrix")
      object <- nullResiduals(object, assay="counts",
                              fam = fam, type = model_type)
      
      pcs_out <- runPCA(assay(object, .id), rank = 10, 
                        BSPARAM = FastAutoParam(), get.rotation = FALSE,
                        BPPARAM = MulticoreParam(workers = num_workers))
      rm(object)
      Rprof(NULL)
      
      profile <- summaryRprof(filename = out_name, chunksize = -1L, 
                              memory = "tseries", diff = FALSE)
      max_mem_realized[i] <- max(rowSums(profile[,1:3]))*0.00000095367432
      file.remove(out_name)
    }

    dat <- data.frame("run_id" = run_id,
                      "data_name" = data_name, 
                      "fam" = fam, 
                      "model_type" = model_type, 
                      "num_workers" = num_workers, 
                      "data_type" = data_type,
                      "n_obs" = n_obs,
                      "max_mem" = max_mem_realized/1000) # GB
    write_csv(dat, 
              path = here("output_mem", 
                          paste0(paste("mem", run_id, data_name, fam, model_type, 
                                       paste0(num_workers, "workers"), data_type, sep="_"), 
                                 ".csv")))
  }

  if(data_type == "ondisk"){
    max_mem_da <- vector(length = length(n_obs))
    
    now <- format(Sys.time(), "%b%d%H%M%OS3")
    for(i in seq_along(n_obs)){
      object <- TENxBrainData()[,seq_len(n_obs[i])]
      
      out_name <- here("output_mem", "tmp_files",
                       paste0(paste("mem", run_id, data_name, fam, model_type, 
                                    paste0(num_workers, "workers"), data_type, 
                                    n_obs[i], now, sep="_"), 
                              ".out"))
      
      Rprof(filename = out_name, append = FALSE, memory.profiling = TRUE)
      object <- nullResiduals(object, assay="counts",
                              fam = fam, type = model_type)

      pcs_out <- runPCA(assay(object, .id), rank = 10, 
                        BSPARAM = FastAutoParam(), get.rotation = FALSE,
                        BPPARAM = MulticoreParam(workers = num_workers))
      Rprof(NULL)
      
      profile <- summaryRprof(filename = out_name, chunksize = -1L, 
                              memory = "tseries", diff = FALSE)
      max_mem_da[i] <- max(rowSums(profile[,1:3]))*0.00000095367432
      file.remove(out_name)
    }
    
    dat <- data.frame("run_id" = run_id,
                      "data_name" = data_name, 
                      "fam" = fam, 
                      "model_type" = model_type, 
                      "num_workers" = num_workers, 
                      "data_type" = data_type,
                      "n_obs" = n_obs,
                      "max_mem" = max_mem_da/1000) # GB
    write_csv(dat, 
              path = here("output_mem", 
                          paste0(paste("mem", run_id, data_name, fam, model_type, 
                                       paste0(num_workers, "workers"), data_type, sep="_"), 
                                 ".csv")))
  }
}
