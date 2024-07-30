geomx_import_fun_HRK <- function(countFile, sampleAnnoFile, featureAnnoFile,
                             rmNegProbe, NegProbeName,
                             colnames.as.rownames,
                             coord.colnames) {
  
  # remove the NegProbe gene from the count matrix and save it in the metadata
  if (rmNegProbe) {
    if(is.data.frame(countFile)){
      countdata <- countFile
    } else {
      countdata <- as.data.frame(readr::read_tsv(countFile), optional = TRUE)
    }
    
    # raw count without negprobes
    # make sure count data have the gene column name as pre-defined, such as TargetName.
    if (!colnames.as.rownames[1] %in% colnames(countdata)) {
      stop("colnames.as.rownames[1] is not in the column names of your count file.")
    }
    # make sure the name of negprobe is in the gene column of count data.
    if (!NegProbeName %in% as.matrix(countdata[, colnames.as.rownames[1]])) {
      stop("NegProbeName is not found in your count file.")
    }
    
    # filter the count data, remove the negprobe.
    countdata_filtered0 <- countdata[countdata[, colnames.as.rownames[1]] != NegProbeName, ]
    countdata_filtered <- countdata_filtered0[, !colnames(countdata_filtered0) %in%
                                                colnames.as.rownames[1]]
    rownames(countdata_filtered) <- as.vector(as.matrix(countdata_filtered0[, colnames.as.rownames[1]]))
    
    
    # gene meta without negprobes
    if (!all(is.na(featureAnnoFile))) {
      if(is.data.frame(featureAnnoFile)){
        genemeta <- featureAnnoFile
      } else {
        genemeta <- as.data.frame(readr::read_tsv(featureAnnoFile), optional = TRUE)
      }
      
      stopifnot(colnames.as.rownames[3] %in% colnames(genemeta)) # make sure column name is there in the gene meta.
      
      genemeta_filtered0 <- as.data.frame(genemeta[genemeta[, colnames.as.rownames[3]] != NegProbeName, ])
      genemeta_filtered <- as.data.frame(genemeta_filtered0[, !colnames(genemeta_filtered0) %in%
                                                colnames.as.rownames[3]])
      rownames(genemeta_filtered) <- genemeta_filtered0[, colnames.as.rownames[3]]
      genemeta_filtered <- genemeta_filtered[rownames(countdata_filtered), ]
      # arrange the gene meta, as the same order as count table.
    } else {
      genemeta_filtered <- data.frame(Type = rep("gene", nrow(countdata_filtered)))
    }
    
    # sample meta
    if(is.data.frame(sampleAnnoFile)){
      samplemeta <- as.data.frame(sampleAnnoFile)
    } else {
      samplemeta <- as.data.frame(readr::read_tsv(sampleAnnoFile), optional = TRUE)
    }
    
    stopifnot(colnames.as.rownames[2] %in% colnames(samplemeta)) # make sure column name is there.
    
    
    samplemeta_filtered <- as.data.frame(samplemeta[, !colnames(samplemeta) %in%
                                        colnames.as.rownames[2]])
    rownames(samplemeta_filtered) <- samplemeta[, colnames.as.rownames[2]]
    samplemeta_filtered <- samplemeta_filtered[colnames(countdata_filtered), ]
    # arrange according to count table.
    
    # negprobe raw count
    negprobecount <- countdata[countdata[, colnames.as.rownames[1]] == 
                                 NegProbeName, ]
    nprobename <- as.vector(as.matrix(negprobecount[,colnames.as.rownames[1]]))
    if(length(nprobename) != length(unique(nprobename))){
      nprobename <- paste0(nprobename,"_",seq(length(nprobename)))
    }
    negprobecount <- negprobecount[, !colnames(negprobecount) %in% 
                                     colnames.as.rownames[1]]
    rownames(negprobecount) <- nprobename
    
    
    # logCPM count
    countdata_filtered_lcpm <- edgeR::cpm(countdata_filtered, log = TRUE)
    
    
    # output spe
    spe <- SpatialExperiment::SpatialExperiment(
      assay = list(
        counts = countdata_filtered,
        logcounts = countdata_filtered_lcpm
      ),
      colData = samplemeta_filtered,
      rowData = genemeta_filtered,
      metadata = list(NegProbes = negprobecount),
      spatialCoords = as.matrix(samplemeta_filtered[, coord.colnames])
    )
  } else {
    # it doesn't remove the NegProbe genes, leave them in the count matrix
    # raw count
    if (is.data.frame(countFile)) {
      countdata0 <- countFile
    }
    else {
      countdata0 <- as.data.frame(readr::read_tsv(countFile), 
                                  optional = TRUE)
    }
    stopifnot(colnames.as.rownames[1] %in% colnames(countdata0))
    
    countdata <- countdata0[, !colnames(countdata0) %in% colnames.as.rownames[1]]
    rownames(countdata) <- as.vector(as.matrix(countdata0[, colnames.as.rownames[1]]))
    
    if (!all(is.na(featureAnnoFile))) {
      if (is.data.frame(featureAnnoFile)) {
        genemeta0 <- featureAnnoFile
      }
      else {
        genemeta0 <- as.data.frame(readr::read_tsv(featureAnnoFile), 
                                   optional = TRUE)
      }
      stopifnot(colnames.as.rownames[3] %in% colnames(genemeta0))
      
      genemeta <- genemeta0[, !colnames(genemeta0) %in% colnames.as.rownames[3]]
      rownames(genemeta) <- genemeta0[, colnames.as.rownames[3]]
      genemeta <- genemeta[rownames(countdata), ]
    }
    else {
      genemeta <- data.frame(Type = rep("gene", nrow(countdata)))
    }
    
    if (is.data.frame(sampleAnnoFile)) {
      samplemeta0 <- sampleAnnoFile
    }
    else {
      samplemeta0 <- as.data.frame(readr::read_tsv(sampleAnnoFile), 
                                   optional = TRUE)
    }
    stopifnot(colnames.as.rownames[2] %in% colnames(samplemeta0))
    
    samplemeta <- samplemeta0[, !colnames(samplemeta0) %in% colnames.as.rownames[2]]
    rownames(samplemeta) <- samplemeta0[, colnames.as.rownames[2]]
    samplemeta <- samplemeta[colnames(countdata), ]
    
    countdata_lcpm <- edgeR::cpm(countdata, log = TRUE)
    spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = countdata, 
                                                              logcounts = countdata_lcpm), colData = samplemeta, 
                                                rowData = genemeta, spatialCoords = as.matrix(samplemeta[, 
                                                                                                         coord.colnames]))
  }
  return(spe)
}

