# Load necessary libraries
library(ezRun)
library(SummarizedExperiment)
library(DESeq2)


# Load input dataset
input_dataset <- "/path/to/input_dataset.tsv"  # Change this to the actual path of your input dataset
dataRoot <- "/srv/gstore/projects"  # Change this to the actual path of your data root
de_output_dataset <- "/path/to/dataset.tsv"  # Change to actual path of output
input_dataset <- EzDataset$new(file = input_dataset, dataRoot = dataRoot)

## change to corresponding sample list for sc
sampleList <- c("vis-1033", "vis-1037", "vis-1083", "vis-1144", "vis-1162", "vis-1190", "vis-1199", "vis-1202", "vis-1222", "vis-1254", "vis-1264", "vis-1274", "vis-1282", "vis-1284", "vis-1298", "vis-1303", "vis-1307", "vis-1349", "vis-1375", "vis-1394", "vis-1424", "vis-1444", "vis-1459", "vis-1468", "vis-1480", "vis-1491", "vis-1493", "vis-1497", "vis-1498", "vis-1562", "vis-1583", "vis-1710", "vis-1713", "vis-1993", "vis-2010", "vis-2081", "vis-2086", "vis-2088", "vis-2091", "vis-2107", "vis-2111", "vis-2115", "vis-2162", "vis-2171", "vis-2195", "vis-2199", "vis-223", "vis-228", "vis-2335", "vis-2438", "vis-2447", "vis-2598", "vis-2714", "vis-2883", "vis-2898", "vis-3036", "vis-314", "vis-3202", "vis-458", "vis-497", "vis-615", "vis-627", "vis-650", "vis-665", "vis-695", "vis-715", "vis-743", "vis-750", "vis-784", "vis-794", "vis-833", "vis-880", "vis-945", "vis-949")

# Load output dataset
output_dataset <- EzDataset$new(file = de_output_dataset, dataRoot = dataRoot)

# Define the parameters
params <- list(
  cores = 4,
  ram = 12,
  scratch = 10,
  partition = "employee",
  process_mode = "DATASET",
  samples = sampleList,
  refBuild = "Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_32-2019-11-05",
  refFeatureFile = "genes.gtf",
  featureLevel = "gene",
  grouping = "Status",
  sampleGroup = "muo",
  refGroup = "mho",
  onlyCompGroupsHeatmap = TRUE,
  grouping2 = "EMR",
  backgroundExpression = 10,
  transcriptTypes = "protein_coding",
  runGO = TRUE,
  pValThreshGO = 0.05,
  log2RatioThreshGO = 0,
  fdrThreshORA = 0.05,
  fdrThreshGSEA = 0.05,
  useRefGroupAsBaseline=FALSE,
  onlyCompGroupsHeatmap=FALSE,
  testMethod="deseq2",
  normMethod="DESeq2_MedianRatio",
  specialOptions = "",
  expressionName = "",
  mail = "",
  Rversion = "Dev/R/4.3.2",
  sushi_app = "DESeq2App",
  comparison = "muo--over--mho",
  name = "muo--over--mho"
)
params <- ezParam(params)

cleanupTwoGroupsInput <- function(input, param) {
  dataset <- input$meta
  if (param$useFactorsAsSampleName) {
    dataset$Name <- rownames(dataset)
    rownames(dataset) <- addReplicate(apply(ezDesignFromDataset(dataset, param), 1, paste, collapse = "_"))
  }
  inputMod <- EzDataset(meta = dataset, dataRoot = param$dataRoot)
  if (!is.null(param$removeOtherGroups) && param$removeOtherGroups) {
    grouping <- inputMod$getColumn(param$grouping)
    keep <- grouping %in% c(param$sampleGroup, param$refGroup)
    inputMod <- inputMod$subset(keep)
  }
  return(inputMod)
}

twoGroupCountComparison <- function(rawData) {
  require(SummarizedExperiment)
  x <- assays(rawData)$counts
  param <- metadata(rawData)$param
  presentFlag <- assays(rawData)$presentFlag
  job <- ezJobStart("twoGroupCountComparison")
  metadata(rawData)$analysis <- "NGS two group analysis"
  if (is.null(param$testMethod)) {
    param$testMethod <- "glm"
    metadata(rawData)$param <- param
  }
  if (param$testMethod != "glm") {
    param$deTest <- NULL
    metadata(rawData)$param <- param
  }
  metadata(rawData)$method <- param$testMethod

  if (ezIsSpecified(param$grouping2)) {
    if (param$testMethod %in% c("glm", "limma", "deseq2")) {
      metadata(rawData)$method <- paste(
        metadata(rawData)$method,
        "using two factors"
      )
    } else {
      return(list(error = paste("Second factor only supported for the test methods glm, sam and deseq2")))
    }
  }

  isSample <- param$grouping == param$sampleGroup
  isRef <- param$grouping == param$refGroup
  isPresent <- ezPresentFlags(x,
    presentFlag = presentFlag, param = param,
    isLog = FALSE
  )
  useProbe <- logical(nrow(x))
  useProbe[rowMeans(isPresent[, isRef, drop = FALSE]) >= 0.5] <- TRUE
  useProbe[rowMeans(isPresent[, isSample, drop = FALSE]) >= 0.5] <- TRUE
  rowData(rawData)$isPresentProbe <- useProbe
  assays(rawData)$isPresent <- isPresent

  res <- switch(param$testMethod,
    deseq2 = runDeseq2(round(x), param$sampleGroup, param$refGroup,
      param$grouping,
      grouping2 = param$grouping2,
      isPresent = useProbe,
      cooksCutoff = ezIsSpecified(param$cooksCutoff) && param$cooksCutoff
    ),
    stop("unsupported testMethod: ", param$testMethod)
  )
  pValue <- res$pval
  pValue[is.na(pValue)] <- 1
  useProbe[is.na(useProbe)] <- FALSE
  fdr <- rep(NA, length(pValue))
  fdr[useProbe] <- p.adjust(pValue[useProbe], method = "fdr")
  
  rowData(rawData)$log2Ratio <- res$log2FoldChange
  colData(rawData)$sf <- res$sf
  rowData(rawData)$pValue <- pValue
  rowData(rawData)$fdr <- fdr
  rowData(rawData)$usedInTest <- useProbe
  metadata(rawData)$nativeResult <- res
  assays(rawData)$xNorm <- ezScaleColumns(x, colData(rawData)$sf)

  ezWriteElapsed(job, status = "done")
  metadata(rawData)$summary <- c(
    "Name" = param$name,
    "Reference Build" = param$refBuild,
    "Feature Level" = metadata(rawData)$featureLevel,
    "Normalization" = param$normMethod
  )

  seqAnno <- data.frame(rowData(rawData),
                        row.names = rownames(rawData), check.names = FALSE)
                        
  if (doGo(param, seqAnno)) {
    metadata(rawData)$enrichInput <- compileEnrichmentInput(param, rawData)
    metadata(rawData)$enrichResult <- ezEnricher(metadata(rawData)$enrichInput, param)
    metadata(rawData)$enrichResultGSEA <- ezGSEA(metadata(rawData)$enrichInput, param)
  }
  return(rawData)
}

runDeseq2 <- function(x, sampleGroup, refGroup, grouping, grouping2 = NULL,
                      isPresent = NULL, cooksCutoff = FALSE) {
  require(DESeq2)
  colData <- data.frame(
    grouping = as.factor(grouping),
    row.names = colnames(x)
  )
  dds <- DESeqDataSetFromMatrix(
    countData = x, colData = colData,
    design = ~grouping
  )
  dds <- estimateSizeFactors(dds, controlGenes = isPresent)
  sf <- 1 / dds@colData$sizeFactor

  isSample <- grouping == sampleGroup
  isRef <- grouping == refGroup
  grouping <- grouping[isSample | isRef]
  x <- x[, isSample | isRef]
  if (ezIsSpecified(grouping2)) {
    grouping2 <- grouping2[isSample | isRef]
  }
  
  if (ezIsSpecified(grouping2)) {
    colData <- data.frame(
      grouping = as.factor(grouping), grouping2 = grouping2,
      row.names = colnames(x)
    )
    colData$grouping <- factor(x = colData$grouping, 
                               levels = c(refGroup, sampleGroup))
    dds <- DESeqDataSetFromMatrix(
      countData = x, colData = colData,
      design = ~ grouping + grouping2
    )
  } else {
    colData <- data.frame(
      grouping = as.factor(grouping), row.names = colnames(x)
      )
    colData$grouping <- factor(x = colData$grouping, 
                               levels = c(refGroup, sampleGroup))
    dds <- DESeqDataSetFromMatrix(
      countData = x, colData = colData,
      design = ~grouping
    )
  }
  dds <- estimateSizeFactors(dds, controlGenes = isPresent)
  dds <- DESeq(dds, quiet = FALSE, minReplicatesForReplace = Inf)
  res <- results(dds,
    contrast = c("grouping", sampleGroup, refGroup),
    cooksCutoff = cooksCutoff
  )
  res <- as.list(res)
  res$sf <- sf
  res$dds = dds
  return(res)
}


# Define the main ezMethodDeseq2 function
ezMethodDeseq2 = function(input=NA, output=NA, param=NA){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))
  stopifnot(param$sampleGroup != param$refGroup)
  
  input = cleanupTwoGroupsInput(input, param)
  param$groupingName <- param$grouping
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1){
    param$grouping2Name <- param$grouping2
    param$grouping2 = input$getColumn(param$grouping2)
  }
  
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport("00index.html", param=param, error=rawData$error)
    return("Error")
  }
  
  deResult = twoGroupCountComparison(rawData)
  if (isError(deResult)){
    writeErrorReport("00index.html", param=param, error=deResult$error)
    return("Error")
  }
  dds = metadata(deResult)$nativeResult$dds
  dataset <- data.frame(colData(deResult), check.names = FALSE)
  dataset <- dataset[rownames(dataset) %in% rownames(dds@colData), ]
  seqAnno <- data.frame(rowData(deResult),
                        row.names = rownames(deResult),
                        check.names = FALSE)
  
  makeRmdReport(output=output, param=param, deResult=deResult, rmdFile="twoGroups.Rmd", reportTitle = param$comparison)
  rmStatus <- file.remove(list.files(pattern="enrichr-.*rds"))
  return("Success")
}

ezMethodDeseq2 = function(input=NA, output=NA, param=NA){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))
  stopifnot(param$sampleGroup != param$refGroup)
  
  input = cleanupTwoGroupsInput(input, param)
  param$groupingName <- param$grouping
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1){
    param$grouping2Name <- param$grouping2
    param$grouping2 = input$getColumn(param$grouping2)
  }
  
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport("00index.html", param=param, error=rawData$error)
    return("Error")
  }
  
  deResult = twoGroupCountComparison(rawData)
  if (isError(deResult)){
    writeErrorReport("00index.html", param=param, error=deResult$error)
    return("Error")
  }
  dds = metadata(deResult)$nativeResult$dds
  dataset <- data.frame(colData(deResult), check.names = FALSE)
  dataset <- dataset[rownames(dataset) %in% rownames(dds@colData), ]
  seqAnno <- data.frame(rowData(deResult),
                        row.names = rownames(deResult),
                        check.names = FALSE)
  
  makeRmdReport(output=output, param=param, deResult=deResult, rmdFile="twoGroups.Rmd", reportTitle = param$comparison)
  rmStatus <- file.remove(list.files(pattern="enrichr-.*rds"))
  return("Success")
}

# Run the DESeq2 analysis
ezMethodDeseq2(input = input_dataset, output = output_dataset, param = params)
