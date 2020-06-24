################################################################################
### R script to import and trasform FCS files
### Simone Puccio
### June 18th, 2020
### designed to be executed with Cytophenograph
### run "Rscript Import_trasformationFCS.R --help" to get some help
################################################################################
.libPaths("/home/spuccio/miniconda3/envs/dev_cyto/lib/R/library")
suppressMessages(require("flowCore"))
suppressMessages(require("ggcyto"))
suppressMessages(require("flowWorkspace"))
suppressMessages(require("optparse"))
suppressMessages(require("Biobase"))
suppressMessages(require("ggplot2movies"))
suppressMessages(require("flowWorkspace"))
#
option_list = list(
    make_option(c("-f", "--fcsDir"), dest="fcsDir", default=NA, type='character',
              help="Path to the directory containing the Flow Cytometry Standard (FCS) version 3.0 files."),
    make_option(c("-p", "--phenoData"),dest="phenoData", default=NA, type='character',
                help="Path to the design/info file."),
    make_option(c("-m", "--markerDescription"),dest="markerDescription", default=NA, type='character',
                help="Path to the ."),
    make_option(c("-t", "--trasformationType"),dest="trasformationType", default='biexponential', type='character',help="Transformation definitions for the fluorescent channels."),
    make_option(c("-n", "--projectName"),dest="projectName", default=NA, type='character',help="Name of the project."),
    make_option(c("-o", "--outputDir"),dest="outputDir", default=NA, type='character',help="Path to the directory where store output files."),
    make_option(c("-w", "--widthBasis"),dest="widthBasis", default=-100, type='integer',help="Applied only for biexponential trasformation.The value for widthBasis will determine the amount of channels to be compressed into linear space around zero.  The space of linear does not change, but rather the number of channels or bins being compressed into the linear space. Default value is setting to -100.")
)
opt = parse_args(OptionParser(option_list=option_list))
# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
                       option_list=option_list, 
                       description="TO DO",
                       epilogue="For comments, bug reports etc... please contact Simone Puccio <simone.puccio@humanitasresearch.it>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options


fcsDir <- opt$fcsDir                         # path to the directory containing fcs files
phenoData <- opt$phenoData                   # path to the design/info file
markerDescription <- opt$markerDescription # transformation definitions for the fluorescent channels
trasformationType <- opt$trasformationType   # transformation definitions for the fluorescent channels
projectName <- opt$projectName               # name of the project
outputDir <- opt$outputDir                   # path to the directory containing raw counts files
widthBasis <- opt$widthBasis
#

print(paste("projectName", projectName))
print(paste("phenoData", phenoData))
print(paste("fcsDir", fcsDir))
print(paste("markerDescription", markerDescription))
print(paste("trasformationType", trasformationType))
print(paste("outputDir", outputDir))
################################################################################
###                             running script                               ###
################################################################################

checkParameters <- function(fcsDir,phenoData,markerDescription,trasformationType,projectName,outputDir){
  problem <- FALSE
  if (!is.character(projectName) | length(projectName)!=1){
    message("projectName must be a character vector of minimum length 3")
    problem <- TRUE
  }
  if (!is.character(phenoData) | length(phenoData)!=1 || !file.exists(phenoData)){
    message("phenoData must be a character vector of length 1 specifying an accessible file")
    problem <- TRUE
  }
  if (!is.character(fcsDir) | length(fcsDir)!=1 || is.na(file.info(fcsDir)[1,"isdir"]) | !file.info(fcsDir)[1,"isdir"]){
    message("fcsDir must be a character vector of length 1 specifying an accessible directory")
    problem <- TRUE
  }  
  if (!is.character(markerDescription) | length(markerDescription)!=1 || !file.exists(markerDescription)){
    message("markerDescription must be a character vector of length 1 specifying an accessible file")
    problem <- TRUE
  }
  if (!is.character(trasformationType) | length(trasformationType)!=1 || !I(trasformationType %in% c("log","arcsinh","biexponential","logicle"))){
    message("trasformationType must be equal to 'log', 'arcsinh','biexponential','logicle'")
    problem <- TRUE
  }
  if (!is.character(outputDir) | length(outputDir)!=1 || !file.exists(outputDir)){
    message("outputDir must be a character vector of length 1 specifying an accessible file")
    problem <- TRUE
  }
  if (!problem){
    print("All the parameters are correct")
  }
  return(invisible(problem))
}
###################################
loadInfoFile <- function(phenoData){
  target <- read.table(phenoData, header=TRUE, sep="\t", na.strings="")
  if (!I("File_name" %in% names(target))) stop(paste("The factor of interest", "File_name", "is not in the target file"))
  if (any(is.na(target[,c("File_name")]))) stop("NA are present in the File_name column")
  if (!I("Experiment_name" %in% names(target))) stop(paste("The factor of interest", "Experiment_name", "is not in the target file"))
  if (any(grepl("[[:punct:]]", as.character(target[,"Experiment_name"])))) stop(paste("The", "Experiment_name", "variable contains punctuation characters, please remove them"))
  if (!I("Conditions" %in% names(target))) stop(paste("The factor of interest", "Conditions", "is not in the target file"))
  if (any(grepl("[[:punct:]]", as.character(target[,"Conditions"])))) stop(paste("The", "Conditions", "variable contains punctuation characters, please remove them"))
  if (!I("Instrument" %in% names(target))) stop(paste("The factor of interest", "Instrument", "is not in the target file"))
  if (any(grepl("[[:punct:]]", as.character(target[,"Instrument"])))) stop(paste("The", "Instrument", "variable contains punctuation characters, please remove them")) 
  if (!I("Dosages" %in% names(target))) stop(paste("The factor of interest", "Dosages", "is not in the target file"))
  if (any(grepl("[[:punct:]]", as.character(target[,"Dosages"])))) stop(paste("The", "Dosages", "variable contains punctuation characters, please remove them"))
  if (!I("Timepoint" %in% names(target))) stop(paste("The factor of interest", "Timepoint", "is not in the target file"))
  if (any(grepl("[[:punct:]]", as.character(target[,"Timepoint"])))) stop(paste("The", "Timepoint", "variable contains punctuation characters, please remove them"))
  if (!I("Individuals" %in% names(target))) stop(paste("The factor of interest", "Individuals", "is not in the target file"))
  if (any(grepl("[[:punct:]]", as.character(target[,"Individuals"])))) stop(paste("The", "Individuals", "variable contains punctuation characters, please remove them"))
  if (!I("Date" %in% names(target))) stop(paste("The factor of interest", "Date", "is not in the target file"))
  if (any(grepl("[[:punct:]]", as.character(target[,"Date"])))) stop(paste("The", "Date", "variable contains punctuation characters, please remove them"))
  if (!I("CompensationMatrix" %in% names(target))) stop(paste("The factor of interest", "CompensationMatrix", "is not in the target file"))
  cat("InfoFile file:\n")
  print(target)
  pheno_data <- AnnotatedDataFrame(read.csv(phenoData, row.names = 1, stringsAsFactors = F,sep = "\t"))
  return(pheno_data)
}
###################################
CreateFolder <- function(outputDir,typeofoutput){
    if(dir.exists(outputDir)){
        dir.create(paste(outputDir,typeofoutput,sep="/"), showWarnings = FALSE)
        return(paste(outputDir,typeofoutput,sep="/"))
    }
    else {stop("Output directory is not valid.") }
}
###################################
CheckChannel <- function(fcsDir,markerDescription){
    df <- read.csv(markerDescription,header = FALSE)
    files <- list.files(fcsDir, pattern = "\\.fcs$", full.names = TRUE)
    if (any(is.na(df[,c("V1","V2")]))) {
        stop("Exit na value found.")
    }
    for (lines in files) {
        x <- read.FCS(lines, transformation=FALSE)
        if (length(unique(df[,c("V2")])) != length(unique(x@parameters@data$name))){
            stop("Exit number of marker is different to number of channels.")
        }
        else { map <- data.frame(alias = df$V1, channels = df$V2 )
        }
    }
return(map)
}
###################################
BioexpTrasformation <- function(fcsDir,phenoData,alias,outputdir,widthBasis){
    setwd(fcsDir)
    pheno_data <- AnnotatedDataFrame(read.csv(file=phenoData, row.names = 1, stringsAsFactors = F,sep = "\t"))
    if (dir.exists(fcsDir)){
        files <- list.files(fcsDir, pattern = "\\.fcs$", full.names = TRUE)
        fs <- read.flowSet(files=files,phenoData=pheno_data,transformation=FALSE) 
        gs <- GatingSet(fs)
        #print("#################")
        #print(head(exprs(gs_cyto_data(gs)[[1]])))
        biexpTrans <- flowjo_biexp_trans(channelRange=4096, maxValue=262144
, pos=4.5,neg=0, widthBasis=widthBasis)
        chnls <- parameters(fs[[1]])
        tf <- transformerList(colnames(exprs(gs_cyto_data(gs)[[1]]))[grep('^*[FS]SC', as.vector(colnames(exprs(gs_cyto_data(gs)[[1]]))),invert = TRUE)], biexpTrans)
        gs <- transform(gs, tf)
        #logTrans <- flowjo_log_trans(decade = 3, offset = 30)
        #print("#################")
        #print(head(exprs(gs_cyto_data(gs)[[1]])))
        tf <- transformerList(colnames(exprs(gs_cyto_data(gs)[[1]]))[grep('^*[FS]SC', as.vector(colnames(exprs(gs_cyto_data(gs)[[1]]))),invert = TRUE)], biexpTrans)
        gs <- transform(gs, tf)
        #print("#################")
        #print(head(exprs(gs_cyto_data(gs)[[1]])))
        write.flowSet(gs_cyto_data(gs), outputdir)
    }
}
###################################
CompensationCheck <- function (fcsDir,phenoData,outputdir,alias){
    if (dir.exists(fcsDir)){
        files <- list.files(fcsDir, pattern = "\\.fcs$", full.names = TRUE)
        for (lines in files) {
            x <- read.FCS(lines, transformation=FALSE,channel_alias = alias)
            if (isFCSfile(lines)) {
                if (basename(lines) %in% row.names(phenoData) && !is.na(pData(phenoData)[basename(lines),]$CompensationMatrix)){
                    comp.mat <- read.table(pData(phenoData)[basename(lines),]$CompensationMatrix, header=TRUE, skip=2,sep="\t", check.names = FALSE)
                    comp <- compensation(comp.mat)
                    temp <- compensate(x, comp)
                    write.FCS(temp,paste(outputdir,basename(lines), sep="/"), what="numeric", delimiter = "\\")
                }
                else if (basename(lines) %in% row.names(phenoData) && is.na(pData(phenoData)[basename(lines),]$CompensationMatrix)){
                    write.FCS(x,paste(outputdir,basename(lines), sep="/"), what="numeric", delimiter = "\\")
                }
                else if (!(basename(lines) %in% row.names(phenoData))){
                    stop("Error File is not listed in the InfoFile.")
                }}
            else {
                stop("File is not FCS file.")
            }
        }
    }
    else {
        stop("Input FCS directory does not exists.") 
    }
}
# checking parameters
problem <- checkParameters(fcsDir=fcsDir,phenoData=phenoData,markerDescription=markerDescription,trasformationType=trasformationType,projectName=projectName,outputDir=outputDir)
if (problem) quit(save="yes")
################ 
# loading info file
target <- loadInfoFile(phenoData=phenoData)
#
alias <- CheckChannel(fcsDir=fcsDir,markerDescription=markerDescription)
# create directory
CompensateDir <- CreateFolder(outputDir=outputDir,typeofoutput="FCScompensated")
# compensate 
CompensationCheck(fcsDir=fcsDir,phenoData=target,outputdir=CompensateDir,alias=alias)
#
TraformedDir <- CreateFolder(outputDir=outputDir,typeofoutput="FCStrasformed")
#
if (trasformationType %in% "biexponential"){
    BioexpTrasformation(fcsDir=outputDir,phenoData=phenoData,alias=alias,outputdir=TraformedDir,widthBasis=widthBasis)
} 



