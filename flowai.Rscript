require(flowAI)
require(flowCore)
options(warn=-1)

args = commandArgs(trailingOnly = TRUE)
fcs = args[1]
outputfolder = args[2]

frames <- read.FCS(fcs)
as(frames, "flowSet")

resQC <- flow_auto_qc(frames,fcs_QC = "_concatenate_after_QC",output = 2,mini_report = FALSE,
                      second_fractionFR=0.01,
                      folder_results =outputfolder) 
