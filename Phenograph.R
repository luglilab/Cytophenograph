options( warn = -1 )
suppressMessages(require("Rphenograph"))
suppressMessages(require("optparse"))
set.seed(123456)
#
option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="inputmatrix"),
  make_option(c("-k", "--kvalue"), action="store", default=NA, type='integer',
              help=""),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help=""),
  make_option(c("-n", "--name"), action="store", default=NA, type='character',
              help="")
)
opt = parse_args(OptionParser(option_list=option_list))
# read input
inputmatrix <- read.csv(opt$input,sep="\t")
# run phenograph
cluster_PhenoGraph <- Rphenograph(as.data.frame(inputmatrix), k=opt$kvalue)
# list with cluster
cluster_PhenoGraph<- membership(cluster_PhenoGraph[[2]])
# append column with cluster
inputmatrix[ , "Phenograph"] <- cluster_PhenoGraph
#
outname <- paste(opt$name,".tsv",sep="")
#
write.table(inputmatrix,paste(opt$output,outname,sep ="/"),sep = "\t",quote = F,row.names = F)

