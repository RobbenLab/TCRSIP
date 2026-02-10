####################################################################
####################################################################
#########                                         #################
#########       Generate Analysis plots for data ##################
#########                                        ##################    
###################################################################
####################################################################


setwd("D://Lab/TCRBinding/")

c('dplyr','tidyr','data.table','ggpubr','ggalign','ggpmisc','raster','lattice','matter','corrplot','ggrepel','gplots','RColorBrewer','pheatmap','ggbreak'
require(corrplot)
require(ggrepel)
library(gplots)
require(RColorBrewer)
require(VennDiagram)
require(pheatmap)
require(rjson)
require(cowplot)
require(ggbreak)
require(ggvenn)
require(bio3d)
require(plotly)
require(FactoMineR)
require(factoextra)
require(umap)
require(caTools)
require(pROC)
require(glmnet)
require(slider)
require(imagefx)
require(ggridges)



require(dplyr)
require(tidyr)
require(data.table)
require(Biostrings)
require(stringdist)
require(pracma)
require(akima)
require(stringr)
require(svMisc)
require(ggpubr)
require(ggalign)
require(ggpmisc)
require(raster)
require(lattice)
require(matter)
require(corrplot)
require(ggrepel)
library(gplots)
require(RColorBrewer)
require(VennDiagram)
require(pheatmap)
require(rjson)
require(cowplot)
require(ggbreak)
require(ggvenn)
require(bio3d)
require(plotly)
require(FactoMineR)
require(factoextra)
require(umap)
require(caTools)
require(pROC)
require(glmnet)
require(slider)
require(imagefx)
require(ggridges)
getSubstring <- function(x,s,n){
  return(unlist(lapply(strsplit(x,s),function(x) x[[n]])))
}