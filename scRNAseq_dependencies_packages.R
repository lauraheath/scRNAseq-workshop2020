#dependencies to run in terminal window before installing packages:
#sudo apt-get update -y
#sudo apt-get install -y libudunits2-dev libgdal-dev

## also set up synapse login configuration in terminal:
#nano .synapseConfig
##in nano window type the following:
#[authentication]
#username = <username>
#password = <password>

BiocManager::install(c("Biobase", "BiocGenerics","DelayedArray","DelayedMatrixStats","limma","S4Vectors","SingleCellExperiment","SummarizedExperiment", "batchelor", "Matrix.utils"))
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")

library(dplyr)
library(ggplot2)
library(monocle3)
