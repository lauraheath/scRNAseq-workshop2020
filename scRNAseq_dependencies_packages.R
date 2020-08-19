#dependencies to run before installing packages:
#apt-get update -y\
#apt-get install -y curl dpkg-dev zlib1g-dev libssl-dev libffi-dev curl nano less libcurl4-openssl-dev gawk build-essential chrpath libssl-dev libxft-dev libfreetype6 libfreetype6-dev libfontconfig1 libfontconfig1-dev libudunits2-dev libgdal-dev libpython-dev libpython3-dev libgsl-dev

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase", "BiocGenerics","DelayedArray","DelayedMatrixStats","limma","S4Vectors","SingleCellExperiment","SummarizedExperiment", "batchelor","edgeR","qvalue", "multtest", "scater", "scran"))
#install.packages("devtools")  #may already be installed, check
install.packages("boot")
install.packages("foreign")
install.packages("MASS")
install.packages("nlme")
install.packages("Rfast")
install.packages("Seurat")
install.packages("readxl")
install.packages("rqdatatable")
install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
devtools::install_url('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.7.tar.gz')
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
