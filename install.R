if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install(c("qvalue", "limma", "impute", "imputeLCMD",
                       "pcaMethods", "biomaRt"))


install.packages(c("sourcetools", "missForest", "VIM", "ggplot2", "ggridges", 
                  "car", "Hmisc", "mice", "ggfortify", "reshape2", "cluster", 
                   "stringr", "dplyr", "factoextra", "FactoMineR", "naniar", 
                   "finalfit", "matrixTests", "ggrepel", "tidyr", 
                   "viridis", "hrbrthemes"))
