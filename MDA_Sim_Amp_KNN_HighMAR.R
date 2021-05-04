#LOAD LIBRARIES FOR ANALYSIS
pkg_str <- c("sourcetools", "qvalue",
             "limma", "missForest",
             "VIM", "ggplot2", "ggridges", 
             "car", "Hmisc", "mice", "ggfortify", 
             "reshape2", "cluster", "stringr",
             "dplyr", "factoextra", "FactoMineR",
             "naniar","finalfit", "impute",
             "pcaMethods", "imputeLCMD", "matrixTests", 
             "ggrepel", "biomaRt", "tidyr", "dplyr", 
             "viridis", "hrbrthemes", "ggpubr", 
             "gridExtra", "ggplotify", "fBasics")

for (i in seq_along(pkg_str)) {
     library(pkg_str[i], character.only = T)
}


#ROUND DECIMAL PLACES
roundUp <- function(x,to=10)
{
     to*(x%/%to + as.logical(x%%to))
}

impute.wrapper.MLE.rnd = function (dataSet.mvs) 
{
     s <- prelim.norm(dataSet.mvs)
     thetahat <- em.norm(s, showits = FALSE)
     rngseed(sample(1:1000,1))
     dataSet.imputed <- imp.norm(s, thetahat, dataSet.mvs)
     return(dataSet.imputed)
}

impute.wrapper.KNN.rnd = function (dataSet.mvs,I)
{
     resultKNN = impute.knn(dataSet.mvs, k = I, rowmax = 0.99,
                            colmax = 0.99, maxp = 1500, rng.seed = sample(1:362436069,1))
     dataSet.imputed = resultKNN[[1]]
     return(dataSet.imputed)
}

impute.MAR.rnd = function (dataSet.mvs, model.selector, method = "MLE",I) 
{
     if (length(which(model.selector[[1]] == 1)) == 0) {
          dataSet.imputed = dataSet.mvs
     }
     else {
          dataSet.MCAR = dataSet.mvs[which(model.selector[[1]] == 
                                                1), ]
          switch(method, MLE = {
               dataSet.MCAR.imputed = impute.wrapper.MLE.rnd(dataSet.MCAR)
          }, SVD = {
               dataSet.MCAR.imputed = impute.wrapper.SVD(dataSet.MCAR, 
                                                         K = I)
          }, KNN = {
               dataSet.MCAR.imputed = impute.wrapper.KNN.rnd(dataSet.MCAR,I)
          })
          dataSet.imputed = dataSet.mvs
          dataSet.imputed[which(model.selector[[1]] == 1), ] = dataSet.MCAR.imputed
     }
     return(dataSet.imputed)
}


#CHANGE DATAFILE TO DIRECTORY AND FILE NAME, ADDITIONAL PARAMETERS TO CHANGE
now <- format(Sys.time(), "%m-%d-%Y_%H-%M-%S")
datafile = "MDA_Sim_Data.csv"

maxiterations = 25
min_num_observations = 0
min_counts = 2^15
mar_imputation = "KNN"
mnar_imputation = "QRILC"
tag = "MDA_Sim_Amp_HighMAR"
tag2 = "MDA_Sim_Amp_KNN_HighMAR"
tag3 = "MDA_Sim_Amp_MLE_HighMAR"
tag4 = "MDA_Sim_Amp_SVD_HighMAR"
tag5 = "MDA_Sim_Amp_MinDet_HighMAR"
tag6 = "MDA_Sim_Amp_MinProb_HighMAR"
tag7 = "MDA_Sim_Amp_QRILC_HighMAR"
tag8 = "MDA_Sim_Amp_Hybrid_HighMAR"
limma_norm = "Quantile" #limma_norm - "Quantile", "none"
pvalue_threshold = 0.05
qvalue_threshold = 0.05
max_rows_to_display = 100
max_columns_to_display = 50

options(repr.matrix.max.cols=max_rows_to_display, repr.matrix.max.rows=max_columns_to_display)


#READ DATA FILE
df = read.csv(datafile,stringsAsFactors = F)
head(df)
colnames(df)


#FILTER DATA TO SELECT PROTEIN ID AND INTENSITIES
df2 <- df[, c(1, 17:22)]
df2


#CHECK COLUMN NAMES - HAS TO MATCH TREATMENT NAMES!!!! 
cat("Use the following columns names for selection in next cell\n")
cat(paste('"',paste(sort(colnames(df2[c(-1)])),collapse='","'),'"',sep=""))


#GROUP DATA BY COLUMN NAMES
MDA_GD_sim = c("MDA_GD1_sim","MDA_GD2_sim","MDA_GD3_sim")
MDA_HG_sim  = c("MDA_HG1_sim","MDA_HG2_sim","MDA_HG3_sim")


#CREATE FACTORS FROM DATA 
fac.all <- as.factor(c(rep("MDA_GD_sim",3), rep("MDA_HG_sim",3)))
data_columns = c(MDA_GD_sim, MDA_HG_sim)
all_columns = c("ProteinID", data_columns)


#SPECIFY APPROPRIATE FACTOR LABELS FOR THE SAMPLE NAMES ABOVE - NAMES NEED TO MATCH 
cat("Factors for data table are:\n", paste(fac.all,collapse=","), "\n")
df2_all = df2[all_columns]
df3 = df2_all[, -1]
row.names(df3) = df2_all[, 1]
df3[is.na(df3)] <- 0


#SIMULATE THE DATA FOR MISSING VALUES FOR MAR PATTERN AND MNAR PATTERN
result <- ampute(data = df3, mech = "MAR")
my_pattern <- result$pattern


result2 <- ampute(data = df3, mech = "MNAR")
my_pattern2 <- result2$pattern


#AMPUTE THE COMPLETE DATA ONCE FOR EVERY MECHANISM
ampdata1 <- ampute(df3, patterns = my_pattern, 
                   prop = 0.8, mech = "MAR")$amp
ampdata2 <- ampute(df3, patterns = my_pattern2, 
                   prop = 0.2, mech = "MNAR")$amp


#CREATE A RANDOM ALLOCATION VECTOR
# USE PROB TO SPECIFY HOW MUCH OF EACH MECH SHOULD BE CREATED
#HERE, MAR = 0.5 & MNAR = 0.5 
indices <- sample(x = c(1, 2), size = nrow(df3), 
                  replace = TRUE, prob = c(0.5, 0.5))


#CREATE AN EMPTY DATA MATRIX
#FILL MATRIX W/VALUES FROM EITHER OF THE TWO AMPUTED DATASETS
ampdata <- matrix(NA, nrow = nrow(df3), 
                  ncol = ncol(df3))
ampdata[indices == 1, ] <- as.matrix(ampdata1[indices == 1, ])
ampdata[indices == 2, ] <- as.matrix(ampdata2[indices == 2, ])
colnames(ampdata) <- colnames(df3)
temp_ampdata <- as.data.frame(ampdata)
rownames(temp_ampdata) <- rownames(df3)
temp_ampdata  <- tibble::rownames_to_column(temp_ampdata, "ProteinID")

#OUTPUT INDICES FOR "1" MAR AND "2" NULL
#CONVERT TO DATAFRAME AND FILTER FOR TYPE OF MISSINGNESS
MV_index <- as.data.frame(indices)
names(MV_index) <- "Type_MV"
MV_output <- cbind(temp_ampdata, MV_index)


MAR_Out <- MV_output %>%
        dplyr:: filter(Type_MV == 1) 
MAR_Out <- MAR_Out %>%
        dplyr:: mutate(Total = sum(is.na(MAR_Out)))


MNAR_Out <- MV_output %>%
        dplyr:: filter(Type_MV == 2)
MNAR_Out <- MNAR_Out %>%
        dplyr:: mutate(Total = sum(is.na(MNAR_Out)))

All_MV <- rbind(MAR_Out, MNAR_Out)
MV_Final <- merge(MV_output, All_MV, by = c("ProteinID", "MDA_HG1_sim", "MDA_HG2_sim", "MDA_HG3_sim", 
                                            "MDA_GD1_sim", "MDA_GD2_sim", "MDA_GD3_sim", "Type_MV"))

write.csv(MV_Final, "AmputedData_MVType_HighMAR.csv")

data2 <- ampdata
rownames(data2) <- rownames(df3)
data2[data2 == 0] <- NA
df4 = data.frame(ProteinID = row.names(df3), data2, 
                 row.names = NULL)
write.csv(df4, "MDA_Sim_Amp_Data_HighMAR.csv", row.names = F)

#PLOT MISSING DATA PATTERN AND ADD UP NUMBER OF MISSING VALUES
pdf("MD_Pattern_Sim_Amp_HighMAR.pdf")
z <- md.pattern(data2, rotate.names = TRUE)
sum(is.na(ampdata))
dev.off()


#PLOT UNNORMALIZED DATA AS INTENSITY VALUES 
pdf(paste0("Density_Plots_", tag, ".pdf"))
ampdata3 = melt(data2)

ggplot(ampdata3, aes(x = value, y = Var2, height = ..density..)) + 
        geom_density_ridges2(stat = "density_ridges", fill = alpha("darkred", 0.7)) +
        xlab("Distribution of Intensity Values") +
        ylab("Treatments") +
        theme_bw() +
        theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))

dev.off()


#MISSING VALUE COMBINED VIOLIN + BOX PLOT
pdf(paste0("Violin_Box_Prenorm", tag, ".pdf"))
ggplot(ampdata3, aes(x=Var2, y=value, fill=Var2)) +
        geom_violin(trim = F) +
        geom_boxplot(width=0.25, fill="white") +
        ggtitle("Pre Imputation\n") +
        xlab("") +
        ylab ("log2 Transformed Intensity") +
        scale_fill_viridis(discrete = TRUE) +
        theme_bw() +
        theme(legend.position="none",
              plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
              axis.text.x = element_text(face = "bold", color = "black", angle = 90),
              axis.text.y = element_text(face = "bold", color = "black"), 
              axis.title.y = element_text(face = "bold", color = "black", hjust = 0.5, size = 12)
                
        ) 
dev.off()


#VISUALIZE MISSING VALUES IN THE DATA 
pdf(paste0("Missing_Data_Exam_", tag, ".pdf"), width = 12)


#NUMBER OF MISSING VALUES BY PROTEIN ID
gg_miss_case(as.data.frame(data2))  + labs(x = "Protein ID", y = "Number of Missing Values") +
        theme_minimal() +
        theme(axis.text = element_text(face = "bold", size = 12, family = "Helvetica", color = 'black'), 
              axis.title = element_text(face = "bold", size = 18, family = "Helvetica")) 


#CUMULATIVE SUM OF MISSING VALUES BY PROTEIN ID
gg_miss_case_cumsum(as.data.frame(data2)) + labs(x = "Protein ID", y = "Missing Cumulative Sum") +
        theme_minimal() +
        theme(axis.text = element_text(face = "bold", size = 12, family = "Helvetica", color = 'black'), 
              axis.title = element_text(face = "bold", size = 18, family = "Helvetica"),
              axis.text.x= element_blank(), 
              axis.ticks.x = element_blank()) 


#NUMBER OF MISSING VALUES BY SAMPLE
gg_miss_var(as.data.frame(data2)) + labs(x = "Sample", y = "Number of Missing Values") +
        theme_minimal() +
        theme(axis.text = element_text(face = "bold", size = 12, family = "Helvetica", color = 'black'), 
              axis.title = element_text(face = "bold", size = 18, family = "Helvetica"))


#CUMULATIVE SUM OF MISSING VALUES BY SAMPLE
gg_miss_var_cumsum(as.data.frame(data2)) + labs(x = "Sample", y = "Missing Cumulative Sum") +
        theme_minimal() +
        theme(axis.text = element_text(face = "bold", size = 12, family = "Helvetica", color = 'black', angle = 45), 
              axis.title = element_text(face = "bold", size = 18, family = "Helvetica"))


#MISSING VALUE BY SAMPLE BY PROTEIN ID
vis_miss(as.data.frame(data2)) +
        coord_flip() +
        labs(y = "") +
        theme(axis.text = element_text(face = "bold", size = 24, family = "sans", color = "black"), 
              axis.title = element_text(face = "bold", size = 24, family = "sans"), 
              plot.title = element_blank(), 
              axis.text.x= element_blank(), 
              axis.ticks.x = element_blank()
        )


dev.off()


#CALCULATE PRINCIPAL COMPONENTS FOR PCA PLOT FOR UNIMPUTED DATA
pdf(paste0("PCA_Scree_Plots_", tag, ".pdf"))
data2_temp = data2
data2_temp[is.na(data2_temp)] <- 0
PC <- prcomp(data2_temp, scale = TRUE)
scores <- PC$rotation[,1:5]  #SCORES FOR FIRST 5 PCs


# K-MEANS CLUSTERING [ASSUMES 2 CLUSTERS]
km     <- kmeans(scores, centers=2, nstart=10)

ggdata <- data.frame(scores, Cluster=km$cluster, Sample=fac.all)

ggplot(ggdata, aes(x= PC1,y = PC2, col = Sample)) +
        geom_point(size=4)+ 
        stat_ellipse(aes(x = PC1,y=PC2), fill = 'white', 
                     geom="polygon",level=0.95, alpha=0.2) +
        scale_color_manual(values = c("purple4","darkred","blue4")) + 
        
        labs(title="LFQ Intensity Values\n",
             y = 'Principal Component 1', x = "Principal Component 2") +
        
        theme_bw() +
        theme(plot.title = element_text(family = "Helvetica",face = "bold", hjust = 0.5, size = 22),
              axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
              axis.text = element_text(family = "Helvetica", face = "bold", size = 14, color = "black"),
              legend.title = element_text(family = "Helvetica",face = "bold", size = 16),
              legend.text = element_text(family = "Helvetica",face = "plain", size = 13),
              axis.ticks.x=element_blank(),
              axis.line = element_line(color = "black"),
              legend.position = 'right'
        )


#EIGEN VALUES, VARIANCE, AND CUMULATIVE VARIANCE FOR SCREEPLOT 
eig.val <- get_eigenvalue(PC)

ggplot(eig.val, aes(x=1:nrow(eig.val), y=variance.percent)) +
        geom_bar(stat = 'identity', col="black", fill= 'darkgreen') +
        geom_text(aes(label=format(round(variance.percent, digits = 1)), vjust = -.5, hjust = 0), angle = 30) +
        labs(title = "Variance Analysis of Factors\n LFQ Intensities",
             x = "Principal Components",y = "Percentage of Variances") +
        ylim(0, 100) +
        theme_bw() +
        theme(plot.title = element_text(family = "Helvetica",face = "bold", hjust = 0.5, size = 22),
              axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
              axis.text = element_text(family = "Helvetica", face = "bold", size = 14, color = "black"),
              legend.title = element_text(family = "Helvetica",face = "bold", size = 16),
              legend.text = element_text(family = "Helvetica",face = "plain", size = 13),
              axis.line = element_line(color = "black"),
              legend.position = 'right'
        )


dev.off()


#ADD NUMBER OF MISSING FOR EACH SAMPLE AND TOTAL 
my_data = data2_temp
my_data <- data.frame(ProteinID = row.names(data2_temp), my_data, 
                      row.names = NULL)


df_combined <-  subset(my_data,select=c(ProteinID))
df_combined$num_missing <- rowSums(my_data == 0)
df_combined$a_missing <-rowSums(my_data[2:4] == 0)
df_combined$b_missing <- rowSums(my_data[5:7] == 0)
df_combined$both_missing <- df_combined$a_missing > 0 & df_combined$b_missing > 0


#REPLACE 0s WITH NAs
sel_data_imp <- data2_temp
sel_data_imp[sel_data_imp == 0 ] <- NA
sum(is.na(sel_data_imp))
nrow(sel_data_imp)


#CONVERT DATA TO MATRIX FOR IMPUTATION 
exprsData.MD <- as.matrix(sel_data_imp)

#RUN MODEL SELECTOR ON DATA 
#"1" INDICATES MAR/MCAR METHODS (MLE, SVD, KNN) AND 
#"0" INDICATES LEFT-CENSORED MNAR METHODS (QRILC, MinDet, MinProb)
m.s = model.Selector(exprsData.MD)
m.s
length(m.s)
str(m.s)
m.s[1] <- list(rep(1, 3807))  #IMPORT LIST OF 1s SO ALL MAR EQUAL TO LIST OF PROTEINS

#KNN IMPUTATION
#CREATE LIST TO CONCATENATE ROWS 
datalist = list() 


#PERFORM MAR IMP
for (i in seq(1, maxiterations))  {
     exprsData.imputed <- impute.MAR.rnd (exprsData.MD, m.s, method = "KNN", I=i)
     
     if (i == 1) {
          df.threshold = data.frame(a_threshold = m.s[[2]])
     }
     else {
          df.threshold[nrow(df.threshold) + 1,] = c(m.s[[2]])
     }
     
     imputed_data <- as.data.frame(exprsData.imputed)
     write.csv(imputed_data, paste0("Imp_Post_", tag2, i,".csv"), row.names = F)
     write.csv(df.threshold, paste0("Model_Selector_", tag2, ".csv"), row.names = F)
     
     
     #VISUALIZE DATA BEFORE AND AFTER IMP (HISTOGRAM)
     pdf(paste0("Histogram_Overlay_", tag2, ".pdf"))
     
     hist(as.matrix(exprsData.MD[,]),
          main  = "RAW Intensity: Sim Amp MDA",
          xlab = "Log2 Peptide Intensities",
          col = "red4")
     box()
     
     hist(exprsData.imputed, 
          main = "MAR KNN Imputation: Sim Amp MDA",
          xlab = "Log2 Peptide Intensities",
          col = "blue4", add = F)
     box()
     
     hist(exprsData.imputed,
          main  = "MAR KNN Imputation: Sim Amp MDA",
          xlab = "Log2 Peptide Intensities",
          col = "blue4")
     
     hist(as.matrix(exprsData.MD[,]), 
          #breaks = 25,
          xlab = "Log2 Peptide Intensities",
          col = "red4",
          add = T)
     legend("topright", 
            c("Raw Intensity", "Post KNN Imputation"), 
            fill= c("red4","blue4"), 
            bty = "n")
     box()
     
     dev.off()
     
     
     #VISUALIZE DATA BEFORE AND AFTER IMP (DENSITY PLOTS)
     df_mis_sel <- sel_data_imp
     df_mis_sel[is.na(df_mis_sel)] <- 0
     plot_raw1_sel = melt(df_mis_sel)
     
     
     #BEFORE IMP
     pdf(paste0("Density_Plot_", tag2, ".pdf"))
     ggplot() + 
          geom_density_ridges2(data = plot_raw1_sel, aes(x = value, y = Var2, height = ..density..), 
                               stat = "density_ridges", fill = alpha("darkred", 0.7)) +
          xlab("Distribution Pre-Imputation: Sim Amp MDA") +
          ylab("Treatments") +
          theme_bw() +
          theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))
     
     
     #PLOT ATER IMP
     plot_raw2_all = melt(exprsData.imputed)
     ggplot() + 
          geom_density_ridges2(data = plot_raw2_all, aes(x = value, y = Var2, height = ..density..), 
                               stat = "density_ridges", fill = alpha("darkblue", 0.7)) +
          xlab("Distribution Post KNN Imputation: Sim Amp MDA") +
          ylab("Treatments") +
          theme_bw() +
          theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))
     
     
     #OVERLAY
     ggplot() + 
          geom_density_ridges2(data = plot_raw1_sel, aes(x = value, y = Var2, height = ..density..), 
                               stat = "density_ridges", fill = alpha("darkred", 0.7)) +
          geom_density_ridges2(data = plot_raw2_all, aes(x = value, y = Var2, height = ..density..), 
                               stat = "density_ridges", fill = alpha("darkblue", 0.5)) +
          xlab("Distribution Overlay KNN Imputation: Sim Amp MDA") +
          ylab("Treatments") +
          theme_bw() +
          theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))
     
     dev.off()
     
     
     #PCA AND SCREE PLOTS FOR DATASET POST-IMPUTATION
     df_merge_temp <- exprsData.imputed
     pdf(paste0("Scree_PCA_Post_Imp",tag2, ".pdf"))
     PC3 <- prcomp(exprsData.imputed, scale = TRUE)
     scores3 <- PC3$rotation[, 1:5] #SCORES FIRST 5 PCs
     
     
     #K-MEANS CLUSTERING [ASSUME 2 CLUSTERS]
     km3     <- kmeans(scores3, centers=2, nstart=10)
     
     ggdata3 <- data.frame(scores3, Cluster=km3$cluster, Sample=fac.all)
     
     ggplot(ggdata3, aes(x= PC1,y = PC2, col = Sample)) +
          geom_point(size=4)+ 
          stat_ellipse(aes(x = PC1,y=PC2), fill = 'white', 
                       geom="polygon",level=0.95, alpha=0.2) +
          scale_color_manual(values = c("darkred","blue4","green4")) +
          labs(title="Post KNN Imputation Intensities:\n Sim Amp MDA",
               y = 'Principal Component 1', x = "Principal Component 2") +
          theme_bw() +
          theme(plot.title = element_text(family = "Helvetica",face = "bold", hjust = 0.5, size = 22),
                axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
                axis.text = element_text(family = "Helvetica", face = "bold", size = 14, color = "black"),
                legend.title = element_text(family = "Helvetica",face = "bold", size = 16),
                legend.text = element_text(family = "Helvetica",face = "plain", size = 13),
                axis.ticks.x=element_blank(),
                axis.line = element_line(color = "black"),
                legend.position = 'right'
          )
     
     
     #EIGEN VALUES, VARIANCE, AND CUMULATIVE VARIANCE FOR SCREEPLOT
     eig.val3 <- get_eigenvalue(PC3)
     ggplot(eig.val3, aes(x=1:nrow(eig.val3), y=variance.percent)) +
          geom_bar(stat = 'identity', col="black", fill= 'darkgreen') +
          geom_text(aes(label=format(variance.percent, digits = 1), vjust = -.5, hjust = 0), angle = 30) +
          labs(title = "Variance Analysis of Factors\nPost KNN Imputation: Sim Amp MDA",
               x = "Principal Components",y = "Percentage of Variances") +
          ylim(0,100) +
          theme_bw() +
          theme(plot.title = element_text(family = "Helvetica",face = "bold", hjust = 0.5, size = 22),
                axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
                axis.text = element_text(family = "Helvetica", face = "bold", size = 14, color = "black"),
                legend.title = element_text(family = "Helvetica",face = "bold", size = 16),
                legend.text = element_text(family = "Helvetica",face = "plain", size = 13),
                axis.line = element_line(color = "black"),
                legend.position = 'right'
          )
     
     dev.off()
     
     
     #VIOLIN & BOXPLOT POST_IMP
     data_box3 <- stack(as.data.frame(exprsData.imputed))
     pdf(paste0("Violin_Box_Post_Imp_", tag2, ".pdf"))
     ggplot(data_box3, aes(x=ind, y=values, fill=ind)) +
          geom_violin(trim = F) +
          geom_boxplot(width=0.25, fill="white") +
          ggtitle("Intensity Distribution Post KNN: Sim Amp MDA") +
          xlab("") +
          ylab ("log2 Transformed Intensity") +
          scale_fill_viridis(discrete = TRUE) +
          theme_bw() +
          theme(
               legend.position="none",
               plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
               axis.text.x = element_text(face = "bold", color = "black", angle = 90),
               axis.text.y = element_text(face = "bold", color = "black"), 
               axis.title.y = element_text(face = "bold", color = "black", hjust = 0.5, size = 12)
               
          ) 
     dev.off()
     
     
     #PERFORM LIMMA TEST FOR DE MDA_WCP
     pdf(paste0("Pvalues_", tag2, ".pdf"))
     cat(paste("Performing", limma_norm, "normalization on data"))
     #PERFORM QUANTILE NORM OF DATA AND LIMMA ANALYSIS FOR DE
     if (limma_norm == "Quantile") {
          exprSet_quantile = normalizeQuantiles(exprsData.imputed)
     } else {
          exprSet_quantile = (exprsData.imputed)
     }
     
     
     fit <- lmFit(exprSet_quantile, design = model.matrix(~ fac.all))
     fit <- eBayes(fit)
     tt <- topTable(fit, coef=2)
     limmares <- data.frame(dm=coef(fit)[,2], p.value = fit$p.value[, 2])
     results <- decideTests(fit)
     pvalues <- fit$p.value
     qobj <- qvalue(p = pvalues[,2])
     qvalues <- qobj$qvalues
     pi0 <- qobj$pi0
     lfdr <- qobj$lfdr
     
     
     #HISTOGRAM OF PVALUES
     hist(qobj)
     dev.off()
     
     
     #PRINT QVALUE PLOTS FOR QC
     pdf(paste0("Qvalues_", tag2, ".pdf"))
     plot(qobj)
     dev.off()
     
     
     #KNIT COLUMNS BACK TOGETHER
     temp_cnames_imp <- colnames(df_merge_temp) 
     temp_cnames_norm <- colnames(exprSet_quantile) 
     colnames(df_merge_temp) <- paste(colnames(df_merge_temp), "imp", sep = "_")
     colnames(exprSet_quantile) <- paste(colnames(exprSet_quantile), "quant_norm", sep = "_")
     output_data = cbind(data2_temp, df_merge_temp, exprSet_quantile, data.frame(topTable(fit, coef=2, sort="none", n=Inf)))
     output_data['qvalue'] <- qvalues
     
     temp_ids = stringr::str_match(row.names(output_data),'^([a-z]*\\|([a-zA-Z0-9]*)\\|)([a-zA-Z0-9]*)_[a-zA-Z0-9]*')
     colnames(temp_ids) <- c("Protein_Description", "UniprotID", "UniprotAC", "Protein_Name")
     output_data = cbind(temp_ids, output_data, protein_IDs=row.names(output_data))
     colnames(exprsData.imputed) <-  temp_cnames_imp 
     colnames(exprSet_quantile) <-  temp_cnames_norm 
     
     
     # OUTPUT RESULTS TO CSV
     write.csv(output_data[order(output_data$qvalue),], 
               paste0("DE_Ordered_", tag2, i, ".csv"), row.names = F)
     
     
     #ADD NUMBER OF MISSING VALUES FOR EACH PROTEIN BY SAMPLE TO DATA
     temp_data1 <- output_data[, c("Protein_Description", "logFC", "P.Value", "qvalue", "AveExpr")]
     datalist[[i]] <- temp_data1
     names(temp_data1) = paste0(names(temp_data1) ,"_",i)
     df_combined <- merge(df_combined, temp_data1, by.x = c("ProteinID"),by.y=c(paste0("Protein_Description" ,"_",i)), 
                          row.names = F)
     
}


#RBIND DATA FOR LONG FORMAT
big_data = do.call(rbind, datalist)
write.csv(big_data, paste0("Long_data_", tag2, ".csv"), row.names = F)


#CALCULATE LOGFC MEAN AND STDDEV, QVALUE MEAN AND STDDEV
df_combined$logFC_mean = rowMeans(df_combined[,grepl("logFC",names(df_combined))])
df_combined$logFC_stdev = apply(df_combined[,grepl("logFC",names(df_combined))],1,sd)
df_combined$P.Value_mean = rowMeans(df_combined[,grepl("P.Value",names(df_combined))])
df_combined$P.Value_stdev = apply(df_combined[,grepl("P.Value",names(df_combined))],1,sd)
df_combined$qvalue_mean = rowMeans(df_combined[,grepl("qvalue",names(df_combined))])
df_combined$qvalue_stdev = apply(df_combined[,grepl("qvalue",names(df_combined))],1,sd)
df_combined$AveExpr_mean = rowMeans(df_combined[,grepl("AveExpr",names(df_combined))])
df_combined$AveExpr_stdev = apply(df_combined[,grepl("AveExpr",names(df_combined))],1,sd)
df_combined$Factor = ifelse(df_combined$AveExpr_mean > m.s[[2]], "Green", "Red")


#OUTPUT RESULTS AS CSV AND BOXPLOTS
write.csv(df_combined, paste0("Combined_logFC_qvalue_pvalue_", tag2, ".csv"), 
          row.names = F)
pdf(paste0("Boxplots_Variance_", tag2, ".pdf"))


p1 <- as.grob(~boxplot(df_combined$logFC_stdev~df_combined$num_missing*df_combined$both_missing, 
                       ylab = "logFC Std Dev",
                       xlab = "Missing Value Combinations"))

p2 <- as.grob(~boxplot(df_combined$P.Value_stdev~df_combined$num_missing*df_combined$both_missing, 
                       ylab = "p-value Std Dev",
                       xlab = "Missing Value Combinations"))

p3 <- as.grob(~boxplot(df_combined$qvalue_stdev~df_combined$num_missing*df_combined$both_missing,
                       ylab = "q-value Std Dev",
                       xlab = "Missing Value Combinations"))

p4 <- as.grob(~boxplot(df_combined$AveExpr_stdev~df_combined$num_missing*df_combined$both_missing,
                       ylab = "Mean Expression Std Dev",
                       xlab = "Missing Value Combinations"))

figure <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, 
                    labels = c("A", "B", "C", "D"))

annotate_figure(figure, top = text_grob("MAR KNN Imputation Sim Amp MDA: HighMAR", face = "bold", size = 18), 
                fig.lab = "", fig.lab.face = "bold")

dev.off()

