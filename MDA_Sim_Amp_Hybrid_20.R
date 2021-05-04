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
tag = "MDA_Sim_Amp_LowMV"
tag2 = "MDA_Sim_Amp_KNN_LowMV"
tag3 = "MDA_Sim_Amp_MLE_LowMV"
tag4 = "MDA_Sim_Amp_SVD_LowMV"
tag5 = "MDA_Sim_Amp_MinDet_LowMV"
tag6 = "MDA_Sim_Amp_MinProb_LowMV"
tag7 = "MDA_Sim_Amp_QRILC_LowMV"
tag8 = "MDA_Sim_Amp_Hybrid_LowMV"
limma_norm = "Quantile" #limma_norm - "Quantile", "none"
pvalue_threshold = 0.05
qvalue_threshold = 0.05
max_rows_to_display = 100
max_columns_to_display = 50

options(repr.matrix.max.cols=max_rows_to_display, repr.matrix.max.rows=max_columns_to_display)


#READ SIM DATA FILE
df <- read.csv("MDA_Sim_Amp_Data_LowMV.csv", stringsAsFactors = F)
head(df)


#CHECK COLUMN NAMES - HAS TO MATCH TREATMENT NAMES!!!! 
cat("Use the following columns names for selection in next cell\n")
cat(paste('"',paste(sort(colnames(df[c(-1)])),collapse='","'),'"',sep=""))


#GROUP DATA BY COLUMN NAMES
MDA_GD_sim = c("MDA_GD1_sim","MDA_GD2_sim","MDA_GD3_sim")
MDA_HG_sim  = c("MDA_HG1_sim","MDA_HG2_sim","MDA_HG3_sim")


#CREATE FACTORS FROM DATA 
fac.all <- as.factor(c(rep("MDA_GD_sim",3), rep("MDA_HG_sim",3)))
data_columns = c(MDA_GD_sim, MDA_HG_sim)
all_columns = c("ProteinID", data_columns)


#SPECIFY APPROPRIATE FACTOR LABELS FOR THE SAMPLE NAMES ABOVE - NAMES NEED TO MATCH 
cat("Factors for data table are:\n", paste(fac.all,collapse=","), "\n")
df2_all = df[all_columns]
colnames(df2_all)


#ADD NUMBER OF MISSING FOR EACH SAMPLE AND TOTAL 
my_data = df2_all
my_data[is.na(my_data)] <- 0
df_combined <-  subset(my_data,select=c(ProteinID))
df_combined$num_missing <- rowSums(my_data == 0)
df_combined$a_missing <-rowSums(my_data[2:4] == 0)
df_combined$b_missing <- rowSums(my_data[5:7] == 0)
df_combined$both_missing <- df_combined$a_missing > 0 & df_combined$b_missing > 0


#REPLACE 0s WITH NAs
df2_all_filt_temp = as.matrix(df2_all[, -1])
rownames(df2_all_filt_temp) <- df2_all[, 1]
sel_data_imp <- df2_all_filt_temp
sel_data_imp[sel_data_imp == 0] <- NA
sum(is.na(sel_data_imp))
nrow(sel_data_imp)


#CONVERT DATA TO MATRIX FOR IMPUTATION 
exprsData.MD <- as.matrix(sel_data_imp)


#Hybrid IMPUTATION
#CREATE LIST TO CONCATENATE ROWS
datalist = list() 


#PERFORM HYBRID IMP
for (i in seq(1, maxiterations))  {
        #HYBRID MODEL SELECTOR (MAR, MNAR) - 
        #CHOOSE FIRST COMPARISON COLUMNS FOR BIOREPS
        sel_data_imp_a = sel_data_imp[,1:3]
        m.s.all.a = model.Selector(sel_data_imp_a)
        
        df.mis.all.ms.a = impute.MAR.MNAR(sel_data_imp_a, m.s.all.a, 
                                          method.MAR = mar_imputation,
                                          method.MNAR = mnar_imputation)
        
        
        #HYBRID MODEL SELECTOR (MAR, MNAR) 
        #CHOOSE SECOND COMPARISON COLUMNS FOR BIOREPS 
        sel_data_imp_b = sel_data_imp[,4:6]
        m.s.all.b = model.Selector(sel_data_imp_b)
        df.mis.all.ms.b = impute.MAR.MNAR(sel_data_imp_b, m.s.all.b, 
                                          method.MAR = mar_imputation,
                                          method.MNAR = mnar_imputation)
        if (i == 1) {
                df.threshold = data.frame(a_threshold = m.s.all.a[[2]],b_threshold=m.s.all.b[[2]])
        }
        else {
                df.threshold[nrow(df.threshold) + 1,] = c(m.s.all.a[[2]],m.s.all.b[[2]])
        }
        
        df_merge = data.frame(df.mis.all.ms.a,df.mis.all.ms.b)
        nrow(df_merge)
        head(df_merge)
        write.csv(df_merge, paste0("Post_Imp_", tag8, i, ".csv"), row.names = F)
        write.csv(df.threshold, paste0("Model_Selector_", tag8, ".csv"), row.names = F)
        exprsData.imputed <- as.matrix(df_merge)
        
        
        
        #VISUALIZE DATA BEFORE AND AFTER IMP (HISTOGRAM)
        pdf(paste0("Histogram_Overlay_", tag8, ".pdf"))
        
        hist(as.matrix(exprsData.MD[,]),
             main  = "RAW Intensity: Sim Amp MDA",
             xlab = "Log2 Peptide Intensities",
             col = "red4")
        box()
        
        hist(exprsData.imputed, 
             main = "Hybrid Imputation: Sim Amp MDA",
             xlab = "Log2 Peptide Intensities",
             col = "blue4", add = F)
        box()
        
        hist(exprsData.imputed,
             main  = "Hybrid Imputation: Sim Amp MDA",
             xlab = "Log2 Peptide Intensities",
             col = "blue4")
        
        hist(as.matrix(exprsData.MD[,]), 
             #breaks = 25,
             xlab = "Log2 Peptide Intensities",
             col = "red4",
             add = T)
        legend("topright", 
               c("Raw Intensity", "Post Hybrid Imputation"), 
               fill= c("red4","blue4"), 
               bty = "n")
        box()
        
        dev.off()
        
        
        #VISUALIZE DATA BEFORE AND AFTER IMP (DENSITY PLOTS)
        df_mis_sel <- sel_data_imp
        df_mis_sel[is.na(df_mis_sel)] <- 0
        plot_raw1_sel = melt(df_mis_sel)
        
        
        #BEFORE IMP
        pdf(paste0("Density_Plot_", tag8, ".pdf"))
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
                xlab("Distribution Post Hybrid Imputation: Sim Amp MDA") +
                ylab("Treatments") +
                theme_bw() +
                theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))
        
        
        #OVERLAY
        ggplot() + 
                geom_density_ridges2(data = plot_raw1_sel, aes(x = value, y = Var2, height = ..density..), 
                                     stat = "density_ridges", fill = alpha("darkred", 0.7)) +
                geom_density_ridges2(data = plot_raw2_all, aes(x = value, y = Var2, height = ..density..), 
                                     stat = "density_ridges", fill = alpha("darkblue", 0.5)) +
                xlab("Distribution Overlay Hybrid Imputation: Sim Amp MDA") +
                ylab("Treatments") +
                theme_bw() +
                theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))
        
        dev.off()
        
        
        #PCA AND SCREE PLOTS FOR DATASET POST-IMPUTATION
        df_merge_temp <- exprsData.imputed
        pdf(paste0("Scree_PCA_Post_Imp_",tag8, ".pdf"))
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
                labs(title="Post Hybrid Imputation Intensities:\n Sim Amp MDA",
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
                labs(title = "Variance Analysis of Factors\nPost Hybrid Imputation: Sim Amp MDA",
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
        pdf(paste0("Violin_Box_Post_Imp_", tag8, ".pdf"))
        ggplot(data_box3, aes(x=ind, y=values, fill=ind)) +
                geom_violin(trim = F) +
                geom_boxplot(width=0.25, fill="white") +
                ggtitle("Intensity Distribution Post Hybrid: Sim Amp MDA") +
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
        pdf(paste0("Pvalues_", tag8, ".pdf"))
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
        pdf(paste0("Qvalues_", tag8, ".pdf"))
        plot(qobj)
        dev.off()
        
        
        #KNIT COLUMNS BACK TOGETHER
        temp_cnames_imp <- colnames(df_merge_temp) 
        temp_cnames_norm <- colnames(exprSet_quantile) 
        colnames(df_merge_temp) <- paste(colnames(df_merge_temp), "imp", sep = "_")
        colnames(exprSet_quantile) <- paste(colnames(exprSet_quantile), "quant_norm", sep = "_")
        output_data = cbind(df2_all, df_merge_temp, exprSet_quantile, data.frame(topTable(fit, coef=2, sort="none", n=Inf)))
        output_data['qvalue'] <- qvalues
        
        temp_ids = stringr::str_match(row.names(output_data),'^([a-z]*\\|([a-zA-Z0-9]*)\\|)([a-zA-Z0-9]*)_[a-zA-Z0-9]*')
        colnames(temp_ids) <- c("Protein_Description", "UniprotID", "UniprotAC", "Protein_Name")
        output_data = cbind(temp_ids, output_data, protein_IDs=row.names(output_data))
        colnames(exprsData.imputed) <-  temp_cnames_imp 
        colnames(exprSet_quantile) <-  temp_cnames_norm 
        
        
        # OUTPUT RESULTS TO CSV
        write.csv(output_data[order(output_data$qvalue),], 
                  paste0("DE_Ordered_", tag8, i, ".csv"), row.names = F)
        
        
        #ADD NUMBER OF MISSING VALUES FOR EACH PROTEIN BY SAMPLE TO DATA
        temp_data <- output_data[, c("Protein_Description", "logFC", "P.Value", "qvalue", "AveExpr")]
        datalist[[i]] <- temp_data
        names(temp_data) = paste0(names(temp_data) ,"_",i)
        df_combined <- merge(df_combined, temp_data, by.x = c("ProteinID"),by.y=c(paste0("Protein_Description" ,"_",i)), 
                             row.names = F)
        
}


#RBIND DATA FOR LONG FORMAT
big_data = do.call(rbind, datalist)
write.csv(big_data, paste0("Long_data_", tag8, ".csv"), row.names = F)


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
write.csv(df_combined, paste0("Combined_logFC_qvalue_pvalue_", tag8, ".csv"), 
          row.names = F)
pdf(paste0("Boxplots_Variance_", tag8, ".pdf"))


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

annotate_figure(figure, top = text_grob("Hybrid Imputation Sim Amp MDA: LowMV", face = "bold", size = 18), 
                fig.lab = "", fig.lab.face = "bold")

dev.off()

