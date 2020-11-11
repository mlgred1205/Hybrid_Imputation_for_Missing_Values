
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
             "gridExtra", "ggplotify")

for (i in seq_along(pkg_str)) {
  library(pkg_str[i], character.only = T)
}


#ROUND DECIMAL PLACES

roundUp <- function(x,to=10)
{
  to*(x%/%to + as.logical(x%%to))
}


#CHANGE DATAFILE TO DIRECTORY AND FILE NAME, ADDITIONAL PARAMETERS TO CHANGE

now <- format(Sys.time(), "%m-%d-%Y_%H-%M-%S")
datafile = "Parsed_Data_MDA_HG_GD.txt"

maxiterations = 25
min_num_observations = 3
min_counts = 2^15
limma_norm = "Quantile" #limma_norm - "Quantile", "none"

mnar_imputation = "QRILC" #mnar_imputation - "MinDet", "MinProb", QRILC" 
mar_imputation = "KNN" # mar_imputation - "KNN", MLE", "SVD"

pvalue_threshold = 0.05
qvalue_threshold = 0.05
max_rows_to_display = 100
max_columns_to_display = 50

options(repr.matrix.max.cols=max_rows_to_display, repr.matrix.max.rows=max_columns_to_display)

#LOAD DATA FILE 

df = read.table(datafile,sep='\t',header=T)
colnames(df) = make.names(colnames(df))
colnames(df)[colnames(df) == 'protein'] <- 'ProteinID'
names(df)[5:10] <- c("MDA_GD3", "MDA_GD2", "MDA_GD1", "MDA_HG1", "MDA_HG3" ,"MDA_HG2")

head(df)
colnames(df)

#FILTER DATA TO SELECT PROTEIN ID AND INTENSITIES
df2 <- df[, c(-2, -3, -4)]
df2

#CHECK COLUMN NAMES - HAS TO MATCH TREATMENT NAMES!!!! 
cat("Use the following columns names for selection in next cell\n")
cat(paste('"',paste(sort(colnames(df2[c(-1)])),collapse='","'),'"',sep=""))

#GROUP DATA BY COLUMN NAMES
MDA_GD = c("MDA_GD1","MDA_GD2","MDA_GD3")
MDA_HG  = c("MDA_HG1","MDA_HG2","MDA_HG3")


#CREATE FACTORS FROM DATA 
fac.all <- as.factor(c(rep("MDA_GD",3), rep("MDA_HG",3)))

data_columns = c(MDA_GD, MDA_HG)
all_columns = c("ProteinID", data_columns)

head(df2)

# SPECIFY APPROPRIATE FACTOR LABELS FOR THE SAMPLE NAMES ABOVE - NAMES NEED TO MATCH 
cat("Factors for data table are:\n", paste(fac.all,collapse=","), "\n")
df2_all = df2[all_columns]

# #FILTER DATA SO THAT REQUIRES INTENSITY THRESHOLD >= 2^15 AND IN 3 SAMPLES) 
cat(paste("Number of rows before filtering:\n",nrow(df2_all),"\n"))
df2_all_filtered = df2_all[rowSums(df2_all[,c(-1)] > min_counts) >= min_num_observations , ]
cat(paste("Number of rows after filtering:\n",nrow(df2_all_filtered)))
head(df2_all_filtered)
tail(df2_all_filtered)

#PLOT UNNORMALIZED DATA AS INTENSITY VALUES 
pdf("Density_Plots_RAW_Intensity_Unimputed_MDA_WCP.pdf")
df_mis_raw_all <- subset(df2_all_filtered,select=-c(ProteinID))
rownames(df_mis_raw_all) <- df2_all_filtered[,1]
plot_raw1_all = melt(df_mis_raw_all)

ggplot(plot_raw1_all, aes(x = value, y = variable, height = ..density..)) + 
  geom_density_ridges2(stat = "density_ridges", fill = alpha("darkred", 0.7)) +
  xlab("Distribution of Intensity Values: MDA-MB-468") +
  ylab("Treatments") +
  theme_bw() +
  theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))

#PLOT UNNORMALIZED DATA AS LOG2 TRANSFORMED WITH MISSING VALUES 

plot_raw2_all = melt(log2(df_mis_raw_all))
plot_raw2_all[plot_raw2_all == -Inf] <- 0

ggplot(plot_raw2_all, aes(x = value, y = variable, height = ..density..)) + 
  geom_density_ridges2(stat = "density_ridges", fill = alpha("darkred", 0.7)) +
  xlab("Distribution of Log2 Transformed Intensity Values: MDA-MB-468") +
  ylab("Treatments") +
  theme_bw() +
  theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))

dev.off()

#MISSING VALUE COMBINED VIOLIN + BOX PLOT
pdf("Violin_Box_Prenormalization_MDA_WCP.pdf")
ggplot(plot_raw2_all, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = F) +
  geom_boxplot(width=0.25, fill="white") +
  ggtitle("Pre Hybrid Imputation:\n MDA-MB-468") +
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

#CONVERT MISSING VALUES (0s) TO NAs AND TRANSFORM DATA (log2) 

df_mis_raw_na_all = df_mis_raw_all
df_mis_raw_na_all[df_mis_raw_na_all == 0] <- NA

head(df_mis_raw_na_all)

df_mis_all = log2(df_mis_raw_na_all)
head(df_mis_all)

# VISUALIZE MISSING VALUES IN THE DATA 

pdf("Missing_Data_Examination_RAW_MDA_WCP.pdf", width = 12)

#NUMBER OF MISSING VALUES BY PROTEIN ID
gg_miss_case(df_mis_all)  + labs(x = "Protein ID", y = "Number of Missing Values") +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 12, family = "Helvetica", color = 'black'), 
        axis.title = element_text(face = "bold", size = 18, family = "Helvetica")) 

#CUMULATIVE SUM OF MISSING VALUES BY PROTEIN ID
gg_miss_case_cumsum(df_mis_all) + labs(x = "Protein ID", y = "Missing Cumulative Sum") +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 12, family = "Helvetica", color = 'black'), 
        axis.title = element_text(face = "bold", size = 18, family = "Helvetica"),
        axis.text.x= element_blank(), 
        axis.ticks.x = element_blank()) 

#NUMBER OF MISSING VALUES BY SAMPLE
gg_miss_var(df_mis_all) + labs(x = "Sample", y = "Number of Missing Values") +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 12, family = "Helvetica", color = 'black'), 
        axis.title = element_text(face = "bold", size = 18, family = "Helvetica"))

#CUMULATIVE SUM OF MISSING VALUES BY SAMPLE
gg_miss_var_cumsum(df_mis_all) + labs(x = "Sample", y = "Missing Cumulative Sum") +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 12, family = "Helvetica", color = 'black', angle = 45), 
        axis.title = element_text(face = "bold", size = 18, family = "Helvetica"))

#MISSING VALUE BY SAMPLE BY PROTEIN ID
vis_miss(df_mis_all) +
  coord_flip() +
  labs(y = "") +
  theme(axis.text = element_text(face = "bold", size = 24, family = "sans", color = "black"), 
        axis.title = element_text(face = "bold", size = 24, family = "sans"), 
        plot.title = element_blank(), 
        axis.text.x= element_blank(), 
        axis.ticks.x = element_blank()
  )


dev.off()



#CONVERT DATA TO MATRIX FOR PCA, SCREE AND BOX PLOTS
df_mis_all_temp <- df_mis_all 
df_mis_all_temp[is.na(df_mis_all_temp)] <- 0
all_data <- as.matrix(sapply(df_mis_all_temp, as.numeric))

#CALCULATE PRINCIPAL COMPONENTS FOR PCA PLOT FOR UNIMPUTED DATA
pdf("PCA_Scree_Plots_RAW_MDA_WCP.pdf")
PC <- prcomp(all_data, scale = TRUE)
scores <- PC$rotation[,1:5]                        # scores for first five PC's

# K-MEANS CLUSTERING [Assume 2 Clusters]
km     <- kmeans(scores, centers=2, nstart=10)

ggdata <- data.frame(scores, Cluster=km$cluster, Sample=fac.all)

ggplot(ggdata, aes(x= PC1,y = PC2, col = Sample)) +
  geom_point(size=4)+ 
  stat_ellipse(aes(x = PC1,y=PC2), fill = 'white', 
               geom="polygon",level=0.95, alpha=0.2) +
  scale_color_manual(values = c("purple4","darkred","blue4")) + 
  
  labs(title="LFQ Intensity Values:\nRAW Intensity MDA-MB-468",
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
  ylim(0, 70) +
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

#VISUALIZE BOXPLOTS OF UNIMPUTED DATA 
pdf("Boxplots_RAW_MDA_WCP.pdf")

data_box <- stack(as.data.frame(all_data))
data_box

box1 <- ggplot(data_box, aes(x = ind, y = values)) +
  geom_boxplot(fill = c(rep("darkred", 3), rep("purple4", 3)), alpha = 0.75) +
  geom_hline(aes(yintercept = mean(data_box$values)), lty = "dashed", color = "chartreuse3", size = 1.25) +
  labs(title = "RAW Intensity Values:\n MDA-MB-468", y = "LFQ Intensity",
       x = "Sample") +
  theme_bw() +
  theme(plot.title = element_text(family = "Helvetica",face = "bold", hjust = 0.5,
                                  size = 22),
        axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
        axis.text = element_text(family = "Helvetica", face = "bold", size = 14, 
                                 color = "black"),
        legend.title = element_text(family = "Helvetica",face = "bold", size = 
                                      16),
        legend.text = element_text(family = "Helvetica",face = "plain", size = 
                                     13),
        axis.text.x = element_text(angle = 90),
        axis.ticks.x=element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = 'right')

box1

dev.off()

#ADD NUMBER OF MISSING FOR EACH SAMPLE AND TOTAL 

df_combined <-  subset(df2_all_filtered,select=c(ProteinID))
df_combined$num_missing <- rowSums(df2_all_filtered == 0)
df_combined$a_missing <-rowSums(df2_all_filtered[2:4] == 0)
df_combined$b_missing <- rowSums(df2_all_filtered[5:7] == 0)
df_combined$both_missing <- df_combined$a_missing > 0 & df_combined$b_missing > 0


#REPLACE 0s WITH NAs and LOG2 TRANSFORM DATA
df2_all_filtered[df2_all_filtered == 0 ] <- NA

df_counts_sel_trans <- log2(df2_all_filtered[, c(-1)])
rownames(df_counts_sel_trans) <- df2_all_filtered[, 1]
df_counts_sel_trans


#CONVERT DATAFRAME TO MATRIX BEFORE IMPUTATION
sel_data_imp <- as.matrix(sapply(df_counts_sel_trans, as.numeric))

#CREATE LIST TO CONCATENATE ROWS 
datalist = list() 

for (i in seq(1, maxiterations))  {
  
  #HYBRID MODEL SELECTOR (MAR, MNAR) - CHOOSE FIRST COMPARISON # OF COLUMNS FOR BIOREPS
  sel_data_imp_a = sel_data_imp[,1:3]
  m.s.all.a = model.Selector(sel_data_imp_a)
  
  df.mis.all.ms.a = impute.MAR.MNAR(sel_data_imp_a, m.s.all.a, 
                                    method.MAR = mar_imputation,
                                    method.MNAR = mnar_imputation)
  
  #HYBRID MODEL SELECTOR (MAR, MNAR) - CHOOSE SECOND COMPARISON # OF COLUMNS FOR BIOREPS 
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
  write.csv(df_merge, paste0("Imp_Post_Hybrid_MDA_WCP_", i, ".csv"), row.names = F)
  write.csv(df.threshold, "Model_Selector_Hybrid_MDA_WCP.csv", row.names = F)
  
  #VISUALIZE DATA BEFORE AND AFTER IMPUTATIONS AS HISTOGRAM
  
  pdf("Intensity_Distribution_Post_Hybrid_MDA_WCP.pdf")
  
  hist(sel_data_imp[, c(1:6)],
       breaks = 25,
       main  = "Raw Data MDA-MB-468",
       xlab = "Log2 Peptide Intensities",
       col = "red4")
  box()
  
  hist(as.matrix(df_merge[, c(1:6)]), 
       breaks = 25,
       main = "Hybrid Imputation: MDA-MB-468",
       xlab = "Log2 Peptide Intensities",
       col = "blue4", add = F)
  box()
  
  hist(as.matrix(df_merge[, c(1:6)]),
       breaks = 25,
       main  = "Hybrid Imputation: MDA-MB-468",
       xlab = "Log2 Peptide Intensities",
       col = "blue4")
  
  hist(sel_data_imp[, c(1:6)], 
       breaks = 25,
       xlab = "Log2 Peptide Intensities",
       col = "red4",
       add = T)
  legend("topright", 
         c("Raw Intensity", "Hybrid Imputation"), 
         fill= c("red4","blue4"), 
         bty = "n")
  box()
  
  dev.off()
  
  #VISUALIZE DATA AS DENSITY PLOTS BEFORE AND AFTER IMPUTATION
  df_mis_sel <- df_counts_sel_trans
  df_mis_sel[is.na(df_mis_sel)] <- 0
  plot_raw1_sel = melt(df_mis_sel)
  
  #BEFORE IMPUTATION
  pdf("Density_Plot_Post_Hybrid_WCP_MDA.pdf")
  ggplot() + 
    geom_density_ridges2(data = plot_raw1_sel, aes(x = value, y = variable, height = ..density..), 
                         stat = "density_ridges", fill = alpha("darkred", 0.7)) +
    xlab("Distribution Pre-Imputation: MDA-MB-468") +
    ylab("Treatments") +
    theme_bw() +
    theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))
  
  #PLOT ATER IMPUTATION
  plot_raw2_all = melt(df_merge)
  
  ggplot() + 
    geom_density_ridges2(data = plot_raw2_all, aes(x = value, y = variable, height = ..density..), 
                         stat = "density_ridges", fill = alpha("darkblue", 0.7)) +
    xlab("Distribution Post Hybrid Imputation: MDA-MB-468") +
    ylab("Treatments") +
    theme_bw() +
    theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))
  
  #OVERLAY THE PLOTS TOGETHER
  
  ggplot() + 
    geom_density_ridges2(data = plot_raw1_sel, aes(x = value, y = variable, height = ..density..), 
                         stat = "density_ridges", fill = alpha("darkred", 0.7)) +
    geom_density_ridges2(data = plot_raw2_all, aes(x = value, y = variable, height = ..density..), 
                         stat = "density_ridges", fill = alpha("darkblue", 0.5)) +
    xlab("Distribution Overlay Hybrid Imputation: MDA-MB-468") +
    ylab("Treatments") +
    theme_bw() +
    theme(axis.title=element_text(size=14,  family="Helvetica", face = "bold"))
  
  dev.off()
  
  
  #PCA AND SCREE PLOTS FOR SELECTED DATASET POST-IMPUTATION
  pdf("PCA_Scree_Post_Hybrid_MDA_WCP.pdf")
  df_merge_temp <- as.matrix(df_merge)
  
  PC3 <- prcomp(df_merge_temp, scale = TRUE)
  scores3 <- PC3$rotation[,1:5]                        # scores for first five PC's
  
  # K-MEANS CLUSTERING [Assume 2 Clusters]
  km3 <- kmeans(scores3, centers=2, nstart=10)
  
  ggdata3 <- data.frame(scores3, Cluster=km3$cluster, Treatment=fac.all)
  
  ggplot(ggdata3, aes(x= PC1,y = PC2, col = Treatment)) +
    geom_point(size=4) + 
    stat_ellipse(aes(x = PC1,y=PC2), fill = 'white', 
                 geom="polygon",level=0.95, alpha=0.2) +
    scale_color_manual(values = c("purple4", "darkred")) + 
    labs(title="LFQ Intensity Values:\nPost Hybrid Imputation MDA-MB-468",
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
    labs(title = "Variance Analysis of Factors\nLFQ Intensities: MDA-MB-468",
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
  
  #BOXPLOT FOR DATA POST-IMPUTATION
  pdf("BoxPlots_Post_Hybrid_MDA_WCP.pdf")
  data_box3 <- stack(as.data.frame(df_merge_temp))
  
  ggplot(data_box3, aes(x = ind, y = values)) +
    geom_boxplot(fill = c(rep("darkred", 3), rep("purple4", 3)), alpha = 0.75) +
    geom_hline(aes(yintercept = mean(df_merge_temp)), lty = "dashed", color = "chartreuse3", size = 1.25) +
    labs(title = "LFQ Intensity Values:\nPost Hybrid Imputation MDA-MB-468", y = "LFQ Intensity",
         x = "Treatment") +
    theme_bw() +
    theme(plot.title = element_text(family = "Helvetica",face = "bold", hjust = 0.5,
                                    size = 22),
          axis.title = element_text(family = "Helvetica", face = "bold", size = 18),
          axis.text = element_text(family = "Helvetica", face = "bold", size = 14, 
                                   color = "black"),
          legend.title = element_text(family = "Helvetica",face = "bold", size = 
                                        16),
          legend.text = element_text(family = "Helvetica",face = "plain", size = 
                                       13),
          axis.text.x = element_text(angle = 90),
          axis.ticks.x=element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = 'right')
  
  dev.off()
  
  #VIOLIN & BOXPLOT POST_IMPUTATION AND NORM
  pdf("Violin_Box_Post_Hybrid_MDA_WCP.pdf")
  ggplot(data_box3, aes(x=ind, y=values, fill=ind)) +
    geom_violin(trim = F) +
    geom_boxplot(width=0.25, fill="white") +
    ggtitle("Post Hybrid Imputation:\n MDA-MB-468") +
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
  
  #PERFORM LIMMA TEST FOR DE
  pdf("Pvalues_Post_Hybrid_MDA_WCP.pdf")
  
  cat(paste("Performing", limma_norm, "normalization on data"))
  # perform quantile normalization of data points and then limma analysis
  if (limma_norm == "Quantile") {
    exprSet_quantile = normalizeQuantiles(as.matrix((df_merge_temp)))
  } else {
    exprSet_quantile = as.matrix((df_merge_temp))
  }
  
  fit <- lmFit(as.matrix(exprSet_quantile), design = model.matrix(~ fac.all))
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=2)
  limmares <- data.frame(dm=coef(fit)[,2], p.value = fit$p.value[, 2])
  results <- decideTests(fit)
  pvalues <- fit$p.value
  qobj <- qvalue(p = pvalues[,2])
  qvalues <- qobj$qvalues
  pi0 <- qobj$pi0
  lfdr <- qobj$lfdr
  
  hist(qobj)
  dev.off()
  
  pdf("QValue_Post_Hybrid_MDA_WCP.pdf")
  plot(qobj)
  dev.off()
  
  
  #KNIT COLUMNS BACK TOGETHER
  temp_cnames_imp <- colnames(df_merge_temp) 
  temp_cnames_norm <- colnames(exprSet_quantile) 
  colnames(df_merge_temp) <- paste(colnames(df_merge_temp), "imp", sep = "_")
  colnames(exprSet_quantile) <- paste(colnames(exprSet_quantile), "quant_norm", sep = "_")
  output_data = cbind(df2_all_filtered, df_merge_temp, exprSet_quantile, data.frame(topTable(fit, coef=2, sort="none", n=Inf)))
  output_data['qvalue'] <- qvalues
  
  temp_ids = stringr::str_match(output_data$ProteinID,'^([a-z]*\\|([a-zA-Z0-9]*)\\|)([a-zA-Z0-9]*)_[a-zA-Z0-9]*')
  colnames(temp_ids) <- c("Protein_Description", "UniprotID", "UniprotAC", "Protein_Name")
  output_data = cbind(temp_ids, output_data[,c(-1)], protein_IDs=output_data[, c(1)] )
  colnames(df_merge) <-  temp_cnames_imp 
  colnames(exprSet_quantile) <-  temp_cnames_norm 
  
  # OUTPUT RESULTS TO CSV
  write.csv(output_data[order(output_data$qvalue),], paste0("DE_Ordered_Hybrid_MDA_WCP_", i,".csv"), 
            row.names = F)
  
  
  #ADD NUMBER OF MISSING VALUES FOR EACH PROTEIN BY SAMPLE TO DATA
  
  temp_data <- output_data[, c("Protein_Description", "logFC", "P.Value", "qvalue", "AveExpr")]
  datalist[[i]] <- temp_data
  names(temp_data) = paste0(names(temp_data) ,"_",i)
  df_combined <- merge(df_combined,temp_data, by.x = c("ProteinID"),
                       by.y=c(paste0("Protein_Description" ,"_",i)))
  
}
#RBIND DATA FOR LONG FORMAT
big_data = do.call(rbind, datalist)
write.csv(big_data, "Long_data_qvalue_logFC_Hybrid_MDA_WCP.csv", row.names = F)

  #CALCULATE LOGFC MEAN AND STDDEV, QVALUE MEAN AND STDDEV
  
df_combined$logFC_mean = rowMeans(df_combined[,grepl("logFC",names(df_combined))])
df_combined$logFC_stdev = apply(df_combined[,grepl("logFC",names(df_combined))],1,sd)
df_combined$P.Value_mean = rowMeans(df_combined[,grepl("P.Value",names(df_combined))])
df_combined$P.Value_stdev = apply(df_combined[,grepl("P.Value",names(df_combined))],1,sd)
df_combined$qvalue_mean = rowMeans(df_combined[,grepl("qvalue",names(df_combined))])
df_combined$qvalue_stdev = apply(df_combined[,grepl("qvalue",names(df_combined))],1,sd)
df_combined$AveExpr_mean = rowMeans(df_combined[,grepl("AveExpr",names(df_combined))])
df_combined$AveExpr_stdev = apply(df_combined[,grepl("AveExpr",names(df_combined))],1,sd)
df_combined$Factor = ifelse(df_combined$AveExpr_mean > m.s.all.a[[2]], "Green", "Red")

write.csv(df_combined, "Combined_logFC_qvalue_pvalue_Hybrid_MDA_WCP.csv", row.names = F)

pdf("Boxplots_Variance_Hybrid_MDA_WCP.pdf")

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

annotate_figure(figure, top = text_grob("Hybrid Imputation MDA-MB-468", face = "bold", size = 18), 
                fig.lab = "Fig. S2", fig.lab.face = "bold")

dev.off()
