
library(tidyverse)

otu = fread("feature-table1000_resampled.tsv",sep="\t",skip = 1)
env <- read.csv("global_phoDenv.csv")
env_rarefied <- filter(env,Sample_name %in% names(otu))

totu <- t(otu)
method = "bray"
threads = 50
bray_rarefied <- parallelDist::parallelDist(as.matrix(totu), 
                                    method = method, 
                                    diag = FALSE, upper = FALSE, 
                                    threads = threads)
str(bray_rarefied)
filtered_bray_rarefied <- as.matrix(bray_rarefied)

saveRDS(filtered_bray_rarefied,"bray_total_rarefied.RDS")

# running NMDS model
mod_NMDS  <- metaMDS(bray_rarefied)
points <- mod_NMDS$points

########### Distance-based redundancy analysis############
sample_ID <- env$Sample_name
sample_id_df <- as.data.frame(sample_ID)
names(sample_id_df) <- "Sample_name" 
# 使用这个布尔向量来选择对应的列
matched_columns <- names(otu) %in% sample_ID
env2_rarefied <- inner_join(sample_id_df,env_rarefied)
rdaenvall1 <- env2_rarefied[,c(24,15:16,18:23)]
rdaenvall1_subset <- as.data.frame(scale(rdaenvall1))
otu3 <- as.data.frame(as.data.frame(otu)[, matched_columns])


bray_rarefied2 <- parallelDist::parallelDist(as.matrix(t(otu3)), 
                                    method = method, 
                                    diag = FALSE, upper = FALSE, 
                                    threads = threads)

bray_rarefied2_dist <- as.matrix(bray_rarefied2)

# 使用 capscale 并传递加速后的距离矩阵
dbrda.res <- capscale(bray_rarefied2_dist ~ ., data = rdaenvall1_subset)

#dbrda.res <- capscale(t(otu3) ~ ., data=rdaenvall1_subset, dist = 'bray')

fit <- envfit(dbrda.res, rdaenvall1_subset, perm = 999, display = "lc")
# Relative contribution of each variable

RDA_variable <- rdacca.hp::rdacca.hp(as.dist(bray_rarefied2_dist),rdaenvall1_subset,method="dbRDA",type="adjR2",var.part = TRUE)

RDA_importance <- as.data.frame(RDA_variable$Hier.part)
names(RDA_importance)[4] <- "PercentImportance"
RDA_importance$variable <- rownames(RDA_importance)
global_phoD$dbRDA$plot_dbRDA_importance <- RDA_importance