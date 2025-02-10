# Ramdom forest models for alpha diversity and beta diversity 
commcomp.rf.res <-mclapply(X = colnames(NMDSscore), mc.cores = 30, FUN = function(x){
  axis.used.in.rf <- NMDSscore[, x]
  scaled_datamain <- as.data.frame(scale(datamain, center = TRUE, scale = TRUE))
  rfs.dat <- data.frame(value = axis.used.in.rf, scaled_datamain)
  trf <- randomForest::tuneRF(rfs.dat[, 2:ncol(rfs.dat)], rfs.dat[, "value"], ntreeTry  = 999, plot = FALSE)
  mt <- trf[which.min(trf[,2]), 1]
  rf.res <- randomForest::randomForest(value ~ ., data = rfs.dat, mtry = mt, importance = TRUE, na.action = na.omit, ntree  = 999)
  rf.Pres <- rfPermute::rfPermute(value ~ ., data=rfs.dat,mtry = mt,importance=T,
                       ntree  = 999)
  
  importance.res <- as.data.frame(importance( rf.Pres, scale = T))
  importance.res.sorted <- importance.res[order(importance.res$`%IncMSE`, decreasing = TRUE),]
  xx <- capture.output(rf.res)[length(capture.output(rf.res))]
  exp.perc <- as.numeric(strsplit(xx, split = ': ', fixed = T)[[1]][2])
  res <- data.frame(axis = x, exp.perc = exp.perc)
  result_list <- list(importance_res_sorted = importance.res.sorted, res = res)
  return(result_list)
})

# Ramdom forest models for the 19,194 domiant phoD-harboring phylotypes

perform_rf_analysis <- function(data, num_cores = 36) {
  rf_results <- mclapply(X = colnames(data), mc.cores = num_cores, FUN = function(x) {
    axis_used_in_rf <- data[, x]
    scaled_data <- as.data.frame(scale(datamain_otu, center = TRUE, scale = TRUE))
    rfs_data <- data.frame(value = axis_used_in_rf, scaled_data)
    trf <- randomForest::tuneRF(rfs_data[, 2:ncol(rfs_data)], rfs_data[, "value"], ntreeTry  = 999, plot = FALSE)
    mt <- trf[which.min(trf[,2]), 1]
    rf_res <- randomForest::randomForest(value ~ ., data = rfs_data, mtry = mt, importance = TRUE, na.action = na.omit, ntree  = 999)
    rf_Pres <- rfPermute::rfPermute(value ~ ., data=rfs_data,mtry = mt,importance=T, ntree  = 999)
    importance_res <- as.data.frame(importance(rf_Pres, scale = T))
    importance_res_sorted <- importance_res[order(importance_res$`%IncMSE`, decreasing = TRUE),]
    xx <- capture.output(rf_res)[length(capture.output(rf_res))]
    exp_perc <- as.numeric(strsplit(xx, split = ': ', fixed = T)[[1]][2])
    res <- data.frame(axis = x, exp_perc = exp_perc)
    result_list <- list(importance_res_sorted = importance_res_sorted, res = res)
    return(result_list)
  })

  return(rf_results)
}


asv.use <- read.csv("Dominat_otu_table.csv")
# Perform analysis
rf.res <- perform_rf_analysis(asv.use)


### Semi-partial correlation 

#### Semi-partial correlation
asv.use.sel <- asv.use[which(rownames(asv.use) %in% unique.asv),]
rownames(datamain_otu) <- env_replaced3$Sample_name
var.use.sel <- datamain_otu[,which(colnames(datamain_otu) %in% unique.var)]


var.use.sel <- datamain_otu
if (length(unique(colnames(asv.use.sel) == rownames(var.use.sel))) > 1 |
    !unique(colnames(asv.use.sel) == rownames(var.use.sel))){
  stop("Sample order error!")
}
library(ppcor)
asv.use.sel.t <- t(asv.use.sel)
asv.use.sel.id <- colnames(asv.use.sel.t)
pc.res <- mclapply(X = 1:length(asv.use.sel.id), mc.cores = 100, FUN = function(x){
  asv.id <- asv.use.sel.id[x]
  as.semi.pcor <- asv.use.sel.t[,which(colnames(asv.use.sel.t) == asv.id)]
  asv.dat <- data.frame(row.names = rownames(asv.use.sel.t), asv = as.semi.pcor)
  colnames(asv.dat) <- asv.id
  dat.main <- as.data.frame(cbind(asv.dat, var.use.sel))
  dat.main <- scale(dat.main,center = TRUE, scale = TRUE)
  cor.rst <- ppcor::spcor(dat.main, method = "spearman")
  cor.res <- as.data.frame(cor.rst$estimate[which(rownames(cor.rst$estimate) == asv.id), which(colnames(cor.rst$estimate) != asv.id)])
  pva.res <- as.data.frame(cor.rst$p.value[which(rownames(cor.rst$p.value) == asv.id), which(colnames(cor.rst$p.value) != asv.id)])
  colnames(cor.res) <- colnames(pva.res) <- asv.id
  res <- list(r = cor.res, p = pva.res)
})
res.r <- res.p <- NULL
for (i in seq(length(pc.res))){
  tmp <- pc.res[[i]]
  if (i == 1){
    res.r <- tmp$r
    res.p <- tmp$p
  }else{
    res.r <- cbind(res.r, tmp$r)
    res.p <- cbind(res.p, tmp$p)
  }
}
res.r <- as.data.frame(res.r)
res.p <- as.data.frame(res.p)
cor.m <- as.matrix(t(res.r))
pva.m <- as.matrix(t(res.p))
if (length(unique(colnames(cor.m) == colnames(pva.m))) > 1 |
    !unique(colnames(cor.m) == colnames(pva.m))){
  stop("Order error!")
}
if (length(unique(rownames(cor.m) == rownames(pva.m))) > 1 |
    !unique(rownames(cor.m) == rownames(pva.m))){
  stop("Order error!")
}

# Pick out significant results of semi-partial correlations
cor.m.sel <- cor.m
for (i in seq(nrow(cor.m.sel))){
  for (j in seq(ncol(cor.m.sel))){
    if (pva.m[i, j] > 0.05 | abs(cor.m[i, j]) < 0.1) { # 保证相关性系数大于等于0.1
      cor.m.sel[i, j] <- 0
    }
  }
}
cor.m.sel.bin <- cor.m.sel; cor.m.sel.bin[cor.m.sel.bin != 0] <- 1
cor.m.sel.bin.sum <- data.frame(asv = rownames(as.data.frame(rowSums(cor.m.sel.bin))),  as.data.frame(rowSums(cor.m.sel.bin)))
valid.asv <- cor.m.sel.bin.sum[which(cor.m.sel.bin.sum[,2] > 0),]$asv
cor.m.sel <- cor.m.sel[which(rownames(cor.m.sel) %in% valid.asv),]
#pheatmap::pheatmap(cor.m.sel, show_rownames = F)

# Identify habitat cluster
# We only consider those clusters with members more than or equal to 30

habitat.clusters <- lapply(X = 1:nrow(cor.m.sel), function(x){
  asv.id <- rownames(cor.m.sel)[x]
  tmp <- data.frame(var = colnames(cor.m.sel), cor.val = cor.m.sel[x,])
  tmp$cor.val <- abs(tmp$cor.val)
  tmp <- tmp[order(tmp$cor.val, decreasing = T),]
  corr.var <- tmp$var[1]
  cor.val <- cor.m.sel[x,which(colnames(cor.m.sel) == corr.var)]
  cor.typ <- ifelse(cor.val > 0, "High", "Low")
  var.count <- nrow(tmp[which(tmp$cor.val > 0),]) # How many variables correlated with this ASV
  res.dat <- data.frame(asv.id = asv.id, corr.var = corr.var, cor.val = cor.val,
                        cor.typ = cor.typ, group = paste(corr.var, cor.typ, sep = '_'), var.count = var.count)
}) %>% do.call("rbind", .)
#habitat.clusters <- habitat.clusters[which(abs(habitat.clusters$cor.val) < 0.1),]
#use.var.count <- length(unique(habitat.clusters$corr.var))
#habitat.clusters <- habitat.clusters[which(habitat.clusters$var.count < round(use.var.count * 0.5, 0)),]
cluster.m.count <- as.data.frame(table(habitat.clusters$group))
valid.clusters <- cluster.m.count[which(cluster.m.count$Freq >= 30),]$Var1
habitat.clusters_valid <- habitat.clusters[which(habitat.clusters$group %in% valid.clusters),]

