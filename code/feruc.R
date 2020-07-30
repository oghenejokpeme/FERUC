library(data.table)
library(ranger)
library(smotefamily)
library(caret)
library(arules)

model_performance <- function(y_true, pred){
  rsquared <- 1 - (sum((y_true - pred)^2) / sum((y_true - mean(y_true))^2))
  mse <- sum((y_true - pred)^2) / length(y_true)
  rmse <- sqrt(sum((y_true - pred)^2) / length(y_true))
  c(rsquared = rsquared, mse = mse, rmse = rmse)
}

random_forests <- function(x, y, ntrees=1000, prob=F){
  set.seed(45322)
  ranger(x=x, y=y, importance="impurity", num.trees = ntrees, 
         probability=prob, verbose=T)
}

read_data <- function(){
  bp <- "../input/datasets/"
  x_train <- data.frame(fread(paste0(bp, "base_fp_train.csv"), header=T), 
                        row.names=1)
  x_test <- data.frame(fread(paste0(bp, "base_fp_test.csv"), header=T), 
                       row.names=1)
  ys_train <- data.frame(fread(paste0(bp, "expression_train.csv"), header=T), 
                         row.names=1)
  ys_test <- data.frame(fread(paste0(bp, "expression_test.csv"), header=T), 
                        row.names=1)
  list(x_train=x_train, x_test=x_test, 
       ys_train=round(ys_train, 3), ys_test=round(ys_test, 3))
}

get_sorted_bins <- function(y, n){
  y <- sort(y)
  split(y, sort(seq_along(y)%%n))
}

get_disc_bins <- function(y, n, type){
  idx <- discretize(y, method=type, breaks=n)
  ybins <- split(y, idx)
  names(ybins) <- 0:(n-1)
  ybins
}

get_random_bins <- function(y, n){
  amnts <- seq_along(y)%%n
  set.seed(237)
  idx <- sample(1:length(y))
  split(y, amnts[idx])
}

print_data_shape <- function(title, df){
  dfs <- dim(df)
  dfs <- paste0(dfs[1], " ", dfs[2])
  print(paste0("       ", title, " -> ", dfs))
}

# Build classifier without considering class imbalance.
build_rg_classifier <- function(x_pos, x_neg){
  x_class <- rbind(x_pos, x_neg)
  y_pos <- rep(1, nrow(x_pos))
  y_neg <- rep(0, nrow(x_neg))
  y_class <- as.factor(c(y_pos, y_neg))
  print_data_shape("rg class", x_class)
  random_forests(x_class, y_class, prob=T)
}

# Build classifier with random down-sampling of majority class.
build_dn_classifier <- function(x_pos, x_neg){
  set.seed(478)
  if (nrow(x_pos) > nrow(x_neg)){
    x_pos <- x_pos[sample(1:nrow(x_pos), nrow(x_neg)), ]
  } else {
    x_neg <- x_neg[sample(1:nrow(x_neg), nrow(x_pos)), ]
  }
  x_class <- rbind(x_pos, x_neg)
  y_pos <- rep(1, nrow(x_pos))
  y_neg <- rep(0, nrow(x_neg))
  y_class <- as.factor(c(y_pos, y_neg))
  print_data_shape("dn class", x_class)
  random_forests(x_class, y_class, prob=T)
}

# Build classifier with upsampling of minority class with SMOTE.
build_up_classifier <- function(x_pos, x_neg, k=5){
  y_pos <- rep(1, nrow(x_pos))
  y_neg <- rep(0, nrow(x_neg))
  y_class <- as.factor(c(y_pos, y_neg))
  x_class <- rbind(x_pos, x_neg)
  if (dim(x_pos)[1] < k){
    k <- dim(x_pos)[1] - 1
  }
  set.seed(375)
  bdata <- SMOTE(x_class, y_class, K=k)
  cidx <- ncol(bdata$data)
  nx_class <- bdata$data[, -cidx]
  ny_class <- as.factor(bdata$data[, cidx])
  print_data_shape("up class", nx_class)
  random_forests(nx_class, ny_class, prob=T)
}

build_models <- function(x, y_bins){
  models <- vector(mode="list", length=length(y_bins))
  for (i in 1:length(y_bins)){
    print(paste0("     model: ", i))
    # Process current bin samples.
    y_bin <- y_bins[[i]]
    bin_samples <- names(y_bin)
    rem_samples <- setdiff(rownames(x), bin_samples)

    x_pos <- x[bin_samples, ]
    x_neg <- x[rem_samples, ]
    
    # Build regressor.
    regressor <- random_forests(x_pos, y_bin)

    # Build classifiers.
    rg_classifier <- build_rg_classifier(x_pos, x_neg)
    dn_classifier <- build_dn_classifier(x_pos, x_neg)
    up_classifier <- build_up_classifier(x_pos, x_neg)

    # Add models to list
    models[[i]]$regressor <- regressor
    models[[i]]$rg_classifier <- rg_classifier
    models[[i]]$dn_classifier <- dn_classifier
    models[[i]]$up_classifier <- up_classifier
  }
  models
}

categorise_test_samples <- function(y_chunks, y_test){
  classes <- vector(mode="list", length=length(y_chunks))
  for (i in 1:length(y_chunks)){
    imin <- min(y_chunks[[i]])
    imax <- max(y_chunks[[i]])
    csamps <- names(which((y_test >= imin) & (y_test <= imax)))
    classes[[i]] <- csamps
  }
  classes
}

get_predictions <- function(models, x_test, nbins){
  rgs <- NULL # Ignore class imbalance predictions.
  dns <- NULL # Down-sampling classifier predictions.
  ups <- NULL # Up-sampling classifier predictions.
  pvals <- NULL

  all_bins <- 1:nbins
  for (i in all_bins){
    rgcls <- predict(models[[i]]$rg_classifier, x_test)$predictions[, "1"]
    dncls <- predict(models[[i]]$dn_classifier, x_test)$predictions[, "1"]
    upcls <- predict(models[[i]]$up_classifier, x_test)$predictions[, "1"]
    pval <- predict(models[[i]]$regressor, x_test)$predictions
    
    rgs <- cbind(rgs, rgcls)
    dns <- cbind(dns, dncls)
    ups <- cbind(ups, upcls)
    pvals <- cbind(pvals, pval)
  }

  rgs <- round(rgs, 3)
  dns <- round(dns, 3)
  ups <- round(ups, 3)
  colnames(rgs) <- all_bins
  colnames(dns) <- all_bins
  colnames(ups) <- all_bins
  colnames(pvals) <- all_bins

  list(rg=rgs, dn=dns, up=ups, pvals=pvals)
}

make_predictions <- function(raw_preds){
  # Simple averaging.
  savg <- rowMeans(raw_preds$pvals)
  # Class imbalance is ignored.
  rgweights <- raw_preds$rg/rowSums(raw_preds$rg)
  rgpreds <- rowSums(raw_preds$pvals*rgweights)  
  # Down-sampling classifier weighting.
  dnweights <- raw_preds$dn/rowSums(raw_preds$dn)
  dnpreds <- rowSums(raw_preds$pvals*dnweights)
  # Up-sampling classifier weighting.
  upweights <- raw_preds$up/rowSums(raw_preds$up)
  uppreds <- rowSums(raw_preds$pvals*upweights)
  
  list(ag=savg, rg=rgpreds, dn=dnpreds, up=uppreds)
}

write_reg_performance <- function(gene, nbins, type, y_test, predictions){
  for (pred in names(predictions)){
    fp <- paste0("../output/performance/reg_", nbins, "_", pred, "_", 
                 type, ".txt")
    nperf <- round(model_performance(y_test, predictions[[pred]]), 3)
    wperf <- c(gene, as.vector(nperf))
    write(wperf, append=T, ncolumns=length(wperf), file=fp)
  }
}

convert_probs_to_class <- function(df, test_samples){
  ndf <- NULL
  for (i in 1:ncol(df)){
    nwcol <- sapply(df[,i], function(x){if(x>0.5){1}else{0}})
    ndf <- cbind(ndf, nwcol)
  }
  colnames(ndf) <- colnames(df)
  rownames(ndf) <- test_samples
  ndf
}

classifier_perf <- function(clsdf, test_cats){
  cperfs <- NULL
  test_samples <- rownames(clsdf)
  for (i in 1:ncol(clsdf)){
    pos_samples <- test_cats[[i]]
    neg_samples <- setdiff(test_samples, pos_samples)
    actual <- c(rep(1, length(pos_samples)), rep(0, length(neg_samples)))
    pred <- clsdf[c(pos_samples, neg_samples), i]
    actual <- factor(actual, levels=c("0", "1"))
    pred   <- factor(pred, levels=c("0", "1"))
    cf <- confusionMatrix(actual, pred, positive="1")
    cperf <- round(c(cf$overall, cf$byClass), 3)
    cperfs <- rbind(cperfs, cperf)
  }
  rownames(cperfs) <- colnames(clsdf)
  cperfs
}

get_classifier_perf <- function(raw_preds, test_samples, test_cats){
  # rgcls
  rg <- convert_probs_to_class(raw_preds$rg, test_samples)
  rg_perf <- classifier_perf(rg, test_cats)
  # dncls
  dn <- convert_probs_to_class(raw_preds$dn, test_samples)
  dn_perf <- classifier_perf(dn, test_cats) 
  # upcls
  up <- convert_probs_to_class(raw_preds$up, test_samples)
  up_perf <- classifier_perf(up, test_cats)

  list(rg=rg_perf, dn=dn_perf, up=up_perf)
}

write_cls_performance <- function(gene, nbins, type, cperfs){
  for (classifier in names(cperfs)){
    cperf <- cperfs[[classifier]]
    afp <- paste0("../output/performance/all_class/", gene, 
                  "_", nbins, "_", classifier, "_", type, ".csv")
    write.csv(cperf, file=afp)

    sfp <- paste0("../output/performance/cls_", nbins, "_", classifier, 
                  "_", type, ".txt")
    if(!file.exists(sfp)){
      wperf <- c("Gene", colnames(cperf))
      write(wperf, append=T, ncolumns=length(wperf), sep="\t", file=sfp)  
    }
    nperf <- as.vector(colMeans(cperf))
    wperf <- c(gene, round(nperf, 3))
    write(wperf, append=T, ncolumns=length(wperf), sep="\t", file=sfp)
  }
}

log_models <- function(models, gene, nbins, type){
  for (i in 1:nbins){
    bp <- paste0("../output/models/", type, "_", nbins, "_", gene, "_", i, "_")
    saveRDS(models[[i]]$regressor, file=paste0(bp, "regressor.RDS"))
    saveRDS(models[[i]]$rg_classifier, file=paste0(bp, "rg.RDS"))
    saveRDS(models[[i]]$dn_classifier, file=paste0(bp, "dn.RDS"))
    saveRDS(models[[i]]$up_classifier, file=paste0(bp, "up.RDS"))
  }
}

feruc <- function(x_train, x_test, y_train_bins, y_test, nbins, gene, type){
  models <- build_models(x_train, y_train_bins)
      
  # Make regression predictions.
  raw_predictions <- get_predictions(models, x_test, nbins)
  predictions <- make_predictions(raw_predictions)
  write_reg_performance(gene, nbins, type, y_test, predictions)
      
  # Test bin classifier accuracy.
  test_cats <- categorise_test_samples(y_train_bins, y_test)
  cperfs <- get_classifier_perf(raw_predictions, names(y_test), test_cats)
  write_cls_performance(gene, nbins, type, cperfs)

  log_models(models, gene, nbins, type)
}

main <- function(){
  datasets <- read_data()
  x_train <- datasets$x_train
  x_test  <- datasets$x_test
  ys_train <- datasets$ys_train
  ys_test  <- datasets$ys_test

  genes <- read.table(file="../input/genes.txt", header=F)[,1]
  for (gene in genes){
    print(gene)
    y_train <- ys_train[, gene]
    y_test  <- ys_test[, gene]
    names(y_test)  <- rownames(ys_test)
    names(y_train) <- rownames(ys_train)

    # Normal case
    print(paste0("  Performing base case.")) 
    rfnorm <- random_forests(x_train, y_train)
    spreds <- predict(rfnorm, x_test)$predictions
    sperf  <- round(model_performance(y_test, spreds), 3)
    wperf <- c(gene, as.vector(sperf))
    saveRDS(rfnorm, file=paste0("../output/models/base_", gene, "_regressor.RDS"))
    write(wperf, append=T, ncolumns=length(wperf),
          file="../output/performance/base.txt") 
    
    # Bin case
    print(paste0("  Performing FERUC.")) 
    for (nbins in 2:5){
      print(paste0("  Bin size: ", nbins))
      print("    Sorted.")    
      y_train_bins <- get_sorted_bins(y_train, nbins)
      feruc(x_train, x_test, y_train_bins, y_test, nbins, gene, "srt")
     
      print("    Frequency.")    
      y_train_bins <- get_disc_bins(y_train, nbins, "frequency")
      feruc(x_train, x_test, y_train_bins, y_test, nbins, gene, "frq")
      
      print("    K-means.")    
      y_train_bins <- get_disc_bins(y_train, nbins, "cluster")
      feruc(x_train, x_test, y_train_bins, y_test, nbins, gene, "kmn")
      
      print("    Random.")    
      y_train_bins <- get_random_bins(y_train, nbins)
      feruc(x_train, x_test, y_train_bins, y_test, nbins, gene, "rdm")
    }
  }
}

main()