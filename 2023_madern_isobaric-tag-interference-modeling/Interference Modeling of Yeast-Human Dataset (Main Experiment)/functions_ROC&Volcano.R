## This R-file contains all the required functions for the evaluation of interference correcton result ("Evaluation_InterferenceCorrectedProteinIntensities.Rmd")


### Write function that calculates ROC statistics via limma ###
ROC_function <- function (matrix_group1,
                          matrix_group2,
                          true_hypothesis = NULL,
                          ROC = TRUE, 
                          limma = TRUE,
                          trend_limma=TRUE){
  # note: true_hypothesis should be a logical vector; with TRUE when Nullhypothesis is true, else FALSE
  # note: the adjusted p-values are calculated via Benjamini Hochberg, i.e. they control the FDR 
  # The ROC-curve uses the adjusted p-values
  
  ## read limma package
  library(limma)
  library(ROCR)
  
  ## check if dimensions are equal
  if( nrow(matrix_group1) != nrow(matrix_group2)){
    stop("dimension of two matrixes has to be equal")
  }
  
  # initiate/create output vectors
  fc <- apply(matrix_group2, MARGIN=1, FUN=mean, na.rm=TRUE) - apply(matrix_group1, MARGIN=1, FUN=mean, na.rm=TRUE)
  p_val <- numeric(nrow(matrix_group1))
  data_plot <- data.frame(fc=fc, true_hypothesis=true_hypothesis)
  
  # do limma linear model DE testing and store p-values
  if(limma){
    design <- model.matrix(~factor(c(rep("group1", times=ncol(matrix_group1)), rep("group2", times=ncol(matrix_group2))), levels=c("group1", "group2")))
    colnames(design) <-  c("group1","group2")
    fit <- lmFit(cbind(matrix_group1, matrix_group2), design)
    fit_ebayes <- eBayes(fit, trend=trend_limma)  
    plotSA(fit_ebayes)
    abline(v=15, col="red")
    abline(v=20, col="violet")
    abline(h=0.2, col="green")
    abline(h=0.4, col="orange")
    LIMMAresults <- topTable(fit_ebayes, number=Inf, coef= "group2", adjust="BH", sort.by="none")
    p_val <- LIMMAresults$P.Value
    adj_p_val <- LIMMAresults$adj.P.Val
    
  # else (if limma = FALSE) do simple t-test  
  } else {
    for (i in 1:nrow(matrix_group1)){
      x <- matrix_group1[i,]
      y <- matrix_group2[i,]
      res_t.test <- t.test(x=x,y=y)
      p_val[i] <- res_t.test$p.value
    }
    adj_p_val <- p.adjust(p_val, method = "BH")
  }
  
  # calculate ROC curves
  if (ROC){
    pval_roc <- p_val
    pred <- prediction(pval_roc, as.numeric(true_hypothesis))  
    perf <- performance(pred,"tpr","fpr")
    par(xpd=FALSE)
    cutoff <- sort(unique(c(10^(-(30:3)), 1.2*10^(-(30:3)), 1.4*10^(-(30:3)), 1.6*10^(-(30:3)), 1.8*10^(-(30:3)), 2.010^(-(30:3)), seq(from=0, to=1, by=0.000001))))
    roc_plot_values <- sapply(cutoff, FUN=function(i){
      tpr_roc_i <- sum(pval_roc[!true_hypothesis] <= i)/(sum(!true_hypothesis))
      fpr_roc_i <- sum(pval_roc[true_hypothesis] <= i)/(sum(true_hypothesis))
      res_i <- cbind(tpr_roc_i, fpr_roc_i)
      return(res_i)
    }) 
    
    # get AUC:
    auc <- performance(pred, measure = "auc")
    auc_value <- auc@y.values[[1]]
    
    # print fp and tp numbers at p-value cutoff 0.05
    positives <- sum(!true_hypothesis)
    tp_pval <- sum(p_val[!true_hypothesis] <= 0.05)
    fn_pval <- sum(p_val[!true_hypothesis] > 0.05)
    fp_pval <- sum(p_val[true_hypothesis] <= 0.05)
    tn_pval <- sum(p_val[true_hypothesis] > 0.05)
    conf_matrix_pval_0.05 <- matrix(c(tp_pval, fn_pval, fp_pval, tn_pval), byrow = FALSE, ncol=2)
    rownames(conf_matrix_pval_0.05) <- c("predicted positive", "predicted negative")
    colnames(conf_matrix_pval_0.05) <- c("positive", "negative")
    # print(conf_matrix_pval_0.05)
    
    # print fp and tp numbers at adj. p-value cutoff 0.05
    tp_adj_pval <- sum(adj_p_val[!true_hypothesis] <= 0.05)
    fn_adj_pval <- sum(adj_p_val[!true_hypothesis] > 0.05)
    fp_adj_pval <- sum(adj_p_val[true_hypothesis] <= 0.05)
    tn_adj_pval <- sum(adj_p_val[true_hypothesis] > 0.05)
    conf_matrix_adj_pval_0.05 <- matrix(c(tp_adj_pval, fn_adj_pval, fp_adj_pval, tn_adj_pval), byrow = FALSE, ncol=2)
    rownames(conf_matrix_adj_pval_0.05) <- c("predicted positive", "predicted negative")
    colnames(conf_matrix_adj_pval_0.05) <- c("positive", "negative")
    # print(conf_matrix_adj_pval_0.05)
  }
  
  return(list(FC=fc,
              true_hypothesis=true_hypothesis,
              p_val=p_val,
              adj_p_val=adj_p_val, 
              ROC_object=perf, 
              auc_value=auc_value,
              conf_matrix_pval=conf_matrix_pval_0.05,
              conf_matrix_adj_pval=conf_matrix_adj_pval_0.05,
              ROC_coordinates = t(roc_plot_values),
              ROC_coordinates_pval_0.05= c(fp_pval/(fp_pval+tn_pval),tp_pval/(tp_pval + fn_pval))))
}



### Write function that produces volcano plot from ROC_function output ###
volcano_plot_from_roc <- function(list_roc,
                                  type_colors=type_colors,
                                  theo_yeastFC=1.333,
                                  x_limit) {
  ## create output list
  list_plots <- vector(mode="list", length=length(list_roc))
  
  ## go over each roc analysis to create ggplot. save in output list
  for (i in 1:length(list_roc)){
    # exctract parameters
    roc_i <- list_roc[[i]]
    fc_i <- roc_i$FC
    pval_i <- roc_i$p_val
    DE_i <- !roc_i$true_hypothesis
    df_plot <- data.frame(fc = fc_i,
                          pval = pval_i,
                          DE = DE_i)
    gg <- ggplot(data=df_plot) +
      geom_vline(xintercept=log2(theo_yeastFC),col="black", linetype=2, size=0.4) +
      geom_point(aes(x=fc, y=-log10(pval), col=DE), alpha=0.2, cex=1) +
      scale_color_manual(values=c("grey", type_colors[i])) +
      theme_bw() +
      xlab("ratio") +
      theme(legend.position = "none") +
      xlim(-x_limit,x_limit)
    list_plots[[i]] <- gg
  }
  ## return output list
  return(list_plots)
}



### Write function that creates ROC curve plots from ROC_function output ###
roc_plot_from_roc <- function(list_roc,
                              type_colors,
                              name){
  
  ## create empty plot window (full size)
  pdf(paste0("ROC_", name, "_all.pdf"), width=5, height=5)
  par(pty="s")
  par(mar=c(4,4,4,5))
  par(mgp=c(2.5,0.8,0))
  plot(1,1,xlim=c(0,1), ylim=c(0,1), cex.lab=1,  font.lab=2, type="n", ylab="True positive rate", xlab="False positive rate", xaxt="n", yaxt="n")
  axis(side=1, cex.axis=1, at=seq(0,1,by=0.2), cex.axis=0.9)
  axis(side=2, cex.axis=1, at=seq(0,1,by=0.2), cex.axis=0.9)
  abline(a=0,b=1, col="grey", lty="dashed")
  abline(v=0, col="grey")
  abline(v=1, col="grey")
  abline(h=0, col="grey")
  abline(h=1, col="grey")
  
  ## go over each ROC analysis to add lines
  N <- length(list_roc)
  for (i in 1:N){
    roc_i <- list_roc[[i]]
    lines(x=roc_i$ROC_coordinates[,2],
          y=roc_i$ROC_coordinates[,1],
          lwd=3, col=type_colors[i])
  }
  dev.off()
  
  ## create empty plot window (reduced size top-left corner)
  pdf(paste0("ROC_", name, "_small.pdf"), width=5, height=5)
  par(pty="s")
  par(mar=c(4,4,4,5))
  par(mgp=c(2.5,0.8,0))
  plot(1,1,xlim=c(0,0.1), ylim=c(0.9,1), cex.lab=1,  font.lab=2, type="n", ylab="True positive rate", xlab="False positive rate", xaxt="n", yaxt="n")
  axis(side=1, cex.axis=1, at=seq(0,1,by=0.05), cex.axis=0.9)
  axis(side=2, cex.axis=1, at=seq(0,1,by=0.05), cex.axis=0.9)
  abline(a=0,b=1, col="grey", lty="dashed")
  abline(v=0, col="grey")
  abline(v=1, col="grey")
  abline(h=0, col="grey")
  abline(h=1, col="grey")
  
  ## go over each ROC analysis to add lines
  for (i in 1:N){
    roc_i <- list_roc[[i]]
    lines(x=roc_i$ROC_coordinates[,2],
          y=roc_i$ROC_coordinates[,1],
          lwd=3, col=type_colors[i])
  }
  dev.off()
}








