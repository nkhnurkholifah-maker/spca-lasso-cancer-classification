library(readr)
library(caret)
library(nnet)
library(glmnet)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
set.seed(42)
#===============================================================
#Section 1 : Load Data
#===============================================================
df <- read_csv("/Users/nurkholifah/Downloads/Doğrusal İstatistiksel Modeller/TCGA-PANCAN-HiSeq-801x20531/data.csv")
labels <- read_csv("/Users/nurkholifah/Downloads/Doğrusal İstatistiksel Modeller/TCGA-PANCAN-HiSeq-801x20531/labels.csv")
all(df$...1 == labels$...1) # cek kecocokan sample
y=factor(labels$Class)
X=df[,-1]
cat('Data dimension(nxp):',dim(X),"\n")
cat("Class distribution:\n")
print(table(y))

#===============================================================
#Section 2 : Train - Test Split (Stratified 70:30)
#===============================================================

cat("\n== Train-Test Split (Stratified 70:30) ==\n")
train_idx=createDataPartition(y,p=0.70,list=FALSE) #library_CARET
X_train=X[train_idx,]
X_test=X[-train_idx,]
y_train=y[train_idx] #hanya ada 1 variabel jadi tanpa (,)
y_test=y[-train_idx]
cat("Train size:",nrow(X_train),"|Test size",nrow(X_test),"\n")
cat("Class Distribution-Train:\n");print(table(y_train))
cat("Class Distribution-Test:\n");print(table(y_test))

#===============================================================
#Section 3 : Data Preprocessing (Standardisation)
#===============================================================
cat("\n== Standardisation ==\n")
preproc_params=preProcess(X_train,method=c("center","scale"))
X_train_sc=predict(preproc_params,X_train)
X_test_sc=predict(preproc_params,X_test) # Test set menggunakan mean & sd yang SAMA dari training set
cat("Standardisation completed. Mean of the first feature (should be close to 0):",round(mean(X_train_sc[[1]]),6),"\n")
cat("Standard deviation of the first feature (should be close to 1):",round(sd(X_train_sc[[1]]),6),"\n")

#===============================================================
#Section 4 : Supervised PCA (SPCA)
#Step 1. Calculate F-Statistic ANOVA for every Gene
#===============================================================

cat("\n== SPCA Procedure ==\n")

cat("\nStep 1: Calculate ANOVA F-Statistic for",ncol(X_train_sc),"gen...\n")
compute_anova_f=function(x_matrix,y_labels){
  x_matrix=as.matrix(x_matrix)
  storage.mode(x_matrix)="numeric"
  y_labels=factor(y_labels)
  classes=levels(y_labels)
  K=length(classes)
  n=nrow(x_matrix)
  p=ncol(x_matrix)
  grand_mean=colMeans(x_matrix)
  MSB=sapply(1:p,function(j){
    ssb=sum(sapply(classes,function(k){
      idx=which(y_labels==k)
      nk=length(idx)
      xbar_k=mean(x_matrix[idx,j])
      nk*(xbar_k-grand_mean[j])^2
    }))
    ssb/(K-1)
  })
  MSW=sapply(1:p,function(j){
    ssw=sum(sapply(classes,function(k){
      idx=which(y_labels==k)
      sum((x_matrix[idx,j]-mean(x_matrix[idx,j]))^2)
    }))
    ssw/(n-K)
  })
  F_stats=MSB/MSW
  names(F_stats)=colnames(x_matrix)
  F_stats[is.na(F_stats)]=0
  F_stats[is.infinite(F_stats)]=0
  return(F_stats)
}
F_stats=compute_anova_f(X_train_sc,y_train)

cat("Calculate F-Statistic. Top 5 gene\n")
print(sort(F_stats,decreasing = TRUE)[1:5])

#Step 2. Choose Top-k Gene (with 5- fold CV)
cat("\nStep 2. Choose Top-k with 5-fold Cross Validation\n")

k_candidates=c(100,500,1000)
spca_cv_accuracy <- function(k, x_mat, y_lab, n_folds = 5) {

  folds <- createFolds(y_lab, k = n_folds, list = TRUE)

  fold_acc <- sapply(folds, function(val_idx) {

    x_cv_train <- x_mat[-val_idx, ]
    x_cv_val   <- x_mat[val_idx, ]
    y_cv_train <- y_lab[-val_idx]
    y_cv_val   <- y_lab[val_idx]

    F_fold <- compute_anova_f(x_cv_train, y_cv_train)
    gene_order_fold <- names(sort(F_fold, decreasing = TRUE))
    top_genes <- gene_order_fold[1:k]

    x_cv_train_sub <- x_cv_train[, top_genes]
    x_cv_val_sub   <- x_cv_val[, top_genes]

    pca_cv <- prcomp(x_cv_train_sub, center = FALSE, scale. = FALSE)

    cum_var <- cumsum(pca_cv$sdev^2) / sum(pca_cv$sdev^2)
    n_pc <- max(2, which(cum_var >= 0.90)[1])

    Z_train <- pca_cv$x[, 1:n_pc, drop = FALSE]
    Z_val <- predict(pca_cv, newdata = x_cv_val_sub)[, 1:n_pc, drop = FALSE]

    df_cv <- data.frame(y = y_cv_train, Z_train)

    fit_cv <- multinom(y ~ ., data = df_cv, trace = FALSE, MaxNWts = 10000)

    pred_cv <- predict(fit_cv, newdata = data.frame(Z_val))

    mean(pred_cv == y_cv_val)
  })

  mean(fold_acc)
}
#Calculate Accuracy of CV for every k
cv_results <- data.frame(k = k_candidates, cv_accuracy = NA)

for(i in seq_along(k_candidates)){
  k_val <- k_candidates[i]
  cat("Calculate CV for k =", k_val, "genes...\n")

  cv_results$cv_accuracy[i] <- spca_cv_accuracy(
    k = k_val,
    x_mat = X_train_sc,
    y_lab = y_train
  )

  cat("k =", k_val,
      "CV accuracy =", round(cv_results$cv_accuracy[i], 4), "\n")
}

cat("The best CV:\n")
print(cv_results)

k_best <- cv_results$k[which.max(cv_results$cv_accuracy)]
cat("\nThe best k:", k_best, "genes\n")

#Step 3. PCA on Top k Selected Genes
cat("\n---Step 3: PCA on Top",k_best,"Selected Genes---\n")
top_k_genes=names(sort(F_stats,decreasing = TRUE))[1:k_best]
X_train_topk=as.matrix(X_train_sc[,top_k_genes])
#matrix n_train x k_best
X_test_topk=as.matrix(X_test_sc[,top_k_genes])
#matrix n_test x k_best

#Fit PCA with SVD
#center=FALSE & scale.=FALSE because data have standardisation in section 3
pca_model=prcomp(X_train_topk,center=FALSE,scale.=FALSE)
pve=pca_model$sdev^2/sum(pca_model$sdev^2)
cum_pve=cumsum(pve)
n_pc_best=max(2,which(cum_pve>=0.90)[1])
#min 2 PC
cat("Number of PCs retained (cumulative PVE>=90%):",n_pc_best,'\n')
cat("Cummulative PVE with",n_pc_best,"PCs:", round(cum_pve[n_pc_best]*100,2),"%\n")

#PC Score for training set
Z_train=pca_model$x[,1:n_pc_best,drop=FALSE]
Z_test=predict(pca_model,newdata=X_test_topk)[,1:n_pc_best,drop=FALSE]
cat("Z_train dimension:",dim(Z_train),"\n") #563xn_pc_best
cat("Z_test dimension:",dim(Z_test),"\n") #238xn_pc_best

#Step 4. Multinomial Logistic Regression on PC Score
cat("\n---Step 4: Multinomial Logistic Regression on PC Score---\n")
df_spca_train=data.frame(y=y_train,Z_train)
df_spca_test=data.frame(Z_test)
#Fit Model-MaxNWTs dinaikkan karena banyak PC x 5 kelas
spca_model=multinom(y~.,data=df_spca_train,trace=FALSE,MaxNWTs=1000000)
#Probability and Prediction Class on Test Set
spca_pred_class=predict(spca_model,newdata=df_spca_test,type="class")
spca_pred_prob=predict(spca_model,newdata=df_spca_test,type="probs")
spca_pred_class=factor(spca_pred_class,levels=levels(y_test))

cat("SPCA prediction done.\n")
cat("Confusion Matrix(Test set):\n")
print(table(Predicted =spca_pred_class,Actual=y_test))

#===============================================================
#Section 5 : LASSO Multinomial Logistic Regression
#===============================================================
cat("\n\n=== LASSO Multinomial Logistic Regression ===\n")
cat("Fitting LASSO with 10-fold CV on full",ncol(X_train_sc),"genes.... (may take a few minutes)\n")
X_train_mat=as.matrix(X_train_sc)  #glmnet butuh matrix bukan dataframe
X_test_mat=as.matrix(X_test_sc)
lasso_cv=cv.glmnet(
  x=X_train_mat,
  y=y_train,
  family="multinomial",
  alpha=1, #remember lasso 1, ridge 0
  nfolds=10,
  type.measure="class"
  )
lambda_min=lasso_cv$lambda.min
cat("Optimal lambda (lambda.min):",round(lambda_min,6),"\n")
#Class Prediction
lasso_pred_class=predict(lasso_cv,newx=X_test_mat,
                         s="lambda.min",
                         type="class")
lasso_pred_class=factor(as.vector(lasso_pred_class),levels=levels(y_test))

#Probability Prediction (AUC)
lasso_pred_prob=predict(lasso_cv,
                        newx=X_test_mat,
                        s="lambda.min",
                        type="response")
lasso_pred_prob=lasso_pred_prob[,,1] #from y(nxKx1) matrix (nxK) array
cat("LASSO prediction done.\n")
cat("Confusion Matrix (Test Set):\n")
print(table(Predicted=lasso_pred_class,Actual=y_test))
#Calculate Genes which choosen by LASSO (koef!= 0  minimal 1 class)
lasso_coef_list=coef(lasso_cv,s="lambda.min")
selected_per_class=lapply(lasso_coef_list,function(coef_k){
  nm=rownames(coef_k)[which(coef_k!=0)]
  nm[nm!="(Intercept)"]

})
selected_genes_lasso=unique(unlist(selected_per_class))
n_genes_lasso=length(selected_genes_lasso)
cat("Number of genes selected by LASSO:",n_genes_lasso,"\n")

#===============================================================
#Section 6 : Model Evaluation
#===============================================================
cat("\n\n=== Model Evaluation on Test Set ===\n")
evaluate_model=function(pred_class,pred_prob,true_class,model_name){
  #Confusion Matrix Detail
  cm=confusionMatrix(pred_class,true_class, mode="everything")
  #Overall accuracy
  accuracy=as.numeric(cm$overall["Accuracy"])
  #Per-class metrics from byClass
  sens=cm$byClass[,"Sensitivity"]
  spec=cm$byClass[,"Specificity"]
  prec=cm$byClass[,"Precision"]
  f1=cm$byClass[,"F1"]
  class_names=gsub("Class:","",rownames(cm$byClass))
  #Macro-averaged F1
  macro_f1=mean(f1,na.rm=TRUE)
  #Multiclass AUC
  colnames(pred_prob)=levels(true_class)
  roc_obj=multiclass.roc(response=true_class,
                         predictor=pred_prob,
                         levels=levels(true_class))
  multi_auc=as.numeric(roc_obj$auc)
  cat("\n",strrep("=",55),"\n",sep="")
  cat("MODEL:",model_name,"\n")
  cat(strrep("=",55),"\n",sep="")
  cat(sprintf("Overall Accuracy:%.4f (%.2f%%)\n",accuracy,accuracy*100))
  cat(sprintf("Macro F1-Score:%.4f\n",macro_f1))
  cat(sprintf("Multiclass AUC:%.4f\n",multi_auc))
  cat("\nPer-Class Metric:\n")
  per_class_df=data.frame(
    Class=class_names,
    Sensitivity=round(sens,4),
    Specificity=round(spec,4),
    Precision=round(prec,4),
    F1_Score=round(f1,4),
    row.names=NULL
  )
  print(per_class_df)
  invisible(list(
    accuracy=accuracy,
    macro_f1=macro_f1,
    auc=multi_auc,
    per_class=per_class_df,
    conf_table=cm$table
  ))
}

results_spca=evaluate_model(spca_pred_class,spca_pred_prob,y_test,
                           "SPCA+Multinomial LogReg")
results_lasso=evaluate_model(lasso_pred_class,lasso_pred_prob,y_test,
                            "LASSO Multinomial LogReg")

#===============================================================
#Section 7 : Summary Comparison & Feature Reduction
#===============================================================

cat("\n",strrep("=",55),"\n",sep="")
cat("Summary Comparison\n")
cat(strrep("=",55),"\n",sep="")
summary_df=data.frame(
  Metric=c("Overall accuracy","Macro F1-Score","Multiclass AUC"),
  SPCA=round(c(results_spca$accuracy,
               results_spca$macro_f1,
               results_spca$auc),4),
  LASSO=round(c(results_lasso$accuracy,
                results_lasso$macro_f1,
                results_lasso$auc),4)
)
print(summary_df)
cat("\n Feature Reduction Summary:\n")
feat_df=data.frame(
  Method=c("SPCA","LASSO"),
  Total_Genes=c(ncol(X_train_sc),ncol(X_train_sc)),
  Selected_Genes=c(k_best,n_genes_lasso),
  Final_Features=c(n_pc_best,n_genes_lasso),
  Reduction_pct=c(
    round((1-n_pc_best/ncol(X_train_sc))*100,2),
    round((1-n_genes_lasso/ncol(X_train_sc))*100,2)
  )
)
print(feat_df)

#===============================================================
#Section 8 : Visualization
#===============================================================
cat("\n\n=== Generating Plots ===\n")
plot_cm=function(cm_table,title_str,fill_high){
  df=as.data.frame(cm_table)
  colnames(df)=c("Predicted","Actual","Freq")
  ggplot(df,aes(x=Predicted,y=Actual,fill=Freq))+
    geom_tile(color="white",linewidth=0.5)+
    geom_text(aes(label=Freq),size=5,fontface="bold")+
    scale_fill_gradient(low="#F7FBFF",high=fill_high)+
    labs(title=title_str,x="Predicted Class",y="Actual Class")+
    theme_minimal(base_size=13)+
    theme(plot.title=element_text(face="bold",hjust=0.5),
          axis.text=element_text(size=11),
          legend.title=element_text(size=10))
}
p_cm_spca=plot_cm(results_spca$conf_table,
                  paste0("Confusion Matrix-SPCA\n Accuracy:",
                         round(results_spca$accuracy*100,2),"%"),
                  "#2171B5")
p_cm_lasso=plot_cm(results_lasso$conf_table,
                  paste0("Confusion Matrix-LASSO\n Accuracy:",
                         round(results_lasso$accuracy*100,2),"%"),
                  "#D94701")
ggsave("confusion_matrix_SPCA.png",p_cm_spca,width=7,height=5.5,dpi=300)
ggsave("confusion_matrix_LASSO.png",p_cm_lasso,width=7,height=5.5,dpi=300)
cat("Saved:confusion_matrix_SPCA.png & confusion_matrix_LASSO.png\n")

#F1-Score perClass: SPCA vs LASSO
cancer_classes=c("BRCA","KIRC","COAD","LUAD","PRAD")
f1_df=data.frame(
  Class=rep(cancer_classes,2),
  Model=rep(c("SPCA","LASSO"),each=5),
  F1=c(results_spca$per_class$F1_Score,
       results_lasso$per_class$F1_Score)
)
p_f1=ggplot(f1_df,aes(x=Class,y=F1,fill=Model))+
  geom_col(position = "dodge",width=0.65, color="white")+
  geom_text(aes(label=sprintf("%.3f",F1)),
            position = position_dodge(0.65),vjust=-0.4,size=3.5)+
  scale_fill_manual(values=c("SPCA"="#2171B5","LASSO"="#D94701"))+
  scale_y_continuous(limits=c(0,1.1),labels=percent_format(accuracy=1))+
  labs(title="F1-Score per Cancer Class:SPCA vs LASSO",
       x="Cancer Type",y="F-Score",fill="Method")+
  theme_minimal(base_size = 13)+
  theme(plot.title = element_text(face="bold",hjust=0.5),
        legend.position = "top")
ggsave("f1_per_class_comparison.png",p_f1,width=8,height=5,dpi=300)
cat("Saved:f1_per_class_comparison.png\n")

#Scree plot PCA
max_pc=min(30,length(pve))
scree_df=data.frame(
  PC=1:max_pc,
  PVE=pve[1:max_pc]*100,
  CumPVE=cum_pve[1:max_pc]*100
)
p_scree=ggplot(scree_df,aes(x=PC))+
  geom_col(aes(y=PVE),fill="#9ECAE1",color="white")+
  geom_line(aes(y=CumPVE),color="#08519C",linewidth=1,group=1)+
  geom_point(aes(y=CumPVE),color="#08519C",size=2)+
  geom_hline(yintercept = 90,linetype="dashed",color="red")+
  annotate("text",x=max_pc*0.8,y=91.5,
           label="90% threshold",color="red",size=4)+
  annotate("point",x=n_pc_best,y=cum_pve[n_pc_best]*100,
           shape=17,color="red",size=4)+
  annotate("text",x=n_pc_best+1,y=cum_pve[n_pc_best]*100-4,
           label=paste0("PC",n_pc_best),color="red",size=3.5)+
  labs(title=paste0("Scree Plot-PCA on Top-",k_best,"Genes"),
       subtitle=paste0(n_pc_best,"PCs explain>=90% variances"),
       x="Principal Component",y="Variance Explained (%)")+
  theme_minimal(base_size = 13)+
  theme(plot.title = element_text(face="bold"))
ggsave("scree_plot.png",p_scree,width=8,height=5,dpi=300)
cat("Saved:scree_plot.png\n")

#PCA Score Plot (PC1 vs PC2)
pca_df=data.frame(
  PC1=Z_train[,1],
  PC2=Z_train[,2],
  Class=y_train
)
p_pca=ggplot(pca_df,aes(x=PC1,y=PC2,color=Class))+
  geom_point(alpha=0.65,size=2)+
  stat_ellipse(level=0.90,linetype="dashed",linewidth=0.7)+
  scale_color_brewer(palette="Set1")+
  labs(title=paste0("PCA Score Plot(Top-",k_best,"Genes)"),
       subtitle=paste0("PC1:", round(pve[1]*100,1),
                       "%PC2:", round(pve[2]*100,1),"%"),
       x=paste0("PC1(",round(pve[1]*100,1),"%)"),
       y=paste0("PC2(",round(pve[2]*100,1),"%)"),
       color="Cancer Type")+
  theme_minimal(base_size = 13)+
  theme(plot.title=element_text(face="bold"),legend.position="right")
ggsave("pca_score_plot.png",p_pca,width=8,height=6,dpi=300)
cat("Saved:pca_score_plot.png\n")

#CV k-selection bar chart
cv_df=cv_results
p_cvk=ggplot(cv_df,aes(x=factor(k),y=cv_accuracy,
                       fill=factor(k==k_best)))+
  geom_col(width=0.5,color="white")+
  geom_text(aes(label=sprintf("%.2f%%",cv_accuracy*100)),
            vjust=-0.4,size=4,fontface="bold")+
  scale_fill_manual(values=c("FALSE"="#9ECAE1","TRUE"="#08519C"),
                    guide="none")+
  scale_y_continuous(limits=c(0,1.15),labels=percent_format(accuracy = 1))+
  labs(title="5-Fold CV Accuracy for Gene Selection (SPCA)",
       subtitle = paste0("Best k=",k_best,"genes(darker bar)"),
       x="Number of Selected Genes (k)",y="CV Accuracy")+
  theme_minimal(base_size = 13)+
  theme(plot.title=element_text(face="bold"))
ggsave("cv_k_selection.png",p_cvk,width=6,height=5,dpi=300)
cat("Saved:cv_k_selection.png\n")

#LASSO CV Lambda Plot
png("lasso_cv_lambda.png",width=900,height=550,res=150)
plot(lasso_cv,main="LASSO 10-Fold CV:Misclassification Error vs log(lambda)")
abline(v=log(lambda_min),col="red",lty=2,lwd=2)
legend("topleft",paste0("lamda.min=",round(lambda_min,5)),
       col="red",lty=2,lwd=2,bty="n")
dev.off()
cat("Saved:lasso_cv_lambda.png\n")


#===============================================================
#Export Tabel CSV
#===============================================================
write.csv(summary_df,"table1_main_metrics.csv",row.names=FALSE)
write.csv(results_spca$per_class,"table2_perclass_SPCA.csv",row.names = FALSE)
write.csv(results_lasso$per_class,"table3_perclass_LASSO.csv",row.names = FALSE)
write.csv(feat_df,"table4_feature_reduction.csv",row.names = FALSE)
cat("\nAll CSV table save.\n")

