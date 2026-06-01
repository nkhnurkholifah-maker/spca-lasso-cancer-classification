# SPCA vs LASSO for Multiclass Cancer Classification

## Overview

This repository contains the R code and results for a comparative study of:

- Supervised Principal Component Analysis (SPCA) + Multinomial Logistic Regression
- LASSO Multinomial Logistic Regression

applied to the TCGA PANCAN RNA-Seq dataset.

## Dataset

TCGA PANCAN RNA-Seq Dataset

- Samples: 801
- Genes: 20,531
- Cancer Types:
  - BRCA
  - COAD
  - KIRC
  - LUAD
  - PRAD

Dataset source:

https://archive.ics.uci.edu/dataset/401/gene+expression+cancer+rna+seq

## Methods

### SPCA Pipeline

1. Standardization
2. ANOVA F-statistic gene screening
3. 5-fold cross-validation for k selection
4. PCA
5. Multinomial Logistic Regression

### LASSO Pipeline

1. Standardization
2. 10-fold cross-validation
3. LASSO feature selection
4. Multinomial Logistic Regression

## Main Results

| Method | Final Model Input | Accuracy |
|----------|----------|----------|
| SPCA | 6 Principal Components | 100% |
| LASSO | 33 Genes | 100% |

