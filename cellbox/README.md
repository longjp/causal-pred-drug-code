## Comparison of Cellbox model and Linear Regression on Melanoma cell line perturbation data.

### data folder

These data are copied from [https://github.com/sanderlab/CellBox](https://github.com/sanderlab/CellBox). See that repository for more information.

```expert.csv```: Protein expression data from RPPA for the protein and phenotype responses. Each row is a condition and each column is a node. Additional columns encode drug applied in condition. 

```random_partition_average_testhat_929.csv```: Predicted responses from CellBox model after randomly split the data 1000 times for training and testing.

### code_script folder

```LR_RF.ipynb```: Fit linear model for random fold validation and make scatterplot comparing predictions to experimental values.

```Cellbox_RF.ipynb```: Regenerate scatterplot of CellBox predictions versus and experimental values across all conditions.

```LODO.ipynb```: Performs LODO validation for linear regression models and makes a dotplot and boxplot for comparing results with Cellbox.
