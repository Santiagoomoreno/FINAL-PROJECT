library(DynamicCancerDriverKM)
library(e1071)
library(caret)
library(dplyr)
library(pROC)
library(tidyverse)
library(class)
library(rpart)
library(glmnet)

view(DynamicCancerDriverKM::BRCA_normal)
view(DynamicCancerDriverKM::BRCA_PT)

load("C:\\Users\\Usuario\\Documents\\GENESCORE\\geneScore.rdata")

npt_combined <- rbind(BRCA_normal, BRCA_PT)
df_modified <- npt_combined[, !(names(npt_combined) %in% c("barcode", "bcr_patient_barcode", "bcr_sample_barcode", "vital_status", "days_to_death", "treatments_radiation_treatment_or_therapy"))]
any_na <- any(is.na(df_modified))

sample_matrix <- as.matrix(df_modified[, -1])
threshold_value <- 0.0002 * max(sample_matrix)
expressed_genes <- sample_matrix > threshold_value
true_per_gene <- colSums(expressed_genes)
threshold_to_keep_column <- nrow(sample_matrix) * 0.2
columns_to_keep <- which(true_per_gene >= threshold_to_keep_column)
filtered_genes <- df_modified[, c(1, columns_to_keep + 1)]

gene_score <- prub$features
genes_mod <- colnames(filtered_genes)[-1]
common_genes <- intersect(gene_score, genes_mod)
common_genes_df <- prub[gene_score %in% common_genes, ]
sorted_genes <- common_genes_df %>% arrange(desc(score))
top_genes_modified <- sorted_genes[1:100, ]
top_100_genes <- top_genes_modified$features

response_variable <- filtered_genes$sample_type
X_modified <- filtered_genes[, top_100_genes]
response_variable <- as.factor(response_variable)

set.seed(200)
trainIndex_modified <- createDataPartition(response_variable, p = 0.8, list = FALSE)
train_data_mod <- X_modified[trainIndex_modified, ]
test_data_mod <- X_modified[-trainIndex_modified, ]
train_labels_mod <- response_variable[trainIndex_modified]
test_labels_mod <- response_variable[-trainIndex_modified]

model_svm <- svm(train_labels_mod ~ ., data = cbind(train_data_mod, train_labels_mod), kernel = "linear")
predictions_svm <- predict(model_svm, newdata = cbind(test_data_mod, test_labels_mod))
conf_matrix_svm <- confusionMatrix(predictions_svm, test_labels_mod)
conf_matrix_svm
roc_curve_svm <- roc(as.numeric(predictions_svm), as.numeric(test_labels_mod))
roc_curve_svm

logistic_model <- cv.glmnet(as.matrix(train_data_mod), train_labels_mod, family = "binomial")
predictions_logistic <- predict(logistic_model, newx = as.matrix(test_data_mod), s = "lambda.1se", type = "response")
predicted_labels_logistic <- as.factor(ifelse(predictions_logistic > 0.5, levels(response_variable)[2], levels(response_variable)[1]))
conf_matrix_logistic <- confusionMatrix(predicted_labels_logistic, test_labels_mod)
conf_matrix_logistic
precision_logistic <- posPredValue(predicted_labels_logistic, test_labels_mod)
paste("Precisión del modelo de regresión logística:", precision_logistic)
roc_curve_logistic <- roc(as.numeric(predicted_labels_logistic), as.numeric(test_labels_mod))
roc_curve_logistic

normalized_train_data_knn <- scale(train_data_mod)
normalized_test_data_knn <- scale(test_data_mod)
knn_model <- knn(train = normalized_train_data_knn, test = normalized_test_data_knn, cl = train_labels_mod, k = 5)
knn_conf_matrix <- confusionMatrix(knn_model, test_labels_mod)
knn_conf_matrix

tree_model <- rpart(train_labels_mod ~ ., data = train_data_mod, method = "class")
plot(tree_model)
text(tree_model, pretty = 0)
tree_predictions <- predict(tree_model, newdata = test_data_mod, type = "class")
tree_conf_matrix <- confusionMatrix(tree_predictions, test_labels_mod)
print(tree_conf_matrix)

