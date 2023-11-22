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
view(DynamicCancerDriverKM::PPI)

normal_pt <- rbind(BRCA_normal, BRCA_PT)

df <- normal_pt[, !(names(normal_pt) %in% c("barcode", "bcr_patient_barcode", "bcr_sample_barcode", "vital_status", "days_to_death", "treatments_radiation_treatment_or_therapy"))]

any_na <- any(is.na(df))

samples_matrix <- as.matrix(df[, -1])
threshold <- 0.0002 * max(samples_matrix)
expressed_genes <- samples_matrix > threshold
true_per_gene <- colSums(expressed_genes)
threshold_column <- nrow(samples_matrix) * 0.2
columns_to_keep <- which(true_per_gene >= threshold_column)
filtered_df <- df[, c(1, columns_to_keep + 1)]

gene_names <- colnames(filtered_df)[-1]
gene_names_conversion <- AMCBGeneUtils::changeGeneId(gene_names, from = "Ensembl.ID")
colnames(filtered_df)[-1] <- gene_names_conversion$HGNC.symbol

ppi_genes <- PPI$`Input-node Gene Symbol`
genes <- colnames(filtered_df)[-1]
common_genes <- intersect(ppi_genes, genes)
common_genes_df <- PPI[ppi_genes %in% common_genes, ]

input_gene_count <- common_genes_df %>%
  group_by(`Input-node Gene Symbol`) %>%
  summarize(count_input = n())

output_gene_count <- common_genes_df %>%
  group_by(`Output-node Gene Symbol`) %>%
  summarize(count_output = n())

counts <- full_join(input_gene_count, output_gene_count, by = c("Input-node Gene Symbol" = "Output-node Gene Symbol"))
counts[is.na(counts)] <- 0
counts$total_count <- rowSums(counts[, c("count_input", "count_output")])
sorted_counts <- counts %>%
  arrange(desc(total_count))
top_genes <- head(sorted_counts, 100)$`Input-node Gene Symbol`

y <- filtered_df$sample_type
X <- filtered_df[, top_genes]
y <- as.factor(y)

set.seed(200)
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
train_data <- X[trainIndex, ]
test_data <- X[-trainIndex, ]
train_labels <- y[trainIndex]
test_labels <- y[-trainIndex]

model_svm <- svm(train_labels ~ ., data = cbind(train_data, train_labels), kernel = "linear")
predictions_svm <- predict(model_svm, newdata = cbind(test_data, test_labels))
confusionMatrix(predictions_svm, test_labels)
roc_curve_svm <- roc(as.numeric(predictions_svm), as.numeric(test_labels))

logistic_model <- cv.glmnet(as.matrix(train_data), train_labels, family = "binomial")
predictions_logistic <- predict(logistic_model, newx = as.matrix(test_data), s = "lambda.1se", type = "response")
predicted_labels_logistic <- as.factor(ifelse(predictions_logistic > 0.5, levels(y)[2], levels(y)[1]))
conf_matrix_logistic <- confusionMatrix(predicted_labels_logistic, test_labels)
precision_logistic <- posPredValue(predicted_labels_logistic, test_labels)

normalized_train_data_knn <- scale(train_data)
normalized_test_data_knn <- scale(test_data)
knn_model <- knn(train = normalized_train_data_knn, test = normalized_test_data_knn, cl = train_labels, k = 5)
knn_conf_matrix <- confusionMatrix(knn_model, test_labels)

tree_model <- rpart(train_labels ~ ., data = train_data, method = "class")
tree_predictions <- predict(tree_model, newdata = test_data, type = "class")
tree_conf_matrix <- confusionMatrix(tree_predictions, test_labels)

plot(tree_model)
text(tree_model, pretty = 0)

print(confusionMatrix(predictions_svm, test_labels))
print(roc_curve_svm)

print(conf_matrix_logistic)
print(precision_logistic)

print(knn_conf_matrix)

print(tree_conf_matrix)

