多组学整合分析：构建完整的生物调控网络
接上期： 我们已经深入分析了转录组数据，识别了关键调控基因。现在让我们将视野扩展到多组学层面，整合转录组、蛋白质组和代谢组数据，揭示完整的生物调控网络！

模块十一：多组学数据整合策略
多组学数据类型与预处理
转录组数据
基因表达矩阵
差异表达分析结果
基因功能注释
蛋白质组数据
蛋白丰度矩阵
磷酸化修饰数据
蛋白互作网络
代谢组数据
代谢物浓度矩阵
代谢通路信息
代谢物鉴定结果
整合分析工具
MOFA+
  mixOmics
iCluster
R代码：多组学数据加载与预处理

# 加载必要的包
library(tidyverse)
library(mixOmics)
library(MOFA2)
library(ComplexHeatmap)
library(circlize)

# 设置随机种子保证可重复性
set.seed(123)

# 模拟多组学数据
# 转录组数据 (100个基因，10个样本)
transcriptome_data <- matrix(rnorm(100*10, mean=10, sd=2), nrow=100, ncol=10)
rownames(transcriptome_data) <- paste0("Gene_", 1:100)
colnames(transcriptome_data) <- paste0("Sample_", 1:10)

# 蛋白质组数据 (80个蛋白，10个样本)
proteome_data <- matrix(rnorm(80*10, mean=8, sd=1.5), nrow=80, ncol=10)
rownames(proteome_data) <- paste0("Protein_", 1:80)
colnames(proteome_data) <- paste0("Sample_", 1:10)

# 代谢组数据 (50个代谢物，10个样本)
metabolome_data <- matrix(rnorm(50*10, mean=12, sd=3), nrow=50, ncol=10)
rownames(metabolome_data) <- paste0("Metabolite_", 1:50)
colnames(metabolome_data) <- paste0("Sample_", 1:10)

# 样本分组信息
sample_groups <- data.frame(
  sample_id = paste0("Sample_", 1:10),
  group = rep(c("Treatment", "Control"), each=5),
  batch = rep(c("Batch1", "Batch2"), 5)
)

print("多组学数据维度:")
print(paste("转录组:", dim(transcriptome_data)[1], "基因 x", dim(transcriptome_data)[2], "样本"))
print(paste("蛋白质组:", dim(proteome_data)[1], "蛋白 x", dim(proteome_data)[2], "样本"))
print(paste("代谢组:", dim(metabolome_data)[1], "代谢物 x", dim(metabolome_data)[2], "样本"))
模块十二：多组学数据整合分析
使用MOFA+进行多组学因子分析
R代码：MOFA+多组学整合分析

# 准备MOFA数据
multiomics_data <- list(
  transcriptome = transcriptome_data,
  proteome = proteome_data,
  metabolome = metabolome_data
)

# 创建MOFA对象
mofa_object <- create_mofa(multiomics_data)

# 设置模型选项
model_options <- get_default_model_options(mofa_object)
model_options$num_factors <- 5  # 设置因子数量

# 设置训练选项
train_options <- get_default_training_options(mofa_object)
train_options$convergence_mode <- "slow"
train_options$seed <- 123

# 准备MOFA对象
mofa_prepared <- prepare_mofa(
  object = mofa_object,
  model_options = model_options,
  training_options = train_options
)

# 运行MOFA分析
mofa_trained <- run_mofa(mofa_prepared)

print("MOFA分析完成!")
print(paste("获得的因子数量:", get_dimensions(mofa_trained)$K))
因子解释与可视化
R代码：MOFA结果解释与可视化

# 查看因子方差解释度
variance_explained <- get_variance_explained(mofa_trained)
print("各因子对多组学数据的方差解释度:")
print(variance_explained$r2_total)

# 绘制方差解释度热图
plot_variance_explained(mofa_trained)

# 样本在因子空间的分布
factors <- get_factors(mofa_trained)[[1]]

# 添加样本分组信息
factors_with_groups <- as.data.frame(factors)
factors_with_groups$group <- sample_groups$group
factors_with_groups$batch <- sample_groups$batch

# 绘制因子散点图
ggplot(factors_with_groups, aes(x = Factor1, y = Factor2, color = group, shape = batch)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "样本在多组学因子空间的分布",
       x = "Factor 1",
       y = "Factor 2") +
  scale_color_manual(values = c("Control" = "blue", "Treatment" = "red"))
模块十三：跨组学关联网络分析
构建多组学关联网络
R代码：跨组学相关性网络分析

# 计算跨组学相关性
library(WGCNA)

# 选择代表性的特征进行关联分析
transcriptome_features <- transcriptome_data[1:30, ]  # 前30个基因
proteome_features <- proteome_data[1:25, ]           # 前25个蛋白
metabolome_features <- metabolome_data[1:20, ]       # 前20个代谢物

# 合并所有特征
all_features <- rbind(transcriptome_features, proteome_features, metabolome_features)

# 计算特征间相关性矩阵
feature_correlation <- cor(t(all_features), method = "spearman")

# 设置相关性阈值
correlation_threshold <- 0.7
significant_correlations <- abs(feature_correlation) > correlation_threshold

# 创建关联网络边列表
network_edges <- which(significant_correlations & upper.tri(significant_correlations), arr.ind = TRUE)
network_edges_df <- data.frame(
  from = rownames(feature_correlation)[network_edges[, 1]],
  to = rownames(feature_correlation)[network_edges[, 2]],
  correlation = feature_correlation[network_edges]
)

print(paste("获得的显著关联数量:", nrow(network_edges_df)))

多组学网络可视化
R代码：多组学网络构建与可视化

# 构建多组学网络
library(igraph)

multiomics_network <- graph_from_data_frame(
  d = network_edges_df,
  directed = FALSE
)

# 设置节点属性（根据组学类型）
V(multiomics_network)$type <- ifelse(
  grepl("Gene", V(multiomics_network)$name), "Transcriptome",
  ifelse(grepl("Protein", V(multiomics_network)$name), "Proteome", "Metabolome")
)

# 设置节点颜色
type_colors <- c("Transcriptome" = "#e74c3c", "Proteome" = "#3498db", "Metabolome" = "#2ecc71")
V(multiomics_network)$color <- type_colors[V(multiomics_network)$type]

# 设置节点大小基于度中心性
V(multiomics_network)$size <- sqrt(degree(multiomics_network)) * 3

# 设置边的颜色和宽度基于相关性
E(multiomics_network)$color <- ifelse(E(multiomics_network)$correlation > 0, "red", "blue")
E(multiomics_network)$width <- abs(E(multiomics_network)$correlation) * 3

# 绘制多组学网络
set.seed(123)
plot(multiomics_network,
     layout = layout_with_fr(multiomics_network),
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.frame.color = "gray",
     main = "多组学关联网络")

# 添加图例
legend("bottomright", 
       legend = names(type_colors),
       fill = type_colors,
       cex = 0.8,
       title = "组学类型")
模块十四：通路水平整合分析
多组学通路富集分析
R代码：整合通路富集分析

# 加载通路数据库
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# 准备多组学特征列表
# 假设我们已经识别了重要的多组学特征
important_features <- c(
  "Gene_1", "Gene_5", "Gene_12",  # 重要基因
  "Protein_3", "Protein_8", "Protein_15",  # 重要蛋白
  "Metabolite_2", "Metabolite_7"  # 重要代谢物
)

# 基因和蛋白的通路富集分析
gene_entrez <- bitr(grep("Gene", important_features, value = TRUE),
                    fromType = "SYMBOL", toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

protein_entrez <- bitr(grep("Protein", important_features, value = TRUE),
                       fromType = "SYMBOL", toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)

# 合并基因和蛋白
all_entrez <- unique(c(gene_entrez$ENTREZID, protein_entrez$ENTREZID))

# Reactome通路富集分析
reactome_enrich <- enrichPathway(
  gene = all_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

# 可视化通路富集结果
dotplot(reactome_enrich, 
        title = "多组学整合通路富集分析",
        showCategory = 15) +
  theme(axis.text.y = element_text(size = 9))
代谢通路可视化
R代码：代谢通路整合分析

# 代谢通路分析
library(pathview)

# 准备代谢物数据（转换为KEGG ID）
# 这里使用示例代谢物，实际分析中需要正确的KEGG ID
example_metabolites <- c("C00031", "C00022", "C00118")  # 葡萄糖, 丙酮酸, 甘油-3-磷酸

# 获取相关代谢通路
metabolite_pathways <- search_kegg_compound(example_metabolites)

# 可视化特定代谢通路（示例：糖酵解通路）
# pathview(gene.data = NULL, 
#          cpd.data = example_metabolites,
#          pathway.id = "hsa00010",  # 糖酵解通路
#          species = "hsa",
#          out.suffix = "multiomics_integration")

print("代谢通路分析完成!")
print(paste("分析的代谢物数量:", length(example_metabolites)))
模块十五：生物学解释与机制推断
构建调控机制假说
R代码：调控网络机制分析


# 整合多组学证据构建调控假说
library(visNetwork)

# 创建机制网络节点
mechanism_nodes <- data.frame(
  id = 1:8,
  label = c("Drug Treatment", "TF Activation", "Gene Expression", 
            "Protein Synthesis", "Enzyme Activity", "Metabolite Production",
            "Pathway Activation", "Phenotype"),
  group = c("Input", "Regulation", "Transcriptome", "Proteome", 
            "Function", "Metabolome", "Pathway", "Output"),
  value = c(10, 8, 6, 6, 7, 6, 8, 10)
)

# 创建机制网络边
mechanism_edges <- data.frame(
  from = c(1, 2, 3, 4, 5, 6, 7),
  to = c(2, 3, 4, 5, 6, 7, 8),
  label = c("induces", "regulates", "translates to", "affects", 
            "produces", "activates", "leads to")
)

# 交互式可视化调控机制
visNetwork(mechanism_nodes, mechanism_edges, 
           main = "多组学调控机制假说", 
           submain = "基于整合分析结果构建") %>%
  visGroups(groupname = "Input", color = "#e74c3c") %>%
  visGroups(groupname = "Regulation", color = "#3498db") %>%
  visGroups(groupname = "Transcriptome", color = "#9b59b6") %>%
  visGroups(groupname = "Proteome", color = "#2ecc71") %>%
  visGroups(groupname = "Function", color = "#f1c40f") %>%
  visGroups(groupname = "Metabolome", color = "#e67e22") %>%
  visGroups(groupname = "Pathway", color = "#1abc9c") %>%
  visGroups(groupname = "Output", color = "#34495e") %>%
  visLegend() %>%
  visPhysics(stabilization = FALSE)
生成综合分析报告
R代码：多组学整合报告生成

# 生成多组学整合分析总结
multiomics_summary <- list(
  datasets_integrated = length(multiomics_data),
  total_features = sum(sapply(multiomics_data, nrow)),
  samples_analyzed = ncol(transcriptome_data),
  factors_identified = get_dimensions(mofa_trained)$K,
  network_nodes = vcount(multiomics_network),
  network_edges = ecount(multiomics_network),
  enriched_pathways = nrow(reactome_enrich@result),
  key_insights = c(
    "转录组-蛋白质组协调变化揭示翻译调控",
    "代谢物积累反映通路活性改变",
    "多组学因子识别样本异质性"
  )
)

print("多组学整合分析总结:")
print(multiomics_summary)

# 保存关键结果
saveRDS(mofa_trained, "mofa_integration_results.rds")
saveRDS(multiomics_network, "multiomics_network.rds")
write.csv(network_edges_df, "cross_omics_correlations.csv", row.names = FALSE)

# 生成分析报告摘要
report_data <- data.frame(
  Analysis_Step = names(multiomics_summary),
  Results = sapply(multiomics_summary, function(x) {
    if(is.character(x)) paste(x, collapse = "; ") else as.character(x)
  })
)

write.csv(report_data, "multiomics_integration_report.csv", row.names = FALSE)
实用技巧与最佳实践
多组学整合分析工作流
数据质控：
确保各組学数据质量一致
数据标准化：
统一数据尺度，消除技术偏差
特征选择：
选择生物学相关的特征进行整合
整合分析：
使用适当算法发现跨组学模式
生物学解释：
结合领域知识解释分析结果
实验验证：
设计实验验证关键发现
多组学整合分析优势：

提供更全面的生物学视角
识别传统单组学分析遗漏的模式
增强生物学发现的可靠性
支持机制性假说的构建
注意事项：

不同组学数据的噪声水平和动态范围不同
需要考虑样本匹配和技术批次效应
生物学解释需要谨慎，避免过度解读
计算资源需求较高，需要合理规划
扩展应用与前沿方向
单细胞多组学：
整合scRNA-seq, scATAC-seq, CITE-seq数据
空间多组学：
结合空间转录组和蛋白质组数据
时间序列多组学：
分析动态生物过程
机器学习整合：
使用深度学习发现复杂模式
临床数据整合：
结合电子病历和影像数据

