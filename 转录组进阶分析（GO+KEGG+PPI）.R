转录组进阶分析（GO+KEGG+PPI）：如何从差异表达基因挖掘生物学意义？
接上期： 我们已经完成了差异表达分析，获得了显著差异表达基因列表。现在让我们深入挖掘这些基因的生物学意义！
模块七：差异表达结果的深度解读

差异基因的功能分类与注释
R代码：差异基因功能分类

# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# 读取上期获得的差异表达基因
# significant_genes <- read.csv("significant_genes.csv")

# 提取基因列表（假设我们已经有了significant_genes数据框）
deg_genes <- significant_genes$gene

# 基因ID转换（从基因名到Entrez ID）
gene_df <- bitr(deg_genes, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)

# GO富集分析 - 更详细的设置
go_bp <- enrichGO(gene = gene_df$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",           # 生物过程
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

go_mf <- enrichGO(gene = gene_df$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", 
                  ont = "MF",           # 分子功能
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

go_cc <- enrichGO(gene = gene_df$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "CC",           # 细胞组分
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

# 查看最显著的GO条目
head(go_bp@result, 10)

富集结果可视化
R代码：富集分析可视化

# 富集结果点图
dotplot(go_bp, showCategory=15, title="GO Biological Process") +
  theme(axis.text.y = element_text(size=10))

# 富集网络图
emapplot(go_bp, showCategory=20)

# 富集结果条形图
barplot(go_bp, showCategory=15, title="GO Biological Process") +
  theme(axis.text.y = element_text(size=10))

# KEGG通路富集分析
kegg_enrich <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = 'hsa',           # 人类
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# KEGG结果可视化
dotplot(kegg_enrich, title="KEGG Pathway Enrichment")
模块八：蛋白质互作网络(PPI)分析

使用STRING数据库获取PPI数据
R代码：PPI网络数据获取

# 安装和加载STRINGdb包
# BiocManager::install("STRINGdb")
library(STRINGdb)

# 创建STRING数据库连接
string_db <- STRINGdb$new(version="11.5", 
                          species=9606,  # 人类
                          score_threshold=400,  # 互作得分阈值
                          input_directory="")

# 映射基因到STRING ID
mapped_genes <- string_db$map(gene_df, "ENTREZID", removeUnmappedRows=TRUE)

# 获取互作网络
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# 查看互作网络基本信息
print(paste("获得互作对数量:", nrow(interactions)))
print(paste("涉及基因数量:", length(unique(c(interactions$from, interactions$to)))))


注意： STRING数据库需要网络连接，首次运行时会自动下载数据库文件。

PPI网络构建与分析
R代码：PPI网络分析与可视化


# 加载网络分析包
library(igraph)

# 创建igraph网络对象
ppi_network <- graph_from_data_frame(
  d = interactions[, c("from", "to", "combined_score")],
  directed = FALSE
)

# 网络基本信息
print(paste("节点数:", vcount(ppi_network)))
print(paste("边数:", ecount(ppi_network)))
print(paste("网络密度:", graph.density(ppi_network)))

# 计算节点中心性指标
degree_centrality <- degree(ppi_network)
betweenness_centrality <- betweenness(ppi_network)
closeness_centrality <- closeness(ppi_network)

# 创建节点属性数据框
node_attributes <- data.frame(
  gene_id = V(ppi_network)$name,
  degree = degree_centrality,
  betweenness = betweenness_centrality,
  closeness = closeness_centrality
)

# 识别Hub基因（度中心性前10）
hub_genes <- node_attributes[order(-node_attributes$degree), ][1:10, ]
print("Top 10 Hub基因:")
print(hub_genes)

PPI网络可视化
R代码：PPI网络可视化

# 设置节点大小和颜色基于度中心性
V(ppi_network)$size <- sqrt(degree_centrality) * 2
V(ppi_network)$color <- ifelse(
  V(ppi_network)$name %in% hub_genes$gene_id[1:5], 
  "red", "lightblue"
)

# 设置边的宽度基于互作得分
E(ppi_network)$width <- E(ppi_network)$combined_score / 200

# 绘制网络图
set.seed(123)  # 保证布局可重复
plot(ppi_network,
     layout = layout_with_fr(ppi_network),
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.frame.color = "gray",
     main = "蛋白质互作网络(PPI) - 差异表达基因")

# 添加图例
legend("bottomright", 
       legend = c("Hub基因", "其他基因"),
       fill = c("red", "lightblue"),
       cex = 0.8)
模块九：关键调控基因识别

整合表达变化与网络中心性
R代码：关键基因识别

# 合并表达变化与网络属性
key_genes_analysis <- merge(
  significant_genes, 
  node_attributes, 
  by.x = "gene", 
  by.y = "gene_id",
  all.x = TRUE
)

# 计算综合评分（考虑表达变化和网络中心性）
key_genes_analysis$composite_score <- 
  abs(key_genes_analysis$log2FoldChange) * 
  key_genes_analysis$degree * 
  key_genes_analysis$betweenness

# 筛选潜在的关键调控基因
potential_regulators <- key_genes_analysis[
  order(-key_genes_analysis$composite_score), 
][1:20, ]

# 添加基因描述信息（需要相应的注释包）
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

gene_descriptions <- select(org.Hs.eg.db, 
                            keys = potential_regulators$gene,
                            columns = c("SYMBOL", "GENENAME"),
                            keytype = "SYMBOL")

potential_regulators <- merge(potential_regulators, 
                              gene_descriptions, 
                              by.x = "gene", 
                              by.y = "SYMBOL")

print("潜在关键调控基因:")
print(potential_regulators[, c("gene", "log2FoldChange", "degree", "composite_score", "GENENAME")])

关键基因的功能验证
R代码：关键基因功能富集

# 对关键基因进行专门的功能富集分析
key_gene_entrez <- bitr(potential_regulators$gene, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)

key_go_enrich <- enrichGO(
  gene = key_gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.01,  # 更严格的阈值
  qvalueCutoff = 0.05,
  readable = TRUE
)

# 可视化关键基因的富集结果
dotplot(key_go_enrich, 
        title = "关键调控基因的GO富集分析",
        showCategory = 15) +
  theme(axis.text.y = element_text(size=9))

模块十：结果整合与生物学解释

创建综合分析报告
R代码：结果整合与报告生成

# 创建结果汇总数据框
final_results <- potential_regulators[, c("gene", "log2FoldChange", "padj", 
                                          "degree", "betweenness", "composite_score", "GENENAME")]

# 添加功能注释
top_go_terms <- head(go_bp@result$Description, 5)
top_kegg_pathways <- head(kegg_enrich@result$Description, 5)

# 生成总结报告
analysis_summary <- list(
  total_deg = nrow(significant_genes),
  up_regulated = sum(significant_genes$log2FoldChange > 0),
  down_regulated = sum(significant_genes$log2FoldChange < 0),
  ppi_network_nodes = vcount(ppi_network),
  ppi_network_edges = ecount(ppi_network),
  top_hub_genes = hub_genes$gene_id[1:5],
  key_regulators = potential_regulators$gene[1:10],
  significant_go_terms = top_go_terms,
  significant_pathways = top_kegg_pathways
)

print("分析总结:")
print(analysis_summary)

# 保存关键结果
write.csv(final_results, "key_regulatory_genes.csv", row.names = FALSE)
write.csv(analysis_summary, "analysis_summary.csv", row.names = FALSE)

生物学意义解读指南
步骤1：识别核心调控网络
关注Hub基因和高度连接的模块，这些通常是信号转导通路的核心组件。

步骤2：验证功能一致性
检查关键基因的富集通路是否与实验处理的预期效应一致。

步骤3：寻找新的调控关系
注意那些表达变化显著且在网络中处于关键位置，但功能未知的基因。

步骤4：构建调控假说
基于网络拓扑和表达模式，提出可能的调控机制并进行实验验证。

实用技巧与注意事项
网络分析最佳实践：

使用合适的互作得分阈值（通常400-700）平衡网络规模和质量
考虑使用Cytoscape进行更复杂的网络可视化和分析
结合其他证据（如文献挖掘）验证关键基因的重要性
注意不同组织或条件下的PPI网络可能有所不同
生物学解释注意事项：

相关性不等于因果关系
转录组变化可能不直接反映蛋白水平变化
需要考虑转录后调控和蛋白修饰的影响
实验验证是确认生物学功能的关键步骤
扩展分析方向
时间序列分析：
研究基因表达的动态变化模式
转录因子调控网络：
整合ChIP-seq数据识别直接调控关系
多组学整合：
结合蛋白质组、代谢组数据获得更全面的认识
机器学习应用：
使用随机森林等算法识别生物标志物


下期预告： 我们将探索多组学整合分析策略，展示如何将转录组数据与蛋白质组、代谢组数据结合，构建完整的生物调控网络！

往期推荐：
转录组测序分析完整指南：从数据到发现

点击蓝字

关注我们

科研猫猫猫

图片
图片
微信号丨x17585577064

B站丨稳妥学长




留言
都在搜:转录组分析全流程指南
写留言