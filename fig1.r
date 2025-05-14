
rm(list=ls())
setwd("/home/zhangxp/testt")
###########lib############
library(limma)
library(Seurat)
library(ggplot2)
library(stringr)
library(tibble)
library(patchwork)
library(clustree)
library(DoubletFinder)
library(harmony)
library(SingleR)
library(celldex)
library(infercnv)
library(RColorBrewer)
#BiocManager::install("celldex")
#BiocManager::install("infercnv")



############read all samples########
samples <- list.files(pattern = "^(CID|GSM)", full.names = F)#^(CID|GSM)
samples <- samples[!grepl("(*_HER2_*|_Epi$|*_Total$)", samples)]
to_remove <- c("CID4463_ER", "CID4471_ER", "CID4530N_ER", "CID4535_ER", "CID3946_TN", 
               "CID4465_TN", "GSM4909272_N_MH0021_Total", "GSM4909286_TN_B1_MH0131")
samples <- setdiff(samples, to_remove)
samples

seurat_list <- list()
# 遍历样本列表
for (sample in samples) {
  print(sample)
  seurat_data <- Read10X(data.dir = sample, gene.column = 1)
  # 使用文件夹名称作为项目名称
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   project = sample,  # 保留完整的 sample 名称
                                   min.features = 200, 
                                   min.cells = 3)
  print(head(rownames(seurat_obj,10)))
  seurat_list <- append(seurat_list, seurat_obj)}

# 合并所有样本
seurat_combined <- merge(seurat_list[[1]],
                         y = seurat_list[-1],
                         add.cell.ids = samples)  # 使用 samples 中的完整名称

pbmc <- JoinLayers(seurat_combined)
# 替换 Seurat 对象中的行名，去掉 "-1"
pbmc@meta.data$orig.ident <- row.names(pbmc@meta.data)
pbmc@meta.data$orig.ident <-  gsub("-1$", "", pbmc@meta.data$orig.ident)
pbmc@meta.data$orig.ident <- gsub("^(.*)_(.*)_.*$", "\\1_\\2", pbmc@meta.data$orig.ident)

dim(pbmc)

table(pbmc@meta.data$orig.ident)
DefaultAssay(pbmc)
dim(pbmc)
#saveRDS(pbmc,"./pbmc_raw.rds")


############QC#####
#pbmc <- readRDS("pbmc_raw.rds")
# 计算线粒体基因比例
mito_genes <- rownames(pbmc)[grep("^MT-", rownames(pbmc))]
#线粒休基因比例
pbmc[["percent.mt"]]<- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 计算核糖体基因比例
ribo_genes <- rownames(pbmc)[grep("^RP[SL]", rownames(pbmc), ignore.case = TRUE)]
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# 计算红血细胞基因比例
hb_genes <- rownames(pbmc)[grep("^HB", rownames(pbmc))]
pbmc[["percent.hb"]]  <- PercentageFeatureSet(pbmc, pattern = "^HB")
head(pbmc@meta.data)

#设置质控标准#基因
minGene=300
maxGene=6500
minUMI=600
pctMT=20#线粒体
pctHB=0.2#血细胞
pctRIB=40#核糖体

#数据质控并绘制小提琴图
pbmc_filter <- subset(pbmc,subset = nFeature_RNA > minGene & nFeature_RNA <maxGene &
                        nCount_RNA > minUMI & percent.mt < pctMT & percent.hb < pctHB & percent.ribo < pctRIB)
table(pbmc_filter@meta.data$orig.ident)
table(pbmc@meta.data$orig.ident)
dim(pbmc_filter)

pbmc <- NormalizeData(object = pbmc_filter, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
VariableFeatures(object = pbmc)

pbmc = RunPCA(pbmc, pc.genes=VariableFeatures(object = pbmc))#, verbose=FALSE)
pcs <- 1:50


pbmc <- RunHarmony(pbmc, group.by.vars="orig.ident", assay.use="RNA", max_iter = 50, theta=3)
pbmc <- FindNeighbors(pbmc, dims = pcs, reduction="harmony", k.param=20) %>% FindClusters(resolution =0.1)
pbmc <- RunUMAP(pbmc, dims = pcs, reduction='harmony') %>%  RunTSNE(dims = pcs, reduction='harmony')

pbmc <- readRDS("/home/zhangxp/testt/pbmc_er_tn.rds")
DimPlot(object = pbmc, reduction = 'umap',raster=FALSE,
        group.by = 'celltype', pt.size = 1, label=FALSE)


##############1a
#######################
library(ggunchull)
plotData <- as.data.frame(pbmc1[["umap"]]@cell.embeddings)
plotData$celltype <-  pbmc1@meta.data$celltype
plotData$sample <-  pbmc1@meta.data$orig.ident


cell_colors <- c("#9c8a86","#c893c6","#a85343","#80a589","#e4db95","#6badd5","#ffcaff","#dccbd6")

plot5 = ggplot(plotData, aes(x = umap_1, y = umap_2, fill =celltype, color = celltype)) +
  stat_unchull(alpha = 0.25, size = 0.25, lty = 2, delta = 0.5) +
  geom_point(size = 0.01, show.legend = FALSE) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(hjust = 0.05, face = "italic", size = 14)
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_manual(values = cell_colors) +
  scale_color_manual(values = cell_colors)
plot5


##############1b

pbmc1 <- pbmc
pbmc1@meta.data$celltype[pbmc1@meta.data$celltype == "Tissue_stem_cells" ] <- "Fibroblasts"
pbmc1@meta.data$celltype[pbmc1@meta.data$celltype == "Monocyte" ] <- "Macrophage"
pbmc1@meta.data$celltype[pbmc1@meta.data$celltype == "MSC" ] <- "Epithelial_cells"

# 定义基因名称向量
gene_names <- c(
  "CD79A", "IGHG4", "IGKC", "JCHAIN", ###B CELL
  "CLDN5", "PECAM1", "VWF", ###Endothelial cell
  "CLDN4", "EPCAM", "KRT18", "KRT8", ####Epithelial cell
  "COL1A1",  "DCN",  "LUM", ###Fibroblast
  "AIF1", "CD163",  "CD68", ###Macrophage
  "CD7", "CD2", "CD3D"####T cell
)


DotPlot(pbmc1, features = gene_names, group.by = "celltype") + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colors = c("white", "darkred"))



##############1c

# 加载 ggplot2 包
library(ggplot2)

# 使用 table 函数计算细胞类型的频率
cell_type_freq <- table(pbmc1@meta.data$celltype)

# 将频率转换为 data.frame
cell_type_freq_df <- data.frame(celltype = names(cell_type_freq), freq = as.numeric(cell_type_freq))

# 指定颜色
color_vector <- c("#9c8a86","#c893c6","#a85343","#80a589","#e4db95","#6badd5","#ffcaff","#dccbd6")

# 画一个横向柱状图
ggplot(cell_type_freq_df, aes(x = freq, y = reorder(celltype, freq), fill = factor(celltype))) + 
  geom_col() + 
  geom_text(aes(label = freq), hjust = 0, size = 4, color = "black") + 
  scale_fill_manual(values = color_vector) + 
  labs(title = "Cell numbers", x = "counts", y = "", fill = "") + 
  theme(
    panel.grid.major = element_blank(),  # 去掉白色网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.background = element_blank(),  # 去掉灰色背景
    axis.line = element_line(colour = "black"),  # 坐标轴为黑色（非红色）
    axis.text.y = element_text(size = 10, hjust = 0),
    legend.position="none"
  )


####################D,E

pbmc <- readRDS("/home/zhangxp/testt/pbmc_er_tn.rds")
DimPlot(object = pbmc, reduction = 'umap',raster=FALSE,
        group.by = 'celltype', pt.size = 1, label=FALSE)

dim(pbmc)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(Seurat) 
human_countexp_Seurat <-sc.metabolism.Seurat(obj = pbmc,####### 固定pbmc
                                            method = "AUCell",
                                            imputation = F, metabolism.type = "KEGG")


print(DimPlot.metabolism(obj = human_countexp_Seurat, pathway = "Glycolysis / Gluconeogenesis ", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.05))

abc <- as.data.frame(human_countexp_Seurat@assays[["METABOLISM"]][["score"]]);abc <- t(abc)
abc <- as.data.frame(abc)
rownames(abc) <-  gsub(".1$", "-1", rownames(abc))

pbmc1 <- AddMetaData(pbmc, metadata = abc)
head(pbmc1@meta.data)


library(ggraph)
# Combine meta data with UMAP embeddings and ensure that the Glycolysis / Gluconeogenesis column is included
df_temp <- data.frame(pbmc1@meta.data, pbmc1@reductions$umap@cell.embeddings)

# Check if the 'Glycolysis / Gluconeogenesis' column is present
colnames(df_temp)[10] <- "Aucell"

# Plotting
ggplot(df_temp, aes(x = umap_1, y = umap_2, color = `Aucell`)) + 
  geom_point(size = 0.05) + 
  scale_color_viridis(option = "A") + 
  theme_light(base_size = 15) +
  labs(title = "Glycolysis / Gluconeogenesis") +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(hjust = 0.5))




# 将 MSC 改为 Epithelial_cells
df_temp$celltype[df_temp$celltype == "MSC"] <- "Epithelial_cells"
df_temp$celltype[df_temp$celltype == "Monocyte"] <- "Macrophage"
df_temp$celltype[df_temp$celltype == "Tissue_stem_cells" ] <- "Fibroblasts"

# 定义自定义的颜色
colors <- c("B_cell" = "#66c2a5", 
            "Endothelial_cells" = "#fc8d62", 
            "Epithelial_cells" = "#8da0cb", 
            "Fibroblasts" = "#e78ac3",
            "Macrophage" = "#83a0b6",
            "T_cells" = "#e09761")

library(ggsignif)
library(RColorBrewer);library(ggsignif)
compare_list <- list(
  c("Epithelial_cells", "B_cell"), c("Epithelial_cells", "Endothelial_cells"),c("Epithelial_cells", "Fibroblasts"),
  c("Epithelial_cells", "Macrophage"),c("Epithelial_cells", "T_cells"))


df_temp$group <- factor(df_temp$group, levels = c("B_cell", "Endothelial_cells", "Epithelial_cells", "Fibroblasts", "Macrophage", "T_cells"))

p <- ggplot(df_temp, aes(x = celltype, y = `Aucell`, fill = celltype)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "",
       x = "",
       y = "Glycolysis / Gluconeogenesis") +
  scale_fill_manual(values = colors) + 
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),  # Set axis text to black
        axis.title = element_text(color = "black"), # Set axis titles to black
        axis.ticks = element_line(color = "black"), # Set tick lines to black
        axis.line = element_line(color = "black"))+  # Set axis lines to black
  geom_signif(comparisons = compare_list,        # Define comparison groups
              map_signif_level = TRUE,          # Add significance levels
              textsize = 4,                     # Size of significance text
              test = "t.test",             # Use Wilcoxon test
              step_increase = 0.15)
p

###############f
pbmc_epi <- subset(pbmc1, celltype =='Epithelial_cells')

pbmc_epi <- NormalizeData(object = pbmc_epi, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_epi <- FindVariableFeatures(object = pbmc_epi, selection.method = "vst", nfeatures = 2000)
pbmc_epi <- ScaleData(pbmc_epi)
VariableFeatures(object = pbmc_epi)

pbmc_epi = RunPCA(pbmc_epi, pc.genes=VariableFeatures(object = pbmc_epi))#, verbose=FALSE)
pcs <- 1:30

pbmc_epi <- RunHarmony(pbmc_epi, group.by.vars="orig.ident", assay.use="RNA", max_iter = 50)#, theta=3)
pbmc_epi <- FindNeighbors(pbmc_epi, dims = pcs, reduction="harmony", k.param=20) %>% FindClusters(resolution =0.1)
pbmc_epi <- RunUMAP(pbmc_epi, dims = pcs, reduction='harmony') %>%  RunTSNE(dims = pcs, reduction='harmony')

#saveRDS(pbmc_epi,"copycakt2/pbmc_er_tn_epi.rds")

pbmc_epi <- readRDS("copycakt2/pbmc_er_tn_epi.rds")

library(RColorBrewer)
set3_colors <- brewer.pal(n = 12, name = "Set3")
print(DimPlot(pbmc_epi, reduction = "umap", label = F, raster = FALSE, cols = set3_colors))


##########g


samples <- list.files(path="/home/zhangxp/testt/copycakt2/", pattern = "*_copykat_prediction.txt", full.names = F)

first_file <- paste0("/home/zhangxp/testt/copycakt2/", samples[1])
common_df <- read.table(first_file, header = TRUE, sep = "\t")


for (i in 2:length(samples)) {
  file_path <- paste0("/home/zhangxp/testt/copycakt2/", samples[i])
  print(samples[i])
  current_df <- read.table(file_path, header = TRUE, sep = "\t")
  print(table(current_df$copykat.pred))
  common_df <- rbind(common_df,current_df)
}

malignant <- data.frame(copykat.pred = common_df$copykat.pred, row.names = common_df$cell.names)

pbmc_epi <- AddMetaData(pbmc_epi, metadata = malignant)
table(pbmc_epi$copykat.pred)

pbmc_epi$copykat.pred[pbmc_epi$copykat.pred == "c1:diploid:low.conf" | pbmc_epi$copykat.pred == "c2:aneuploid:low.conf"] <- "not.defined"

# 查看更新后的结果
table(pbmc_epi$copykat.pred)

##############H,I
library(AUCell) 
library(clusterProfiler)

exprMat <- pbmc_epi@assays$RNA@layers$data;
rownames(exprMat) <- rownames(pbmc_epi)
colnames(pbmc_epi) <- colnames(pbmc_epi)

cells_rankings <- AUCell_buildRankings(exprMat,nCores= 40) 
Hallmarker <- read.gmt("/home/zhangxp/R/x86_64-pc-linux-gnu-library/4.3/scMetabolism/data/KEGG_metabolism_nc.gmt") 
geneSets <- lapply(unique(Hallmarker$term), function(x){print(x);Hallmarker$gene[Hallmarker$term == x]})
names(geneSets) <- unique(Hallmarker$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)


##set gene set of interest here for plotting
geneSet <- "Glycolysis / Gluconeogenesis"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
pbmc_epi$AUC  <- aucs

library(ggraph)
ggplot(data.frame(pbmc_epi@meta.data, pbmc_epi@reductions$umap@cell.embeddings), aes(umap_1, umap_2, color=AUC)) + 
  geom_point( size=0.1) + scale_color_viridis(option="A")  + 
  theme_light(base_size = 15)+labs(title = "Glycolysis / Gluconeogenesis")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))






df_temp <- data.frame(pbmc_epi@meta.data);
df_temp <- subset(df_temp, df_temp$copykat.pred != "not.defined")


df_temp <- df_temp%>%
  mutate(group = str_split(orig.ident, "_", simplify = TRUE)[, 2]) 


df_temp$group <- paste(df_temp$copykat.pred, "_", df_temp$group)
table(df_temp$group)


# 定义自定义的颜色
colors <- c("aneuploid _ ER" = "#66c2a5", 
            "aneuploid _ TN" = "#fc8d62", 
            "diploid _ ER" = "#8da0cb", 
            "diploid _ TN" = "#e78ac3")

library(ggsignif)
library(RColorBrewer);library(ggsignif)
compare_list <- list(
  c("diploid _ ER", "aneuploid _ ER"), 
  c("diploid _ TN", "aneuploid _ TN"))


df_temp$group <- factor(df_temp$group, levels = c("diploid _ ER", "aneuploid _ ER",  "diploid _ TN", "aneuploid _ TN" ))

p <- ggplot(df_temp, aes(x = group, y = `AUC`, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "",
       x = "",
       y = "Glycolysis / Gluconeogenesis") +
  scale_fill_manual(values = colors) + 
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),  # Set axis text to black
        axis.title = element_text(color = "black"), # Set axis titles to black
        axis.ticks = element_line(color = "black"), # Set tick lines to black
        axis.line = element_line(color = "black"))+  # Set axis lines to black
  geom_signif(comparisons = compare_list,        # Define comparison groups
              map_signif_level = TRUE,          # Add significance levels
              textsize = 4,                     # Size of significance text
              test = "t.test",             # Use Wilcoxon test
              step_increase = 0.15)
p



# 提取 pbmc_epi 中的细胞名称
cells_epi <- colnames(pbmc_epi)

# 使用相同的细胞名过滤 pbmc 对象
pbmc_subset <- pbmc[, cells_epi]

# 检查维度是否匹配
dim(pbmc_subset)
dim(pbmc_epi)



#################J,K
###################################epi_diff################
meta_data <- pbmc_epi@meta.data
meta_data <- meta_data %>%
  mutate(minitype = str_extract(orig.ident, "(?<=_)[^_]+(?=_)"))
pbmc_epi@meta.data <- meta_data
pbmc_epi@meta.data$copy_type <- paste(pbmc_epi@meta.data$copykat.pred,"_",pbmc_epi@meta.data$minitype) 

head(pbmc_epi@meta.data)


Idents(pbmc_epi) <- pbmc_epi$copy_type
pbmc_epi.markers <- FindAllMarkers(pbmc_epi, only.pos = FALSE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.5)


# Identify differentially expressed genes between TN and N groups
TN_vs_N_markers <- FindMarkers(pbmc_epi, ident.1 = "aneuploid _ TN", ident.2 = "diploid _ TN",
                               only.pos = FALSE, 
                               min.pct = 0.5, 
                               logfc.threshold = 0)
#write.csv(TN_vs_N_markers,"TN_vs_N_markers.csv")

TN_vs_N_markers$cluster <- c(rep("TN_vs_N",nrow(TN_vs_N_markers)))
TN_vs_N_markers$gene <- rownames(TN_vs_N_markers)

ER_vs_N_markers <- FindMarkers(pbmc_epi, ident.1 = "aneuploid _ ER", ident.2 = "diploid _ ER",
                               only.pos = FALSE, 
                               min.pct = 0.5, 
                               logfc.threshold = 0)
#write.csv(ER_vs_N_markers,"ER_vs_N_markers.csv")
ER_vs_N_markers$cluster <- c(rep("ER_vs_N",nrow(ER_vs_N_markers)))
ER_vs_N_markers$gene <- rownames(ER_vs_N_markers)


pbmc_epi.markers1 <- rbind(TN_vs_N_markers,ER_vs_N_markers)
library(dplyr)

# 统计ER_vs_N中avg_log2FC > 0 和 avg_log2FC < 0 的个数
ER_vs_N_counts <- pbmc_epi.markers1 %>%
  filter(cluster == "ER_vs_N") %>%
  summarise(
    greater_than_0 = sum(avg_log2FC > 0),
    less_than_0 = sum(avg_log2FC < 0)
  )

# 统计TN_vs_N中avg_log2FC > 0 和 avg_log2FC < 0 的个数
TN_vs_N_counts <- pbmc_epi.markers1 %>%
  filter(cluster == "TN_vs_N") %>%
  summarise(
    greater_than_0 = sum(avg_log2FC > 0),
    less_than_0 = sum(avg_log2FC < 0)
  )

# 输出结果
print(ER_vs_N_counts)
print(TN_vs_N_counts)



library(scRNAtoolVis)
jjVolcano(diffData = pbmc_epi.markers1,
          tile.col = corrplot::COL2('RdBu', 15))[4:5])

#write.csv(pbmc_epi.markers, "pbmc.markers222333.csv",quote = F)


