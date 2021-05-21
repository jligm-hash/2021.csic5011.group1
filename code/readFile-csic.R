
library(data.table)
readTab = fread('~/Documents/1_course/CSIC5011/final_projects/GSE104276_all_pfc_2394_UMI_TPM_NOERCC.xls',
               sep = '\t')


readDf = readTab
row.names(readDf) = readTab$V1
readDf$V1 = NULL
library(pheatmap)
# pheatmap(readd)



splitFirstOne = function(x){
  return(strsplit(x, "_")[[1]][1])
}

labelDf = data.frame("colName" = colnames(readDf))
labelDf$GW = apply(labelDf[1], 1, splitFirstOne)

# write.csv(labelDf,'~/Documents/1_course/CSIC5011/final_projects/label.csv')



library(slingshot, quietly = FALSE)
countMat = readDf[,]
dim(countMat)
# deng_SCE = SingleCellExperiment(list(counts=readDf))
sce = SingleCellExperiment(list(counts=countMat))

pca_data <- prcomp(t(countMat))
library(gmodels)
pca_data <- fast.prcomp(t(countMat))

# library(Rtsne)
# set.seed(5252)
# tsne_data <- Rtsne(pca_data$x[,1:50], pca = FALSE)

reducedDims(sce) <- list(PCA=pca_data$x[,1:10])
sce <- slingshot(sce, reducedDim = "PCA")
t = slingPseudotime(sce)

library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
slingshot(sce, reducedDim = 'PCA')

summary(sce$slingPseudotime_1)
write.csv(t, '~/Documents/1_course/CSIC5011/final_projects/time2.csv', quote = F)

reducedDims(deng_SCE)
reducedDim(deng_SCE, "PCA")

table(deng_SCE$cell_type1)
t = slingshot(sce, reducedDim = "PCA")

example("SingleCellExperiment")
reducedDim(sce, "PCA")


cellTime = fread('/Users/jiabaoli/Documents/1_course/CSIC5011/final_projects/label.csv')
timeSat = data.frame(table(cellTime$week))
colnames(timeSat) = c('GW', 'cellNumber')
ggplot(timeSat, mapping = aes(x=reorder(as.factor(GW),-cellNumber), 
                              y=cellNumber,fill=GW)) +
  geom_bar(stat = 'identity') +
  theme_test() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,860)) +
  geom_text(aes(label=cellNumber), vjust=-1) +
  labs( y = '# of cells',
        x = 'Gestational weeks') +
  scale_fill_brewer(palette="Set3") +
  guides(fill = F)+
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12)) 

summaryReadCounts = apply(readDf, 2, summary)
summaryReadCounts = data.frame(summaryReadCounts)

rownames(readDf)
apply(readDf, 1, mean)

geneSta = data.frame('gene' = rownames(readDf),
                     'averageTPM' = apply(readDf, 1, mean))

geneStaInp = head(geneSta[order(geneSta$averageTPM, decreasing=TRUE), ], 30)

# ----- This section prepare a dataframe for labels ---- #
# Get the name and the y position of each label
label_data <- geneStaInp
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
label_data$id = 1:number_of_bar
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar    
# I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)
# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)

library("ggsci")
library(viridis)
library(pals)
library(jcolors)
ggplot(geneStaInp, mapping = aes(x = reorder(gene, -averageTPM), 
                                 y = as.factor(averageTPM),
                                 fill=reorder(gene, -averageTPM))) +
  geom_bar(stat="identity") + #alpha("blue", 0.3)
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
    # plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) +
  guides(fill=F)  +
  # scale_fill_viridis(discrete = T, direction =-1) +
  # scale_fill_brewer(palette = 'pal10') +
  # scale_fill_jcolors(palette = "pal12") +
  labs(title = 'Top 30 genes in the dataset') +
  geom_text(data=label_data, aes(x=id, y=as.factor(averageTPM), label=gene, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0)

