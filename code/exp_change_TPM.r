
setwd("/Users/songdong/Dropbox/Dropbox/CSIC5011/single_Cell_SongD/")
library(ggplot2)
library(DESeq2)
library(stringr)
library("pheatmap")
library("RColorBrewer")
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library("SC3")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

#raw data preprocessing
Sample_features <- read.delim("./SampleInfo.txt",sep = "\t",header = T)

samples <- data.frame(Sample_features)
rownames(samples) <- Sample_features$ID


##############################################33
## create DESeq2 object
newExp <- read.csv("./markers_exprs.csv",header = T)

genes<-newExp$X
rownames(newExp) <- genes
newExp<-newExp[,2:2310]

## create DESeq2 object
pesudoT <- read.delim("./zongchao.txt",sep = "\t",header = T)


Sample_features <- read.delim("./meta.data.Serat.Cluster.txt",sep = "\t",header = T)

samples <- data.frame(Sample_features)




newInfo <- merge(samples,pesudoT,1,1)

genelist <- c("PAX6","SOX2","HOPX","EOMES","OLIG1","AQP4","GAD1","NEUROD2","PDGFRA","SFRP1","PTPRC","P2RY12",
              "COL20A1","PMP2","RBFOX1","GAD1","PDE4DIP","GFAP","SLCO1C1")
selectedG <- t(as.data.frame(newExp[(rownames(newExp) %in% genelist),]))
tmp1 <- data.frame(cbind(rownames(selectedG),selectedG))

myTable <- merge(newInfo,tmp1,1,1)
summary(log10(as.numeric(as.character(myTable$OLIG1))+1))


cleanTable <- myTable[which(myTable$traj.coord != "Inf" ),]
##for all cells
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(OLIG1))+1)), se = TRUE,span = 0.5,color="black")
new.plot <- new.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(OLIG1))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of OLIG1')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(OLIG1))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(OLIG1))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of OLIG1')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="OLIG1_exp.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)



######PAX6
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PAX6))+1)), se = TRUE,span = 0.5,color="black")
new.plot <- new.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PAX6))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of PAX6')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PAX6))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PAX6))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of PAX6')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="PAX6_exp.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)



######PDGFRA
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PDGFRA))+1)), se = TRUE,span = 0.5,color="black")
new.plot <- new.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PDGFRA))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of PDGFRA')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PDGFRA))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PDGFRA))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of PDGFRA')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="PDGFRA_exp.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)

######GFAP
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(GFAP))+1)), se = TRUE,span = 0.5,color="black")
new.plot <- new.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(GFAP))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of GFAP')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(GFAP))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(GFAP))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of GFAP')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="GFAP_exp.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)


######GAD1
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(GAD1))+1)), se = TRUE,span = 0.5,color="black")
new.plot <- new.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(GAD1))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of GAD1')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(GAD1))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(GAD1))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of GAD1')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="GAD1_exp.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)




######AQP4
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(AQP4))+1)), se = TRUE,span = 0.5,color="black")
new.plot <- new.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(AQP4))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of AQP4')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(AQP4))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(AQP4))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of GAD1')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="AQP4_exp.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)




######NEUROD2
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1)), se = TRUE,span = 0.5,color="black")
new.plot <- new.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of NEUROD2')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of NEUROD2')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="NEUROD2_exp.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)




tableNeuron <-cleanTable[which(cleanTable$CellType == "NPCs" | cleanTable$CellType == "ExcitatoryNeurons"),]
tableGlial <- cleanTable[which(cleanTable$CellType == "OPCs" | cleanTable$CellType == "Astrocytes"),]
######NEUROD2
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = tableNeuron, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1)), se = F,span = 5,color=gg_color_hue(6)[6])
new.plot <- new.plot+   geom_smooth(data = tableGlial, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1)), se = F,span = 5,color=gg_color_hue(2)[2])
new.plot <- new.plot+  geom_point(data = tableNeuron, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=23) 
new.plot <- new.plot+  geom_point(data = tableGlial, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=22) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,2,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of NEUROD2')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of NEUROD2')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="NEUROD2_exp2.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)




tableNeuron <-cleanTable[which(cleanTable$CellType == "Interneurons" | cleanTable$CellType == "ExcitatoryNeurons"),]
tableGlial <- cleanTable[which(cleanTable$CellType == "OPCs" | cleanTable$CellType == "Astrocytes"),]
######PAX6
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = tableNeuron, aes(x = traj.coord, y = log10(as.numeric(as.character(HOPX))+1)), se = F,span = 5,color=gg_color_hue(2)[1])
new.plot <- new.plot+   geom_smooth(data = tableGlial, aes(x = traj.coord, y = log10(as.numeric(as.character(HOPX))+1)), se = F,span = 5,color=gg_color_hue(2)[2])
new.plot <- new.plot+  geom_point(data = tableNeuron, aes(x = traj.coord, y = log10(as.numeric(as.character(HOPX))+1),fill='Neuron'),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
new.plot <- new.plot+  geom_point(data = tableGlial, aes(x = traj.coord, y = log10(as.numeric(as.character(HOPX))+1),fill='Glial'),alpha=0.9,stroke=0.1,size=2.5,shape=22) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,2,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of PAX6')+scale_fill_manual(name=NULL,values=gg_color_hue(2))#values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(HOPX))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(HOPX))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of HOPX')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="HOPX_exp2.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)



tableNeuron <-cleanTable[which(cleanTable$CellType == "Interneurons" | cleanTable$CellType == "ExcitatoryNeurons"),]
tableGlial <- cleanTable[which(cleanTable$CellType == "OPCs" | cleanTable$CellType == "Astrocytes"),]
######PAX6
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = tableNeuron, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1)), se = F,span = 5,color=gg_color_hue(2)[1])
new.plot <- new.plot+   geom_smooth(data = tableGlial, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1)), se = F,span = 5,color=gg_color_hue(2)[2])
new.plot <- new.plot+  geom_point(data = tableNeuron, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1),fill='Neuron'),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
new.plot <- new.plot+  geom_point(data = tableGlial, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1),fill='Glial'),alpha=0.9,stroke=0.1,size=2.5,shape=22) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,2,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of PAX6')+scale_fill_manual(name=NULL,values=gg_color_hue(2))#values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of NEUROD2')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))


new2.plot<-ggplot() 
new2.plot <- new2.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1)), se = TRUE,span = 0.5,color="black")
new2.plot <- new2.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(NEUROD2))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new2.plot<- new2.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new2.plot<- new2.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of NEUROD2')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new2.plot<- new2.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new2.plot




plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new2.plot),ggplotGrob(new.plot),size="first")

ggsave(file="NEUROD2_exp3.pdf", plot=plot7,bg = 'white', width = 43, height = 15, units = 'cm', dpi = 600)









######PTPRC
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PTPRC))+1)), se = TRUE,span = 0.5,color="black")
new.plot <- new.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PTPRC))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of PTPRC')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PTPRC))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PTPRC))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of PTPRC')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="PTPRC_exp.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)


######SOX2
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(SOX2))+1)), se = TRUE,span = 0.5,color="black")
new.plot <- new.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(SOX2))+1),fill=CellType),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,0.2),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_blank())

new.plot<- new.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of SOX2')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

##for all cells
new1.plot<-ggplot() 
new1.plot <- new1.plot+   geom_smooth(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(SOX2))+1)), se = TRUE,span = 0.5,color="black")
new1.plot <- new1.plot+  geom_point(data = cleanTable, aes(x = traj.coord, y = log10(as.numeric(as.character(SOX2))+1),fill=week),alpha=0.9,stroke=0.1,size=2.5,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new1.plot<- new1.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                            legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new1.plot<- new1.plot+ggtitle(NULL)+xlab("Pseudo-time")+ylab('Expression of SOX2')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new1.plot<- new1.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.3),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(-0.5,27),breaks=seq(0,50,5))
new1.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new1.plot),ggplotGrob(new.plot),size="first")

ggsave(file="SOX2_exp.pdf", plot=plot7,bg = 'white', width = 30, height = 15, units = 'cm', dpi = 600)



#Astrocytes="Astrocytes" ,Excitatory Neurons="Excitatory Neurons" ,Internuerons ="Internuerons,Microglia="Microglia" ,NCPs="NPCs" ,OPCs="OPCs")



















##for NPCs
table1 <- myTable[which(myTable$CellType=="NCPs"),]
table2 <- myTable[which(myTable$CellType=="OPCs"),]
table3 <- myTable[which(myTable$CellType=="Internuerons"),]
table4 <- myTable[which(myTable$CellType=="Microglia"),]
table5 <- myTable[which(myTable$CellType=="Excitatory Neurons"),]
table6 <- myTable[which(myTable$CellType=="Astrocytes"),]
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = myTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PAX6))+1)), se = TRUE,span = 0.1,color="black")
new.plot <- new.plot+  geom_point(data = myTable, aes(x = traj.coord, y = log10(as.numeric(as.character(PAX6))+1),fill=CellType),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of PAX6')#+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_PAX6.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


#GAD1






###OLIG1
new.plot <-ggplot()
new.plot <- new.plot+   geom_smooth(data = table5, aes(x = traj.coord, y = as.numeric(as.character(OLIG1))), se = TRUE,span = 0.1,color="black")
new.plot <- new.plot+  geom_point(data = table5, aes(x = traj.coord, y = as.numeric(as.character(OLIG1)),fill=week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of PAX6')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,20),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_OLIG1.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


###OLIG1
new.plot <-ggplot()
new.plot <- new.plot+   geom_smooth(data = table4, aes(x = traj.coord, y = 2^as.numeric(as.character(OLIG1))), se = TRUE,span = 0.1,color="black")
new.plot <- new.plot+  geom_point(data = table4, aes(x = traj.coord, y = 2^as.numeric(as.character(OLIG1)),fill=week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of PAX6')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,20),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_OLIG1_2.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)

###PTPRC
new.plot <-ggplot()
new.plot <- new.plot+   geom_smooth(data = table4, aes(x = traj.coord, y = log10(as.numeric(as.character(PTPRC))+1)), se = TRUE,span = 0.1,color="black")
new.plot <- new.plot+  geom_point(data = table4, aes(x = traj.coord, y = log10(as.numeric(as.character(PTPRC))+1),fill=week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of PTPRC')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,20),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_PTPRC.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)





##############################################################################################################################
#________________DE_analysis_____________#

library(tibble)
library(dplyr)
library(tidyr)

newExp$rowSD <- rowVars(log10(as.matrix(newExp)+1))
ranks <-  newExp$rowSD
names(ranks) <- rownames(newExp)
#names(ranks) <- filtterData$geneName
#colnames(ranks) <- c("gene","rank")
ranks <- na.omit(ranks)

pathways.hallmark <- gmtPathways("/Users/songdong/Dropbox/Dropbox/germ_cell_tumor/analysis/expression/lib/msigdb_v7.2_GMTs/h.all.v7.2.symbols.gmt")
head(pathways.hallmark)

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) # order by normalized enrichment score (NES)
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "sig", "c_ns")
fgseaResTidy$adjPvalue[fgseaResTidy$NES>0 & fgseaResTidy$padj <= 0.05] <- "a_up"
fgseaResTidy$adjPvalue[fgseaResTidy$NES<0 & fgseaResTidy$padj <= 0.05] <- "b_down"
fgseaResTidy$pathway <- tolower(fgseaResTidy$pathway)

plot1 <- ggplot() +
  geom_col(data=fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue),alpha=0.7) +   
  coord_flip() +  labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways Enrichment Score from GSEA")+
  scale_fill_manual(NULL,values = c(non_significant = "grey", High_variance = "red",Low_variance = "blue")) + theme_classic() +
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(2,1,0.5,1),'lines'),plot.title=element_text(size=20,vjust=0.5,hjust=0.5,face='plain'),
        text=element_text(size=14,face='bold'),legend.key.width=unit(1,'cm'),legend.key.height=unit(1,'cm'),legend.position='right',
        legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=14,face='plain'),axis.text.y=element_text(size=12,face='plain',color='black'),
        axis.text.x=element_text(size=12,face='plain',color='black'),axis.title.x=element_text(size=16,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

figure_1<-rbind(ggplotGrob(plot1),size="last")

ggsave(file="Hallmarks1_rowSD.pdf", plot=figure_1,bg = 'white', width = 30, height = 36, units = 'cm', dpi = 600)


plotEnrichment(pathway = pathways.hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, ranks)

fgseaRes_pos <- fgseaRes[(fgseaRes$NES>0 & fgseaRes$padj < 0.01),]
temp <- order(-fgseaRes_pos$NES,fgseaRes_pos$padj)
fgseaRes_pos <- fgseaRes_pos[temp,]
plotGseaTable(pathways.hallmark[fgseaRes_pos$pathway[fgseaRes_pos$padj < 0.01]], ranks, fgseaRes_pos,  gseaParam=0.5)

fgseaRes_neg <- fgseaRes[(fgseaRes$NES<0 & fgseaRes$padj < 0.01),]
temp <- order(fgseaRes_neg$NES,fgseaRes_neg$padj)
fgseaRes_neg <- fgseaRes_neg[temp,]
plotGseaTable(pathways.hallmark[fgseaRes_neg$pathway[fgseaRes_neg$padj < 0.01]], ranks, fgseaRes_neg,  gseaParam=0.5)
















###############################################################################################################################





      

paper <- read.delim("./Pseudotime_paper.txt",sep = "\t",header = T)

colnames(paper) <- c("ID","PseudoTime","State","CellGroup","Week","RealTime")

new.plot <-ggplot()
new.plot <- new.plot+   geom_smooth(data = paper, aes(x = RealTime, y = as.numeric(as.character(PseudoTime))), se = TRUE,span = 0.1,color="black")
new.plot <- new.plot+  geom_point(data = paper, aes(x = RealTime, y = as.numeric(as.character(PseudoTime)),fill=Week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("Week")+ylab('PseudoTime')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,20),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_time.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


our <- data.frame(myTable$ID,myTable$traj.coord)
colnames(our) <- c("ID","Monocle")
compare <- merge(paper,our,1,1)
no_Inf <- compare[which(compare$Monocle!="Inf"),]
new.plot <-ggplot()
new.plot <- new.plot+   geom_smooth(data = no_Infe, aes(x = RealTime, y = as.numeric(as.character(Monocle))), se = TRUE,span = 0.1,color="black")
new.plot <- new.plot+  geom_point(data = no_Inf, aes(x = RealTime, y = as.numeric(as.character(Monocle)),fill=Week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("Week")+ylab('PseudoTime')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,20),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_time_our.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)



new.plot <-ggplot()
new.plot <- new.plot+   geom_smooth(data = no_Inf, aes(x = Monocle, y = as.numeric(as.character(PseudoTime))), se = TRUE,span = 0.1,color="black")
new.plot <- new.plot+  geom_point(data = no_Inf, aes(x = Monocle, y = as.numeric(as.character(PseudoTime)),fill=Week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("Monocle_PseudoTime")+ylab('PseudoTime predicted in paper')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,20),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_time_vs.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)




theTable <- merge(no_Inf,tmp1,1,1)
##for NPCs
tableA <- theTable[which(theTable$CellGroup=="Astrocytes"),]
tableB <- theTable[which(theTable$CellGroup=="ExcitatoryNeurons"),]
tableC <- theTable[which(theTable$CellGroup=="NPCs"),]
tableD <- theTable[which(theTable$CellGroup=="OPCs"),]
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(PAX6))+1)), se = TRUE,span = 3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(PAX6))+1),fill=CellGroup),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of PAX6')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_PAX6_allCells.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(PAX6))+1)), se = TRUE,span = 3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(PAX6))+1),fill=Week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of PAX6')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))

#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_PAX6_allweeks.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)




new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = tableC, aes(x = Monocle, y = log10(as.numeric(as.character(PAX6))+1)), se = TRUE,span = 2,color="black")
new.plot <- new.plot+  geom_point(data = tableC, aes(x = Monocle, y = log10(as.numeric(as.character(PAX6))+1),fill=Week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of PAX6')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,23),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_PAX6_our_NPCs_week.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)

###PTPRC
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(PTPRC))+1)), se = TRUE,span = 3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(PTPRC))+1),fill=CellGroup),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of PTPRC')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_PTPRC_allCells.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(PTPRC))+1)), se = TRUE,span = 3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(PTPRC))+1),fill=Week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of PTPRC')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))

#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_PTPRC_allweeks.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


###P2RY12
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(P2RY12))+1)), se = TRUE,span = 3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(P2RY12))+1),fill=CellGroup),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of P2RY12')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_P2RY12_allCells.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(P2RY12))+1)), se = TRUE,span = 3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(P2RY12))+1),fill=Week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of P2RY12')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))

#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_P2RY12_allweeks.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)



###OLIG1
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1)), se = TRUE,span = 3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1),fill=CellGroup),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of OLIG1')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_OLIG1_allCells.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1)), se = TRUE,span = 3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1),fill=Week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of OLIG1')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))

#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_OLIG1_allweeks.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


###OLIG1
new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1)), se = TRUE,span = 0.3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1),fill=CellGroup),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of OLIG1')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_OLIG1_allCells.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)

temp<- theTable[which(theTable$CellGroup!="ExcitatoryNeurons" & theTable$CellGroup!="NPCs"),]

new.plot<-ggplot() 
#new.plot <- new.plot+   geom_smooth(data = temp, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1)), se = TRUE,span = 0.3,color="black")
new.plot <- new.plot+  geom_point(data = temp, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1),fill=CellGroup),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of OLIG1')+scale_fill_manual(name=NULL,values=c(Interneurons="#BB9F00",ExcitatoryNeurons="#FF6C67",Microglia="#FF52E9",OPCs="#509BFF",Astrocytes="#00C2C6",NPCs="#00BF0D"))#,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))
new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.1,4.2),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_OLIG1_glial.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)


new.plot<-ggplot() 
new.plot <- new.plot+   geom_smooth(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1)), se = TRUE,span = 3,color="black")
new.plot <- new.plot+  geom_point(data = theTable, aes(x = Monocle, y = log10(as.numeric(as.character(OLIG1))+1),fill=Week),alpha=0.7,size=3,shape=21) 
#new.plot <- new.plot+  geom_point(aes(x=tmpsyn2$Age[which(tmpsyn2$ID == 'F007')],y=tmpsyn2$MB[which(tmpsyn2$ID == 'F007')]),color=gg_color_hue(3)[3],alpha=0.7,size=4) 
new.plot<- new.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
                          legend.position='top',legend.text=element_text(size=16,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),
                          axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                          axis.title.x=element_text(size=18,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=18,face='plain',color='black'))

new.plot<- new.plot+ggtitle(NULL)+xlab("PesudoTime")+ylab('Expression of OLIG1')+scale_fill_manual(name=NULL,values=c(GW08="#f7fcfd",GW09="#e0ecf4",GW10="#bfd3e6",GW12="#9ebcda" ,GW13 ="#8c96c6",GW16="#8c6bb1" ,GW19="#88419d" ,GW23="#810f7c" ,GW26="#4d004b"))

#new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits = c(-0.5,4),breaks=seq(0,500,2))+scale_x_continuous(expand=c(0,0),limits = c(0,30),breaks=seq(0,50,5))
new.plot
#Fra.plot<-Fra.plot+scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10))+scale_x_discrete(labels=c(Germline='Germline\nmutation',Somatic='Somatic\nmutation',F007='Patient\nF007',Others='Other\npatients'))

plot7<-cbind(ggplotGrob(new.plot),size="first")

ggsave(file="regression_OLIG1_allweeks.pdf", plot=plot7,bg = 'white', width = 18, height = 15, units = 'cm', dpi = 600)

