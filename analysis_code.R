###################
# Terminal
# ################# Prep whitelist
# umi_tools whitelist --stdin rawdata/AK-WTA_S14_R1_001.fastq.gz \
# --method=reads \
# --extract-method=regex \
# --bc-pattern="(?P<cell_1>.{9})(?P<discard_1>.{12})(?P<cell_2>.{9})(?P<discard_2>.{13})(?P<cell_3>.{9})(?P<umi_1>.{8}).*" \
# --set-cell-number=10000 \
# --ed-above-threshold=correct \
# --error-correct-cell  \
# --log2stderr> whitelist.txt
# 
# whitelist = read.table("whitelist.txt",sep = "\t",header = F)
# plot(sort(log10(whitelist$V3),decreasing = T))
# 
# ################# Add cell barcode and UMI to read ID
# umi_tools extract  --stdin rawdata/AK-WTA_S14_R1_001.fastq.gz \
# --read2-in rawdata/AK-WTA_S14_R2_001.fastq.gz \
# --read2-stdout  \
# --extract-method=regex \
# --bc-pattern="(?P<cell_1>.{9})(?P<discard_1>.{12})(?P<cell_2>.{9})(?P<discard_2>.{13})(?P<cell_3>.{9})(?P<umi_1>.{8}).*" \
# --error-correct-cell \
# --whitelist=whitelist.txt \
# --log2stderr \
# > umi_tools/AK-WTA_S14_R2_001extracted.fastq.gz
# 
# ################# Trimming read2
# fastp --trim_poly_x -w 20 -o trimmed/AK-WTA_S14_R2_001extracted_fastp.fastq.gz  --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -i umi_tools/AK-WTA_S14_R2_001extracted.fastq.gz -l 31 -h trimmed/fastp.html
# 
# ################# Mapping
# bwa_mem.rb -i list -d /mnt/hdd/hdd/data/species/D.japonica/DjTrascriptome/DjTrascriptome.fasta -a 64

# ################# Remove redundancy
# umi_tools dedup -I bwa_mem/AK-WTA_S14_R2_001extracted_fastp.fastq.gz --out-sam  --log2stderr --per-contig --per-gene --umi-separator="_" -S bwa_mem/AK-WTA_S14_R2_001extracted.uniq.sam
# ################# Split sam into each cell sam
# sam_demutiplexing.rb -i bwa_mem/AK-WTA_S14_R2_001extracted.uniq.sam
# ################# Convert sam to bam
# sam2bam.rb -i list -f /mnt/hdd/hdd/data/species/D.japonica/DjTrascriptome/DjTrascriptome.fasta

# ################# Quantification
# salmon.rb -i list -a 60 -f /mnt/hdd/hdd/data/species/D.japonica/DjTrascriptome/DjTrascriptome.fasta

###################
# R
ggColorHue <- function(n, l=65) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}
################# Load count data
dirs = list.dirs("salmon/",recursive = F) 
i = 1
tmp = read.table(sprintf("%s/quant.sf",dirs[i]),header = T)
rawcnt = matrix(0,ncol = length(dirs),nrow = nrow(tmp))
rawcnt[,i] = tmp$NumReads
for(i in 2:length(dirs)){
  cat(i,"\n")
  tmp = read.table(sprintf("%s/quant.sf",dirs[i]),header = T)
  rawcnt[,i] = tmp$NumReads
}
rownames(rawcnt) = tmp$Name
colnames(rawcnt) = gsub("salmon//","", dirs)

library(Seurat)
load("CLS-barcodes")
barcodes = barcodes[is.element(barcodes,colnames(rawcnt))]
seurat  = CreateSeuratObject(rawcnt[,barcodes],project = "intact-TAS")
seurat = NormalizeData(seurat)
plot(sort(log10(seurat@meta.data$nCount_RNA),decreasing = T))
abline(h = log10(200))
seurat = subset(seurat,subset = nCount_RNA >= 200)
median(seurat.intactTAS$nCount_RNA)

# Iterative Clustering
Idents(seurat) = factor(rep("0", length(seurat@active.ident)), levels = c("0"))
system("rm -r clustering")
dir.create("clustering")
markers.2nd = NULL
min.cell.in.cluster = 10
min.gene.in.cluster = 5
start = 0.1
step = 0.1
th = 1.5
library(ggplot2)
set.seed(1)
fin = NULL
for(i in 1:100){
  seurat.2 = seurat
  ok = 0
  pdf = 0
  cat(max(table(seurat.2@active.ident)),":")
  for(c in unique(seurat.2@active.ident)){
    # for(c in names(table(seurat.2@active.ident))[table(seurat.2@active.ident)==max(table(seurat.2@active.ident))]){
    if(sum(!is.element(c,fin))>0){
      r = start
      round.1 = 0
      cat(i,"\t",c,"\t","Round.1","\t")
      while(round.1 == 0 & r<=20){
        seurat.43 = subset(seurat.2, idents = c(as.character(c)))
        seurat.43 = ScaleData(seurat.43, verbose = F,model.use = "poisson")
        seurat.43 = FindVariableFeatures(seurat.43,verbose = F)
        if((length(seurat.43@active.ident)-1)>50){
          npc = 50
        } else {
          npc = length(seurat.43@active.ident)-1
        }
        seurat.43 = suppressWarnings(RunPCA(seurat.43,npcs = npc,verbose = F))
        seurat.43 = FindNeighbors(seurat.43,verbose = F,dims = 1:npc)
        seurat.43 = try(FindClusters(seurat.43,resolution = r,verbose = F),silent = T)
        if (class(seurat.43) != "try-error") {
          num.cluster = length(unique(seurat.43@active.ident))
          if(num.cluster>1){
            round.1 = 1
            cat(r,"\t")
          } else {
            r = r+1
          }
        }else {
          r = r+1
        }
      }
      cat("\t","Round.3","\t")
      round.3 = 0
      if(r == 20.1){
        round.3 = 1
        fin = c(fin,c)
      }
      if(r>1){
        r = r-1
      }
      while (round.3 == 0 & r<=20.1) {
        if(r >= 20){
          fin = c(fin,c)
        }
        cat(r,"\t")
        seurat.43 = subset(seurat.2, idents = c(as.character(c)))
        seurat.43 = ScaleData(seurat.43, verbose = F,model.use = 'poisson')
        seurat.43 = FindVariableFeatures(seurat.43,verbose = F)
        if((length(seurat.43@active.ident)-1)>50){
          npc = 50
        } else {
          npc = length(seurat.43@active.ident)-1
        }
        seurat.43 = suppressWarnings(RunPCA(seurat.43,npcs = npc,verbose = F))
        seurat.43 = FindNeighbors(seurat.43,verbose = F,dims = 1:npc)
        seurat.43 = try(FindClusters(seurat.43,resolution = r,verbose = F,group.singletons = F),silent = T)
        if (class(seurat.43) != "try-error" & sum(seurat.43@active.ident=="singleton")==0) {
          num.cluster = length(unique(seurat.43@active.ident))
          min.cluster = min(table(seurat.43@active.ident))
          if(num.cluster>1){
            markers.43  = FindAllMarkers(seurat.43,verbose = F, logfc.threshold = th,only.pos = T,min.pct = 0.5)
            num.marker = table(markers.43$cluster)
            pass = names(num.marker)[num.marker>=min.gene.in.cluster]
            pass = pass[table(seurat.43@active.ident)[pass]>=min.cell.in.cluster]
            ng = unique(seurat.43@active.ident)[!is.element(unique(seurat.43@active.ident),pass)]
            if(length(pass)>0){
              for(c2 in ng){
                cells = names(seurat.43@active.ident)[seurat.43@active.ident==as.character(c2)]
                Idents(seurat.43,cells) <- min(as.numeric(as.character(ng)))
              }
              markers.43  = FindAllMarkers(seurat.43,verbose = F, logfc.threshold = th,only.pos = T,min.pct = 0.5)
              num.marker = table(markers.43$cluster)
              pass = names(num.marker)[num.marker>=min.gene.in.cluster]
              pass = pass[table(seurat.43@active.ident)[pass]>=min.cell.in.cluster]
              pass = c(pass, min(as.numeric(as.character(unique(seurat.43@active.ident)))))
              ng = as.character(unique(seurat.43@active.ident))[!is.element(as.character(unique(seurat.43@active.ident)), pass)]
              if(length(ng)==0){
                ok = ok + 1
                cat(num.cluster)
                markers.43 = cbind(markers.43, data.frame("parent-cluster" = as.character(c), "round" = i,"r" = r))
                markers.43 = markers.43[is.element(markers.43$cluster,pass),]
                markers.2nd = rbind(markers.2nd, markers.43)
                g = DoHeatmap(seurat.43,features = unique(markers.43$gene))+ggtitle(c)+scale_fill_gradientn(colors = c("black", "yellow"))
                if(pdf==0){
                  pdf(sprintf("clustering/%s-clustering.pdf",i),width = 14)
                  pdf = 1
                }
                print(g)
                for(c2 in unique(as.character(seurat.43@active.ident))){
                  cells = names(seurat.43@active.ident)[seurat.43@active.ident==as.character(c2)]
                  Idents(seurat.2,cells) <- sprintf("%s:%s",c,c2)
                }
                round.3 = 1
              } else {
                r = r+0.1
              }
            } else {
              r = r+0.1
            } 
          } else {
            r = r+0.1
          }
        }else {
          r = r+0.1
        }
      }
      cat("\n")
    }
  }
  if(pdf == 1){
    dev.off()
  }
  seurat = seurat.2
  
  # save(seurat, file = sprintf("clustering/seurat-%s",i))8773
}


markers.intactTAS = markers.2nd
seurat.intactTAS = seurat
ctl = sort(as.character(unique(seurat.intactTAS@active.ident)))
seurat.intactTAS@active.ident = factor(seurat.intactTAS@active.ident,levels = ctl)
seurat.intactTAS = ScaleData(seurat.intactTAS,model.use = 'poisson',features = rownames(seurat.intactTAS@assays$RNA@data))

# Convert Cell type ID and select markers
tmp2 = NULL
c = 0
for(cell in levels(seurat.intactTAS@active.ident)){
  cat(cell, "\n")
  c = c + 1
  lineage = unlist(strsplit(cell,":"))
  for (i in 2:length(lineage)) {
    parent = lineage[i-1]
    cluster = as.numeric(lineage[i])
    tmp = markers.intactTAS[markers.intactTAS$parent.cluster==parent & markers.intactTAS$cluster==cluster,]
    if(nrow(tmp)>0){
      out = data.frame("Original.Cell.Type" = cell, "Cell.Type" = c, "gene" = tmp$gene, "is.Positive" = T)
      tmp2 = rbind(tmp2,out)
    }
    tmp = markers.intactTAS[markers.intactTAS$parent.cluster==parent & markers.intactTAS$cluster!=cluster,]
    if(nrow(tmp)>0){
      out = data.frame("Original.Cell.Type" = cell, "Cell.Type" = c, "gene" = tmp$gene, "is.Positive" = F)
      tmp2 = rbind(tmp2,out)
    }
  }
}
Idents(seurat.intactTAS) = factor(as.numeric(seurat.intactTAS@active.ident), levels = 1:length(unique(seurat.intactTAS@active.ident))) 
markers.intactTAS.original = markers.intactTAS
markers.intactTAS = FindAllMarkers(seurat.intactTAS,logfc.threshold = 1,only.pos = T,min.pct = 0.5)
save(seurat.intactTAS,markers.intactTAS,markers.intactTAS.original,file = "seurat_intactTAS_200")
length(unique(markers.intactTAS$cluster))
# one-step clustering performance
seurat = seurat.intactTAS
seurat = FindVariableFeatures(seurat)
seurat = RunPCA(seurat)
seurat = FindNeighbors(seurat)
all.pre = 1
res.one.clustering = NULL
for(r in seq(0.1,20,0.1)){
  cat(r,"\t")
  seurat =FindClusters(seurat, resolution = r, group.singletons = F, verbose = F)
  all = length(unique(seurat@active.ident))
  if(all > all.pre){
    m = FindAllMarkers(seurat,logfc.threshold = 1,only.pos = T,min.pct = 0.5,verbose = F)
    ID = length(unique(m$cluster))
    all.pre = all
    res.one.clustering = rbind(res.one.clustering, data.frame("r" = r, "all" = all, "IDed" = ID))
    cat(all,"\t",ID)
  } 
  cat("\n")
}

res.one.clustering = cbind(res.one.clustering, data.frame("noID" = res.one.clustering$all - res.one.clustering$IDed))

library(reshape2)
tmp = melt(res.one.clustering[res.one.clustering$all<100,-1],id.vars = "all")
tmp$variable = factor(tmp$variable, levels = c("noID","IDed"))
g1 = ggplot(tmp, aes(x = all,y = value, fill = variable))+
  geom_bar(stat = "identity")+
  scale_fill_grey()+
  scale_x_continuous(breaks = seq(0,100,5))+
  scale_y_continuous(breaks = seq(0,100,5),limits = c(0,100))+
  theme_bw()+
  NoLegend()
tmp2 = melt(data.frame("all" = length(levels(Idents(seurat.intactTAS)))-length(unique(markers.intactTAS$cluster)), "value" = length(unique(markers.intactTAS$cluster)), "variable" = "IDed"))
tmp2$variable = factor(c("noID","IDed"),levels = c("noID","IDed"))
g2 = ggplot(tmp2[,-2], aes(x = "1",y = value, fill = variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(breaks = seq(0,100,5),limits = c(0,100))+
  scale_fill_grey()+
  theme_bw()+
  NoLegend()

# Evaluation of clustering
pdf("figure/FigS1AC.pdf",width = 10,heigh = 3)
library(gridExtra)
grid.arrange(g1,g2,ncol = 2, layout_matrix = matrix(c(rep(1,10),2),nrow = 1),top = "FindAllMarkers(seurat,logfc.threshold = 1,only.pos = T,min.pct = 0.5,verbose = F)")
dev.off()


# Heatmap for signature genes
library(ggplot2)
pdf("figure/Fig2A.pdf",width = 12)
DoHeatmap(seurat.intactTAS,features = unique(markers.intactTAS$gene),lines.width = 3)+
scale_fill_gradientn(colours = c("black","yellow"))+
NoLegend()
dev.off()

pdf("figure/Fig2A-legend.pdf",width = 12)
DoHeatmap(seurat.intactTAS,features = unique(markers.intactTAS$gene),lines.width = 3)+
  scale_fill_gradientn(colours = c("black","yellow"))
dev.off()


# Calusulating cell cycling score
library(openxlsx)
cell.cycle = read.xlsx("/mnt/hdd/hdd/data/species/D.japonica/DjTrascriptome/known_marker_X-ray.xlsx",sheet = 2)
cell.cycle$DjID = gsub("_","-",cell.cycle$DjID)
cell.cycle = cell.cycle[!is.na(cell.cycle$DjID),]
rownames(cell.cycle) = cell.cycle$DjID
seurat.intactTAS <- CellCycleScoring(seurat.intactTAS, s.features = cell.cycle$DjID[cell.cycle$lineage =="S-phase"], g2m.features = cell.cycle$DjID[cell.cycle$lineage =="M-phase"], set.ident = F)
VlnPlot(seurat.intactTAS, features = c("G2M.Score","S.Score"),ncol = 1,pt.size = 0)
data = data.frame("Cell.type" = seurat.intactTAS@active.ident, "Phase" = seurat.intactTAS@meta.data$Phase)

# Compariosn Dj and Smed scRNA-Seq 
rawcnt.smed = read.table("/mnt/hdd/hdd/data/species/S.mediterranea/Fincher2018/GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt")
at.smed = read.xlsx("/mnt/hdd/hdd/data/species/S.mediterranea/Fincher2018/aaq1736_Table-S1.xlsx")
markers.major.smed = read.xlsx("/mnt/hdd/hdd/data/species/S.mediterranea/Fincher2018/aaq1736_Table-S2.xlsx",sheet = 1,startRow = 6)
rownames(at.smed) = at.smed$Cell.ID
at.smed = at.smed[colnames(rawcnt.smed),]
seurat.smed = CreateSeuratObject(rawcnt.smed,project = "Smed_Fincher")
# Idents(seurat.smed@active.ident) = factor(at.smed$Subcluster.ID,levels = unique(at.smed$Subcluster.ID))
seurat.smed@meta.data$section = at.smed$Section
seurat.smed@meta.data$FACS = at.smed$FACS_Gate
seurat.smed@meta.data$Major.cluster = gsub(" ","",at.smed$Major.cluster.description)
seurat.smed@meta.data$Sub.cluster = gsub(" ","",at.smed$Subcluster.ID)
seurat.smed = NormalizeData(seurat.smed)
seurat.smed = ScaleData(seurat.smed)
seurat.smed = FindVariableFeatures(seurat.smed)
seurat.smed = RunPCA(seurat.smed)
Idents(seurat.smed) = sprintf("%s@%s",seurat.smed@meta.data$Major.cluster,seurat.smed@meta.data$Sub.cluster)
ElbowPlot(seurat.smed)
seurat.smed = RunUMAP(seurat.smed,dims = 1:50)
DimPlot(seurat.smed)+NoLegend()
# VlnPlot(seurat.smed,features = "nCount_RNA",group.by = "Major.cluster")+NoLegend()
# VlnPlot(seurat.smed,features = "nFeature_RNA",group.by = "Major.cluster")+NoLegend()
save(seurat.smed, file = "/mnt/hdd/hdd/data/species/S.mediterranea/Fincher2018/seurat.smed")

library(dplyr)
library(reshape2)
library(ggplot2)
Dj.dd4 = read.table("/mnt/hdd/hdd/data/species/D.japonica/DjTrascriptome/DjTrascriptome-Sned-scRNA-dd4.txt",sep = "\t")
Dj.dd4$V1 = gsub("_","-",Dj.dd4$V1)
Dj.dd4$V2 = gsub("_","-",Dj.dd4$V2)
Dj.dd4.2 = as.data.frame(Dj.dd4 %>%  dplyr::distinct(V1,V2,.keep_all = T))[,1:2]
rownames(Dj.dd4.2) = Dj.dd4.2$V1

Dj.dd4 = read.table("/mnt/hdd/hdd/data/species/D.japonica/DjTrascriptome/DjTrascriptome-Smed-scRNA-dd4.txt",sep = "\t")
# Dj.dd4 = rbind(Dj.dd4[which(Dj.dd4$V1=="dd_Smed_v4_659_0_1"),],Dj.dd4)
Dj.dd4$V1 = gsub("_","-",Dj.dd4$V1)
Dj.dd4$V2 = gsub("_","-",Dj.dd4$V2)
Dj.dd4 = as.data.frame(Dj.dd4 %>%  dplyr::distinct(V1,V2,.keep_all = T))[,1:2]

Dj.smed.ortholog = NULL
for(i in 1:nrow(Dj.dd4)){
  tmp = Dj.dd4[i,]
  tmp2 =Dj.dd4.2[tmp$V2,]
  if(!is.na(tmp2$V1)){
    if(tmp[1]==tmp2[2]){
      Dj.smed.ortholog = rbind(Dj.smed.ortholog,tmp2)
    }
  }
}
library(openxlsx)
write.xlsx(Dj.smed.ortholog,file = "TableS3.xlsx")
Dj.dd4 = Dj.smed.ortholog

seurat.intactTAS = ScaleData(seurat.intactTAS,features = rownames(seurat.intactTAS@assays$RNA@data))
seurat.smed = ScaleData(seurat.smed,features = rownames(seurat.smed@assays$RNA@data))

Ave.ex.Dj = AverageExpression(seurat.intactTAS,slot = "data",features = Dj.dd4$V2)
Ave.ex.Dj$RNA = t(apply(Ave.ex.Dj$RNA, 1, scale))
rownames(Ave.ex.Dj$RNA) = Dj.dd4$V2
Ave.ex.Dj$RNA[is.na(Ave.ex.Dj$RNA)] = 0

Ave.ex.smed = AverageExpression(seurat.smed,slot = "data",features = Dj.dd4$V1[is.element(Dj.dd4$V1, rownames(seurat.smed@assays$RNA@counts))])
Ave.ex.smed$RNA = t(apply(Ave.ex.smed$RNA, 1, scale))
rownames(Ave.ex.smed$RNA) = Dj.dd4$V1[is.element(Dj.dd4$V1, rownames(seurat.smed@assays$RNA@counts))]
Ave.ex.smed$RNA[is.na(Ave.ex.smed$RNA)] = 0

##### Major cell type estimation based on Smed marker genes
#  Comparison the expressions of Smed marker gene in Smed and Dja
Idents(seurat.smed) = seurat.smed@meta.data$Major.cluster
Ave.ex.smed = AverageExpression(seurat.smed,slot = "data")
Ave.ex.smed$RNA = t(apply(Ave.ex.smed$RNA, 1, scale))
rownames(Ave.ex.smed$RNA) = rownames(seurat.smed@assays$RNA@counts)
colnames(Ave.ex.smed$RNA) = levels(seurat.smed)
Ave.ex.smed$RNA[is.na(Ave.ex.smed$RNA)] = 0
Dj.homolog = Ave.ex.Dj$RNA[Dj.dd4$V2[is.element(Dj.dd4$V2,rownames(Ave.ex.Dj$RNA))],]
Smed.homolog = Ave.ex.smed$RNA[Dj.dd4$V1[is.element(Dj.dd4$V1,rownames(Ave.ex.smed$RNA))],]
colnames(Smed.homolog) = levels(Idents(seurat.smed))
common = rownames(Dj.homolog)[is.element(rownames(Dj.homolog),Dj.dd4$V2[is.element(Dj.dd4$V1,rownames(Smed.homolog))])]
Dj.homolog = Dj.homolog[common,]
Smed.homolog = Smed.homolog[Dj.dd4[is.element(Dj.dd4$V2,common),]$V1,]
cors.Dj.Smed.cell.type = NULL
for(i in 1:ncol(Dj.homolog)){
  tmp = apply(Smed.homolog, 2, FUN = function(x){
    return(cor(Dj.homolog[,i],x))
  })
  cors.Dj.Smed.cell.type = rbind(cors.Dj.Smed.cell.type,tmp)
}
rownames(cors.Dj.Smed.cell.type) = colnames(Dj.homolog)
cors.Dj.Smed.cell.type[cors.Dj.Smed.cell.type<0] = 0
tmp = apply(cors.Dj.Smed.cell.type, 1, FUN = function(x){
  max = max(x)
  tmp = rep(0, length(x))
  tmp[which(x==max)]=1
  return(tmp)
})
tmp = t(tmp)
rownames(tmp) = levels(Idents(seurat.intactTAS))
colnames(tmp) = colnames(Smed.homolog)
major.cell.type = NULL
for(i in 1:nrow(tmp)){
  major.cell.type = c(major.cell.type,names(which(tmp[i,]==1)))
}
tmp2 = NULL
for(i in 1:length(seurat.intactTAS@active.ident)){
  tmp2 = c(tmp2,major.cell.type[seurat.intactTAS@active.ident[i]])
}
seurat.intactTAS@meta.data$major.cell.type = tmp2

write.xlsx(data.frame("Minor.cell.type" = levels(Idents(seurat.intactTAS)),"Major.cell.type"=major.cell.type,"Num.of.Cell" = table(seurat.intactTAS@active.ident)),file = "TableS4.xlsx",overwrite = T)

tmp = NULL
list = list("Neoblast"=c("piwiA", "DjGI005146-001"),
            "Epidermal"=c("prog1/2","DjGI012450-001"),
            "Neural" = c("syt","DjGI007010-001"),
            "Muscle" = c("col4A6A","DjGI001767-002"),
            "Pharynx" = c("dd1071","DjGI008213-001"),
            "Intestine" = c("smed_dd75","DjGI000320-002"),
            "Protonephridia" = c("egfr5", "DjGI009250-001"), 
            "Parapharyngeal" = c("ano7","DjGI005169-001"),
            "Cathepsin+cells"=c("ctsl2","DjGI009130-001"))
g0 = NULL
library(patchwork)
for(i in 1:length(list)){
  g = VlnPlot(seurat.intactTAS2,features = list[[i]][2])+NoLegend()+ggtitle(label = sprintf("%s (%s marker)",list[[i]][1],names(list)[i] ),subtitle = list[[i]][2])+xlab(NULL)
  g0 = c(g0,list(g))
}
pdf("figure/Fig2B.pdf",width = 15,height = 7)
wrap_plots(g0)
dev.off()

pdf("figure/FigS2.pdf",height = 5, width = 10)
VlnPlot(seurat.intactTAS2,features = c("nCount_RNA","nFeature_RNA"),group.by = "major.cell.type",pt.size = 0)
dev.off()

markers.major.cell.type = FindAllMarkers(seurat.intactTAS2,logfc.threshold = 2,only.pos = T)

Table.S1 = data.frame("CellID" = names(seurat.intactTAS@active.ident), "Major.cell.type" = seurat.intactTAS@meta.data$major.cell.type,"Minor.cell.type" = seurat.intactTAS@active.ident,"S.Score" = seurat.intactTAS$S.Score,"G2M.Score" = seurat.intactTAS$G2M.Score)
write.xlsx(Table.S1,file = "TableS1.xlsx",overwrite = T)
load("/mnt/hdd/hdd/data/species/D.japonica/DjTrascriptome/DjTranscriptome_des")
rownames(des) = gsub("_","-",rownames(des) )
Table.S2 = cbind(markers.intactTAS, des[gsub("-","_",markers.intactTAS$gene),])
write.xlsx(Table.S2,file = "TableS2.xlsx",overwrite = T)

# Cell type composition
a = sort(table(seurat.intactTAS@meta.data$major.cell.type))
data = data.frame("Major.Cell.Type" =seurat.intactTAS$major.cell.type)
data$Major.Cell.Type = factor(data$Major.Cell.Type,levels = names(a))
ggplot(data = data, aes(fill = Major.Cell.Type)) +
  geom_bar(stat = "count",
           aes(x = 1),
           position = position_stack(), color = "black") +
  geom_text(stat = "count",
            aes(x = 1, y = ..count..,
                label = ..count..),
            position = position_stack(vjust = 0.5)) +
  geom_text(stat = "count",
            aes(x = 1.55, y = ..count..,
                label = Major.Cell.Type),
            position = position_stack(vjust = 0.5)) +
  theme_bw()+
  theme(legend.position = "none") +
  coord_polar(theta = "y")
# Fig.3
library(openxlsx)
known.gene = read.xlsx("/mnt/hdd/hdd/data/species/D.japonica/DjTrascriptome/known_marker_X-ray.xlsx",sheet = 1)
known.gene = known.gene[!is.na(known.gene$DjID),]
known.gene$DjID = gsub("_","-",known.gene$DjID)
rownames(known.gene) = known.gene$GeneName
VlnPlot(seurat.intactTAS,features = known.gene[c("Djp2x-A"),]$DjID)+NoLegend()
cell.posi = c("7","28","54","55","56","63")
col = ggColorHue(length(levels(Idents(seurat.intactTAS))))
col[(1:63)[!is.element(1:63,as.integer(cell.posi))]] = "white"
pdf("figure/Fig4.pdf",width = 18, height = 14)
g4 = VlnPlot(seurat.intactTAS,features = c(known.gene[c("piwiA"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("piwiA")+scale_fill_manual(values = col)+NoLegend()
g1 = VlnPlot(seurat.intactTAS,features = c(known.gene[c("Djp2x-A"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("p2x-A")
g6 = VlnPlot(seurat.intactTAS,features = c(known.gene[c("tgs1"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("tgs1")
g0 = VlnPlot(seurat.intactTAS,features = c(known.gene[c("Djcalu"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("calu")+scale_fill_manual(values = col)
g5 = VlnPlot(seurat.intactTAS,features = c("DjGI002326-001"),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("cdh1")+scale_fill_manual(values = col)
g2 = VlnPlot(seurat.intactTAS,features = "S.Score",ncol = 1,pt.size = 0)+scale_fill_manual(values = col)+NoLegend()
g3 = VlnPlot(seurat.intactTAS,features = "G2M.Score",ncol = 1,pt.size = 0)+scale_fill_manual(values = col)+NoLegend()
grid.arrange(g1,g4,g6,g0,g5,g2,g3,ncol = 1)
dev.off()

gnl.fig3 = rbind(known.gene["Djp2x-A",],known.gene[c(1:3,5:33),])
Ave = AverageExpression(seurat.intactTAS,slot = "data",features = gnl.fig3$DjID)
Ave$RNA = t(apply(Ave$RNA, 1, scale))
rownames(Ave$RNA) = gnl.fig3$DjID
colnames(Ave$RNA) = levels(seurat.intactTAS@active.ident)

tmp = Ave$RNA[gnl.fig3$DjID,]
tmp = melt(tmp)
# tmp$Var2 = factor(tmp$Var2,levels = c("Neoblast","Epidermal","Neural","Muscle","Pharynx","Intestine","Protonephridia","Parapharyngeal","Cathepsin+cells"))
g4 = ggplot(tmp,aes(x = as.factor(Var2), y = rev(Var1), fill = value))+
  geom_tile()+
  scale_fill_gradientn(colours = c("black","yellow","orange"))+
  scale_y_discrete(label = rev(gnl.fig3$GeneName))+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black"))
pdf("figure/Fig3B.pdf",width = 12,height = 7)
print(g4)
dev.off()

markers.63 = markers.intactTAS[markers.intactTAS$cluster=="63",]

# Fig.4
VlnPlot(seurat.intactTAS,features = known.gene[c("MTA-B"),]$DjID)+NoLegend()
cell.posi.B = as.character(c(9,14:17,28:30,34,35,46,48,50:52,54:56,58,60:63))
cell.posi.A = as.character(c(9,15:17,35,51:52,54:56,60,61,63))
col = rep("white",length(levels(Idents(seurat.intactTAS))))
names(col) = 1:length(levels(Idents(seurat.intactTAS)))
col[cell.posi.A] = "cyan"
col[cell.posi.B[!is.element(cell.posi.B,cell.posi.A)]] = "magenta"
g4 = VlnPlot(seurat.intactTAS,features = c(known.gene[c("piwiA"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("piwiA")+scale_fill_manual(values = col)+NoLegend()
g1 = VlnPlot(seurat.intactTAS,features = c(known.gene[c("MTA-A"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("MTA-A")+scale_fill_manual(values = col)+NoLegend()
g6 = VlnPlot(seurat.intactTAS,features = c(known.gene[c("MTA-B"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("MTA-B")+scale_fill_manual(values = col)+NoLegend()
g0 = VlnPlot(seurat.intactTAS,features = c(known.gene[c("Djcalu"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("calu")+scale_fill_manual(values = col)
g5 = VlnPlot(seurat.intactTAS,features = c(known.gene[c("Djsyt"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("syt")+scale_fill_manual(values = col)
g2 = VlnPlot(seurat.intactTAS,features = c("DjGI002326-001"),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("cdh1")+scale_fill_manual(values = col)
g3 = VlnPlot(seurat.intactTAS,features = c("DjGI009461-001"),ncol = 1,pt.size = 0)+NoLegend()+ggtitle("cdc42")+scale_fill_manual(values = col)
pdf("figure/Fig3.pdf",width = 18, height = 14)
grid.arrange(g1,g6,g4,g5,g0,g2,g3,ncol = 1)
dev.off()

# X-ray
load("/mnt/hdd/hdd/data/species/D.japonica/Seurat/seurat.X.ray")
seurat.X.ray = NormalizeData(seurat.X.ray,features = rownames(seurat.X.ray@assays$RNA@data))
seurat.X.ray = ScaleData(seurat.X.ray,do.scale = T,features = rownames(seurat.X.ray@assays$RNA@data))
Idents(seurat.X.ray) = seurat.X.ray@meta.data$day
markers.x = FindAllMarkers(seurat.X.ray,only.pos = T)
markers.x = markers.x %>%  arrange(-avg_log2FC, ) %>% arrange(cluster)
library(openxlsx)
ave.X = AverageExpression(seurat.X.ray,group.by = "day",slot = "data")
data = data.frame("Contig" = DEG, "Day.after.irradiatuion" = ave.X$RNA[DEG,],"blast.result" = des[gsub("-","_",DEG),])
write.xlsx(data, file = "TableS5.xlsx",sep = "\t")


load("/mnt/hdd/hdd/data/species/D.japonica/DjTrascriptome/DjTranscriptome_des")
rownames(des) = gsub("_","-",rownames(des))
library(gridExtra)
library(patchwork)

gnl = unique(markers.x$gene)
gnl = gnl[is.element(gnl, rownames(seurat.intactTAS@assays$RNA@scale.data))]
gnl = rev(gnl)
Ave.sc = AverageExpression(seurat.intactTAS,features = gnl,slot = "data")
Ave.sc$RNA = t(apply(Ave.sc$RNA, 1, scale))
rownames(Ave.sc$RNA) = gnl
Ave.bulk = AverageExpression(seurat.X.ray,features = gnl,slot = "data",group.by = "day")
Ave.bulk$RNA = t(apply(Ave.bulk$RNA, 1, scale))
rownames(Ave.bulk$RNA) = gnl
tmp = Ave.sc$RNA[gnl,]
tmp = melt(tmp)
g4 = ggplot(tmp,aes(x = as.factor(Var2), y = Var1, fill = value))+
  geom_tile()+
  scale_fill_gradientn(colours = c("black","yellow"))+
  scale_y_discrete(label = NULL)+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black"))+
  ggtitle("scRNA-Seq")
print(g4)

tmp = Ave.bulk$RNA[gnl,]
tmp = melt(tmp)
g5 = ggplot(tmp,aes(x = as.factor(Var2), y = Var1, fill = value))+
  geom_tile()+
  scale_fill_gradientn(colours = c("magenta","black","yellow"))+
  scale_y_discrete(label = NULL)+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black"))+
  ggtitle("bulkRNA-Seq")
print(g5)

pdf("figure/Fig5A.pdf", heigh = 7,width = 15)
print(g5|g4)
dev.off()

library(dplyr)
gnl = markers.x[!is.na(des[markers.x$gene,]$SwissProt_name),] %>% group_by(cluster) %>% top_n(n = 5,wt = avg_log2FC)
gnl = unique(gnl$gene)
pdf("figure/Fig5B.pdf",width = 14)
for(gn in gnl){
  g1 = VlnPlot(seurat.intactTAS,features = gn)+NoLegend()+ggtitle(NULL)+xlab(NULL)
  data = data.frame("Expression" = as.numeric(seurat.X.ray@assays$RNA@data[gn,]),"dpi" =seurat.X.ray@meta.data$day)
  a = data %>% group_by(dpi) %>% summarise(mean(Expression))
  b = splinefun(a$dpi,a$`mean(Expression)`)
  g2 = ggplot(data, aes(x = dpi, y = Expression))+
    geom_point()+
    geom_line(data = data.frame("Expression" = b, (seq(0,7,0.01)),"dpi" = seq(0,7,0.01)), col = "red")+
    scale_x_continuous(breaks=unique(seurat.X.ray@meta.data$day))+
    xlab("Days post irradiation")+
    theme_bw()
  g = {g1/g2}+plot_annotation(title = gn, subtitle = des[gn,]$SwissProt_name)& theme(text = element_text(size = 10))
  print(g)
  
}
dev.off()

pdf("figure/FigS4A.pdf")
for(gn in c("piwiA","zfp-1","p53","prog-12")){
  gid = known.gene[gn,]$DjID
  data = data.frame("Expression" = as.numeric(seurat.X.ray@assays$RNA@data[gid,]),"dpi" =seurat.X.ray@meta.data$day)
  ex = data$Expression
  a = data %>% group_by(dpi) %>% summarise(mean(Expression))
  g = ggplot(data, aes(x = dpi, y = Expression))+
    geom_point()+
    geom_line(data = data.frame("Expression" = a$`mean(Expression)`,"dpi" = unique(seurat.X.ray@meta.data$day)), col = "black")+
    scale_x_continuous(breaks=unique(seurat.X.ray@meta.data$day))+
    xlab("Days post irradiation")+
    ggtitle(gn)+
    theme_bw()
  print(g)
  
}
dev.off()

library(openxlsx)
DEG = markers.x$gene
data = data.frame("Contig" = DEG, "Day.after.irradiation" = ave.X$RNA[DEG,],"Group" = markers.x$cluster, "blast.result" = des[DEG,])
data$Group = c(rep(1,307),rep(1,20),rep(3,48), rep(4,7), rep(5, 15), rep(6,154))
write.xlsx(data, file = "TbaleS4.xlsx",sep = "\t")

# progenitor makers
Ave.piwiA = AverageExpression(seurat.intactTAS,features = known.gene[c("piwiA"),]$DjID)
Ave.piwiA = Ave.piwiA$RNA[1,]
plot(Ave.piwiA)
plot(sort(Ave.piwiA,decreasing = T))
abline(v = seq(0,65,5))
abline(h = 0.5)
sum(Ave.piwiA[Ave.piwiA>0.5])
# VlnPlot(seurat.intactTAS,features = c(known.gene[c("piwiA"),]$DjID),ncol = 1,pt.size = 0)+NoLegend()
cells = names(seurat.intactTAS@active.ident)[is.element(seurat.intactTAS@active.ident,names(Ave.piwiA[Ave.piwiA>0.5]))]
tmp0 = subset(seurat.intactTAS,cells = cells)
VlnPlot(tmp,features = "nCount_RNA")
DoHeatmap(tmp, features = known.gene$DjID)
markers8 = NULL
seurat.8 = subset(seurat.intactTAS,idents = "8")
for(i in unique(tmp0@active.ident)[unique(tmp0@active.ident)!=8]){
  tmp = subset(tmp0,idents = i)
  tmp = merge(seurat.8,tmp)
  tmp = FindAllMarkers(tmp,logfc.threshold = 2,only.pos = T)
  markers8 = rbind(markers8,tmp)
}
length(unique(markers8$gene))
DoHeatmap(tmp0,features = unique(markers8$gene))
negative.markers = unique(markers8$gene)
save(negative.markers,file = "/mnt/hdd/hdd/data/species/D.japonica/Seurat/markers.intactTAS.nNb")
