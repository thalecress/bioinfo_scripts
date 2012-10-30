library (ShortRead)
#get the fasta reads from sly and spe
reads.sly <- readFasta("SlITAG1.9_filtered_plus_Heinz.fa") 
sread<- as.vector(sread (reads.sly))
id<- as.vector(id (reads.sly))  #gene name
reads.sly.csv <- cbind(id, sread) #sequence
colnames (reads.sly.csv)<- c("ID","SLY")
colnames (reads.sly.csv)
dim (reads.sly.csv)

reads.spe <- readFasta("SpITAG1.9_filtered_plus_Heinz.fa") 
sread<- as.vector(sread (reads.spe))
id<- as.vector(id (reads.spe))
reads.spe.csv <- cbind(id, sread)
colnames (reads.spe.csv)<- c("ID","SPE")
colnames (reads.spe.csv)
dim (reads.spe.csv)

#merge the reads
reads.merged <- merge (reads.sly.csv, reads.spe.csv, by = 1,1)
dim (reads.merged)
head (reads.merged)
summary (reads.merged)

#get gene coordinates table .tsv
gene_coor<- read.table ("gene_coordinates_table.tsv", sep="\t")
head (gene_coor)
gene_coor$gene.order <- row.names (gene_coor)
plot (gene_coor$V3, gene_coor$gene.order)
colnames (gene_coor)[1]<- "ITAG"
colnames (gene_coor)[2]<- "Chr"
colnames (gene_coor)[3]<- "Start"
colnames (gene_coor)[4]<- "Stop"
head(gene_coor)
write.csv(gene_coor, "gene.order.forDan.csv")
#Reorder merged genes by the proper order based on gene_coor$number
gene_ord<- gene_coor[,c(1,5)]
reads.merged.ord <- merge(reads.merged, gene_ord, by=1,1)
names (reads.merged.ord)
reads.merged.ord$number <- as.numeric (reads.merged.ord$number)
plot (reads.merged.ord$number)  #Some of the genes are not in the right place.  Correct the order for these.
reads.merged.ord<- reads.merged.ord[order(reads.merged.ord$number),]
plot (reads.merged.ord$number)  #order looks correct now
rownames (reads.merged.ord)<- reads.merged.ord$ID
rownames (reads.merged.ord)

#How to see if the swap has worked
Sly_len<-nchar (as.character(reads.merged.ord [,2]))  #Sly lengths
Spe_len<-nchar (as.character(reads.merged.ord [,3]))  #Spe lengths
plot (Sly_len,Spe_len )
reads.merged.ord<- cbind (reads.merged.ord, Sly_len, Spe_len)
names (reads.merged.ord)
summary (reads.merged.ord)

#get the IL_coordinates
system ("mkdir ./Swapped_refs")
#bounds <- read.csv ("il_and_genes.csv", header = TRUE, row.names=1)
bounds <- read.csv ("il_and_genes_new_all.csv", header = TRUE, row.names=1) # 20120820, see new swap steps 
names (bounds)
#Make an empty dataframe with the following columns: IL, Swap outcome, Genes in swap region
d <-  data.frame (IL="", Swap_outcome="",No_swapped_genes=0,  
                 stringsAsFactors=FALSE) # you don't know levels yet
d
n=0
library (seqinr) #needed for the write fasta command
for (i in rownames (bounds)){
  n=n+1
  d[n,1]<- i
  reads.merged.ord.tmp<- reads.merged.ord
  reads.merged.ord.tmp$SLY<- as.character (reads.merged.ord$SLY)
  reads.merged.ord.tmp$SPE<- as.character (reads.merged.ord$SPE)
  #if (grep ("a",i) ne integer(0)){
   # next
 # }
 # if (grep ("b",i) >1){
#    next
#  }
#}
  Start_ITAG<- bounds[i,4]
  End_ITAG<- bounds[i,6]
  Start_ITAG_row<- grep (Start_ITAG, rownames(reads.merged.ord.tmp))
  End_ITAG_row<- grep (End_ITAG, rownames(reads.merged.ord.tmp))
  # For_swap_SPE<- IL_SPE[c(Start_ITAG_row:End_ITAG_row),]
  reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),2] <-  reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),3] #Swap!
  newdf <-reads.merged.ord.tmp[,c(1:2)]
  newdf$SLY<- strsplit(newdf$SLY,"")
  listdf<- as.list (newdf)
  write.fasta( names = listdf$ID, sequences= listdf$SLY,paste ("./Swapped_refs/", i, ".swap.fa", sep=""))
  #write.fasta( names = reads.merged.ord.tmp$ID, sequences= reads.merged.ord.tmp$SLY,paste ("./Swapped_refs/", i, ".swap.fa", sep=""))
  #now check if the swap worked.  Find colsum of the difference of  gene lengths of Spe_len column and the swapped region 
  #if it is 0, print i: "All systems A OK"
  #If it is ne 0 then print i: "Error!  Error!"
  lentest<-sum(nchar (as.character(reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),2])) - reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),6])
  if (lentest == 0){
    d[n,2]<- paste ("lentest is ", lentest,". All systems A OK!", sep="")
  }
  if (lentest != 0){
    d[n,2]<-paste ("Error!  Error! Recalculating", sep="")
  }
#Print the number of genes in the region
  swapgenes<- nrow (reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),])
  d[n,3]<- paste (swapgenes)

}
#write the summary  
#write.csv (d,"swap_summary.csv")
write.csv (d,"swap_summary_all.csv")



#NOW for the ILs that have a and b in it.
#make the batch script into a function (in swap.function.r)
#Open IL2.1.1a
#cd into /Swapped_refs
system ("cd .//Swapped_refs")
source ("./swap.function.r")
setwd ("/mydata/4.Reference_files/1.9/Swapped_refs/")
#swap_a_b ("IL12.2 (a).swap.fa", "IL12.2 (b).swap.fa", "IL12.2 (b)", "IL12.2")
#swap_a_b ("IL9.2.6 (a).swap.fa", "IL9.2.6 (b).swap.fa", "IL9.2.6 (b)", "IL9.2.6")
#swap_a_b ("IL2.1.1 (a).swap.fa", "IL2.1.1 (b).swap.fa", "IL2.1.1 (b)", "IL2.1.1")
#swap_a_b ("IL3.3 (a).swap.fa", "IL3.3 (b).swap.fa", "IL3.3 (b)", "IL3.3")
swap_a_b ("IL2.2(a).swap.fa", "IL2.2(b).swap.fa", "IL2.2(b)", "IL2.2")
swap_a_b ("IL9.3.1(chr09).swap.fa", "IL9.3.1(**chr12**).swap.fa", "IL9.3.1(**chr12**)", "IL9.3.1")
swap_a_b ("IL2.3(a).swap.fa", "IL2.3(b).swap.fa", "IL2.3(b)", "IL2.3")
swap_a_b ("IL2.3.swap.fa", "IL2.3(c).swap.fa", "IL2.3(c)", "IL2.3")
swap_a_b ("IL3.2(a).swap.fa", "IL3.2(b).swap.fa", "IL3.2(b)", "IL3.2")
swap_a_b ("IL3.2.swap.fa", "IL3.2(c).swap.fa", "IL3.2(c)", "IL3.2")

#####END#####

#Get the genes in each IL:
#Make an empty dataframe with the following columns: IL, Swap outcome, Genes in swap region
system ("mkdir IL_Genes")
n=0
library (seqinr) #needed for the write fasta command
for (i in rownames (bounds)){
  reads.merged.ord.tmp<- reads.merged.ord
  reads.merged.ord.tmp$SLY<- as.character (reads.merged.ord$SLY)
  reads.merged.ord.tmp$SPE<- as.character (reads.merged.ord$SPE)
  Start_ITAG<- bounds[i,4]
  End_ITAG<- bounds[i,6]
  Start_ITAG_row<- grep (Start_ITAG, rownames(reads.merged.ord.tmp))
  End_ITAG_row<- grep (End_ITAG, rownames(reads.merged.ord.tmp))
  reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),2] <-  reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),3] #Swap!
  lentest<-sum(nchar (as.character(reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),2])) - reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),6])
  if (lentest == 0){
    print (paste ("lentest is ", lentest,". All systems A OK!", sep=""))
  }
  if (lentest != 0){
    print(paste ("Error!  Error! Recalculating", sep=""))
  }
  #Print the number of genes in the region
  swapgenes<- rownames (reads.merged.ord.tmp[c(Start_ITAG_row:End_ITAG_row),])
  #print(swapgenes)
  write.csv(swapgenes, paste ("./IL_Genes/", i, "_genelist.txt", sep=""))
}

#Make a merged file
All_genes<- as.data.frame(rownames (reads.merged.ord))
names (All_genes)[1]<- "genes"
InputFiles<- as.character(list.files("./IL_Genes/"))
n=1
for (i in InputFiles) {
  n=n+1
  IL_name<- sub ("_.*","",i)  
  gene_list<- read.csv(paste("./IL_Genes/", i, sep=""), row.names=1)
  gene_list [c(1:nrow(gene_list)), 2]<- as.character(IL_name)
  if (exists ("merged_file")){
    merged_file<-merge (merged_file, gene_list, by = 1,1, all.x = T, all.y = T)
    names (merged_file)[n]<- i 
  }
  else{
    merged_file<-merge (All_genes, gene_list, by = 1,1, all.x = T, all.y = T)
    names (merged_file)[n]<- i 
  }
}
names (merged_file) <- sub ("_.*","",colnames(merged_file))  
row.names(merged_file)<- merged_file$genes
merged_file$genes <- NULL
ILs_and_genes<- as.data.frame(apply (merged_file, 1, function(x)paste(x[!is.na(x)], collapse=",")))#paste, collapse=",")
names (ILs_and_genes)[1] <- "IL"
summary (ILs_and_genes)
write.csv(ILs_and_genes, "ILs_and_genes.csv")
getwd()
______________
ncol (merged_file)

merge_rows<- function (merged_file){
for (n in rownames(merged_file)){
  for (j in 2:87){
    if (is.na(merged_file[n,j]) == FALSE){
      print (merged_file[n,j])
    }
  }
}
}
cl <- makeCluster(16, type = "SOCK")
clusterCall(cl, merge_rows(merged_file)) 
stopCluster (cl)

library (multicore)

parallel(merge_rows(merged_file))

merged_file[n,88]<- IL_combined
}
}

