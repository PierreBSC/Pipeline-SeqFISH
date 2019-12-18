args <- commandArgs(trailingOnly = T)

file_directory = as.character(args[1])
output_directory = as.character(args[2])


library(pheatmap)

#Loading and processing the general cell properties 

Cells_QC = read.delim(paste(file_directory,"/Cell_properties_table.txt",sep = ""),row.names = 1)
Cells_QC$Position = factor(Cells_QC$Position)
Cells_QC$Cell_size_Log10 = log10(Cells_QC$Cell_size)

##Plotting the information on a new pdf :

path_QC_cells_pdf = paste(output_directory,"QC_cells.pdf",sep="/")
pdf(path_QC_cells_pdf,width = 9,height = 9)
par(las=1,mfrow=c(2,2))

#Plotting the relation ship between cell size and total RNA count 

plot(log10(Cells_QC$Cell_size),log10(Cells_QC$Total_RNA_molecules),xaxs="i",yaxs="i",
     ylim= c(min(log10(Cells_QC$Total_RNA_molecules))*0.9,max(log10(Cells_QC$Total_RNA_molecules))*1.1),
     xlim= c(min(log10(Cells_QC$Cell_size))*0.9,max(log10(Cells_QC$Cell_size))*1.1),
     xlab="Cell size (Log10 pixel)",ylab="Total RNA spots (Log10)",pch=21,bg="orange")
m = lm(log10(Cells_QC$Total_RNA_molecules)~log10(Cells_QC$Cell_size))
abline(coef(m),lwd=2,lty=2,col="grey")
R = sqrt(summary(m)$r.squared)
legend("topleft",legend = paste("R = ",round(R,2)),bty='n')

N_cell_position = table(Cells_QC$Position)
barplot(N_cell_position,xlab="Position",ylim=c(0,1.2*max(N_cell_position)),cex.lab=1.4)
title(main="Number of cells per position",adj=0,line=-1,cex=2)


Cells_expression = read.delim("Desktop/Output/RNA_expression_table.txt",row.names = 1)
Cells_expression_normalised = Cells_expression/rowSums(Cells_expression)
pheatmap(t(Cells_expression_normalised),clustering_method = "ward",show_colnames = F,
         annotation_col = Cells_QC[,c(1,7)],annotation_legend = F)
dev.off()
