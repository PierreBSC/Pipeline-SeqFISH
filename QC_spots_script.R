args <- commandArgs(trailingOnly = T)

file_directory = as.character(args[1])
output_directory = as.character(args[2])
color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("grey","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.999,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}

string.to.colors = function (string, colors = NULL) {
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}

#Loading the data into the R environment

Data_table = read.delim(paste(file_directory,"Concatenated_spot_information.txt",sep = "/"))
colnames(Data_table) = c("Position_X","Position_Y","Position_Z",'Intensity','Passed_Filtering','Noise','Round','Position','Channel')

Data_table$Individual_gene = paste(Data_table$Round,Data_table$Channel,sep = "_")
List_individual_genes = unique(Data_table$Individual_gene)
List_individual_genes = paste("Gene ",as.numeric(factor(List_individual_genes)),sep="")


Data_table$SNR = Data_table$Intensity/Data_table$Noise

General_information = read.delim(paste(file_directory,"General_experiment_information.txt",sep = "/"))

N_rounds = General_information[1,1]
N_channels = General_information[1,2]
N_positions = General_information[1,3]

###First plot : for each position we show the number of raw spots,
##filtered spots, filtering ratio, and SNR for each gene 

path_QC_spot_pdf = paste(output_directory,"QC_spots.pdf",sep="/")
pdf(path_QC_spot_pdf,width = 9,height = 12)
for (P in 1:N_positions) {
  
  ##Computing the different values for the Position of Interest
  N_spot_raw_temps = table(Data_table$Individual_gene[Data_table$Position == P])
  N_spot_filtered_temps = table(Data_table$Individual_gene[Data_table$Position == P & Data_table$Passed_Filtering==1])
  
  Mean_SNR_raw_temp = aggregate(Data_table$SNR[Data_table$Position == P],by=list(Data_table$Individual_gene[Data_table$Position == P]),FUN=mean)$x
  Sd_SNR_raw_temp = aggregate(Data_table$SNR[Data_table$Position == P],by=list(Data_table$Individual_gene[Data_table$Position == P]),FUN=sd)$x
  
  Mean_SNR_filtered_temp = aggregate(Data_table$SNR[Data_table$Position == P & Data_table$Passed_Filtering==1],by=list(Data_table$Individual_gene[Data_table$Position == P & Data_table$Passed_Filtering==1]),FUN=mean)$x
  Sd_SNR_filtered_temp = aggregate(Data_table$SNR[Data_table$Position == P & Data_table$Passed_Filtering==1],by=list(Data_table$Individual_gene[Data_table$Position == P & Data_table$Passed_Filtering==1]),FUN=sd)$x
  
  Ratio_filtered = rbind(N_spot_filtered_temps/N_spot_raw_temps,1-N_spot_filtered_temps/N_spot_raw_temps)
  Ratio_filtered = 100 *Ratio_filtered
  
  Max_mean_SNR = max(max(Mean_SNR_raw_temp),max(Mean_SNR_filtered_temp))
  Max_sd_SNR = max(max(Sd_SNR_raw_temp),max(Sd_SNR_filtered_temp))
  Max_SNR = (Max_mean_SNR + Max_sd_SNR) * 1.2
  
  ##Plotting the number of spots
  
  N_max_spot = 1.2*max(N_spot_raw_temps)
  
  par(las=1,mfrow=c(6,1),mar=c(3,5,2,2),adj=0)
  plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
  text(x = 0.1,y=0.5,labels = paste("QC analysis of spots from Position",P),cex = 3)
  barplot(N_spot_raw_temps,ylim=c(0,N_max_spot),names.arg = List_individual_genes,
          main="Total number of spots identified",cex.names = 1.4,cex.axis=1.3)
  barplot(N_spot_filtered_temps,ylim=c(0,N_max_spot),names.arg = List_individual_genes,
          main="Total number of filtered spots identified",cex.names = 1.4,cex.axis=1.3)
  
  ##Checking for the effect of filtering
  barplot(Ratio_filtered,main="Percentage of spots passing the filtering",
          names.arg = List_individual_genes,cex.names = 1.4,cex.axis=1.3)
  
  ##Plotting the mean SNR for each gene 
  
  x = barplot(Mean_SNR_raw_temp,main="SNR of all spots",cex.axis=1.3,
          names.arg = List_individual_genes,ylim=c(0,Max_SNR),cex.names = 1.4)
  segments(x,Mean_SNR_raw_temp-Sd_SNR_raw_temp,x,Mean_SNR_raw_temp+Sd_SNR_raw_temp,lwd = 2)
  segments(x-0.2,Mean_SNR_raw_temp-Sd_SNR_raw_temp,x+0.2,Mean_SNR_raw_temp-Sd_SNR_raw_temp,lwd=2)
  segments(x-0.2,Mean_SNR_raw_temp+Sd_SNR_raw_temp,x+0.2,Mean_SNR_raw_temp+Sd_SNR_raw_temp,lwd=2)
  
  
  x=barplot(Mean_SNR_filtered_temp,main="SNR of filtered spots",cex.axis=1.3,
          names.arg = List_individual_genes,ylim=c(0,Max_SNR),cex.names = 1.4)
  segments(x,Mean_SNR_filtered_temp-Sd_SNR_filtered_temp,x,Mean_SNR_filtered_temp+Sd_SNR_filtered_temp,lwd = 2)
  segments(x-0.2,Mean_SNR_filtered_temp-Sd_SNR_filtered_temp,x+0.2,Mean_SNR_filtered_temp-Sd_SNR_filtered_temp,lwd=2)
  segments(x-0.2,Mean_SNR_filtered_temp+Sd_SNR_filtered_temp,x+0.2,Mean_SNR_filtered_temp+Sd_SNR_filtered_temp,lwd=2)
  
  ##Plotting the mean SNR for each gene 
  
}

dev.off()

