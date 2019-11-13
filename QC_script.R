library(spatstat)

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

Path = 'seqFISH_Temporary_file'
Data_table = read.delim(paste(Path,"Concatenated_spot_information.txt",sep = "/"))
colnames(Data_table) = c("Position_X","Position_Y","Position_Z",'Intensity','Passed_Filtering','Round','Position','Channel')
Data_table$Individual_gene = paste(Data_table$Round,Data_table$Channel,sep = "_")

Noise_matrix = read.delim(paste(Path,"Noise_value.txt",sep = "/"))
colnames(Noise_matrix) = c("Noise","Round","Channel","Position")

General_information = read.delim(paste(Path,"General_experiment_information.txt",sep = "/"))

N_rounds = General_information[1,1]
N_channels = General_information[1,2]
N_positions = General_information[1,3]

#First analysis : looking at each channel across all Positions.
#Looking at spot counting and quality
Distribution_unfiltered = table(Data_table$Individual_gene)
Distribution_filtered = table(Data_table$Individual_gene[Data_table$Passed_Filtering==1])

par(las=1,mfrow=c(N_rounds,1))
for (k in 1:N_rounds) {
  #Spot counts
  Distribution_unfiltered = table(factor(Data_table$Channel[Data_table$Round==k],levels = 1:N_channels))
  Distribution_filtered = table(factor(Data_table$Channel[Data_table$Passed_Filtering==1 & Data_table$Round==k],levels = 1:N_channels))
  barplot(Distribution_unfiltered,main="Distribution of all spots",ylab="Number of spots",xlab="Channel",cex.lab=1.3,cex.names =1.3)
  barplot(Distribution_filtered,main="Distribution of filtered spots",ylab="Number of spots",xlab="Channel",cex.lab=1.3,cex.names=1.3)
  
  #Looking at spot filtering
  Ratio_filtered = table(factor(Data_table$Channel[Data_table$Round==k],levels = 1:N_channels),factor(Data_table$Passed_Filtering[Data_table$Round==k],levels = 1:0))
  Ratio_filtered = Ratio_filtered/rowSums(Ratio_filtered)
  Ratio_filtered = t(100*Ratio_filtered)
  barplot(Ratio_filtered,main="Distribution of filtered spots",ylab="Proportion of high-quality spots (%)",xlab="Channel",cex.lab=1.3,cex.names=1.3,col=c("firebrick4","cornsilk"))
}

##Second analysis : study Signal to noise ratio

SNR_matrix = Noise_matrix[,2:4]
SNR_matrix$SNR_raw = 0
SNR_matrix$SNR_raw_sd = 0
SNR_matrix$SNR_filtered = 0
SNR_matrix$SNR_filtered_sd = 0
SNR_matrix$Merged_name = paste("R_",SNR_matrix$Round," ","Channel_",SNR_matrix$Channel,sep="")

for (k in 1:nrow(SNR_matrix)) {
  
  Selected_Channel = SNR_matrix$Channel[k]
  Selected_Round = SNR_matrix$Round[k]
  Selected_Position = SNR_matrix$Position[k]
  
  Noise_background = Noise_matrix[Noise_matrix$Round ==Selected_Round & Noise_matrix$Position ==Selected_Position & Noise_matrix$Channel ==Selected_Channel  ,1]
  Intensity = Data_table$Intensity[Data_table$Round ==Selected_Round & Data_table$Position ==Selected_Position & Data_table$Channel ==Selected_Channel]
  Kept = Data_table$Passed_Filtering[Data_table$Round ==Selected_Round & Data_table$Position ==Selected_Position & Data_table$Channel ==Selected_Channel]==1
  
  SNR_matrix$SNR_raw[k] = 10*log10(mean(Intensity)/Noise_background)
  SNR_matrix$SNR_filtered[k] = 10*log10(mean(Intensity[Kept])/Noise_background)
  
  SNR_matrix$SNR_raw_sd[k] = sd(10*log10((Intensity[Intensity>0])/Noise_background),na.rm = T)
  SNR_matrix$SNR_filtered_sd[k] = sd(10*log10((Intensity[Kept & Intensity>0])/Noise_background),na.rm = T)
}

Max_SNR = max(c(SNR_matrix$SNR_raw,SNR_matrix$SNR_filtered),na.rm = T) + max(c(SNR_matrix$SNR_raw_sd,SNR_matrix$SNR_filtered_sd),na.rm = T)

par(las=1)
barplot_center = barplot(SNR_matrix$SNR_raw,main="SNR of all spots",ylab="SNR (dB)",ylim = c(0,Max_SNR),cex.main=1.3,
        xlab="Channel",cex.lab=1.3,cex.names=1.3,col=c("grey"),names.arg = SNR_matrix$Merged_name,cex.axis=1.3)
arrows(barplot_center, SNR_matrix$SNR_raw-SNR_matrix$SNR_raw_sd, lwd=2,
       barplot_center, SNR_matrix$SNR_raw+SNR_matrix$SNR_raw_sd,angle=90,code=3)

barplot_center=barplot(SNR_matrix$SNR_filtered,main="SNR of filtered spots",ylab="SNR (dB)",ylim = c(0,Max_SNR),
        xlab="Channel",cex.lab=1.3,cex.names=1.3,col=c("orangered3"),names.arg = SNR_matrix$Merged_name,)
arrows(barplot_center, SNR_matrix$SNR_filtered-SNR_matrix$SNR_filtered_sd, lwd=2,
       barplot_center, SNR_matrix$SNR_filtered+SNR_matrix$SNR_filtered_sd,angle=90,code=3)


##Studiyng spatial information : one position by position

#First measure :  L function for the unfiltered and filtered point position by position 
for (P in 1:N_positions) {
  
  window_temp = owin(c(0,max(Data_table[Data_table$Position==P,1])),c(0,max(Data_table[Data_table$Position==P,2])))
  rmax_temp = min(window_temp$xrange[2]/4,window_temp$yrange[2]/4)
  
  ##Generating L plot for each gene individualy
  for (Selected_channels in unique(Data_table$Individual_gene)) {
    x= strsplit(Selected_channels,split = "_")
    R = x[[1]][1]
    Channel = x[[1]][2]
    
    #Selecting the corresponding spots 
    Selected_points = Data_table$Position==P & Data_table$Individual_gene==Selected_channels
    
    X =  Data_table[Selected_points,1:2]
    Passed_filtering = Data_table[Selected_points,5]==1
    Intensity = Data_table[Selected_points,4]
    
    
    #Creating the ppp object
    position_unfiltered = ppp(X[,1],(X[,2]),window=window_temp)
    position_filtered = ppp(X[Passed_filtering,1],(X[Passed_filtering,2]),window=window_temp)
    
    Lest_unfiltered = Lest(position_unfiltered,correction="border",rmax = rmax_temp)
    Lest_filtered = Lest(position_filtered,correction="border",rmax = rmax_temp)
    par(las=1)
    plot(NULL,xaxs="i",yaxs="i",xlim=c(0,min(max(Lest_filtered$r),max(Lest_unfiltered$r))),ylim = c(0,rmax_temp*1.5),
         xlab="r",ylab="L(r)",cex.lab=1.4,cex.axis=1.4,main = paste("L function for Round",R,"Channel",Channel))
    points(Lest_unfiltered$r,Lest_unfiltered$border,type="l",lwd=2,lty=2,col="blue")
    points(Lest_filtered$r,Lest_filtered$border,type="l",lwd=2,lty=4,col="red")
    points(Lest_filtered$r,Lest_filtered$r,type="l",lwd=2,lty=2,col="grey")
    
  }
  ##Generating cross-type L function plot
  Selected_points = Data_table$Position==P
  X =  Data_table[Selected_points,1:2]
  Source = Data_table$Individual_gene[Selected_points]
  Passed_filtering = Data_table[Selected_points,5]==1
  
  position_unfiltered = ppp(X[,1],(X[,2]),window=window_temp,marks = data.frame(Source=Source))
  position_filtered = ppp(X[Passed_filtering,1],(X[Passed_filtering,2]),window=window_temp,marks = data.frame(Source=Source[Passed_filtering]))
  List_L_cross_unfiltered = NULL
  List_L_cross_filtered = NULL
  
  n=1
  for (i in unique(Data_table$Individual_gene)) {
    for (j in unique(unique(Data_table$Individual_gene))) {
      
      L_cross_unfiltered_temp = Lcross(position_unfiltered,i = i,j=j,rmax = rmax_temp,correction="border")
      L_cross_filtered_temp = Lcross(position_filtered,i = i,j=j,rmax = rmax_temp,correction="border")
      List_L_cross_unfiltered[[n]] = L_cross_unfiltered_temp
      List_L_cross_filtered[[n]] = L_cross_filtered_temp
      
      n = n+1
    }
  }
  
}


