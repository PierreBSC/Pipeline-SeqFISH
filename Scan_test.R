suppressMessages(library(spatstat))


args <- commandArgs(trailingOnly = T)

file_directory = as.character(args[1])
Estimated_cell_size = as.numeric(args[2])


n_channel = length(list.files(file_directory)) / 2


##Loading the data


position_list = c()

intensity_list = c()

print(1)
for (k in 1:n_channel) {

  

  u = read.table(paste(file_directory,"Spot_position_channel_",k,".txt",sep = ""),header = T,sep="\t")

  colnames(u) = c("X","Y","Z")


    

  v = read.table(paste(file_directory,"Spot_intensity_channel_",k,".txt",sep = ""),header = T,sep="\t")

  v = v[,1]

  

  checked_points = (complete.cases(v) & complete.cases(u))

  

  intensity_list[[k]] = v[checked_points]

  position_list[[k]] = u[checked_points,]

  

}
print(2)


##Performing the statistical analysis : spatial scan test from Kulldorf's paper (1997)


analysis_granularity = 30;


LR_table = matrix(0,nrow = analysis_granularity,ncol = n_channel)

list_threshold = c()

List_autolag =  c()

for (k in 1:n_channel) {

  
  #Creating the  intensity threshold list
  intensity = intensity_list[[k]]

  max_intensity = quantile(intensity,probs = 0.999)
  min_intensity = quantile(intensity,probs = 0.01)
  
  intensity_range = seq(min_intensity,max_intensity,length.out = analysis_granularity)

  position = position_list[[k]]
  position_base = ppp(position[,1],position[,2],window = owin(c(min(position[,1]),max(position[,1])),c(min(position[,2]),max(position[,]))))

  if (nrow(position)>15) {
    

  for (i in 1:analysis_granularity) {

    if (length(position$X[intensity>=intensity_range[i]])>1) {

        X = ppp(position[,1],position[,2],window = owin(c(min(position[,1]),max(position[,1])),c(min(position[,2]),max(position[,]))),
                marks = factor(ifelse(intensity>=intensity_range[i],"Positive","Negative")))

          M = scanLRTS(X,r = Estimated_cell_size,alternative = "greater",method = "binomial",case = "Positive")
        
        LR_table[i,k] = -max(M)

    }
    if (length(position$X[intensity>=intensity_range[i]])<=1) {
        LR_table[i,k] = 0;

    }

  }


  }
  
  x =  LR_table[,k]
  Autolog = (cor(x[-1],x[-length(x)]))^2  
  List_autolag = c(List_autolag,Autolog)
  
  list_threshold[[k]] = data.frame(Intensity = intensity_range, Score = LR_table[,k] )
}

print(3)
for (k in 1:n_channel)  {

  u = list_threshold[[k]]
  write.table(x = u,file = paste(file_directory,"Table_channel",k,".txt",sep = ""),sep = "\t",quote = F,row.names = F)
}

write.table(x = List_autolag,file = paste(file_directory,"Table_autolag_score.txt",sep = ""),sep = "\t",quote = F,row.names = F,col.names = F)
