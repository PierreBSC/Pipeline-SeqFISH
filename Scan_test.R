library(spatstat)


args <- commandArgs(trailingOnly = T)

file_directory = as.character(args[1])


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


for (k in 1:n_channel) {

  
  #Creating the  intensity threshold list
  intensity = intensity_list[[k]]
  intensity_range = quantile(intensity,probs = seq(from = 0,to = 0.99,length.out =analysis_granularity ))

  position = position_list[[k]]
  position_base = ppp(position[,1],position[,2],window = owin(c(0,2000),c(0,2000)))
  #position_base = ppp(position[,1],position[,2])

  if (nrow(position)>15) {
    
  base_density = density.ppp(position_base,positive=TRUE,diggle =  T)

  

  for (i in 1:analysis_granularity) {

    if (length(position$X[intensity>=intensity_range[i]])>1) {

        X <- spatstat::ppp(position$X[intensity>=intensity_range[i]], position$Y[intensity>=intensity_range[i]], window=spatstat::owin(c(0,2000),c(0,2000)))

        test = spatstat::scan.test(X,r = c(50,100,150),alternative = "greater",verbose = F ,method = "poisson",baseline = base_density,nsim = 2)

        LR_table[i,k] = test$statistic
    }
    if (length(position$X[intensity>=intensity_range[i]])<=1) {
        LR_table[i,k] = 0;

    }


    end
  }


  }





  

  list_threshold[[k]] = data.frame(Intensity = intensity_range, Score = LR_table[,k] )

}

print(3)
for (k in 1:n_channel)  {

  u = list_threshold[[k]]

  write.table(x = u,file = paste(file_directory,"Table_channel",k,".txt",sep = ""),sep = "\t",quote = F,row.names = F)

}