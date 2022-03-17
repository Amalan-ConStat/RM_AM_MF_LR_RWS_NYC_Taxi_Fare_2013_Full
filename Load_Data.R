load("NYC_Full_Data.RData")

NYC_fare_data_Clean$fare_T<-scales::rescale(NYC_fare_data_Clean$fare_T,to=c(0,1))

Original_Data <- NYC_fare_data_Clean[,-6]

no_of_Variables<-ncol(Original_Data[,-1])

view(dfSummary(Original_Data))

unique(NYC_fare_data_Clean[,c(2:5)])

Original_Data<-cbind(Original_Data[,1],1,Original_Data[,-1],Original_Data[,-c(1:5)]^2)
colnames(Original_Data)<-c("Y","X0",paste0("X",1:no_of_Variables),
                           paste0("X",5:no_of_Variables,"^2"))

glm(Y~.-1,data=Original_Data,family = "poisson")->Full_Model

stepAIC(Full_Model)

glm(Y~.-1,data=Original_Data,family = "quasipoisson")->quasi_Full_Model

summary(quasi_Full_Model)

glm(Y~X0+X1+X2+X3+X4+X5+X6-1,data=Original_Data,family = "poisson")->Covariate_Model

stepAIC(Covariate_Model)

glm(Y~X0+X1+X2+X3+X4+X5+X6-1,data=Original_Data,family = "quasipoisson")->quasi_Covariate_Model

Full_Model$aic
Covariate_Model$aic

corrr::correlate(Original_Data[,c(7:10)])

Squared_Terms<-paste0("X",5:no_of_Variables,"^2")
term_no <- 2

All_Models <- list(c("X0",paste0("X",1:no_of_Variables)))

for (i in 1:no_of_Variables)
  {
  x <- as.vector(combn(Squared_Terms,i,simplify = FALSE))
  for(j in 1:length(x))
    {
    All_Models[[term_no]] <- c("X0",paste0("X",1:no_of_Variables),x[[j]])
    term_no <- term_no+1
    }
  }

Subsample_Size<-2000;
r0<-Nc_size<-100; 
Replicates <- 1000; 
Choice <-100*seq(2,20,1)

save(list = c("Original_Data","Subsample_Size","Replicates","Nc_size","r0","All_Models","Choice"),
     file=here("Identical_r0","Generate_Big_Data","Scaled.RData"))

save(list = c("Original_Data","Subsample_Size","Replicates","Nc_size","r0","All_Models","Choice"),
     file=here("Non_Identical_r0","Generate_Big_Data","Scaled.RData"))

rm(list = ls())

