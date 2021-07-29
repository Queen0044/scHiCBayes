#' performance
#' 
#' This function summarizes the accuracy of simulated data.
#'
#' @param observed Observed observed cells matrix, with each column being the upper triangular of a observed cell.
#' @param expected Underline true count of the simulated data.
#' @param imputed Imputed matrix from MCMCImpute.
#'
#' @return A list of accuracy measurements with mean and standard error.
#' @export
#'
#' @examples
#' data("K562_1_true")
#' options(digits = 2)
#' performance(observed=K562_T1_7k, expected=K562_1_true, imputed=T1_7k_imp)

performance <- function(observed, expected, imputed){

  PTSZ=list()
  PTDO=list()
  #MSE=list()
  CIEZ=list()
  CIEA=list()
  AEOZ=data.frame()
  AEOA=data.frame()
  
  for(j in 1:ncol(observed)){
    index0 <- observed[,j]==0
    index1 <- observed[,j]==1
    index2 <- observed[,j]>=2
    
    ## calibration
    m1 <- lm( imputed[,j][observed[,j]>0]~observed[,j][observed[,j]>0])
    imputeduted <- imputed[,j][imputed[,j]>0 & observed[,j]==0]
    obs2 <- (imputeduted-m1$coefficients[1])/(m1$coefficients[2])
    imputed[,j][imputed[,j]>0 & observed[,j]==0] <- obs2
    
    indexnotrue0=(expected[,j]!=0)
    
    predictv=imputed[,j][index0]
    Truevalue=expected[,j][index0]
    
    PTSZ[[j]]=sum(Truevalue==0 & predictv==0)/sum(Truevalue==0)
    PTDO[[j]]=sum(Truevalue>0 & predictv>0)/sum(Truevalue>0)
    
    #MSE[[j]]=mean((imputed[,j]-expected[,j])^2)
    
    AEOZ=rbind(AEOZ, data.frame(AE=(abs(predictv-Truevalue)),Sample=j))
    AEOA=rbind(AEOA, data.frame(AE=abs(imputed[,j]-expected[,j]),Sample=j))
    
    CIEZ[[j]]=cor(Truevalue,predictv)
    CIEA[[j]]=cor(imputed[,j],observed[,j])
  }
  
  PTSZ=unlist(PTSZ)
  PTDO=unlist(PTDO)
  #MSE=unlist(MSE)
  CIEZ=unlist(CIEZ)
  CIEA=unlist(CIEA)
  
  summa_mean=data.frame(PTSZ=mean(PTSZ),PTDO=mean(PTDO),
                        AEOZ=mean(AEOZ$AE),AEOA=mean(AEOA$AE),
                        CIEZ=mean(CIEZ),CIEA=mean(CIEA))
  summa_se=data.frame(PTSZ=sd(PTSZ),PTDO=sd(PTDO),
                      AEOZ=sd(AEOZ$AE),AEOA=sd(AEOA$AE),
                      CIEZ=sd(CIEZ),CIEA=sd(CIEA))
  rownames(summa_mean)=NULL
  rownames(summa_se)=NULL
  
  mylist <- list(summa_mean,summa_se)
  names(mylist) <- c("summary_mean","summary_se")
  return(mylist)
}



