# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

pmh = function(data = j_sight,lambda=0){

  if(1==2){
    lambda = 0
  }

  log0 = function(x){
    if(is.infinite(log(x))){
      res = 0
    }else{
      res = log(x)
    }
    return(res)
  }

  p_ij = dat/sum(dat)
  n.row = n.col = nrow(p_ij)
  dp_ij = matrix(NA,n.row,n.col)


  A = (2^lambda)/((2^lambda)-1)
  A2 = 1/log0(2)
  p_id = c()
  p_di = c()
  p1rat = p2rat = c()
  him = enm = lenm = pistar = c()

  alln = sum(dat)

  #尺度＝０
  #p_ij = cbind(c(1/18,1/18,1/18),c(1/18,3/18,1/18),c(1/18,2/18,7/18))
  #尺度＝１
  #p_ij = cbind(c(0,7/18,0),c(0,0,0),c(0,11/18,0))
  #rmultinom(1,100,c(3/9,2/9,1/18,1/18,1/18,1/18,1/9,1/18,1/18))

  toriihim = toriienm = 1
  for(i in 1:n.row){
    p_id[i] = sum(p_ij[i,])
    p_di[i] = sum(p_ij[,i])
    p1rat[i] = p_id[i]/(p_id[i]+p_di[i])
    p2rat[i] = p_di[i]/(p_id[i]+p_di[i])
    him[i] = 1 - A*(1-((p1rat[i])^(lambda+1))-((p2rat[i])^(lambda+1)))
    enm[i] = 1 + A2*(p1rat[i]*log0(p1rat[i])+p2rat[i]*log0(p2rat[i]))
    pistar[i] = (p_id[i]+p_di[i])/2
    toriihim = toriihim * (him[i]^pistar[i])
    toriienm = toriienm * (enm[i]^pistar[i])
  }

  #toriihimはlambda!=0の提案尺度に入れたやつ！
  #toriienmはlambda=0の提案尺度に入れたやつ！


  #微分したやつに入れてる！
  if(lambda != 0){
    for(i in 1:n.row){
      for(j in 1:n.col){
        if(i == j){
          dp_ij[i,j] = toriihim * (log0(him[i])+A*(lambda+1)*(p_id[i]^lambda-p_di[i]^lambda)*(p_di[i]-p_id[i])/(2*him[i]*(p_id[i]+p_di[i])^(lambda+1)))
        }else{
          dp_ij[i,j]= (1/2) * toriihim * (log0(him[i])+A*(lambda+1)*(p_id[i]^lambda-p_di[i]^lambda)*p_di[i]/(him[i]*(p_id[i]+p_di[i])^(lambda+1)) + log0(him[j])+A*(lambda+1)*(p_di[j]^lambda-p_id[j]^lambda)*p_id[j]/(him[j]*(p_id[j]+p_di[j])^(lambda+1)))
        }
      }
    }

    Sigma1 = Sigma22 = 0
    for(i in 1:n.row){
      for(j in 1:n.col){
        Sigma1 = Sigma1 + p_ij[i,j]*(dp_ij[i,j]^2)
        Sigma22 = Sigma22 + p_ij[i,j]*dp_ij[i,j]
      }
    }
    Sigma2 = Sigma22^2

    Sigma = Sigma1 - Sigma2

    CI.low = toriihim - qnorm(1-0.05/2)*sqrt(Sigma/alln)
    CI.upper = toriihim + qnorm(1-0.05/2)*sqrt(Sigma/alln)

  }


  if(lambda == 0){
    for(i in 1:n.row){
      for(j in 1:n.col){
        if(i == j){
          dp_ij[i,j] = toriienm * (log0(enm[i])+A2*(p_di[i]-p_id[i])*log0(p_id[i]/p_di[i])/(2*(p_id[i]+p_di[i])*enm[i]))
        }else{
          dp_ij[i,j] = toriienm *(1/2) * (log0(enm[i])+A2*p_di[i]*log0(p_id[i]/p_di[i])/((p_id[i]+p_di[i])*enm[i]) + log0(enm[j])-A2*p_id[j]*log0(p_id[j]/p_di[j])/((p_id[j]+p_di[j])*enm[j]))
        }
      }
    }


    Sigma1 = Sigma22 = 0
    for(i in 1:n.row){
      for(j in 1:n.col){
        Sigma1 = Sigma1 + p_ij[i,j]*(dp_ij[i,j]^2)
        Sigma22 = Sigma22 + p_ij[i,j]*dp_ij[i,j]
      }
    }

    Sigma2 = Sigma22^2

    Sigma = Sigma1 - Sigma2


    CI.low = toriienm - qnorm(1-0.05/2)*sqrt(Sigma/alln)
    CI.upper = toriienm + qnorm(1-0.05/2)*sqrt(Sigma/alln)

    toriihim = toriienm
  }


  #print(paste(c("lambda",lambda),collapse="="))
  #print("隔たり")
  #print(toriihim)
  #print("推定分散")
  #print(Sigma)
  #print("信頼区間")
  #print(c(CI.low ,CI.upper))

  return(list(dep = toriihim,vari = Sigma,CI = c(CI.low ,CI.upper)))

}

