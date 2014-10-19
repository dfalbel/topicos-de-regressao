#------------------------------------------------------------#
# Para rodar este programa deixe no objeto fit.model a saída      
# do ajuste da regressão com erros gama e ligação log.  Deixe     
# também os dados disponíveis através do comando attach(...).    
# Depois use o comando source(...) no R ou S-Plus para        
# executar o programa. A sequência de comandos é a seguinte:               
#
#       fit.model <- ajuste
#       attach(dados)
#       source("envel_gama")                                      
#
# A saída será o gráfico de envelope para o resíduo componente 
# do desvio padronizado. Para colocar um título no gráfico use
# o comando title("..."). Para  usar  outras  ligações mudar 
# no programa abaixo o termo family=Gamma(link=log) para 
# family=Gamma no caso de ligação  recíproca ou por 
# family= Gamma(link=identity) no caso de ligação identidade.
#------------------------------------------------------------#
envel_gama <- function(fit.model, dados){
  X <- model.matrix(fit.model, data = dados)
  n <- nrow(X)
  p <- ncol(X)
  w <- fit.model$weights
  W <- diag(w)
  H <- solve(t(X)%*%W%*%X)
  H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h <- diag(H)
  ro <- resid(fit.model,type="response")
  fi <- (n-p)/sum((ro/(fitted(fit.model)))^ 2)
  td <- resid(fit.model,type="deviance")*sqrt(fi/(1-h))
  #
  e <- matrix(0,n,100)
  #
  for(i in 1:100){
    resp <- rgamma(n,fi)
    resp <- (fitted(fit.model)/fi)*resp
    fit <- glm(resp ~ X, family=Gamma(link=log))
    w <- fit$weights
    W <- diag(w)
    H <- solve(t(X)%*%W%*%X)
    H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
    h <- diag(H)
    ro <- resid(fit,type="response")
    phi <- (n-p)/sum((ro/(fitted(fit)))^ 2)
    e[,i] <- sort(resid(fit,type="deviance")*sqrt(phi/(1-h)))}
  #
  e1 <- numeric(n)
  e2 <- numeric(n)
  #
  for(i in 1:n){
    eo <- sort(e[i,])
    e1[i] <- (eo[2]+eo[3])/2
    e2[i] <- (eo[97]+eo[98])/2}
  #
  med <- apply(e,1,mean)
  faixa <- range(td,e1,e2)
  df <- data.frame(li = e1, ls = e2, ye = med, res = sort(td), q = qnorm(seq(0.01,0.99, length.out = length(e1))))
  
  g <- ggplot(df, aes(x = q)) + geom_ribbon(aes(ymin = li, ymax = ls), alpha = 0.4) + geom_line(aes(y = ye), linetype = "dashed") + geom_point(aes(y = res))
  
  return(g)
}
#------------------------------------------------------------#                      