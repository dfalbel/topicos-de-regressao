#------------------------------------------------------------#
# Para rodar este programa  deixe no objeto fit.model a saída 
# do  ajuste  da  regressão com  erros gama.   Deixe  os dados 
# disponíveis  através do comando attach(...). Depois  use  o 
# comando source(...) no S-Plus ou R para executar o programa. 
# A sequência de comandos é a seguinte:
#
#        > fit.model <- ajuste
#        > attach(dados)
#        > source("diag_gama")
#
# A saída terá quatro gráficos: de pontos de alacanca, 
# de pontos influentes  e  dois de resíduos. Para 
# identificar os pontos que  mais  se destacam usar o 
# comando identify(...). Se por exemplo se destacam
# três pontos no plot(fitted(fit.model),h,...), 
# após esse comando coloque
#     
#        > identify(fitted(fit.model),h,n=3)
#
# O mesmo pode ser feito nos demais gráficos. Nos gráficos de 
# resíduos foram colocados os limites ylim=c(a-1,b+1), 
# em que a é o menor valor e b o maior valor para o resíduo. 
# Este programa usa a library MASS para estimar o parâmetro
# fi da gama que  estará guardado no objeto fi.
#------------------------------------------------------------#

diag_gama <- function(fit.model, dados){
  
  X <- model.matrix(fit.model, data = dados)
  n <- nrow(X)
  p <- ncol(X)
  w <- fit.model$weights
  W <- diag(w)
  H <- solve(t(X)%*%W%*%X)
  H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h <- diag(H)
  library(MASS)
  fi <- gamma.shape(fit.model)$alpha
  ts <- resid(fit.model,type="pearson")*sqrt(fi/(1-h))
  td <- resid(fit.model,type="deviance")*sqrt(fi/(1-h))
  di <- (h/(1-h))*(ts^2)
  par(mfrow=c(2,2))
  a <- max(td)
  b <- min(td)
  indice <- 1:nrow(dados)
  
  p1 <- qplot(x = fitted(fit.model), y = h, geom = "point") + 
    xlab("Valor Ajustado") + ylab("Medida h")
  #plot(fitted(fit.model),h,xlab="Valor Ajustado", ylab="Medida h", pch=16)
  #identify(fitted(fit.model), h, n=1)
  #
  p2 <- qplot(x = indice, y = di, geom = "point") + 
    xlab("Índice") + ylab("Distância de Cook")
  #plot(di,xlab="Índice", ylab="Distância de Cook", pch=16)
  #identify(di, n=2)
  #
  p3 <- qplot(x = fitted(fit.model), y = td, geom = "point") + 
    xlab("Valor Ajustado") + ylab("Componente do Desvio") +
    ylim(b-1, a-1) + geom_hline(aes(yintercept = c(-2,2)))
  #plot(fitted(fit.model),td,xlab="Valor Ajustado", ylab="Componente do Desvio",
  #     ylim=c(b-1,a+1),pch=16)
  #abline(2,0,lty=2)
  #abline(-2,0,lty=2)
  #identify(fitted(fit.model),td, n=1)
  #
  w <- fit.model$weights
  eta <- predict(fit.model)
  z <- eta + resid(fit.model, type="pearson")/sqrt(w)
  
  p4 <- qplot(x = predict(fit.model), y = z, geom = "point") +
    xlab("Preditor Linear") + ylab("Variável z")
  
  #plot(predict(fit.model),z,xlab="Preditor Linear", 
  #     ylab="Variável z", pch=16)
  #par(mfrow=c(1,1))
  
  p <- grid.arrange(p1,p2,p3,p4, nrow = 2)
}
#------------------------------------------------------------#
