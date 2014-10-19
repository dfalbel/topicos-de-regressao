#---------------------------------------------------------------#
# Para rodar este programa  deixe no objeto fit.model a saída 
# do  ajuste  da  regressão com  erro normal. Deixe  os dados 
# disponíveis  através do comando attach(...). Depois  use  o 
# comando source(...) no S-Plus ou R para executar o programa. 
# A sequência de comandos é a seguinte:
#
#        > fit.model <- ajuste
#        > attach(dados)
#        > source("diag_norm")
#
# A saída terá quatro gráficos: de pontos de alavanca, de pontos
# influentes  e  dois de resíduos. Para identificar os pontos
# que  mais  se destacam usar o comando identify(...). Se por
# exemplo se destacam três pontos no 
# plot(fitted(fit.model),h,...), após esse comando coloque
#     
#        > identify(fitted(fit.model),h,n=3)
#
# O mesmo pode ser feito nos demais gráficos. Nos gráficos de 
# resíduos foram traçados os limites ylim=c(a-1,b+1), onde a
# é o menor valor e b o maior valor para o resíduo.Mude esses 
# limites  se  necessário.
#---------------------------------------------------------------#

diag_norm <- function(fit.model, dados){
  
  X <- model.matrix(fit.model)
  n <- nrow(X)
  p <- ncol(X)
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  h <- diag(H)
  lms <- summary(fit.model)
  s <- lms$sigma
  r <- resid(lms)
  ts <- r/(s*sqrt(1-h))
  di <- (1/p)*(h/(1-h))*(ts^2)
  si <- lm.influence(fit.model)$sigma
  tsi <- r/(si*sqrt(1-h))
  a <- max(tsi)
  b <- min(tsi)
  indice <- 1:length(h)
  #par(mfrow=c(2,2))
  cut = 2*p/n
  p1 <- qplot(y=h, x=indice, geom = "point") + geom_hline(aes(yintercept = 2*p/n)) + ylim(0,1) + geom_text(hjust=0, vjust=0, aes(x = indice[h>cut], y = h[h>cut], label = indice[h>cut])) + xlab("Índice")
  #identify(h, n=1)
  #title(sub="(a)")
  #
  p2 <- qplot(x = indice, y = di, geom = "point") + xlab("Índice") + ylab("Distância de Cook")
  #identify(di, n=2)
  #
  
  p3 <- qplot(x = indice, y = tsi, geom = "point") + geom_hline(aes(yintercept = c(-2,2))) + xlab("Índice") + ylab("Resíduo Padronizado") #+ geom_text(hjust = 0, vjust = 0, aes(x = indice[tsi <(-2) | tsi >2], y = tsi[tsi>2 | tsi < -2]))
  
  
  p4 <- qplot(x = fitted(fit.model), y = tsi, geom = "point") + xlab("Valor Ajustado") + ylab("Resíduo Padronizado") + geom_hline(aes(yintercept = c(-2,2)))
  
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  
}


