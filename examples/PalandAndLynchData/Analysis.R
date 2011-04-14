# Created on 12/04/2011 by jdutheil

########################################
# Profiling:
free.aic<-read.table("PalandAndLynch.model_free_AIC.log.txt", header=TRUE, sep="\t")
join.aic<-read.table("PalandAndLynch.model_join_AIC.log.txt", header=TRUE, sep="\t")


read.profile<-function(file) {
  path<-file;
  i<-1
  time<-0;
  while (file.exists(path))
  {
    cat(path, "\n")
    tbl<-read.table(path, header=TRUE, sep="\t")
    tbl$Iter<-i
    tbl$Time<-tbl$Time+time
    time<-tbl$Time[length(tbl$Time)]
    #merge:
    if (i == 1)
      prf<-tbl
    else
      prf<-merge(prf, tbl, all=TRUE)

    #Get the next one:
    path<-paste(file, i, sep="")
    i<-i+1
  }
  return(prf[order(prf$Time),])
}

join.aic.prf<-read.profile("PalandAndLynch.join_AIC.profile")
free.aic.prf<-read.profile("PalandAndLynch.free_AIC.profile")
join.bic.prf<-read.profile("PalandAndLynch.join_BIC.profile")
free.bic.prf<-read.profile("PalandAndLynch.free_BIC.profile")


plot.profile<-function(prf, iterations, omegas) {
  plot(YN98.omega~Time, prf, type="l", ylim=c(0,1))
  for (it in iterations) {
    for (omega in omegas) {
      lines(prf[prf$Iter == it, paste("YN98.omega", omega, sep="_")]~prf[prf$Iter == it, "Time"], col=it+1)
    }
  }
}

plot.profile(join.aic.prf, 1:7, 1:10)
plot.profile(free.aic.prf, 1:6, 1:6)


