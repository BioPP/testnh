# Created on 12/04/2011 by jdutheil

########################################
# Profiling:
free.aic<-read.table("MantellidDataset.model_free_AIC.log.txt", header=TRUE, sep="\t")
join.aic<-read.table("MantellidDataset.model_join_AIC.log.txt", header=TRUE, sep="\t")


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

free.aic.prf<-read.profile("MantellidDataset.free_AIC.profile")
join.aic.prf<-read.profile("MantellidDataset.join_AIC.profile")

#Read partition tables:
brlen<-read.table("MantellidDataset.ml_h.brlen.txt", sep="\t", header=TRUE, row.names="Id")

free.aic.tbl<-read.table("MantellidDataset.partitions_record_free_AIC.txt", sep="\t")
join.aic.tbl<-read.table("MantellidDataset.partitions_record_join_AIC.txt", sep="\t")


plot.profile<-function(prf, iterations, omegas, brlen, tbl, ics) {
  bgcol<-rep(c(grey(0.8), grey(0.9)), length.out=length(iterations))
  layout(matrix(1:2, nrow=2), heights=c(5,3))
  par(mar=c(0,4.1,0.1,4.1), xaxt="n", bg="white")
  prf<-subset(prf, Iter %in% iterations)
  plot.new()
  plot.window(xlim=range(prf$Time), ylim=c(0,0.4))
  axis(2)
  mtext(expression(omega), 2, 3)
  totwd<-NA
  for (it in iterations) {
    if( !is.null(tbl)) {
      s<-split(tbl[, 1], tbl[, it])
      totalLength<-sapply(s, function(x) sum(brlen[x, "Branch.length"]))
    }
    time<-prf[prf$Iter == it, "Time"]
    rect(time[1], -0.5, time[length(time)], 1.5, col=bgcol[it], border=NA)
    for (omega in omegas) {
      lines(prf[prf$Iter == it, paste("YN98.omega", omega, sep="_")]~time, lwd=ifelse(is.null(tbl), 1, 1+totalLength[as.character(omega)]*0.5))
    }
    if (it == 2) {
      totwd<-ifelse(is.null(tbl), 1, 1+totalLength[as.character(omegas[1])]*0.5)
    }
    #abline(v=time[1], lty="dotted")
    text(x=mean(time), y=1.1, as.character(it))
  }
  lines(prf$Time, prf$YN98.omega, lwd=totwd)

  #Now plot information criteria:
  draw.local.min<-function(x, y, txt, col) {
    arrows(x0=x, y0=0.9, y1=y, col=col, lty="dotted")
    text(x=x, y=0.9, txt, col=col, pos=3)
  }
  draw.global.min<-function(x, y, txt, col) {
    arrows(x0=x, y0=0.9, y1=y, col=col, lwd=2)
    text(x=x, y=0.9, txt, col=col, pos=3)
  }
  scale<-function(x, r) { return ((x-r[1])/(r[2] - r[1])) }
  r<-function(x, r) c(floor(min(x)/r)*r, ceiling(max(x)/r)*r)
  r.aic=r(ics$AIC, 100)
  r.bic=r(ics$BIC, 100)
  sAIC=scale(ics$AIC, r.aic)
  sBIC=scale(ics$BIC, r.bic)
  par(mar=c(4.1,4.1,0.1,4.1), xaxt="s", bg="white")
  plot.new()
  plot.window(xlim=range(prf$Time), ylim=c(0, 1.1))
  axis(1)
  mtext("Execution time (seconds)", 1, 3)
  lAIC<-NULL
  lBIC<-NULL
  for (it in iterations) {
    x<-range(subset(prf, Iter == it, "Time"))
    lines(x=x, y=rep(sAIC[it],2), lwd=2, col=grey(0.5))
    lines(x=x, y=rep(sBIC[it],2), lwd=2)
    if (it > 1 && sAIC[it] < sAIC[it - 1] && sAIC[it] < sAIC[it + 1]) {
      draw.local.min(mean(x),  sAIC[it], ics$NbPartitions[it], grey(0.5))
      if (is.null(lAIC))
        lAIC<-list(x=mean(x), y=sAIC[it], n=ics$NbPartitions[it])
      else
        if (sAIC[it] < lAIC$y)
          lAIC<-list(x=mean(x), y=sAIC[it], n=ics$NbPartitions[it])
    }
    if (it > 1 && sBIC[it] < sBIC[it - 1] && sBIC[it] < sBIC[it + 1]) {
      draw.local.min(mean(x),  sBIC[it], ics$NbPartitions[it], "black")
      if (is.null(lBIC))
        lBIC<-list(x=mean(x), y=sBIC[it], n=ics$NbPartitions[it])
      else
        if (sBIC[it] < lBIC$y)
          lBIC<-list(x=mean(x), y=sBIC[it], n=ics$NbPartitions[it])
    }
  }
  draw.global.min(lAIC$x, lAIC$y, lAIC$n, grey(0.5))
  draw.global.min(lBIC$x, lBIC$y, lBIC$n, "black")
  
  labels<-seq(r.aic[1], r.aic[2], by=100)
  axis(2, at=(0:(length(labels)-1))/(length(labels)-1), labels, col=grey(0.5), col.axis=grey(0.5))
  mtext(side=2, "AIC", col=grey(0.5), line=3)
  labels<-seq(r.bic[1], r.bic[2], by=100)
  axis(4, at=(0:(length(labels)-1))/(length(labels)-1), labels, col="black", col.axis="black")
  mtext(side=4, "BIC", col="black", line=3)
}

plot.profile(free.aic.prf, 1:20, 1:21, brlen, free.aic.tbl, free.aic)
dev.print(tiff, file="Figures/FreeMantellid.tif", width=1000, height=500)

plot.profile(join.aic.prf, 1:22, 1:31, brlen, join.aic.tbl, join.aic)
dev.print(tiff, file="Figures/JoinMantellid.tif", width=1000, height=500)


