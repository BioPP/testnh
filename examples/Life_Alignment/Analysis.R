# Created on 21/02/2011 by jdutheil

#Patterns of AIC/BIC:

free.aic<-read.table("Life_Alignment.model_free_AIC.log.txt", header=TRUE, sep="\t")
#free.bic<-read.table("Life_Alignment.model_free_BIC.log.txt", header=TRUE, sep="\t")
join.aic<-read.table("Life_Alignment.model_join_AIC.log.txt", header=TRUE, sep="\t")
#join.bic<-read.table("Life_Alignment.model_join_BIC.log.txt", header=TRUE, sep="\t")

plot.ic<-function(data, main) {
  draw.local.min<-function(x, col) {
    m<-x[1]
    p<-1
    mi<-0
    for (i in 2:length(x)) {
      if (x[i] < m) {
        m<-x[i]
        mi<-i
        p<-0
      } else { 
        p<-p+1
        if (p == 3) {
          arrows(x0=data[mi, "NbPartitions"], y0=0.9, y1=x[mi], col=col)
          text(x=data[mi, "NbPartitions"], y=0.9, data[mi, "NbPartitions"], col=col, pos=3)
        }
      }
    }
  }
  par(mar=c(4,4,4,4)+0.1)
  scale<-function(x, r) { return ((x-r[1])/(r[2] - r[1])) }
  r<-function(x, r) c(floor(min(x)/r)*r, ceiling(max(x)/r)*r)
  r.aic=r(data$AIC, 100)
  r.bic=r(data$BIC, 100)
  plot(scale(AIC, r.aic)~NbPartitions, data, type="b", pch=19, col="blue", xlab="Number of partitions", main=main, yaxt="n", ylab=NA, ylim=c(0,1))
  labels<-seq(r.aic[1], r.aic[2], by=100)
  axis(2, at=(0:(length(labels)-1))/(length(labels)-1), labels, col="blue", col.axis="blue")
  mtext(side=2, "AIC", col="blue", line=3)
  draw.local.min(scale(data$AIC, r.aic), "blue")
  lines(scale(BIC, r.bic)~NbPartitions, data, type="b", pch=19, col="orange")
  labels<-seq(r.bic[1], r.bic[2], by=100)
  axis(4, at=(0:(length(labels)-1))/(length(labels)-1), labels, col="orange", col.axis="orange")
  mtext(side=4, "BIC", col="orange", line=3)
  draw.local.min(scale(data$BIC, r.bic), "orange")
}


layout(matrix(1:2, nrow=2, byrow=TRUE))
plot.ic(free.aic, "Free")
plot.ic(join.aic, "Join")

dev.print(pdf, file="AIC_BIC.pdf", width=8, height=12)

#Parameter estimates:
plot.gc<-function(data, ...) {
  data<-data[-nrow(data),] #Remove last line
  l<-data$BrLen2
  q<-quantile(l, probs=(0:10)/10, na.rm=TRUE)
  q[1]<-q[1]-1
  b<-cut(l, breaks=q)
  col<-grey(1-as.numeric(b)/10)
  plot(T92.theta~NHML, data, ylim=c(0,1), xlim=c(0,1), xlab="NHML", ylab="Partitions", pch=19, col=col, ...)
  abline(0,1, col="red")
}

join.bic.gc<-read.table("Life_Alignment.ml_nh_join_BIC.parameters.with_nhml.csv", header=TRUE, sep=",")
join.aic.gc<-read.table("Life_Alignment.ml_nh_join_AIC.parameters.with_nhml.csv", header=TRUE, sep=",")
free.bic.gc<-read.table("Life_Alignment.ml_nh_free_BIC.parameters.with_nhml.csv", header=TRUE, sep=",")
free.aic.gc<-read.table("Life_Alignment.ml_nh_free_AIC.parameters.with_nhml.csv", header=TRUE, sep=",")

layout(matrix(1:4, nrow=2, byrow=TRUE))
plot.gc(join.bic.gc, main="Join, BIC")
plot.gc(join.aic.gc, main="Join, AIC")
plot.gc(free.bic.gc, main="Free, BIC")
plot.gc(free.aic.gc, main="Free, AIC")

dev.print(pdf, file="GC.pdf", width=12, height=12)

########################################
# Profiling:

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

join.aic.prf<-read.profile("Life_Alignment.join_AIC.profile")
free.aic.prf<-read.profile("Life_Alignment.free_AIC.profile")
join.bic.prf<-read.profile("Life_Alignment.join_BIC.profile")
free.bic.prf<-read.profile("Life_Alignment.free_BIC.profile")

#Read partition tables:
brlen<-read.table("Life_Alignment.ml_h.rooted.brlen.txt", sep="\t", header=TRUE, row.names="Id")

join.aic.tbl<-read.table("Life_Alignment.partitions_record_join_AIC.txt", sep="\t")
free.aic.tbl<-read.table("Life_Alignment.partitions_record_free_AIC.txt", sep="\t")
join.bic.tbl<-read.table("Life_Alignment.partitions_record_join_BIC.txt", sep="\t")
free.bic.tbl<-read.table("Life_Alignment.partitions_record_free_BIC.txt", sep="\t")

plot.profile<-function(prf, iterations, thetas, brlen, tbl, ics) {
  layout(matrix(1:2, nrow=2), heights=c(5,3))
  par(mar=c(0,4.1,0.1,4.1), xaxt="n")
  plot(T92.theta~Time, prf, type="l", ylim=c(0,1))
  for (it in iterations) {
    if( !is.null(tbl)) {
      s<-split(tbl[, 1], tbl[, it])
      totalLength<-sapply(s, function(x) sum(brlen[x, "Branch.length"]))
    }
    for (theta in thetas) {
      lines(prf[prf$Iter == it, paste("T92.theta", theta, sep="_")]~prf[prf$Iter == it, "Time"], lwd=ifelse(is.null(tbl), 1, 1+totalLength[as.character(theta)]*2))
    }
    abline(v=prf[prf$Iter == it, "Time"][1], lty="dotted")
    text(prf[prf$Iter == it, "Time"][1], 1, as.character(it))
  }

  #Now plot information criteria:
  draw.local.min<-function(x, txt, col) {
    arrows(x0=x, y0=0.9, y1=x, col=col)
    text(x=x, y=0.9, txt, col=col, pos=3)
  }
  scale<-function(x, r) { return ((x-r[1])/(r[2] - r[1])) }
  r<-function(x, r) c(floor(min(x)/r)*r, ceiling(max(x)/r)*r)
  r.aic=r(ics$AIC, 100)
  r.bic=r(ics$BIC, 100)
  sAIC=scale(ics$AIC, r.aic)
  sBIC=scale(ics$BIC, r.bic)
  par(mar=c(4.1,4.1,1,4.1), xaxt="s")
  plot.new()
  plot.window(xlim=range(prf$Time), ylim=c(0, 1))
  axis(1)
  for (it in iterations) {
    x<-range(subset(prf, Iter == it, "Time"))
    lines(x=x, y=rep(sAIC[it],2), lwd=2, col=grey(0.5))
    lines(x=x, y=rep(sBIC[it],2), lwd=2)
#    if (it > 1 & sAIC[it] < sAIC[it - 1] & sAIC[it] < sAIC[it + 1])
#      draw.local.min(sAIC[it], ics$NbPartitions[it], grey(0.5))
#    if (it > 1 & sBIC[it] < sBIC[it - 1] & sBIC[it] < sBIC[it + 1])
#      draw.local.min(sBIC[it], ics$NbPartitions[it], "black")


  }
  labels<-seq(r.aic[1], r.aic[2], by=100)
  axis(2, at=(0:(length(labels)-1))/(length(labels)-1), labels, col=grey(0.5), col.axis=grey(0.5))
  mtext(side=2, "AIC", col=grey(0.5), line=3)
  labels<-seq(r.bic[1], r.bic[2], by=100)
  axis(4, at=(0:(length(labels)-1))/(length(labels)-1), labels, col="black", col.axis="black")
  mtext(side=4, "BIC", col="black", line=3)
}

plot.profile(free.aic.prf, 1:30, 1:35, brlen, free.aic.tbl, free.aic)

  
plot.profile(join.aic.prf, 1:20, 1:38, brlen, join.aic.tbl)
plot.profile(free.aic.prf, 1:10, 1:10, brlen, free.aic.tbl)

