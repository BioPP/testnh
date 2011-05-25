# Created on 24/03/11 by jdutheil
# Last modif 24/03/11 by jdutheil
# Read dn and ds mapping and merge data with life history traits.

dn<-read.table("PL_dn.csv", header=T, sep=",", stringsAsFactors=FALSE)
ds<-read.table("PL_ds.csv", header=T, sep=",", stringsAsFactors=FALSE)

dn$dn<-dn$Branch.length
dn$Branch.length<-NULL
ds$ds<-ds$Branch.length
ds$Branch.length<-NULL
d<-merge(dn, ds, all=TRUE)
d$Repro<-substr(d$Name,1,1)

write.table(d, "PL_dn_ds.csv", sep=",", row.names=FALSE, quote=FALSE)

boxplot(dn/ds~Repro,d, subset=ds > 0)


