# Created on 09/03/11 by jdutheil
# Last modif 09/03/11 by jdutheil
# Read dn and ds mapping and merge data with life history traits and taxonomy.

dn<-read.table("mammals_dn.csv", header=T, sep=",")
ds<-read.table("mammals_ds.csv", header=T, sep=",")

dn$dn<-dn$Branch.length
dn$Branch.length<-NULL
ds$ds<-ds$Branch.length
ds$Branch.length<-NULL
d<-merge(dn, ds, all=TRUE)
 
#Missing some species...
#t<-read.table("LHV.csv", header=TRUE)
#t$Name<-paste(t$Genus, t$species, sep="_")
##Prepare data for import in phyview:
#write.table(t, file="mammals_taxo_traits.csv", sep=",", row.names=FALSE, quote=FALSE)

m<-read.table("mammals_mass_taxo.csv", header=TRUE, sep=",")
dm<-merge(d, m, all=TRUE)
write.table(dm, "mammals_mass_taxo_dn_ds.csv", sep=",", row.names=FALSE, quote=FALSE)

plot(dn/ds~log(mass),dm)
cor.test(~I(dn/ds)+log(mass),dm,method="kendal")

dn.fam<-sapply(split(dm$dn, dm$Family),sum, na.rm=TRUE)
ds.fam<-sapply(split(dm$ds, dm$Family),sum, na.rm=TRUE)
dnds.fam<-sapply(split(dm$dn/dm$ds, dm$Family),mean, na.rm=TRUE)
mass.fam<-sapply(split(dm$mass, dm$Family),mean, na.rm=TRUE)
brle.fam<-sapply(split(dm$Branch.length, dm$Family),sum, na.rm=TRUE)
plot(dnds.fam~log(mass.fam))
plot(dn.fam/ds.fam~log(mass.fam))


dn.ord<-sapply(split(dm$dn, dm$Order),sum, na.rm=TRUE)
ds.ord<-sapply(split(dm$ds, dm$Order),sum, na.rm=TRUE)
dnds.ord<-sapply(split(dm$dn/dm$ds, dm$Order),mean, na.rm=TRUE)
mass.ord<-sapply(split(dm$mass, dm$Order),mean, na.rm=TRUE)
brle.ord<-sapply(split(dm$Branch.length, dm$Order),sum, na.rm=TRUE)
plot(dnds.ord~log(mass.ord))
plot(dn.ord/ds.ord~log(mass.ord))


