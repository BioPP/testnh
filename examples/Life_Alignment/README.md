1) Perform a Bowker test.
=========================

1.1) Fit a homogeneous model:
-----------------------------

```bash     
bppml --noninteractive=yes param=ML.bpp > bppml_h.out
```
Log likelihood:  -14122.120990606
    
This will store optimized parameter in file `Life_Alignment.params_h.bpp`
and corresponding tree in file              `Life_Alignment.ml_h.dnd.`
A T92+Gamma substitution model is used.

1.2) Compute the Bowker statistic and compare to simulations under the estimated model
--------------------------------------------------------------------------------------

(See Dutheil and Boussau 2008)

```bash
testnh param=TestNH_homogeneous.bpp
```

The test is significant, sign of non-homogeneity.
All simulations details are output in file `TestNH_homogeneous.null.tsv`, and can
be analyzed and/or plotted using R for instance.

1.3) Compare with the full model:
---------------------------------

### 1.3.1) bppML

```bash
bppml --noninteractive=yes param=ML_NH_Full.bpp > bppml_nh_full.out
```

Log likelihood: -13825.2273607301
    
### 1.3.2) NHML

(requires `eval_nh` program from Nicolas Galtier
```bash     
/usr/bin/time -v eval_nhg Life_Alignment.mase Life_Alignment.dnd nhml.opt >& nhml.out
```



2) Try to find a better non-homogeneous model to describe the data.
===================================================================

2.1) Cluster nodes using substitution mapping.
----------------------------------------------

We will map two types of substitutions: AT->GC and GC->AT, and also tha AT->AT and GC->GC cases,
in order to be used with the "equilibrium" mapping. We use the homogeneous tree as input, after
having rooted it with the midpoint method (bppPhyView).

```bash
mapnh --noninteractive=yes param=MapNH.bpp > mapnh.out 
```

We try to cluster nodes freely:

```bash
clustnh --noninteractive=yes param=ClustNH.bpp \
      test.branch.neighbor=no \
      output.cluster_tree.file=Life_Alignment.cluster_free.dnd > clustnh_free.out
```

And also adding the constraint to cluster only adjacent nodes:

```bash
clustnh --noninteractive=yes param=ClustNH.bpp \
      test.branch.neighbor=yes \
      output.cluster_tree.file=Life_Alignment.cluster_join.dnd > clustnh_join.out
```

This will first create a tree in Nhx format, identical to the input one, but with
nodes ID. It then creates a cluster tree in Newick format, with leaves being nodes ID.

2.2) Use the clustering tree to define partitions, and test them using model testing.
-------------------------------------------------------------------------------------

We test both the clustering with nodes constrained to be adjacent or not in clustering,
and use two model selection criteria (AIC and BIC):

```bash     
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=Life_Alignment.cluster_free.dnd \
       partition.test=BIC \
       METHOD=free_BIC > partnh_free_BIC.out

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=Life_Alignment.cluster_join.dnd \
       partition.test=BIC \
       METHOD=join_BIC > partnh_join_BIC.out

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=Life_Alignment.cluster_free.dnd \
       partition.test=AIC \
       METHOD=free_AIC > partnh_free_AIC.out

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=Life_Alignment.cluster_join.dnd \
       partition.test=AIC \
       METHOD=join_AIC > partnh_join_AIC.out
```

3) Test if the resulting non-homogeneous models take into account the underlying non-homogeneity
================================================================================================

Using the Bowker test against the alternative hypothesis (see Dutheil and Boussau 2008):

```bash
testnh --noninteractive=yes param=TestNH_nonhomogeneous.bpp,Life_Alignment.model_free_BIC.bpp \
       bootstrap.dist_file=TestNH_nonhomogeneous.null_free_BIC.txt \
       input.tree.file=Life_Alignment.ml_nh_free_BIC.nhx > testnh_nh_free_BIC.out

testnh --noninteractive=yes param=TestNH_nonhomogeneous.bpp,Life_Alignment.model_join_BIC.bpp\
       bootstrap.dist_file=TestNH_nonhomogeneous.null_join_BIC.txt\
       input.tree.file=Life_Alignment.ml_nh_join_BIC.nhx > testnh_nh_join_BIC.out

testnh --noninteractive=yes param=TestNH_nonhomogeneous.bpp,Life_Alignment.model_free_AIC.bpp\
       bootstrap.dist_file=TestNH_nonhomogeneous.null_free_AIC.txt\
       input.tree.file=Life_Alignment.ml_nh_free_AIC.nhx > testnh_nh_free_AIC.out

testnh --noninteractive=yes param=TestNH_nonhomogeneous.bpp,Life_Alignment.model_join_AIC.bpp\
       bootstrap.dist_file=TestNH_nonhomogeneous.null_join_AIC.txt\
       input.tree.file=Life_Alignment.ml_nh_join_AIC.nhx > testnh_nh_join_AIC.out
```

All tests are now not significant!

4) Assess the robustness of substitution mapping:
=================================================

We rerun the substitution mapping procedure and clustering of branches using the model selected using BIC, for both the free and join approach:

```bash
mapnh --noninteractive=yes \
      param=MapNH.bpp,Life_Alignment.model_free_BIC.bpp\
      input.tree.file=Life_Alignment.ml_nh_free_BIC.nhx\
      input.tree.format=NHX\
      "output.counts=PerBranchPerType(perBranchLength=false, file=Life_Alignment.mapping_post, format=tsv, splitNorm=true)" > mapnh_post.out

clustnh --noninteractive=yes param=ClustNH.bpp \
      input.counts.file=Life_Alignment.mapping_post.counts.tsv \
      input.norms.file=Life_Alignment.mapping_post.norms.tsv \
      test.branch.neighbor=no \
      output.cluster_tree.file=Life_Alignment.cluster_free_post.dnd > clustnh_free_post.out

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=Life_Alignment.cluster_free_post.dnd \
       partition.test=BIC \
       partition.test.stop_condition=3 \
       METHOD=free_BIC_post > partnh_free_BIC_post.out

clustnh --noninteractive=yes param=ClustNH.bpp \
      input.counts.file=Life_Alignment.mapping_post.counts.tsv \
      input.norms.file=Life_Alignment.mapping_post.norms.tsv \
      test.branch.neighbor=yes \
      output.cluster_tree.file=Life_Alignment.cluster_join_post.dnd > clustnh_join_post.out

partnh --noninteractive=yes \
       param=PartNH.bpp \
       input.cluster_tree.file=Life_Alignment.cluster_join_post.dnd \
       partition.test=BIC \
       partition.test.stop_condition=3 \
       METHOD=join_BIC_post > partnh_join_BIC_post.out
```

# Compare values with R:

```R
f1<-read.table("Life_Alignment.ml_nh_free_BIC.parameters.csv", header=TRUE, row.names="NodeId")
f2<-read.table("Life_Alignment.ml_nh_free_BIC_post.parameters.csv", header=TRUE, row.names="NodeId")
j1<-read.table("Life_Alignment.ml_nh_join_BIC.parameters.csv", header=TRUE, row.names="NodeId")
j2<-read.table("Life_Alignment.ml_nh_join_BIC_post.parameters.csv", header=TRUE, row.names="NodeId")

# Merge by node ids:
f<-merge(f1, f2, by="row.names")
j<-merge(j1, j2, by="row.names")

# Plot results:

layout(matrix(1:2, nrow=1))
plot(jitter(T92.theta.y, amount=0.02)~jitter(T92.theta.x, amount=0.02), f, xlab="First step", ylab="Second step", main="Free model, BIC criterion", xlim=c(0,1), ylim=c(0, 1)); abline(0, 1)
plot(jitter(T92.theta.y, amount=0.02)~jitter(T92.theta.x, amount=0.02), j, xlab="First step", ylab="Second step", main="Join model, BIC criterion", xlim=c(0,1), ylim=c(0, 1)); abline(0, 1)

dev.print(pdf, "Figures/StepsComparison.pdf", width=10, height=5)
```

