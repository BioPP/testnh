1) Dn/Ds analysis on Paland and Lynch data set.
===============================================

1.1) Fit a homogeneous model, with a single Dn/Ds for all branches:
-------------------------------------------------------------------
     
```bash
bppml --noninteractive=yes param=ML.bpp > bppml_h.out &
```
logL = -19467.60877771

This will store optimized parameters in file `PalandAndLynch28.params_h.bpp`
and corresponding tree in file               `PalandAndLynch28.ml_h.dnd`

A YN98+F3X4 codon substitution model is used.

For comparison with PAML:
```bash
cd PAML_H
codeml
```
PAML finds logL = -19467.750099

1.2) Same thing as above, but using the general syntax to specify non-homogeneous models:
-----------------------------------------------------------------------------------------

```bash
bppml param=ML_NH_1DNDS.bpp
```

This will store optimized parameters in file `PalandAndLynch28.params_NH_1DNDS.bpp`
and corresponding tree in file               `PalandAndLynch28.ml_NH_1DNDS.dnd`
The resulting log-likelihood is -19420.97


1.3) Same model as Lynch and Paland tried, using 2 models, one for internal branches, one for external branches:
----------------------------------------------------------------------------------------------------------------
     
```bash
bppml param=ML_NH_2DNDS.bpp
```

The resulting log-likelihood is -19394.8

1.4) Full model:
----------------

```bash
bppml --noninteractive=yes param=ML_NH_Full.bpp > bppml_nh_full.out &
```

logL = -19432.3262250668


```bash
cd PAML_NH
codeml
```
logL = -19432.429427

2) Try to find a better non-homogeneous model to describe the data.
===================================================================

2.1) Cluster nodes using substitution mapping.
----------------------------------------------

We will map synonymous and non-synonymous substitutions.

We try to cluster nodes freely:
```bash
mapnh --noninteractive=yes param=MapNH.bpp test.branch.neighbor=no output.cluster_tree.file=PalandAndLynch.cluster_free.dnd >& mapnh_free.out &
```

And also adding the constraint to cluster only adjacent nodes:
```bash
mapnh --noninteractive=yes param=MapNH.bpp test.branch.neighbor=yes output.cluster_tree.file=PalandAndLynch.cluster_join.dnd >& mapnh_join.out &
```

This will first create a tree in Nhx format, identical to the input one, but with node IDs. It then creates a cluster tree in Newick format, with leaves corresponding to node IDs.

2.2) Use the clustering tree to define partitions, and test them using model testing.
-------------------------------------------------------------------------------------

We test both the clustering with nodes constrained to be adjacent or not in clustering, and use two model selection criteria (AIC and BIC):
     
```bash
partnh param=PartNH.bpp --noninteractive=yes\
       input.cluster_tree.file=PalandAndLynch.cluster_free.dnd\
       partition.test=BIC\
       METHOD=free_BIC > partnh_free_BIC.out &
```

The resulting likelihood is -19441.35326020947832148522, for 3 clusters

```bash
partnh param=PartNH.bpp --noninteractive=yes\
       input.cluster_tree.file=PalandAndLynch.cluster_join.dnd\
       partition.test=BIC\
            METHOD=join_BIC > partnh_join_BIC.out &
```

The resulting likelihood is -19462.23612320328174973838, for 2 cluster

```bash
partnh param=PartNH.bpp --noninteractive=yes\
       input.cluster_tree.file=PalandAndLynch.cluster_free.dnd\
       partition.test=AIC\
       METHOD=free_AIC > partnh_free_AIC.out &
```
	    
The resulting likelihood is -19436.67162289058614987880, for 5 clusters

```bash
partnh param=PartNH.bpp --noninteractive=yes\
       input.cluster_tree.file=PalandAndLynch.cluster_join.dnd\
       partition.test=AIC\
       METHOD=join_AIC > partnh_join_AIC.out &
```

The resulting likelihood is -19445.28772075131928431801, for 11 clusters


3) Assess the robustness of substitution mapping:
=================================================

```bash
mapnh param=MapNH.bpp --noninteractive=yes\
      param=PalandAndLynch28.model_free_BIC.bpp\
      input.tree.file=PalandAndLynch28.ml_nh_free_BIC.nhx\
      input.tree.format=NHX\
      test.branch.neighbor=no\
      output.counts.tree.prefix=PalandAndLynch.counts_post\
      output.cluster_tree.file=PalandAndLynch.cluster_equilibrium_free_post.dnd > mapnh_free_post.out &

partnh param=PartNH.bpp --noninteractive=yes \
      input.cluster_tree.file=PalandAndLynch.cluster_equilibrium_free_post.dnd\
      partition.test=BIC\
      partition.test.stop_condition=3\
      METHOD=free_BIC_post > partnh_free_BIC_post.out &


mapnh param=MapNH.bpp --noninteractive=yes\
      param=PalandAndLynch28.model_join_BIC.bpp\
      input.tree.file=PalandAndLynch28.ml_nh_join_BIC.nhx\
      input.tree.format=NHX\
      test.branch.neighbor=yes\
      output.counts.tree.prefix=PalandAndLynch.counts_post\
      output.cluster_tree.file=PalandAndLynch.cluster_equilibrium_join_post.dnd > mapnh_join_post.out &

partnh param=PartNH.bpp --noninteractive=yes \
      input.cluster_tree.file=PalandAndLynch.cluster_equilibrium_join_post.dnd\
      partition.test=BIC\
      partition.test.stop_condition=3\
      METHOD=join_BIC_post > partnh_join_BIC_post.out &
```


