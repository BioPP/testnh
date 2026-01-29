1) Dn/Ds analysis:
==================

1.1) Fit a homogeneous model, with a single Dn/Ds for all branches:
-------------------------------------------------------------------

```{bash}     
     bppml --noninteractive=yes param=ML.bpp > bppml_h.out &
     
     This will store optimized parameter in file *.params_h.bpp
     and corresponding tree in file              *.ml_h.dnd
```

A YN98+F3X4 codon substitution model is used.

1.2) Fit a complete NH model, with one dN/dS per branch.
--------------------------------------------------------

### 1.2.1) With PAML:

```{bash}
codeml
```

### 1.2.2) With BppML:

Note: this takes a lot of time given the number of parameters! (count several days!)
```{bash}
bppml --noninteractive=yes param=ML_NH.bpp > bppml_nh.out &
```

2) Try to find a better non-homogeneous model to describe the data.
===================================================================

2.1) Cluster nodes using substitution mapping.
----------------------------------------------

We will map non-synonymous substitutions.

We try to cluster nodes freely:
```{bash}
mapnh --noninteractive=yes param=MapNH.bpp test.branch.neighbor=no output.cluster_tree.file=mito_mammals.all.cluster_free.dnd >& mapnh_free.out &
```

And also adding the constraint to cluster only adjacent nodes:

```{bash}
mapnh --noninteractive=yes param=MapNH.bpp test.branch.neighbor=yes output.cluster_tree.file=mito_mammals.all.cluster_join.dnd >& mapnh_join.out &
```

This will first create a tree in Nhx format, identical to the input one, but with
node IDs. It then creates a cluster tree in Newick format, with leaves corresponding to node IDs.

2.2) Use the clustering tree to define partitions, and test them using model testing.
-------------------------------------------------------------------------------------

We test both the clustering with nodes constrained to be adjacent or not in clustering,
and use two model selection criteria (AIC and BIC):

```{bash}     
partnh param=PartNH.bpp --noninteractive=yes\
       input.cluster_tree.file=mito_mammals.all.cluster_free.dnd\
       partition.test=BIC\
       METHOD=_free_BIC >& partnh_free_BIC.out &
```

The resulting likelihood is -82741.7, for 7 clusters

```{bash}
partnh param=PartNH.bpp --noninteractive=yes\
       input.cluster_tree.file=mito_mammals.all.cluster_join.dnd\
       partition.test=BIC\
       METHOD=_join_BIC >& partnh_join_BIC.out &
```

(Runs very long...)

```{bash}
partnh param=PartNH.bpp --noninteractive=yes\
       input.cluster_tree.file=mito_mammals.all.cluster_free.dnd\
       partition.test=AIC\
       METHOD=_free_AIC >& partnh_free_AIC.out &
```

     The resulting likelihood is -19408.2, for 7 clusters

```{bash}
partnh param=PartNH.bpp --noninteractive=yes\
       input.cluster_tree.file=mito_mammals.all.cluster_join.dnd\
       partition.test=AIC\
       METHOD=_join_AIC >& partnh_join_AIC.out &
```

     The resulting likelihood is -19417.4, for 14 clusters


