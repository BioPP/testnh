Note: the input tree `lysozymeSmall.dnd` was manually edited to remove the human branch in order to compare with the original publication.
The resulting tree is named `lysozymeSmall2.dnd`.

# Fit a homogeneous model

```bash    
bppml --noninteractive=yes param=ML.bpp >& bppml_h.out
```

This will store optimized parameter in file `lysozymeLarge.mh_h.params.txt`
and corresponding tree in file `lysozymeLarge.ml_h.txt`.
A YN98 substitution model is used.

# Try to find a better non-homogeneous model to describe the data

## Cluster nodes using substitution mapping.

We will map two types of substitutions: N and S.

```bash
mapnh --noninteractive=yes param=MapNH.bpp > mapnh.out 
```

We try to cluster nodes freely:

```bash
clustnh --noninteractive=yes param=ClustNH.bpp \
      test.branch.neighbor=no \
      output.cluster_tree.file=lysozymeSmall.cluster_free.dnd > clustnh_free.out
```

And also adding the constraint to cluster only adjacent nodes:

```bash
clustnh --noninteractive=yes param=ClustNH.bpp \
      test.branch.neighbor=yes \
      output.cluster_tree.file=lysozymeSmall.cluster_join.dnd > clustnh_join.out
```

This will first create a tree in Nhx format, identical to the input one, but with
node IDs. It then creates a cluster tree in Newick format, with leaves being node IDs.

## Use the clustering tree to define partitions, and test them using model testing.

We test both the clustering with nodes constrained to be adjacent or not in clustering,
and use two model selection criteria (AIC and BIC):

```bash
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=lysozymeSmall.cluster_free.dnd \
       partition.test=BIC \
       METHOD=free_BIC >& partnh_free_BIC.out

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=lysozymeSmall.cluster_join.dnd \
       partition.test=BIC \
       METHOD=join_BIC >& partnh_join_BIC.out

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=lysozymeSmall.cluster_free.dnd \
       partition.test=AIC \
       METHOD=free_AIC >& partnh_free_AIC.out

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=lysozymeSmall.cluster_join.dnd \
       partition.test=AIC \
       METHOD=join_AIC >& partnh_join_AIC.out
```

# Compare with Zhang results:

```bash
bppml --noninteractive=yes param=ML_PNAS.bpp >& bppml_nh_PNAS.out
``

