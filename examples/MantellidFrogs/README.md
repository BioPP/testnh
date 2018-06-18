1) Fit a homogeneous model:
===========================

```bash     
bppml --noninteractive=yes param=ML.bpp >& bppml_h.out &
```

logL = -44799.2881756184

This will store optimized parameter in file `MantellidDataset.mh_h.params.txt`
and corresponding tree in file              `MantellidDataset.ml_h.tt`.
A YN98 substitution model is used.

For comparison, here is the full model:
```bash
bppml --noninteractive=yes param=ML_NH_Full.bpp > bppml_nh_full.out &
```

logL = -44715.4752201262

And with PAML:
```bash
cd PAML_NH
codeml
```
logL = -44715.472820


2) Cluster nodes using substitution mapping.
============================================

We will map non-synonymous substitutions.

We try to cluster nodes freely:
```bash
mapnh --noninteractive=yes param=MapNH.bpp \
      test.branch.neighbor=no \
      output.cluster_tree.file=Mantellid.cluster_free.dnd >& mapnh_free.out &
```

And also adding the constraint to cluster only adjacent nodes:
```bash
mapnh --noninteractive=yes param=MapNH.bpp \
      test.branch.neighbor=yes \
      output.cluster_tree.file=Mantellid.cluster_join.dnd >& mapnh_join.out &
```

This will first create a tree in Nhx format, identical to the input one, but with
node IDs. It then creates a cluster tree in Newick format, with leaves corresponding to node IDs.

3) Use the clustering tree to define partitions, and test them using model testing.
===================================================================================

We test both the clustering with nodes constrained to be adjacent or not in clustering,
and use two model selection criteria (AIC and BIC):

```bash     
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=Mantellid.cluster_free.dnd \
       partition.test=BIC \
       METHOD=free_BIC > partnh_free_BIC.out &
```

The resulting log likelihood is -45298.25456129185477038845, for 2 clusters

```
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=Mantellid.cluster_join.dnd\
       partition.test=BIC\
       METHOD=join_BIC > partnh_join_BIC.out &
```
CONTINUE HERE
The resulting log likelihood is -45297.52684668595611583441, for 5 clusters

```bash
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=Mantellid.cluster_free.dnd\
       partition.test=AIC\
       METHOD=free_AIC > partnh_free_AIC.out &
```

The resulting loge likelihood is -45274.84135715108277508989
 for 11 clusters

```bash
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=Mantellid.cluster_join.dnd \
       partition.test=AIC \
       METHOD=join_AIC > partnh_join_AIC.out &
```

The resulting log likelihood is -45280.47620085658127209172 for 12 clusters

4) Assess the robustness of substitution mapping:
=================================================

```bash
mapnh --noninteractive=yes \
      param=MapNH.bpp \
      param=MantellidDataset.model_free_BIC.bpp \
      input.tree.file=MantellidDataset.ml_nh_free_BIC.nhx \
      input.tree.format=Nhx \
      test.branch.neighbor=no \
      output.counts.tree.prefix=MantellidDataset.counts_post \
      output.cluster_tree.file=MantellidDataset.cluster_free_post.dnd > mapnh_free_post.out &

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=MantellidDataset.cluster_free_post.dnd \
       partition.test=BIC \
       partition.test.stop_condition=3 \
       METHOD=free_BIC_post > partnh_free_BIC_post.out &
```
-45298.25456145822681719437, 2 clusters

```bash
mapnh --noninteractive=yes \
      param=MapNH.bpp \
      param=MantellidDataset.model_join_BIC.bpp \
      input.tree.file=MantellidDataset.ml_nh_join_BIC.nhx \
      input.tree.format=Nhx \
      test.branch.neighbor=yes \
      output.counts.tree.prefix=MantellidDataset.counts_post \
      output.cluster_tree.file=MantellidDataset.cluster_join_post.dnd > mapnh_join_post.out &

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=MantellidDataset.cluster_join_post.dnd \
       partition.test=BIC \
       partition.test.stop_condition=3 \
       METHOD=join_BIC_post > partnh_join_BIC_post.out &
```
-45297.52684665376727934927, 5 clusters


     
