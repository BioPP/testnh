Note: after update of the Bio++ libraries, the results in this example are slightly different from the original publication.
The likelihood values are higher for equivalent models, and the clustering is different.

# Fit a homogeneous model:

```bash     
bppml --noninteractive=yes param=ML.bpp >& bppml_h.out
```
logL = -44782.2080388088

This will store optimized parameter in file `MantellidDataset.mh_h.params.txt`
and corresponding tree in file              `MantellidDataset.ml_h.tt`.
A YN98 substitution model is used.

For comparison, here is the full model:

```bash
bppml --noninteractive=yes param=ML_NH_Full.bpp > bppml_nh_full.out
```

logL = -44697.5169182945

And with PAML:
```bash
cd PAML_NH
codeml
```
logL = -44715.472820


# Cluster nodes using substitution mapping.

We will map two types of substitutions: N and S.

```bash
mapnh --noninteractive=yes param=MapNH.bpp > mapnh.out 
```

We try to cluster nodes freely:

```bash
clustnh --noninteractive=yes param=ClustNH.bpp \
      test.branch.neighbor=no \
      output.cluster_tree.file=MantellidDataset.cluster_free.dnd >& clustnh_free.out
```

And also adding the constraint to cluster only adjacent nodes:

```bash
clustnh --noninteractive=yes param=ClustNH.bpp \
      test.branch.neighbor=yes \
      output.cluster_tree.file=MantellidDataset.cluster_join.dnd >& clustnh_join.out
```

This will first create a tree in Nhx format, identical to the input one, but with
node IDs. It then creates a cluster tree in Newick format, with leaves being node IDs.

# Use the clustering tree to define partitions, and test them using model testing.

We test both the clustering with nodes constrained to be adjacent or not in clustering,
and use two model selection criteria (AIC and BIC):

```bash
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=MantellidDataset.cluster_free.dnd \
       partition.test=BIC \
       METHOD=free_BIC >& partnh_free_BIC.out
```

The resulting log likelihood is -45283.21110696586401900277, for 3 clusters

```
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=MantellidDataset.cluster_join.dnd\
       partition.test=BIC\
       METHOD=join_BIC > partnh_join_BIC.out
```

The resulting log likelihood is -45301.12723197422019438818, for 4 clusters

```bash
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=MantellidDataset.cluster_free.dnd\
       partition.test=AIC\
       METHOD=free_AIC > partnh_free_AIC.out
```

The resulting log likelihood is -45280.04409340591519139707 for 4 clusters

```bash
partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=MantellidDataset.cluster_join.dnd \
       partition.test=AIC \
       METHOD=join_AIC > partnh_join_AIC.out
```

The resulting log likelihood is -45277.62673461444501299411 for 12 clusters

# Assess the robustness of substitution mapping:

```bash
mapnh --noninteractive=yes \
      param=MapNH.bpp,MantellidDataset.model_free_BIC.bpp \
      input.tree.file=MantellidDataset.ml_nh_free_BIC.nhx \
      input.tree.format=NHX \
      "output.counts=PerBranchPerType(perBranchLength=false, file=MantellidDataset.mapping_post, format=tsv, splitNorm=true)" > mapnh_post.out

clustnh --noninteractive=yes param=ClustNH.bpp \
      input.counts.file=MantellidDataset.mapping_post.tsv \
      test.branch.neighbor=no \
      output.cluster_tree.file=MantellidDataset.cluster_free_post.dnd > clustnh_free_post.out

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=MantellidDataset.cluster_free_post.dnd \
       partition.test=BIC \
       partition.test.stop_condition=3 \
       METHOD=free_BIC_post > partnh_free_BIC_post.out
```
-45281.66154597837157780305, 3 clusters

```bash
clustnh --noninteractive=yes param=ClustNH.bpp \
      input.counts.file=MantellidDataset.mapping_post.tsv \
      test.branch.neighbor=yes \
      output.cluster_tree.file=MantellidDataset.cluster_join_post.dnd > clustnh_join_post.out

partnh --noninteractive=yes param=PartNH.bpp \
       input.cluster_tree.file=MantellidDataset.cluster_join_post.dnd \
       partition.test=BIC \
       partition.test.stop_condition=3 \
       METHOD=join_BIC_post > partnh_join_BIC_post.out
```
-45293.71063173310540150851, 6 clusters


     
