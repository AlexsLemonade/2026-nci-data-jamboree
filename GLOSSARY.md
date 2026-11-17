# Glossary

## Deconvolution

**Reference**
Labelled single-cell data, summarised per cell type (and optionally cell state), used to estimate the composition of bulk samples.
_Avoid_: signature matrix (a reference may be richer than one), atlas

**Reference source**
A published labelled single-cell dataset, often an integrated atlas, from which one or more references are built. Several references (e.g. different subsets or label choices) can come from one reference source.
_Avoid_: atlas (unqualified), reference (for the source data itself)

**Reference sample**
A single-cell sample whose cells contribute to a reference; several sequencing libraries from the same sample count as one reference sample. Distinct from a bulk sample.
_Avoid_: donor, subject, library, sample (unqualified)

**Cell type**
The label level at which proportions are primarily reported (e.g. T cell, monocyte, malignant).
_Avoid_: class, cluster

**Cell state**
An optional finer label nested within exactly one cell type (e.g. a GMP-like malignant state within the malignant cell type).
_Avoid_: subtype, substate

**Malignant cell type**
The single cell type whose expression profile may be re-estimated per bulk sample, because tumour cells vary between patients more than normal cells do.
_Avoid_: leukemic, tumor, key

**Bulk sample**
A bulk RNA-seq expression profile whose composition is being estimated; the target of deconvolution.
_Avoid_: target, mixture

**Proportion**
The estimated fraction of a bulk sample attributed to a cell type or cell state; proportions at one level sum to 1 per bulk sample. A proportion is either an **mRNA share** (fraction of the bulk sample's transcripts from that cell type) or a **cell share** (fraction of its cells); the two differ when cell types differ in RNA content, and each method estimates one or the other.
_Avoid_: fraction, abundance, theta
