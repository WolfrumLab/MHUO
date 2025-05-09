# Unveiling Adipose Populations Linked to Metabolic Health in Obesity
**Isabel Reinisch, Adhideb Ghosh, Falko Noé, Wenfei Sun, Hua Dong, Peter Leary, Arne Dietrich, Anne Hoffmann, Matthias Blüher, Christian Wolfrum**

This respository contains code and files related to our study: [Unveiling adipose populations linked to metabolic health in obesity. *Cell Metabolism*](https://doi.org/10.1016/j.cmet.2024.11.006).

## Abstract
Precision medicine is still not considered as a standard of care in obesity treatment, despite a large heterogeneity in the metabolic phenotype of individuals with obesity. One of the strongest factors influencing the variability in metabolic disease risk is adipose tissue (AT) dysfunction; however, there is little understanding of the link between distinct cell populations, cell-type-specific transcriptional programs, and disease severity. Here, we generated a comprehensive cellular map of subcutaneous and visceral AT of individuals with metabolically healthy and unhealthy obesity. By combining single-nucleus RNA-sequencing data with bulk transcriptomics and clinical parameters, we identified that mesothelial cells, adipocytes, and adipocyte-progenitor cells exhibit the strongest correlation with metabolic disease. Furthermore, we uncovered cell-specific transcriptional programs, such as the transitioning of mesothelial cells to a mesenchymal phenotype, that are involved in uncoupling obesity from metabolic disease. Together, these findings provide valuable insights by revealing biological drivers of clinical endpoints.

![Graphical Abstract](/images/graphical_abstract.png)

## Interactive web apps to explore data
&emsp;[App to explore correlations between bulk gene expression and clinical parameters](https://shiny-public.fgcz.uzh.ch/app/tnb_ethz_exploreMHUO) <p>
&emsp;[App to explore snRNAseq data from subcutaneous AT](https://shiny-public.fgcz.uzh.ch/app/tnb_ethz_snMHUO_scAT) <p>
&emsp;[App to explore snRNAseq data from visceral AT](https://shiny-public.fgcz.uzh.ch/app/tnb_ethz_snMHUO_visAT) <p>

## Content

* `DEGs`: Contains tissue and cell type specific DEGs for subcutaneous / visceral AT
* `MarkerGenes`: Contains cell type and subpopulation specific top marker genes for subcutaneous / visceral AT
* `SNPdemux`: Contains code to run SNP demultiplexing using cellSNP and vireo
* `Rscripts`: Contains scripts used to analyze data
  * `bulkRNA`: DE analysis, Bisque deconvolution, Clinical correlations
  * `snRNA`: Pre-processing, Sample integration, Multi-cellular factor analysis, Cell type re-clustering


## Contact
For questions regarding data analysis, please write to adhideb.ghosh[AT]hest.ethz.ch. <p>
For questions regarding web applications, please write to falnoe[AT]ethz.ch. <p>
Corresponding authors: Matthias Blüher (matthias.blueher[AT]medizin.uni-leipzig.de) and Christian Wolfrum (christian-wolfrum[AT]ethz.ch).
