# Unveiling Adipose Populations Linked to Metabolic Health in Obesity
**Isabel Reinisch, Adhideb Ghosh, Falko Noé, Wenfei Sun, Hua Dong, Peter Leary, Arne Dietrich, Anne Hoffmann, Matthias Blüher, Christian Wolfrum**

This respository contains code and files related to our study: [Unveiling adipose populations linked to metabolic health in obesity. *Cell Metabolism*](https://doi.org/10.1016/j.cmet.2024.11.006).

## Abstract
Precision medicine is still not considered as a standard of care in obesity treatment, despite a large heterogeneity in the metabolic phenotype of individuals with obesity. One of the strongest factors influencing the variability in metabolic disease risk is adipose tissue (AT) dysfunction; however, there is little understanding of the link between distinct cell populations, cell-type-specific transcriptional programs, and disease severity. Here, we generated a comprehensive cellular map of subcutaneous and visceral AT of individuals with metabolically healthy and unhealthy obesity. By combining single-nucleus RNA-sequencing data with bulk transcriptomics and clinical parameters, we identified that mesothelial cells, adipocytes, and adipocyte-progenitor cells exhibit the strongest correlation with metabolic disease. Furthermore, we uncovered cell-specific transcriptional programs, such as the transitioning of mesothelial cells to a mesenchymal phenotype, that are involved in uncoupling obesity from metabolic disease. Together, these findings provide valuable insights by revealing biological drivers of clinical endpoints.

![Graphical Abstract](/images/graphical_abstract.png)

## Interactive web apps to explore data
&emsp;[App to explore correlations between bulk gene expression and clinical parameters](https://fgcz-shiny.uzh.ch/tnb_ethz_exploreMHUO)  <p>
&emsp;[App to explore snRNAseq data from subcutaneous AT](https://fgcz-shiny.uzh.ch/tnb_ethz_snMHUO_scAT) <p>
&emsp;[App to explore snRNAseq data from visceral AT](https://fgcz-shiny.uzh.ch/tnb_ethz_snMHUO_visAT) <p>

## Content

* `CellTypist`: CellTypist models trained on integrated data
* `FACS`: Raw data for included FACS analysis  
* `Source Data`: Source data for figure panels in the manuscript
* `Rscripts`: contains scripts used to analyze data
  * `bulk deconvolution`: depots specific and clinical RNA seq deconvolution
  * `Cell Types`: Sublcustering
  * `Integration`: Tested integration methods
  * `Network`: Marker gene comparison before integration
  * `Spatial_Deconvolution`: Code to replicate spatial deconvolution



## Contact
For questions regarding data analysis, please write to adhideb.ghosh[AT]hest.ethz.ch. <p>
Corresponding authors: Matthias Blüher (matthias.blueher[AT]medizin.uni-leipzig.de) and Christian Wolfrum (christian-wolfrum[AT]ethz.ch).
