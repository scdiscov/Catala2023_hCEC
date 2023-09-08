## A single-cell RNA-seq analysis unravels the heterogeneity of primary cultured human corneal endothelial cells

This repository contains the scripts used to analyze single cell RNA data published in [Català et al. 2023](https://www.nature.com/articles/s41598-023-36567-6).

### Abstract


The cornea is a transparent and avascular tissue located in front of the eye. Its inner surface is lined by a monolayer of corneal endothelial cells (CECs), which maintain the cornea transparency. CECs remain arrested in a non-proliferative state and damage to these cells can compromise their function leading to corneal opacity. The primary culture of donor-derived CECs is a promising cell therapy. It confers the potential to treat multiple patients from a single donor, alleviating the global donor shortage. Nevertheless, this approach has limitations preventing its adoption, particularly culture protocols allow limited expansion of CECs and there is a lack of clear parameters to identify therapy-grade CECs. To address this limitation, a better understanding of the molecular changes arising from the primary culture of CECs is required. Using single-cell RNA sequencing on primary cultured CECs, we identify their variable transcriptomic fingerprint at the single cell level, provide a pseudo-temporal reconstruction of the changes arising from primary culture, and suggest markers to assess the quality of primary CEC cultures. This research depicts a deep transcriptomic understanding of the cellular heterogeneity arising from the primary expansion of CECs and sets the basis for further improvement of culture protocols and therapies.


### Single-cell RNA sequencing Method
scRNAseq of primary cultured CECs was performed at Single Cell Discoveries (Utrecht, the Netherlands) following standard 10 × Genomics 3′ V3.1 chemistry protocol. Cells were rehydrated and loaded on the 10 × Chromium controller as follows. Approximately 10,000 cells were loaded per each sample specified in Table 2. The resulting sequencing libraries were prepared following a standard 10 × Genomics protocol and sequenced with an Illumina NovaSeq 6000 platform; read length: 150 bp, paired-end.

The following steps and scripts are used to process the raw fastq files and perform downstream analyses: 

### Scripts
#### 1. Mapping of raw data
The mapping of raw fastq files is done using 10X genomics cellranger 3’GEX pipeline: [Single-Library Analysis with Cell Ranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count).
The BCL files resulting from sequencing were transformed to FASTQ files with 10 × Genomics Cell Ranger mkfastq following its mapping with Cell Ranger count. During sequencing, Read 1 was assigned 28 bp, and were used for identification of the Illumina library barcode, cell barcode and unique molecular identifier (UMI). R2 was used to map the human reference genome GRCh38. Filtering of empty barcodes was done in Cell Ranger. Library names used in the paper g1 to g9 correspond to g011 - g020 (same order).

#### 2. Demultiplexing (2_souporcell.sh)
For each 10X sample, the mixed-genotype scRNAseq results are demultiplexed by individual. Check the details and requirements in the souporcell.sh script and [wheaton5/souporcell: Clustering scRNAseq by genotypes (github.com)](https://github.com/wheaton5/souporcell).
Briefly, the bam file and barcodes of each library were used as input together with the reference genome GRCh38. Besides the default parameters, the number of clusters was set to the number of multiplexed samples per library. The demultiplexing information for each cell was added to the metadata object in Seurat, as obtained by the next step.


#### 3. Processing all single cells (3_processing_all.R)
Processing filtering and clustering of the dataset containing all samples. These steps are done using Seurat’s pipeline in R. For information on the process and tutorials check [Satija Lab](https://satijalab.org/).
For each library a UMI cutoff was used to filter out low quality cells because of the differences between the libraries (i.e. g1-500, g2-3000, g3-500, g4-1313, g5-4000, g6-4000, g7-4000, g8-4000, g9-4000). Additionally, cells with less than 10% mitochondrial gene content were retained for analysis. The data of all 10 × libraries were merged and processed together. The merged dataset was normalized for sequencing depth per cell and log-transformed using a scaling factor of 10,000. The patient and library effect was corrected using Harmony43, as implemented in Seurat and used for dimensionality reduction and clustering of all cells. Cells were clustered using graph-based clustering and the original Louvain algorithm was utilized for modularity optimization.
Putative doublets were computationally identified using scDblFinder but did not compose a separate cluster and therefore were not removed from the dataset.

#### 4. Processing confluency cells (4_processing_confluency.R)
Selecting, processing and clustering confluency samples, using Seurat's pipeline as detailed above. 
    
Pseudotime and trajectory analysis is done using monocle3.
    
Integration with endothelial cell atlas is done using Seurat’s integration and cell label transfer. The necessary files to reconstruct the seurat object of the endothelial atlas can be downloaded here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186433
