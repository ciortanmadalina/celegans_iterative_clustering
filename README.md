## Clustering c-elegans data

This project attemps to cluster a particularly sparse c-elegans dataset by iteratively applying the Seurat method until obtaining clusters of expected dimensions.   
These operations are applied on the original input data but also on an imuted version which have been calculated using the ALRA algorithm (Credit to https://github.com/pavlin-policar/ALRA  and https://github.com/KlugerLab/ALRA).

This project outputs the clusters and also the marker genes for each cluster.

## Project structure

- [seurat-celegans.ipynb](https://github.com/ciortanmadalina/celegans_iterative_clustering/blob/master/seurat-celegans.ipynb) : starts from the original scanpy Seurat algorithm and puts together all steps as part of method 'run' which can be then applied on the original dataset or on a selection (allowing to iteratively recluster the largerst groups)

- [ALRA_imputation.ipynb](https://github.com/ciortanmadalina/celegans_iterative_clustering/blob/master/ALRA_imputation.ipynb) : implements imputation

- figures folder contains all plots generated while running the Seurat pipeline

- Report.docx puts together all figures as final report of this exercise
