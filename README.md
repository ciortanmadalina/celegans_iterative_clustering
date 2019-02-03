## Clustering c-elegans data

This project attemps to cluster a particularly sparse c-elegans dataset by iteratively applying the Seurat method until obtaining clusters of expected dimensions.   
These operations are applied on the original input data but also on an imuted version which have been calculated using the ALRA algorithm (Credit to https://github.com/pavlin-policar/ALRA  and https://github.com/KlugerLab/ALRA).

This project outputs the clusters and also the marker genes for each cluster.