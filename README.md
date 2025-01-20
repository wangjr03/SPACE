# SPACE: Improved spatial domain detection by Spatially Variable Gene Clustering Adjusting for Cell Type Effect 

SPACE is a Gaussian process-based method that initially identifies SVGs and then establishes a dependency map among these SVGs using an intuitive approach. This map links each SVG with other SVGs exhibiting similar spatial patterns, ultimately clustering similarly expressed SVGs together. The resulting SVG-clusters from the cSVG algorithm can serve as inputs for further downstream analysis, for example spatial domain detection. 
Please find below a schematic diagram which shows the steps of SPACE for spatial domain detection:


<img width="1189" alt="Github_README" src="https://github.com/wangjr03/SPACE/assets/73495177/1dfe3d5d-996f-498d-9f45-ded17fd486c5">

## Input data

1. The main function fn_cSVG_par can be utilized as follow:

``` r
fn_cSVG_par(data_mat,
            loc_mat,
            method_step1 = "MargcorTest",
            thres_step1 = "standard",
            control = FALSE,
            ncores=7)
```

Arguments:

data_mat: normalized gene expression matrix with rows as genes and columns as spots. data_mat should be a mxN matrix where m is the number of genes and N is the number of spots.

loc_mat: matrix containing the spatial locations of the spots in data_mat. spatial locations are typically represented as two-dimensional coordinates, so loc_mat is typically a Nx2 matrix. 

method_step1: method to select related genes. The default method is "MargcorTest" or marginal correlation test. The other available options are "SIS" or sure independence screening, "rrcs" or Robust rank correlation based screening, "Enet" or Elastic net, "SIS+Enet".  

thres_step1: threshold for any method specified in method_step1. For method_step1= "MargcorTest", this argument is not required. The default is "standard", which calculates standard threshold for other methods. Look at the "screening" Github page for more information.

control = Specifies the intended use of the function. If control=FALSE, the function fn_cSVG_par performs step 1 of SPACE. If control=TRUE, it performs step 2 of SPACE. Please refer to the paper and look at the output section for more details about step 1 and step 2 results.

nodes: no of nodes used in parallel computation.

2. The function fn_cluster_genes can be utilized as follow:

``` r
fn_cluster_genes(SE_genes,
                 list1,
                 DEC_genes)
```

Arguments:

SE_genes: The names of SVGs as a string. 

list1: list of related genes(an output from step 2 of fn_cSVG_par, see below) 

DEC_genes: The names of significant genes from the output matrix of step 2 of fn_cSVG_par(see below).  


## Output

1. The output of the function fn_cSVG_par for step 1 (when argument control=FALSE):

A mx11 matrix where the rows are genes and each row contains p-values corresponding to 10 different kernels (5 Cosine and 5 Gaussian kernels with varying parameter values) and the combined p-value. The genes with significant adjusted(after multiplicity correction) combined p-value are selected as SVG. 

The output of the function fn_cSVG_par for step 2 (when argument control=TRUE):

Step 2 runs on the subset of SVGs found from step 1 result. If there are m1 SVGs, the outputs of this step are:
1. A m1x11 matrix where p-values corresponding 10 different kernels (5 Cosine and 5 Gaussian kernels with varying parameter values) and the combined p-value are given for each of m1 SVGs.
2. A list of length m1 containing information about the related genes for each m1 SVGs.

2. The output of fn_cluster_genes function:

The cluster IDs corresponding to all SVGs.

## Example 

The Simulation/simulated_datasets folder contains a simulated dataset. See [here] for more details. It demonstrates how the SPACE workflow works with step-by-step code, key functions, and helpful visualizations.

[here]: https://github.com/wangjr03/SPACE/blob/main/Simulation/Simulation.md

## Overview of Github Folders 

These folders demonstrate different features of SPACE:

1. Simulation: Demonstrates the effectiveness of the SVG clustering performance by SPACE.
2. Simulation_domainDetection: Assesses whether the SVG clusters identified by SPACE improve domain detection accuracy.
3. Data analysis: Showcases the effectiveness of the SPACE method using two publicly available datasets and a new dataset.
            [Analysis of breast cancer data]
   
            [Analysis of Human DLPFC 10x Genomics Visium dataset]
   
            [Analysis of HER2 breast tumor data]

   [Analysis of breast cancer data]: https://github.com/wangjr03/SPACE/tree/main/data%20analysis%20code/Breast%20Cancer%20data%20analysis
   [Analysis of Human DLPFC 10x Genomics Visium dataset]: https://github.com/wangjr03/SPACE/tree/main/data%20analysis%20code/DLPFC%20data%20analysis
   [Analysis of HER2 breast tumor data]: https://github.com/wangjr03/SPACE/tree/main/data%20analysis%20code/HER2%20data%20analysis
   
   
5. Data: Provides information/ source of datasets analyzed.
