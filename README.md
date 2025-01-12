# SPACE: Improved spatial domain detection by Spatially Variable Gene Clustering Adjusting for Cell Type Effect 

SPACE is a Gaussian process-based method that initially identifies SVGs and then establishes a dependency map among these SVGs using an intuitive approach. This map links each SVG with other SVGs exhibiting similar spatial patterns, ultimately clustering similarly expressed SVGs together. The resulting SVG-clusters from the cSVG algorithm can serve as inputs for further downstream analysis, for example spatial domain detection. 
Please find below a schematic diagram which shows the steps of SPACE for spatial domain detection:


<img width="1189" alt="Github_README" src="https://github.com/wangjr03/SPACE/assets/73495177/1dfe3d5d-996f-498d-9f45-ded17fd486c5">

## Input data

The main function to implement SPACE can be utilized as follow:

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

## Output




