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




