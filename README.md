# SPACE: Improved spatial domain detection by Spatially Variable Gene Clustering Adjusting for Cell Type Effect 

SPACE is a Gaussian process-based method that initially identifies SVGs and then establishes a dependency map among these SVGs using an intuitive approach. This map links each SVG with other SVGs exhibiting similar spatial patterns, ultimately clustering similarly expressed SVGs together. The resulting SVG-clusters from the cSVG algorithm can serve as inputs for further downstream analysis, for example spatial domain detection. 
Please find below a schematic diagram which shows the steps of SPACE for spatial domain detection:

<img width="744" alt="Github_README" src="https://github.com/wangjr03/SPACE/assets/73495177/1ca86c59-e258-47b9-a32e-232268112dda">
