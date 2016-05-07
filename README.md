#Predicting Gene Regulation in Diverse Global Populations

######By Alexa Badalamenti, Jeffrey Ng, Virginia Saulnier, Shyam Shah, and Dr. Heather E. Wheeler

**Introduction**

  Most studies focus on populations of European ancestry. While genomic research has broadened knowledge on human genetics, work is needed to expand this knowledge to more diverse groups, otherwise an entire breadth of useful information is missing.
  
**Objective**

  We are working to expand genetic predictors of gene expression in additional world populations using SNP data alongside gene expression levels from the third phase of the International [*HapMap*](http://hapmap.ncbi.nlm.nih.gov/index.html.en)  Project.
  
- Populations include East African, West African, East Asian, and Mexican ancestry.
- Elastic net modeling is used to select genotypes and weights to best predict expression of each gene using the [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) package for R.

**Methods**

  A large portion of phenotypic variability in disease risk is due to regulatory variants which regulate gene expression levels. [*PrediXcan*](http://www.nature.com/ng/journal/v47/n9/full/ng.3367.html) is a gene-based association method, testing the mediating effects of gene expression levels by quantifying association between genetically regulated expression levels (GReX) and the phenotypic trait of interest. Gene expression can be decomposed into three basic components: what is genetically determined (GReX), what is altered by the phenotypic trait of interest, and remaining factors (including environment)

**Future Applications**

  We hope to interpret results to answer:
  
- Are predictors similar among diverse populations, or unique?
- Can we better predict gene expression when samples from diverse populations are combined rather than modeled singly?
- When testing predicted expression for associated traits, are new genes implicated, or were the genes previously found in European cohorts?

The hope is to advance biological knowledge of the underlying mechanisms of disease risk not assessed in GWAS studies alone. PrediXcan provides direction of effect, which may yield opportunities for therapeutic development and drug targets.

###***Instructions***

  To run our updated version of the PrediXcan pipeline for these diverse populations, you will need
  
  Software Requirements:  
  -Linux or Mac OS  
  -R  
    -Specific R packages:  
      glmnet  
      data.table  
      plotly  
      dply  
      ggplot2  
      GGally  
      RMarkdown  
      Plus all dependencies of aforementioned packages.  
    
  Input Files:   
    -Expression data matrix in the form of:  
    Columns are the snp ids  
    each row is the illumina/gene id  
    with our data we also had an extra row at the beginning  
    In our functionality/example code, we used the file:  
    expression Data/GIH_p3_expression.txt  
    please use this as a reference  
    - Gene annotation file with the following columns  
    id entrez ensembl ref chromosome start stop    
    In our code we used the file:  
    output.csv  
    please use this as a reference  
  
  -The bvm and fam files, please use the following files as guides:  
    hapmap3_r2_b36_fwd.consensus.qc.poly.bim  
    hapmap3_r2_b36_fwd.consensus.qc.poly.fam  
  
  -The SNP files in raw format, please use the following file as a guide:  
    'Chromosomes/Chr21.raw'  
  
  Since the snp files we used were not properly annotated, we had to impute missing values which may not be necessary for fully fitted data.
  Remember to change the directory, with the dir variable at the top
  
  Example Run:  
  -You can do an example run using the Functionality Test.R code
  -The files needed will be in the example run code file
  
  Output Files & Analysis:  
  - We have included the full results of our code on our github, under the results folder    
  - The visualization code is done in rmarkdown, and is the rmd file in the Results/Graphs (HTML included)  
  - Some of the plots in our powerpoint/paper were done in the Plotly browser due to some glitches with R 
