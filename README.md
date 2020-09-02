# GapClust-rge-paper
Materials (datasets & R code) for reproduction of published results in the paper 

## Folder structure
```
|
|____Rfunction        R functions for the existing algorithms
|____Rproj                 R projects including the datasets and results for reproduction of the results in the paper
|____Rscript              R codes to reproduce the results in the paper
```


## Projects introduction
```
|
|____10X_subsample                      Simulation experiments including 99 datasets, each of which comprises 500 CD56+ NK cells, 
|                                                                 500 CD19+ B cells and various number of CD14 + monocyte ranging from 2 to 100  
|
|____jurkat_two_species               Simulation experiments including 10 datasets, each of which comprises 1540 293T cells 
|                                                                  and various number of Jurkat cells ranging from 8 to 82  
|
|____Doublet                                      Simulation experiments including 140 datasets, each of which comprises CD56+ NK cells,  
|                                                                 CD19+ B cells and 2 CD14 + monocyte to evaluate doublet detection of five algorithms
|
|____10X_DE_sensitivity                Simulation experiments to evaluate sensitivity of algorithms to different number of 
|                                                                 differentail expressed (DE) genes.
|
|____Runtime                                     Simulation experiments including 16 datasets subsampled from 68 k PBMC dataset with 
|                                                                 1000 - ~68000 cells to track the computation time and memory usage
|
|____Intestine                                     GapClust applied to the intestine dataset to identify rare cells
|
|
|____10X_full                                      GapClust applied to the 68 k PBMC dataset to identify rare cells

```

## Notes when running codes
* Make sure the folder structure of the whole project  is not changed.
* Modify the homedir variable to change the working directory, aiming to make sure the code can run properly on your machine. 
