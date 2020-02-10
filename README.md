# synal
Synteny analysis of yeast open reading frames.

This package takes syntenic blocks of an ORF and align-find best orthologus ORF-compare those ORFs and creates a dataframe with comparison data between sequences.

This package uses `tidyverse` and `Bioconductor` packages extensively.

This can be installed by:
```
devtools::install_github('oacar/synal')
```

See the `vignettes/vignette.Rmd` for a use case. Don't hesitate to contact or open an issue.

This is used in Vakirlis et al 2020, De novo emergence of adaptive membrane proteins from thymine-rich genomic sequences, Nature Comm paper.

Slightly improved python version can be found in github.com/oacar/SynORFan


