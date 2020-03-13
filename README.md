# dummy_project

Version 0.1.0

This is a dummy project generated to practice using git and github. There is only one script in the root directory that can be executed using python: reproduce_Spall2011.py. It estimates the hydrographic properties of a marginal sea that is subject to buoyancy loss (like the Labrador Sea or the Greenland Sea), based on a publication by Spall: https://doi.org/10.1175/2011JCLI4130.1 . Running the script will output figures in the /results/figures/ directory. 

## Output
- gamma_sensitivity.pdf 

shows sensitivity of the temperature of the outflow and the interior of the basin to changes in atmospheric forcing.

## Basic dependency
Python 3

## Project organization

```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── requirements.txt
├── bin                <- Compiled and external code, ignored by git (PG)
│   └── external       <- Any external source code, ignored by git (RO)
├── config             <- Configuration files (HW)
├── data               <- All project data, ignored by git
│   ├── processed      <- The final, canonical data sets for modeling. (PG)
│   ├── raw            <- The original, immutable data dump. (RO)
│   └── temp           <- Intermediate data that has been transformed. (PG)
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── src                <- Source code for this project (HW)

```


## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)

## Citation

Please [cite this project as described here](/CITATION.md).
