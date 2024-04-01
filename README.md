
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Depth to bedrock prediction

Code repository for the PLOS ONE article ‘Evaluating spatially enabled
machine learning approaches for depth to bedrock mapping’.
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0296881

## Overview

The depth to bedrock prediction model for Alberta uses an R project
workflow based on the tidyverse and tidymodels suite of packages. The
workflow is re-executed via the following steps:

1.  Data preparation (downloading publicly available remote sensing
    datasets and terrain analysis using the [Rsagacmd
    package](https://cran.r-project.org/package=Rsagacmd)).

2.  Water well litholog augmentation by classifying lithological
    descriptions into surficial and bedrock units based on a statistical
    natural language processing approach and text pattern matching.

3.  Model evaluation and selection based on cross validation and quality
    of the predicted DTB maps in several physiographically varying
    sub-regions.

4.  Final DTB prediction.

## Requirements

This R project uses a reproducible environment based on the [renv
package](https://rstudio.github.io/renv/articles/renv.html). This
environment uses R package versions that were installed under Ubuntu
22.04 and is not guaranteed to cleanly install on other operating
systems, particularly Microsoft Windows. To install the required
packages, run the following code in the R console. The packages and
versions that will be installed are specified in the `renv.lock` file.

``` r
renv::restore()
```

This project has an external dependency - SAGA-GIS (\>= 7.3) where the
`saga_cmd` binary needs to be available in the system path. The
[Rsagacmd package](https://cran.r-project.org/package=Rsagacmd) package
is used to run the SAGA-GIS algorithms from R. SAGA-GIS can be installed
in Ubuntu 22.04 using the following commands:

``` bash
sudo apt-get install saga
```

For Windows users, the installation has to be performed manually by
downloading the SAGA-GIS binary from
[Sourceforge](https://sourceforge.net/projects/saga-gis/) and then
adding the path to the SAGA installation directory to PATH.

In addition, the project automatically downloads the required remote
sensing datasets. To download the MODIS data, a free NASA Earthdata
login is required. The login credentials need to be set in a `.Renviron`
file (a plain text file with no extension) in the project root directory
as ‘EARTHDATA_USER’=<username> and ‘EARTHDATA_KEY’=<password>, or set as
environment variables in the R session. This can be performed by:

``` r
Sys.setenv("EARTHDATA_USER" = "username")
Sys.setenv("EARTHDATA_KEY" = "password")
```

For the statistical natural language model prediction, XGBoost using GPU
is used. This requires a CUDA (Compute Unified Device Architecture)
enabled GPU and the CUDA toolkit to be installed. Alternatively, the
model can be trained on CPU by changing `tree_method` = ‘hist’ in the
`02-nlp.qmd` script, but this will result in a significant increase in
training time.

## Folder structure

The folder structure is organized as follows:

- `R`: R scripts used to prepare data and run the model.
- `projdata`: Data required by the model and visualizations.
- `outputs`: Model results (created by the scripts).
- `models`: Trained models (created by the scripts).
- `data/raw`: Raw data downloaded from remote sources (created by the
  scripts).
- `data/processed`: Processed data used in the model (created by the
  scripts).

## Data

The data includes the following:

- a geological pick dataset (‘picks.rds’) as a R data format (RDS)
  binary object. This data is also available from the Alberta
  Geological Survey (<https://ags.aer.ca/data-maps-models/digital-data>)
- water well data (‘lithologs.rds’), created on Nov 14, 2023 from data
  that is available in the Alberta Water Well Information Database.
- polygon features (‘physio-pettapiece.gpkg’) delineating the
  physiographic regions of Alberta, derived from [Physiographic regions
  of
  Alberta](https://open.alberta.ca/publications/physiographic-subdivisions-of-alberta).
- two polyline features, ‘cross-section-line-rainbow-lake.geojson’ and
  ‘cross-section-line-wcab.geojson’, used to create cross-sections of
  the model results shown in the PLOS ONE article.

## Code

Functions used by the scripts are stored in the `R` folder. The scripts
are organized into the following files:

- `01-grids.qmd`: Quarto notebook to download and prepare remote sensing
  data and terrain analysis.
- `02-nlp.qmd`: Quarto notebook to prepare lithological descriptions for
  machine learning.
- `03-training-data.R`: Prepare the training dataset that is used in the
  depth to bedrock machine learning model, i.e., prepare the feature
  matrix.
- `04-experiments.R`: Model evaluation and selection.
- `05-idw.R`: IDW interpolation of bedrock elevation.
- `06-kriging.R`: Kriging interpolation of DTB.
- `07-provincial-model.R`: Final DTB prediction model.
- `08-analysis-results.qmd`: Generate figures and tables that are used
  in the published paper.

Once these scripts are run, additional directories and model outputs
will be created within the project directory. The final directory
structure will look like:

``` r
fs::dir_tree(recurse = 2)
#> .
#> ├── LICENSE
#> ├── R
#> │   ├── 01-grids.qmd
#> │   ├── 02-nlp.qmd
#> │   ├── 03-training-data.R
#> │   ├── 04-experiments.R
#> │   ├── 05-idw.R
#> │   ├── 06-kriging.R
#> │   ├── 07-provincial-model.R
#> │   ├── 08-analysis-results.qmd
#> │   ├── functions
#> │   │   ├── dtb.R
#> │   │   ├── nlp.R
#> │   │   ├── plots.R
#> │   │   ├── predictors.R
#> │   │   └── resampling.R
#> │   ├── plos-one.csl
#> │   ├── plos2015.bst
#> │   └── zotero.bib
#> ├── README.Rmd
#> ├── README.md
#> ├── _dependencies.R
#> ├── data
#> │   ├── processed
#> │   │   ├── picks-nlp.rds
#> │   │   ├── predictors.tif
#> │   │   └── training-data.rds
#> │   └── raw
#> │       ├── alos-dem.tif
#> │       └── modis.tif
#> ├── depth-to-bedrock-plos-one.Rproj
#> ├── models
#> │   ├── cross-validation-dtb-prov.rds
#> │   ├── cross-validation-idw-prov.rds
#> │   ├── cross-validation-kriging-prov.rds
#> │   ├── experiments-cross-validation.rds
#> │   ├── experiments-dtb.rds
#> │   ├── experiments-importances.rds
#> │   ├── experiments-models.rds
#> │   ├── model-dtb-prov.rds
#> │   ├── model-nlp.rds
#> │   └── resamples-nlp.rds
#> ├── outputs
#> │   ├── dtb-rf-prov-pred-int.tif
#> │   ├── dtb-rf-prov.tif
#> │   ├── idw-prov.tif
#> │   ├── kriging-prov.tif
#> │   ├── picked-combined.gpkg
#> │   ├── picks-nlp-cv.csv
#> │   ├── picks-nlp.csv
#> │   └── picks.csv
#> ├── projdata
#> │   ├── cross-section-line-rainbow-lake.geojson
#> │   ├── cross-section-line-wcab.geojson
#> │   ├── lithologs.rds
#> │   ├── physio-pettapiece.gpkg
#> │   └── picks.rds
#> ├── renv
#> │   ├── activate.R
#> │   ├── library
#> │   │   └── R-4.3
#> │   ├── settings.json
#> │   └── staging
#> └── renv.lock
```
