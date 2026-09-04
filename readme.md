## Global GDE outlook repository

This is the code repository accompanying the manuscript: **"Groundwater-dependent ecosystems are missing from global research and policy agendas"**. Huggins, X., Rohde, M. M., Reinecke, R., Gnann, S., Saccò, M., Hose, G. C., Stella, J. C. & Kløve, B. *In review.*

This repository contains all scripts used to preprocess input data, run analyses, and generate the figures presented in the manuscript. The code performs raster harmonisation and spatial summary statistic calculations across several global datasets.

### Repository structure

- **`on_button.R`** calls `here()` to set the project root, then sources `0-setup/`
- **`0-setup/`** loads required packages, sets session parameters (e.g., seed, temp directory), and defines a small set of custom functions
- **`1-preprocessing/`** harmonises input datasets
- **`2-analysis-and-plots/`** generates summary statistics and figures reported and used in the manuscript

Scripts within `1-preprocessing/` and `2-analysis-and-plots/` are numbered in the order they should be run. Some scripts depend on outputs written in scripts earlier in their numbered order.

### System specifications used in code development

**Operating system:** Windows 11.\
**Language version:** R 4.6.1.\
**Non-standard hardware:** none, all scripts can run on a standard desktop/laptop.\
**Package dependencies and versions:** see [session_info.txt](./session_info.txt)

### Install guide

1.  Clone or download this repository
2.  Download all source data listed in Supplementary Table 1
3.  Install packages listed in `0-setup/00-packages-in.R` ( **install time** should only be a few minutes if R/RStudio are already installed.)
4.  Open any script in your IDE. `here::here()` sets the root automatically.

### Use guide

1.  Run `on_button.R` at the start of every session, which sources all scripts in `0-setup/`
2.  Run the scripts in `1-preprocessing/` then `2-analysis-and-plots/` in their numbered order

### License

This repository is under the MIT License. See [LICENSE](./LICENSE) for details.

### Contact

Xander Huggins — <https://orcid.org/0000-0002-6313-8299>\
email: [xander.huggins@glasgow.ac.uk](mailto:xander.huggins@glasgow.ac.uk)
