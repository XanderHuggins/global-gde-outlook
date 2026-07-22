# Global GDE outlook repository

This is the code repository accomapying the manuscript: 
"Groundwater-dependent ecosystems are missing from global research and policy agendas" - Huggins, X., Rohde, M. M., Reinecke, R., Gnann, S., Saccò, M., Hose, G. C., Stella, J. C. & Kløve, B. *In review*.

This repository contains all scripts used to preprocess input data, run analyses, and generate the figures presented in the manuscript.

### Repository structure
- `on_button.R` — calls the `here()` function (sets project root) and sources `00-setup/`
- `00-setup/` — loads required packages and sets global options
- `0-functions/` — custom functions called throughout the analysis
- `1-preprocessing/` — harmonises input datasets to a common grid and resolution
- `2-analysis-and-plots/` — generates figures and summary statistics

### Contact
Xander Huggins - [https://orcid.org/0000-0002-6313-8299](https://orcid.org/0000-0002-6313-8299)
