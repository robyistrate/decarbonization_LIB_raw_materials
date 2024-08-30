# Decarbonizing lithium-ion battery primary raw materials supply chain: Available strategies, mitigation potential and challenges

## Overview
Repository to share the data and code associated with the scientific article **Istrate et al. Decarbonizing lithium-ion battery primary raw materials supply chain: Available strategies, mitigation potential and challenges. Joule (2024)**. The repository contains data files and code to import the life cycle inventories (LCIs), reproduce the results, and generate the figures presented in the article.

## Repository structure
The data folder includes:
- `lci_LIB_raw_materials.xlsx` contains LCI datasets for lithium-ion battery (LIB) raw materials formatted for use with [Brightway](https://github.com/brightway-lca).
- `IEA_Minerals demand.xlsx` contain data on minerals demand based on IEA scenarios, used to plot Figure 1 in the manuscript and to calculate the decoupling level.
- `results` folder within data contains csv files with the results. which are imported into `03_visualization.ipynb` for analysis and visualization purposes.

The notebooks folder includes:
- `00_project_setup.ipynb` sets up a new Brightway project and imports the ecoinvent database and the LCIs for lithium-ion battery raw materials production.
- `01_ghg_breakdown.ipynb` calculates life cycle GHG emissions and performs the breakdown analysis.
- `02_mitigation_potential.ipynb` creates LCI databases implementing decarbonization strategies for LIB raw materials and calculates the mitigation potential and trade-offs (including Monte Carlo uncertainty analysis).
- `03_visualization.ipynb` imports all results and generates the figures presented in the scientific article.
- `supporting_functions.py` contains functions required to run the notebooks

## How to use
To ensure the replication of the results presented in the article, please follow this steps:

**1. Set Up the Environment:**

Using Anaconda, build the environment using `environment.yaml`:
```
conda env create -f environment.yml
```

**2. Run the notebooks:**

Activate the new environment and run the notebooks in the specified order.

Feel free to reach out if you encounter any issues—I'm happy to help!
