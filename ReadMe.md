# Electrode Pooling Data and Code

MM 5/11/2021

This repo accompanies our preprint

**Kyu Hyun Lee, Yu-Li Ni, Jennifer Colonell, Bill Karsh, Jan Putzeys, Marius Pachitariu, Timothy D. Harris, and Markus Meister (2021) Electrode pooling: boosting the yield of extracellular recordings with switchable silicon probes**

It contains all the data and code needed to reproduce the published analysis.

## Contents of the repo
- Jupyter notebooks (Python) and m-files (Matlab): `Theory`, `Saline`, `Invivo`, `Simulation`, `Simulation_Fig`, `accuracyQC.m`. These notebooks develop the various topics of analysis, starting from raw data, producing figure panels and numerical results for the article along the way. They contain a good number of comments and mathematical sections to guide the user. 

- `code/`: Contains routines accessed from the notebooks.

- `data/`: Data files, both input and output. 

- `figs/`: Image files that make up the figure panels in the article.

## How to reproduce all the analysis starting from raw data

0. Read our paper. A version from May 3, 2021 is included in the repo.   
1. Empty the `figs/` directory.
3. Run the notebooks and m-files in this order: 
	* `Theory.ipynb` - Fig 2
	* `Saline.ipynb` - Fig 4, Fig 8
	* `Invivo.m` - Fig 5
	* `Simulation.m` - Fig 6
	* `Simulation_Fig.ipynb` - Fig 6
	* `accuracyQC.m` - Fig 9
4. Now the `figs/` directory should contain all the figure panels. 

## How to find code for a specific figure panel
- The names of all the figure panels (as numbered in the preprint) appear as headings in the Jupyter notebooks. Look through these to find your figure of interest. Or...
- In the `figs/` directory find the name of the figure file of interest, and search for that name in the Jupyter notebooks or m-files.

## Simulation details (optional) 
The m-file `code/fig_6B_code/main_genSimulation.m` creates samples of voltage traces for simulated pooling, with parameters of Amp, Firing rate, and Noise.

## Compatibility 
Matlab code was tested on Matlab version `R2020b`.

Jupyter Notebooks were run on MacOS using Python 3.8.


