# Electrode Pooling Data and Code

MM 1/21/2021

This repo accompanies our preprint

**Kyu Hyun Lee, Yu-Li Ni, Jennifer Colonell, Bill Karsh, Jan Putzeys, Marius Pachitariu, Timothy D. Harris, and Markus Meister (2021) Electrode pooling: boosting the yield of extracellular recordings with switchable silicon probes**

It contains all the data and code needed to reproduce the published analysis.

## Contents of the repo
- Jupyter notebooks and Matlab files: `Theory`,`Saline`,`Invivo`,`main_Sim_fig6_data`,`Simulation`. These notebooks develop the various topics of analysis, starting from raw data, producing figure panels and numerical results for the article along the way. They contain a good number of comments and mathematical sections to guide the user. 

- `code/`: Contains routines accessed from the notebooks.

- `data/`: Data files, both input and output. 

- `figs/`: Image files that make up the figure panels in the article.

## How to reproduce all the analysis starting from raw data

0. Read our paper. A version from Dec 2020 is included in the repo.   
1. Empty the `figs/` directory.
3. Run the notebooks and Matlab files `Theory`,`Saline`,`Invivo`,`main_Sim_fig6_data`, and `Simulation`.
4. Now the `figs/` directory should contain all the figure panels. 
5. (optional) Matlab file `code/fig_6B_code/main_genSimulation.m` creates example pooling voltage traces, with parameters of Amp, Firing rate
and Noise

## How to find code for a specific figure panel
- The names of all the figure panels (as numbered in the preprint of Dec 2020) appear as level-3 headings in the Jupyter notebooks. Look through these to find your figure of interest. Or...
- In the `figs/` directory find the name of the figure file of interest, and search for that name in the notebooks.
