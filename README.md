# Code for double cone transcriptomic analysis

## Table of Contents
1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Getting Started](#getting-started)
4. [Usage](#usage)
5. [Cite](#cite)

## Project Overview
This repository contains analyses for **Transcriptomic insights into the evolutionary origin of tetrapod double cone**. 

The paper can be accessed [here](https://www.cell.com/current-biology/fulltext/S0960-9822(25)00376-8).

## Repository Structure
- **src/AnalysisAndFigures_Part1.Rmd**: R Markdown notebook where most analysis, figures, and supplementary materials were generated.
- **src/AnalysisAndFigures_Part2.Rmd**: R Markdown notebook where most analysis, figures, and supplementary materials were generated.
- **src/SamapAnalysis_v3_5species.ipynb**: Jupyter notebook with step-by-step SAMap integration of 5 species
- **src/SamapAnalysis_v3_swap.ipynb**: Jupyter notebook with step-by-step SAMap integration of chicken, lizard, and zebrafish
- **src/RunCellpose.ipynb**: Cellpose nuclei segmentation for *in situ* hybridization analysis 
- **src/CountCellPuncta-RedOpsin.ipynb**: Nuclei counting for *in situ* hybridization analysis
- **src/utils**: Utility functions used in the notebooks. 

## Getting Started
To run the analyses, clone this repository and ensure required dependencies are installed. 

```bash
git clone https://github.com/shekharlab/DoubleCones.git
cd DoubleCones/src
```

1. Run through the analysis in **AnalysisAndFigures.Rmd** until export to h5ad
2. Run through the **SamapAnalysis.ipynb**
3. Finish the rest of the **AnalysisAndFigures.Rmd**
4. Generate segmentation masks in **RunCellpose.ipynb**
5. Count the positive nuclei using **CountCellPuncta-RedOpsin.ipynb**

## Usage
Due to the size of the files, the data directory is empty by default. For data files needed to run analyses, please email [Dario Tommasini](mailto:dtommasini@berkeley.edu) or [Karthik Shekhar](mailto:kshekhar@berkeley.edu). 

## Cite
If you find our code, analysis, or results useful and use them in your publications, please cite us using the following citation: 

Tommasini, Yoshimatsu, Baden, and Shekhar. Transcriptomic insights into the evolutionary origin of the tetrapod double cone. *Current Biology*. 2025. 

