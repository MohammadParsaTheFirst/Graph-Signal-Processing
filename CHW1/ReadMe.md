# Graph Signal Processing - Assignment Results

## Overview

This repository contains results from a Graph Signal Processing assignment. The project focuses on implementing various graph-based operations using the **GSPBox** library, with MATLAB, covering topics like graph product computations, graph signals, spectral analysis, and community detection through stochastic block models.

## Prerequisites

- **MATLAB** with **mingw-w64** support
- **GSPBox** library, available [here](https://epfl-lts2.github.io/gspbox-html/download.html)

## Problem Descriptions

### Problem 1: Graph Products and Signal Processing

1. **Graph Construction**:  
   - Defined two graphs, **G1** and **G2**, using the `gsp_comet` and `gsp_ring` functions.
   - Computed the **Kronecker** product and **Cartesian** product of these graphs, yielding **Gt** and **Gs** respectively.

   ![Graph G1 and G2 with Products](path/to/your/image1.png)

2. **Random Graph Signal**:  
   - Generated and visualized a random graph signal on **Gt**.
   
   ![Random Signal on Graph H](path/to/your/image2.png)

3. **Spectral Analysis of Graph Signal**:  
   - Computed eigenvalues and eigenvectors of the Laplacian matrix for **Gt** and visualized its spectrum.
   
   ![Spectrum of Graph Laplacian](path/to/your/image3.png)

4. **Eigenvectors as Signals**:  
   - Plotted the four eigenvectors associated with the two largest and two smallest eigenvalues to analyze smoothness on the graph.
   
   ![Eigenvectors as Signals on Graph](path/to/your/image4.png)

### Problem 2: Community Detection in Stochastic Block Models

1. **Stochastic Block Model (SBM) Setup**:  
   - Defined a stochastic block model **SBM(n, p, q)** to generate communities within a graph. Each node has a probability of connection with nodes in its community (p) or in the other community (q).
   
   ![SBM Graph Clustering with Unnormalized Laplacian](path/to/your/image5.png)

2. **Community Detection via Laplacian Eigenvectors**:  
   - Applied spectral clustering by positioning nodes based on the second and third smallest eigenvectors of the Laplacian matrix, achieving community separation.

   ![SBM Graph Clustering with Normalized Laplacian](path/to/your/image6.png)

3. **Accuracy vs. Criterion Plot**:  
   - Evaluated clustering accuracy as a function of **|√α - √β|** to assess the effect of edge density on classification accuracy.
   
   ![Accuracy vs Criterion Plot](path/to/your/image7.png)

## How to Use

1. **Run MATLAB Scripts**: Ensure all prerequisites are met and run the provided MATLAB scripts for graph signal processing tasks.
2. **Plotting Results**: Modify the file paths in the script if necessary to view your plots in MATLAB.

## Results Directory

All computed matrices, including adjacency and weight matrices for graphs, are saved in the `results` directory for further analysis.

## Additional Information

- **Graph Product Computation**: The adjacency matrices of products **Gt** and **Gs** were derived based on Kronecker and Cartesian product rules.
- **Stochastic Block Model Clustering**: Classification was performed using the sign of the second smallest eigenvector of the Laplacian matrix to determine community structure.

---

## Credits

This assignment was completed as part of the **Graph Signal Processing** coursework.
