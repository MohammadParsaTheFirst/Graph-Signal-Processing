# Graph Signal Processing - Community Detection in Stochastic Block Models

## Overview

This project explores community detection in graphs using the **Stochastic Block Model (SBM)**. The work involves generating graphs with known community structures, applying spectral clustering techniques, and analyzing clustering accuracy based on graph parameters.

## Problem Description

### Community Detection with Stochastic Block Models

1. **Stochastic Block Model (SBM) Setup**:  
   - Constructed an SBM graph with two communities where each node pair within the same community has a probability **p** of being connected, and nodes between communities have a connection probability **q**.
   
   ![SBM Graph Structure](path/to/your/image1.png)

2. **Spectral Clustering via Laplacian Eigenvectors**:  
   - Computed the unnormalized and normalized Laplacian matrices of the graph.
   - Used the second and third smallest eigenvectors of the Laplacian for spectral clustering, effectively visualizing community separation.

   ![Spectral Clustering on SBM Graph](path/to/your/image2.png)

3. **Accuracy vs. Criterion Analysis**:  
   - Assessed the accuracy of community detection as a function of **|√α - √β|**, reflecting the impact of within-community and between-community connection probabilities on classification success.
   
   ![Accuracy vs Criterion Plot](path/to/your/image3.png)

## Results Directory

All relevant matrices (adjacency, Laplacian, etc.) are stored in the `results` directory for further analysis and reproducibility.

## Additional Information

- **Spectral Clustering Insight**: The spectral clustering approach utilized eigenvectors of the Laplacian matrix, with classification achieved by examining the sign of specific eigenvalues associated with graph structure.
- **Effect of Edge Density on Community Detection**: This experiment demonstrates how the probability parameters (p and q) influence the ability of spectral clustering to accurately separate communities.

---

## Credits

This community detection project is part of the **Graph Signal Processing** coursework.
