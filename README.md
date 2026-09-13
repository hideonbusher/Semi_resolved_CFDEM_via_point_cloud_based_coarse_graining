# Semi-resolved-CFDEM-via-point-cloud-based-coarse-graining
This repository contains the development code for the point-cloud-based coarse-graining method in the semi-resolved CFDEM framework.

For each particle, a Fibonacci point cloud is first generated in the CGvoidFractionModel. The generated point clouds are then used in the force model for particle–fluid coupling calculations. 

Two force-model variants are provided: CGlucy12Dgidaspow is intended for small-scale simulations (more stable), while CGlucy15Dgidaspow is designed for large-scale parallel simulations. Select only one according to the computational scale; both must be used with CGVoidFraction.

The two-way coupling between CFD and DEM is controlled through the cfdemCloud class and the cfdemSolverGravity (a modified solver including fluid gravity).

Cite: Liu, Y., Jing, L., Fu, X., & Shi, H. (2026). Enhancing semi-resolved CFD-DEM for dilute to dense particle-fluid systems: A point cloud based, two-step mapping strategy via coarse graining. Journal of Computational Physics, 115227.
