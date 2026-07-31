# Semi-resolved-CFDEM-via-point-cloud-based-coarse-graining
This repository contains the development code for the point-cloud-based coarse-graining method in the semi-resolved CFDEM framework.

For each particle, a Fibonacci point cloud is first generated in the voidFractionModel. The generated point clouds are then used in both the force model and the averaging-field model for particle–fluid coupling calculations. The two-way coupling between CFD and DEM is controlled through the cfdemCloud class and the cfdemSolverGravity solver.

Cite: Liu, Y., Jing, L., Fu, X., & Shi, H. (2026). Enhancing semi-resolved CFD-DEM for dilute to dense particle-fluid systems: A point cloud based, two-step mapping strategy via coarse graining. Journal of Computational Physics, 115227.
