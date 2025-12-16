**Dynamically Evolving Bodies Resolved in Interacting Space**

D.E.B.R.I.S is a high-performance, open-source C library for precise orbital propagation and risk assessment in Low Earth Orbit (LEO) environments. Designed for massive-scale simulations of satellites, constellations, and debris clouds, it delivers unprecedented efficiency and accuracy through innovative numerical methods and multi-resolution modeling.

### Key Features

- **Adaptive Multi-Integrator Framework**: Automatic regime detection and switching between Verlet, DOPRI5(4), Yoshida-8 (with Kahan-Babuška-Lang compensated summation), and Gauss-Jackson 12th-order for optimal stability and precision.
- **Spectral Perturbation Baking**: FFT-based harmonic decomposition and frequency-domain resurrection for 300x+ speedups in long-term evolutions (e.g., multi-year constellation decay).
- **Deterministic Debris Cloud Necromancy**: 100-byte seeds regenerate thousands of breakup fragments exactly; hybrid EM/GA/DE/NSGA-II GMM fitting for multi-resolution representations (covariance → GMM → full particles on-demand).
- **Lazy Hierarchical Conjunction Pipeline**: Octree/BVH culling, Öpik/Alfano screening, and Monte Carlo refinement for catalog-scale risk assessment in minutes.
- **LEO-Specific Optimizations**: Drag blending, HAMR affine fallbacks, SRP with eclipse, J2-J6, third-body caching, and GPU-accelerated particle ensembles.

MIT licensed • Actively developed • Contributions welcome
