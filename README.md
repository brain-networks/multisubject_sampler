# multisubject_sampler
multilayer method for sampling community structure from multiple subjects simultaneously

Code for detecting community structure in multiple subjects simultaneously. It works by creating a multi-layer, multi-subject modularity matrix and using the [Generalized Louvain code](http://netwiki.amath.unc.edu/GenLouvain/GenLouvain) to optimize a multi-layer modularity.

The advantage over existing methods is that, by emebdding all subjects in a single matrix and clustering that matrix directly, it becomes trivial to map communities across layers (subjects). In our paper, we used this property to explore inter-subject community variability.

If you use this code, please cite:
> Betzel, R. F., Bertolero, M. A., Gordon, E. M., Gratton, C., Dosenbach, N. U., & Bassett, D. S. (2018). The community structure of functional brain networks exhibits scale-specific patterns of variability across individuals and time. bioRxiv, 413278.
