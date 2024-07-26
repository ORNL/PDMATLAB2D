# PDMATLAB2D
## A simple, educational 2D MATLAB implementation of peridynamics

The PDMATLAB2D code is a meshfree peridynamics implementation in MATLAB suitable for simulation of two-dimensional fracture problems. The current version implements a bond-based brittle elastic peridynamic model and a critical stretch criterion for bond breaking. PDMATLAB2D provides an entry-level peridynamics computational tool for educational and training purposes. It also serves as an accessible and easily modifiable computational tool for peridynamics researchers who would like to adapt the code for a multitude of peridynamics simulations.

## Run

Examples can be run by passing an input file to the top-level function:
```
PDMATLAB2D('WavePropagation')
PDMATLAB2D('CrackBranching')
```

Tests are located in the `Tests/` folder and can be run individually to check individual components.

Individual plots can be created (as in the paper below) from the `PlottingExamples/` folder. Plots can be easily generated within the folder via the function `PlotPaperFigures'.

## Citing PDMATLAB2D

If you use PDMATLAB2D in your work, please cite the following paper:
```
@article{Seleson2024,
  author = {Seleson, Pablo and Pasetto, Marco and John, Yohan and Trageser, Jeremy and Reeve, Samuel Temple},
  title = {PDMATLAB2D: A Peridynamics MATLAB Two-dimensional Code},
  journal = {Journal of Peridynamics and Nonlocal Modeling},
  volume = {6},
  number = {1},
  pages = {149-205},
  year = {2024},
  url = {https://doi.org/10.1007/s42102-023-00104-w}
}
```

If you would like to cite the software itself, cite the current release from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.7348667).

## Contributing

We encourage you to contribute to PDMATLAB2D!

When contributing, please first discuss the changes you wish to make via issue (or directly with the PDMATLAB2D developers). Changes should be made through opening a pull request with `main` as the destination branch. Review from at least one developer is required.

PDMATLAB2D is distributed under the terms of the BSD 3-Clause license. All new contributions must be made under this license.

## License

PDMATLAB2D is distributed under a [BSD 3-Clause license](LICENSE).
