# PDMATLAB2D
## A simple, educational 2D MATLAB implementation of peridynamics

The PDMATLAB2D code is a meshfree peridynamics implementation in MATLAB suitable for simulation of two-dimensional fracture problems. The current version implements bond-based brittle elastic peridynamics models and a critical stretch criterion for bond breaking. PDMATLAB2D provides an entry-level peridynamics computational tool for educational and training purposes. It also serves as an accessible and easily modifiable computational tool for peridynamics researchers who would like to adapt the code for a multitude of peridynamics simulations.

## Run

Examples can be run by passing an input file to the top-level function:
```
PDMATLAB2D('WavePropagation')
PDMATLAB2D('CrackBranching')
```

Tests are located in the `Tests/` folder and can be run individually to check individual components.

Individual plots can be created (as in the manuscript below) from the `PlottingExamples/` folder. Plots can be easily generated within the folder via the function `PlotPaperFigures'.

## Citing PDMATLAB2D

If you use PDMATLAB2D in your work, please cite the following paper:
```
@article{seleson2022,
  title = {PDMATLAB2D: A Peridynamics MATLAB Two-Dimensional Code},
  year = {2022},
  author = {Pablo Seleson and Marco Pasetto and Yohan John and Jeremy Trageser and Samuel Temple Reeve},
  notes = {Under review}
}
```

If you would like to cite the software itself, cite the current release or version used from Zenodo.

## Contributing

We encourage you to contribute to PDMATLAB2D!

When contributing, please first discuss the changes you wish to make via issue (or directly with the PDMATLAB2D developers). Changes should be made through opening a pull request with `main` as the destination branch. Review from at least one developer is required.

PDMATLAB2D is distributed under the terms of the BSD 3-Clause license. All new contributions must be made under this license.

## License

PDMATLAB2D is distributed under a [BSD 3-Clause license](LICENSE).
