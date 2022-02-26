# [Re] The Evolution of Virulence in Pathogens with Vertical and Horizontal Transmission

## 2021 Replication of Lipsitch et al (1996)

This project is a Julia 1.6.0 replication of Lipsitch et al (1996), a paper in evolutionary modelling that challenged notions on the dynamics of viral transmission. Full reference to the original article:

> M. Lipsitch, S. Siller, and M. A. Nowak. “THE EVOLUTION OF VIRULENCE IN PATHOGENS WITH VERTICAL ANDHORIZONTAL TRANSMISSION.” In:Evolution50.5 (1996), pp. 1729–1741. eprint:https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1558-5646.1996.tb03560.x.

The goal of this project was to replicate the epidemiological model for the evolution of virulence presenting trade-offs between horizontal and vertical transmission of pathogens. Our results can be viewed in [ReScience C](https://rescience.github.io/bibliography/Dhane_2021.html) or on [Zenodo](https://zenodo.org/record/5722853#.YZ62ZS3b3-Y). The open review is available [here](https://github.com/ReScience/submissions/issues/56#issuecomment-977556264).

Julia Lang scripts defining model parameters (as found in the original study) are available in the folder ```Model_scripts```.

All figures produced in the replication study can be found in the ```Article/figures``` subfolder.

To run our replication code, clone our repository to your computer, checkout to the master branch, and run the following lines:

``` julia
] activate .
] instantiate
```

Finally, run the ```00run_project.jl``` file in the ```Model_scripts/``` folder.
