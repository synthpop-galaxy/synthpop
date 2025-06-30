---
title: 'SynthPop: A New Framework for Synthetic Milky Way Population Generation'
tags:
  - Python
  - astronomy
  - milky way
  - population synthesis
authors:
  - name: Jonas Kl\"{u}ter
    orcid: 0000-0002-3469-5133
    equal-contrib: true
    affiliation: 1
  - name: Macy J. Huston
    orcid: 0000-0003-4591-3201
    equal-contrib: true
    affiliation: 2
  - name: Abigail Aronica
    affiliation: 3
  - name: Samson A. Johnson
    orcid: 0000-0001-9397-4768
    affiliation: "3, 4"
  - name: Matthew Penny
    orcid: 0000-0001-7506-5640
    affiliation: 1
  - name: Marz Newman
    orcid: 0009-0002-1973-5229
    affiliation: 1
  - name: Farzaneh Zohrabi
    orcid: 0000-0003-2872-9883
    affiliation: 1
  - name: Alison L. Crisp
    orcid: 0000-0003-4310-3440
    affiliation: 3
  - name: Allison Chevis
    orcid: 0000-0003-2558-1748
    affiliation: 5
affiliations:
  - name: Department of Physics \& Astronomy, Louisiana State University, Baton Rouge, LA 70803, USA
    index: 1
  - name: Astronomy Department, University of California, Berkeley, CA, 94720, USA
    index: 2
  - name: Department of Astronomy, The Ohio State University, Columbus, OH 43210, USA
    index: 3
  - name: NASA Jet Propulsion Laboratory, Pasadena, CA, 91109, USA
    index: 4
  - name: Department of Physics, Pennsylvania State University, State College, PA 16802, USA
    index: 5
date: 02 December 2024
bibliography: joss-paper.bib

aas-doi: xx.xxxx/xxxxx
aas-journal: xxxxxxx
---

# Summary

SynthPop is a new open source, modular population synthesis Galactic modeling software to simulate catalogs of Milky Way stars along any sightline outward from the Sun. 
Motivated by a lack flexibility in existing Galactic models, SynthPop is coded entirely in python, can be run standalone or as an imported module, and is configured by json files that allow different model components to be switched out as desired.

# Statement of need

We have built SynthPop to address the need we perceived for a publicly available population synthesis code that can be easily modified, provides an interface that can be run via script, and ideally can allow comparisons between different models within a single framework. To achieve these goals, we built the package in python, which is likely the most widely used language in astronomy. This choice results in slower performance than a compiled language, but a minimal bar to code modification, and an easy way to import the package into other code. We have designed the package so that it can be modified either through the adjustment of parameter files or by adding new code modules to achieve results that are not enabled by existing tools.

The primary driving use case for SynthPop is to provide synthetic stellar population catalogs for microlensing simulations (@Penny2013, @Penny2019, @Johnson2020).

# Ongoing use

Huston et al. (in prep.) will illustrate the development of a SynthPop Model version that matches well to data, and apply the model to the Roman Galactic Bulge Time Domain Survey to explore the Galactic distributions of anticipated microlensing event lenses and sources. 
This model version is also being used as the Galactic model input for gulls (@Penny2013, @Penny2019) and PyLIMASS (@Bachelet2024) by the Roman Galactic Exoplanet Survey Project Infrastructure Team for updated exoplanet mass yield estimates  and field optimization (Terry et al., in prep.; Zohrabi et al., in prep.).

# Acknowledgements

We appreciate conversations with N. Koshimoto in implementing the (@Koshimoto2021) model and (@Sormani2022) nuclear stellar disk into SynthPop.
We also thank the Roman Galactic Exoplanet Survey Project Infrastructure team members and others who helped improve the software through discussions, use during development, and feedback.

J.K. acknowledges support from NASA award NNX16AC62G and the Gordon and Betty Moore Foundation award GBMF10467.

M.J.H. acknowledges support from the Heising-Simons Foundation under grant No. 2022-3542. 
Early work on this software by M.J.H. and A.A. was supported by the Ohio State University Summer Undergraduate Research Program.

S.A.J.’s work was supported by NASA Grant 80NSSC24M0022 and an appointment to the NASA Postdoctoral Program at the NASA Jet Propulsion Laboratory, administered by Oak Ridge Associated Universities under contract with NASA. 

M.T.P. acknowledges support from NASA awards NNX14AF63G, NNX16AC62G, 80NSSC24M0022, 80NSSC24K0881, and Louisiana Board of Regents Support Fund (RCS Award Contract Simple: LEQSF(2020-23)-RD-A-10). 
Early work by M.T.P. was performed in part under contract with the Jet Propulsion Laboratory (JPL) funded by NASA through the Sagan Fellowship Program executed by the NASA Exoplanet Science Institute.

M.N. acknowledges support from NASA award 80NSSC24K0881.
A.L.C. and F.Z. acknowledge support from NASA grant 80NSSC24M0022. 
A.C.'s work was supported by the National Science Foundation under Grant Number 1852454.
Portions of this research were conducted with high performance computational resources provided by Louisiana State University ([http://www.hpc.lsu.edu](http://www.hpc.lsu.edu)).

# References
