---
title: "ebtelplusplus: Efficient Hydrodynamic Modeling of Coronal Loops"
tags:
- Python
- hydrodynamics
- astronomy
- astrophysics
- solar physics
authors:
- name: Will T. Barnes
  orcid: 0000-0001-9642-6089
  affiliation: "1, 2"
affiliations:
- name: Department of Physics, American University, Washington, D.C. 20016, USA
  index: 1
- name: NASA Goddard Space Flight Center, Greenbelt, MD 20771, USA
  index: 2
date: 18 October 2025
bibliography: paper.bib
---

# Summary

# Statement of Need

# The EBTEL Model

![Spatially-averaged hydrodynamic evolution of the temperature (top right), density (bottom left), and temperature-density phase space (bottom right) of a coronal loop with half-length $L=40$ Mm as modeled with `ebtelplusplus` for five different configurations. In all cases, the heating input (top left panel) is the same: a triangular function of total duration 200 s, starting at $t=250$ s. In the nominal case, (blue) the electron and ion populations are kept in equilibrium, the cross-sectional area of the loop does not expand, and the radiative losses are determined by a power-law function. If the electrons (solid) and ions (dashed) are allowed to evolve separately, heating only the electrons (orange) results in the ions taking about 250 s to fully equilibrate with the electrons while heating only the ions (green) results in the ions becoming over three times hotter than the electrons due to the relative inefficiency of ion thermal conduction. Incorporating area expansion through the corona (red) leads to a higher peak temperature and more delayed peak in the density while calculating the radiative losses using a time-varying abundance (purple) leads to a temperature evolution similar to the nominal case, but a higher peak density.](ebtel_example_figure.pdf)

# Comparison to Other Implementations

# References
