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
- name: Jeffrey W. Reep
  orcid:
  affiliation: "3"
- name: Nabil Freij
  orcid:
  affiliation: "4,5"
- name: Sam Schonfeld
  orcid:
  affiliation: "6"
- name: Rowan Collazzo
  orcid:
  affiliation: "1"
affiliations:
- name: Department of Physics, American University, Washington, D.C. 20016, USA
  index: 1
- name: NASA Goddard Space Flight Center, Greenbelt, MD 20771, USA
  index: 2
- name: Institute for Astronomy, University of Hawai'i at Manoa, Pukalani, HI 96768, USA
  index: 3
- name: SETI Institute, Mountain View, CA 94043, USA
  index: 4
- name: Lockheed Martin Solar and Astrophysics Laboratory, Palo Alto, CA 94304, USA
  index: 5
- name: Air Force Research Laboratory, Space Vehicles Directorate, Kirtland AFB, NM 87117, USA
  index: 6
date: 18 October 2025
bibliography: paper.bib
---

# Summary

# The EBTEL Model

Understanding the dynamics of the solar corona, the outermost layer of the Sun's atmosphere, requires detailed numerical modeling of the coronal plasma and magnetic field.
As the solar corona is a magnetized plasma, its dynamics are governed by a set of *magnetohydrodynamic* equations that describe the evolution of the magnetic field and mass, momentum, and energy of the plasma in time and three spatial dimensions.
In general, solving this set of equations is extremely challenging and requires high-performance computing resources and many CPU hours to model even small regions of the solar atomosphere.
So-called field-aligned hydrodynamic models [e.g. HYDRAD; @bradshaw_self-consistent_2003] exploit the fact that magnetic pressure in the corona dominates the gas pressure such that plasma dynamics in this region are largely confined to the direction of the field.
As such, the explicit dependence on the magnetic field can be neglected and the relevant mass, momentum, and energy equations can be reduced to a single-dimension in space, the coordinate along a magnetic field line, as well as time.
This allows these models to achieve much finer spatial resolution than three-dimensional magnetohydrodynamics [e.g. @bradshaw_influence_2013] at relatively little computational cost even on a laptop or desktop computer.
However, due to the large range of spatial and temporal scales that are necessary to resolve this system, solving even this reduced set of equations is computationally expensive enough to make large parameter space explorations difficult.

The enthalpy-based thermal evolution of loops (EBTEL) model [@klimchuk_highly_2008;@cargill_enthalpy-based_2012] mitigates this difficulty by computing spatial integrals of the aforementioned hydrodynamic equations, yielding a set of ordinary differential equations that model the time-dependent behavior of the temperature, density, and pressure spatially-averaged over the corona.
EBTEL accomplishes this by equating an enthalpy flux with a balance between the heat flux out of the corona and the energy lost due to radiation in the transition region.
If transition region radiation cannot balance the downward heat flux, this drives an upflow of material into the corona and if the transition region is radiating away more energy than the coronal heat flux can supply this drives a downflow.
This approximation is valid for short loops and bulk velocities below the local sound speed [@klimchuk_highly_2008].
Because of its relative simplicity and computational efficiency, EBTEL has been widely used since its initial development [e.g. @raftery_multi-wavelength_2009;@qiu_heating_2012;@ugarte-urra_determining_2014].
Comparisons to spatially-averaged results from field-aligned hydrodynamic models show very good agreement [@cargill_enthalpy-based_2012-1].

# Statement of Need

The EBTEL model was originally developed by @klimchuk_highly_2008 and @cargill_enthalpy-based_2012 in the proprietary Interactive Data Language (IDL).
@barnes_inference_2016 extended the EBTEL model to treat electrons and ions separately and developed the initial version of `ebtelplusplus` in the C++ programming language.


`ebtelplusplus` solves the following equations for the spatially-averaged electron pressure ($p_e$), ion pressure ($p_i$), and density ($n$) of a semi-circular coronal loop of half-length $L$,

\begin{eqnarray}
\frac{1}{\gamma-1}\frac{dp_e}{dt} &=& Q_e + \frac{\psi_c}{L_*}\left(1+\frac{A_{TR}\psi_{TR}}{A_c\psi_c}\right) - \frac{R_c}{L_*}\left(1+c_1\frac{A_{TR}}{A_c}\right), \\
\frac{1}{\gamma-1}\frac{dp_i}{dt} &=& Q_i - \frac{\psi_c}{L_*}\left(1 + \frac{A_{TR}\psi_{TR}}{A_c\psi_c}\right), \\
\frac{dn}{dt} &=& -\frac{(\gamma - 1)\xi c_2}{(\xi + 1)\gamma c_3 k_B L_c T_e}\left(\frac{A_{TR}L_c}{A_cL_*}R_c\left(c_1 - \frac{L_{TR}}{L_c}\right) + \frac{A_0}{A_c}(F_{e,0} + F_{i,0})\right),
\end{eqnarray}

where $Q$ is the *ad hoc* heating term,  $A_{c,TR,0}$ is the cross-sectional area averaged over the corona and transition region and at the transition region/corona boundary, $R_c$ is the energy lost to radiation in the corona, $F_0$ is the conductive heat flux at the transition region/corona boundary, and $L_*=L_c + (A_{TR}/A_c)L_{TR}$ with $L=L_c+L_{TR}$.
The remaining terms are explained more fully in the aforementioned publications.
Note that this set of equations is closed by an equation of state for the electrons and ions: $p_e=k_BnT_e,p_i=k_BnT_i$.

![Spatially-averaged hydrodynamic evolution of the temperature (top right), density (bottom left), and temperature-density phase space (bottom right) of a coronal loop with half-length $L=40$ Mm as modeled with `ebtelplusplus` for five different configurations. In all cases, the heating input (top left panel) is the same: a triangular function of total duration 200 s, starting at $t=250$ s. In the nominal case, (blue) the electron and ion populations are kept in equilibrium, the cross-sectional area of the loop does not expand, and the radiative losses are determined by a power-law function. If the electrons (solid) and ions (dashed) are allowed to evolve separately, heating only the electrons (orange) results in the ions taking about 250 s to fully equilibrate with the electrons while heating only the ions (green) results in the ions becoming over three times hotter than the electrons due to the relative inefficiency of ion thermal conduction. Incorporating area expansion through the corona (red) leads to a higher peak temperature and more delayed peak in the density while calculating the radiative losses using a time-varying abundance (purple) leads to a temperature evolution similar to the nominal case, but a higher peak density.](figure.pdf)

- Unifies several different features implemented separately
- Provides a version without the need for proprietary software
- Provides a canonical version of EBTEL
- Provides adaptive time-stepping
- Provides a fast version of the software (written in C++, wrapped in Python).

# Other Implementations

# References
