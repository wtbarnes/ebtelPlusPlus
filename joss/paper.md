---
title: "ebtelplusplus: Efficient Hydrodynamic Modeling of Coronal Loops"
tags:
- Python
- hydrodynamics
- astronomy
- astrophysics
- solar physics
authors:
- name: W. T. Barnes
  orcid: 0000-0001-9642-6089
  affiliation: "1,2"
- name: J. W. Reep
  orcid: 0000-0003-4739-1152
  affiliation: "3"
- name: N. Freij
  orcid: 0000-0002-6253-082X
  affiliation: "4,5"
- name: S. J. Schonfeld
  orcid: 0000-0002-5476-2794
  affiliation: "6"
- name: R. Collazzo
  orcid:
  affiliation: "1"
- name: P. J. Cargill
  orcid:
  affiliation: "7,8"
- name: S. J. Bradshaw
  orcid: 0000-0002-3300-6041
  affiliation: "9"
- name: J. A. Klimchuk
  orcid: 0000-0003-2255-0305
  affiliation: "2"
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
- name: Space and Atmospheric Physics, The Blackett Laboratory, Imperial College, London SW7 2BW, UK
  index: 7
- name: School of Mathematics and Statistics, University of St Andrews, St Andrews KY16 9SS, UK
  index: 8
- name: Department of Physics and Astronomy, Rice University, Houston, TX 77005, USA
  index: 9
date: 18 October 2025
bibliography: paper.bib
---

# Summary

Understanding the dynamics of the solar corona, the outermost layer of the Sun's atmosphere, requires detailed numerical modeling of the coronal plasma and magnetic field.
While solving the full set of three-dimensional *magnetohydrodynamic* (MHD) equations is extremely challenging for even small regions of the solar atmosphere, field-aligned hydrodynamic models [e.g. HYDRAD @bradshaw_self-consistent_2003] exploit the fact that the magnetic pressure in the corona dominates the gas pressure, meaning plasma dynamics in this region are largely confined to the direction of the magnetic field.
As such, the explicit dependence on the magnetic field can be neglected and the relevant mass, momentum, and energy equations can be reduced to a single-dimension in space, the coordinate along a magnetic field line, as well as time.
However, the large range of spatial and temporal scales necessary to resolve this system still makes these model computationally expensive enough to make large parameter space explorations difficult.
The enthalpy-based thermal evolution of loops (EBTEL) model [@klimchuk_highly_2008;@cargill_enthalpy-based_2012] mitigates this difficulty by computing spatial integrals of the aforementioned hydrodynamic equations.
Because of its relative simplicity and computational efficiency, EBTEL has been widely used since its initial development [e.g. @qiu_heating_2012;@raftery_multi-wavelength_2009;@ugarte-urra_determining_2014].
Comparisons to spatially-averaged results from field-aligned hydrodynamic models show very good agreement [@cargill_enthalpy-based_2012-1].

# Statement of Need

EBTEL consists of a set of ordinary differential equations that model the time-dependent behavior of the relevant thermodynamic quantities spatially-averaged over the corona.
To approximate mass transport within a corona loop, EBTEL equates an enthalpy flux with a balance between the heat flux out of the corona and the energy lost due to radiation in the transition region (TR), the thin layer of the solar atmosphere just below the corona.
If TR radiation cannot balance the downward heat flux, this drives an upflow of material into the corona and if the TR is radiating away more energy than the coronal heat flux can supply this drives a downflow.
This approximation is valid for short loops and bulk velocities below the local sound speed [@klimchuk_highly_2008].

The original software implementation of the EBTEL model of @klimchuk_highly_2008 was in the proprietary Interactive Data Language (IDL).
Subsequent improvements to the gravitational stratification and radiative losses by @cargill_enthalpy-based_2012 gave better agreement with field-aligned hydrodynamic models.
These improvements were also implemented in IDL[^ebtel2].
@barnes_inference_2016 extended the EBTEL model to relax the single-fluid assumption and treat electrons and ions separately, allowing for differential heating between the two species.
They implemented this extension of the model in C++.
@cargill_static_2021 later extended EBTEL to include effects due to cross-sectional area expansion and implemented this in IDL.
Additionally, @reep_modeling_2024 added the ability to vary the abundance model for the radiative losses as a function of time.
These modifications were made to the aforementioned C++ model.

`ebtelplusplus` unifies all of the aforementioned features into a single set of equations and provides an implementation of this model in C++ and Python.
In particular, `ebtelplusplus` solves the following equations for the spatially-averaged electron pressure ($p_e$), ion pressure ($p_i$), and density ($n$) of a semi-circular coronal loop of half-length $L$,

\begin{eqnarray}
\frac{1}{\gamma-1}\frac{dp_e}{dt} &=& Q_e + \frac{\psi_c}{L_*}\left(1+\frac{A_{TR}\psi_{TR}}{A_c\psi_c}\right) - \frac{R_c}{L_*}\left(1+c_1\frac{A_{TR}}{A_c}\right), \\
\frac{1}{\gamma-1}\frac{dp_i}{dt} &=& Q_i - \frac{\psi_c}{L_*}\left(1 + \frac{A_{TR}\psi_{TR}}{A_c\psi_c}\right), \\
\frac{dn}{dt} &=& -\frac{(\gamma - 1)\xi c_2}{(\xi + 1)\gamma c_3 k_B L_c T_e}\left(\frac{A_{TR}L_c}{A_cL_*}R_c\left(c_1 - \frac{L_{TR}}{L_c}\right) + \frac{A_0}{A_c}(F_{e,0} + F_{i,0})\right),
\end{eqnarray}

where $Q$ is the user-specified heating term, $A_{c,TR,0}$ is the cross-sectional area averaged over the corona and TR and at the TR/corona boundary, $R_c$ is the energy lost to radiation in the corona, $F_0$ is the conductive heat flux at the TR/corona boundary, and $L_*=L_c + (A_{TR}/A_c)L_{TR}$ with $L=L_c+L_{TR}$.
The remaining terms are explained more fully in the aforementioned publications and the `ebtelplusplus` documentation[^ebteldocsderivation].
Note that this set of equations is closed by equations of state for the electrons and ions: $p_e=k_BnT_e,p_i=k_BnT_i$.
The figure below shows example output from `ebtelplusplus` with different model parameters for the same time-dependent heating function $Q$.

![Spatially-averaged hydrodynamic evolution of the temperature (top right), density (bottom left), and temperature-density phase space (bottom right) of a coronal loop with half-length $L=40$ Mm for five different cases with the same heating input (top left panel). In the nominal case, (blue) the electron and ion populations are kept in equilibrium, the cross-sectional area of the loop does not expand, and the radiative losses are determined by a power-law function. If the electrons (solid) and ions (dashed) are allowed to evolve separately, heating only the electrons (orange) results in the ions taking about 250 s to fully equilibrate with the electrons while heating only the ions (green) results in the ions becoming over three times hotter than the electrons due to the relative inefficiency of ion thermal conduction. Incorporating area expansion through the corona (red) leads to a higher peak temperature and more delayed peak in the density while calculating the radiative losses using a time-varying abundance (purple) leads to a temperature evolution similar to the nominal case, but a higher peak density.](figure.pdf)

`ebtelplusplus` solves the above equations using a Runge-Kutta Cash-Karp integration method [see section 16.2 of @press_numerical_1992] and an (optional) adaptive time-stepping scheme to ensure the principal physical timescales are resolved at each phase of the loop evolution[^boost].
It is implemented in C++ for computational efficiency and is wrapped in Python using `pybind11` [@wenzel_jakob_2025_16929811] to enable easier installation and a user-friendly API.
Precompiled binary wheels are provided for all major operating systems and the versions of Python recommended by SPEC 0[^spec0].
Where appropriate, all inputs and outputs are expressed as `astropy.units.Quantity` objects [@astropy_collaboration_astropy_2022].
Additionally, `ebtelplusplus` is nearly two orders of magnitude faster than the previous IDL implementation because it is implemented in C++ and due to its use of an adaptive time-stepping scheme.

# Other Implementations

The aforementioned IDL implementation, which includes features described in @cargill_enthalpy-based_2012 and @cargill_static_2021, is referred to as `EBTEL-IDL` [@cargill_2024_13351770].
@rajhans_flows_2022 relaxed the assumption of subsonic flows in EBTEL such that the Mach numbers and velocities produced are in better agreement with field-aligned hydrodynamic simulations.
The IDL software implementation of this model is referred to as [`EBTEL3-IDL`](https://github.com/rice-solar-physics/EBTEL3).
`EBTEL3-IDL` uses an adaptive time grid to ensure the appropriate timescales are resolved in the impulsive phase.
As such, there are currently three separate though slightly different implementations of the EBTEL model.
The table below summarizes the features included in each implementation.

| Feature                     | Citation               | `EBTEL-IDL` | `EBTEL3-IDL` | `ebtelplusplus` |
|:----------------------------|:-----------------------|:-----------:|:------------:|:---------------:|
| Decouple electrons and ions | @barnes_inference_2016 | no          | no           | yes             |
| Adaptive time step          | @barnes_inference_2016 | no          | yes          | yes             |
| Area expansion              | @cargill_static_2021   | yes         | no           | yes             |
| Supersonic flows            | @rajhans_flows_2022    | no          | yes          | no              |
| Time-variable abundances    | @reep_modeling_2024    | no          | no           | yes             |

# References

[^ebteldocsderivation]: A full derivation of the set of EBTEL equations presented here can be found in the [`ebtelplusplus`](https://ebtelplusplus.readthedocs.io/en/stable/topic_guides/derivation.html) documentation.
[^boost]: The Runge-Kutta Cash-Karp integrator, including the adaptive time-stepping, is provided by the [Boost Odeint library](https://www.boost.org/library/latest/numericodeint/).
[^spec0]: [Scientific Python Ecosystem Coordination (SPEC) 0](https://scientific-python.org/specs/spec-0000/) recommends a set of minimum supported dependencies for packages commonly used across the scientific Python ecosystem, including Python.
[^ebtel2]: This version is sometimes referred to as “EBTEL2”.
