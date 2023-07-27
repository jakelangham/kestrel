---
title: 'The `Kestrel` software for simulations of morphodynamic Earth-surface flows'
tags:
  - Fortran
  - geophysics
  - shallow water
  - debris flow
  - morphodynamics
authors:
  - name: Jake Langham
    orcid: 0000-0001-9857-7016
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Mark J. Woodhouse
    orcid: 0000-0002-2198-6791
    equal-contrib: true
    affiliation: 2
affiliations:
 - name: School of Mathematics, Fry Building, University of Bristol, Bristol, BS8 1UG, UK
   index: 1
 - name: School of Earth Sciences, Wills Memorial Building, University of Bristol, Bristol, BS8 1RJ, UK
   index: 2
date: \today
bibliography: paper.bib

---

# Summary

`Kestrel` is a program for simulating flows composed of a
mixture of fluid and sediment. It includes the facility to model material
exchange with the topography over which the flow propagates, by incorporating
sediment entrainment and deposition. These physical processes, which mutually
couple the flow with its underlying bed, are sometimes collectively termed
'morphodynamics'.  Simulations may be initiated either on simple surfaces or on
more realistic terrains via a user-specified digital elevation model (DEM). The
latter option enables computations on topographies measured to approximate the
Earth's surface, so that real world events may be reconstructed and potential
future scenarios may be modelled. `Kestrel` has been primarily developed for
Earth sciences research into natural hazards, including volcanic mudflows, flash
floods and landslides. However, it may also be useful for modelling flows of
interest to engineers, applied mathematicians, geophysicists and industry
scientists, see e.g.\ @Capart:1998, @Cao:2004, @Iverson:2015 and
@Langham:2021. The versatility of the code is a deliberate design choice. As
discussed below, many of the key physical processes are implemented in a modular
way, allowing the user to choose between different options, depending on the
problem. Furthermore, for expert users it should be relatively straightforward
to extend the code to support alternative modelling terms that suit individual
needs.

`Kestrel` is predominantly written in `Fortran`, with some `C++` for handling
geospatial data via external libraries.  While expertise in the scientific
background is required to set up simulations and correctly interpret their
results, the program is otherwise intended to be straightfoward to use for a
broad range of scientists.  It has relatively few dependencies (`GDAL`, `PROJ`
and optionally, `NetCDF`), making it easy to build on modern Unix-like
platforms.  After installation, simulations are prepared by writing an input
file specifying suitable parameter choices and run on the command line.
Solution fields are saved at regular intervals, together with spatial maps of
their maximums over the whole simulation. Output is via text file or `NetCDF`
(preferred).  In the latter case, we provide an extension to the open-source
`QGIS` software, which imports `Kestrel` solutions at their georeferenced
coordinates and prepares each of the data fields for visualisation.  This
provides a particularly convenient workflow for geoscientists.  An example of
this output is given in \autoref{fig:example}.

![An example of georeferenced `Kestrel` output, postprocessed in `QGIS`.
The panels show (left-to-right) flow depth, flow speed, solids fraction,
and topographic elevation change of a small volume flow released onto steep
topography, here using `SRTM` data.
\label{fig:example}](example.png)

# Statement of need

All fluid flows that propagate over the Earth's surface transport sediment to
some degree. The presence of sediment at sufficiently high concentrations
substantially complicates the physics of these flows, by modifying their
density, rheology and their ability to entrain or deposit further volumes
of sediment [@Iverson:1997;@Iverson2015]. In many regions of the world, local conditions can trigger
destructive flowing fluid--sediment mixtures that travel over distance of up to tens of
kilometres [@Pierson:1990;@Scott:2005;@Jakob:2005]. Driven by the need to understand their fundamental physics and
ultimately to create predictive tools that can mitigate hazards, the development
of mathematical models for these flows is an active research area and there is a
corresponding need for flexible research codes that can numerically solve these
models. 

There are many existing codes available for simulating these systems. The
programs most closely related to ours are the open source codes `D-CLAW`
[@Iverson:2014;@George:2014;@DCLAW:2023], `TITAN2D`
[@Patra:2005;@Simakov:2019;@Patra:2020], `r.avaflow`
[@Pudasaini:2019;@ravaflow:2023], `IMEX_SfloW2D`
[@IMEX_SfloW2D:2023;@Vitturi:2023a;@Vitturi:2023b] and the debris flow
components of the proprietary software packages `RAMMS`
[@Christen:2010;@Meyrat:2022] and `Flo-2D` [@Flo2D:2023].  Each of these codes
uses a slightly different description of the flow physics and underlying
mathematical framework. In some cases, it can be difficult to discern from the
available documentation exactly which model and assumptions are used in the
latest software version and how the program operates `under the hood' to produce
its results. 

Like all the above software, `Kestrel` numerically approximates solutions to an
underlying system of partial differential equations for the flow, whose
derivation uses the fact that the flow's lateral extent typically greatly
exceeds its thickness, to reduce the spatial dimension by 1, thereby rendering
simulations computationally tractable at geophysical scales. `Kestrel` supports
simulations in either one or two orthogonal coordinate directions, perpendicular
to gravity. It keeps track of the following
observables, which depend on space and time coordinates $\mathbf{x}$ and $t$: the flow thickness $H(\mathbf{x},t)$, velocity
$\mathbf{u}(\mathbf{x},t)$, volumetric solids concentration $\psi(\mathbf{x},t)$, density
$\rho(\mathbf{x},t)$ and the bed elevation $b(\mathbf{x},t)$. 
In two spatial dimensions, the governing equations are
\begin{gather}\frac{\partial H}{\partial t} +
\nabla\cdot(H\mathbf{u}) = \mathcal{E} - \mathcal{D} + \mathcal{Q}_H,\label{eq:governing eqs 1}\\
\frac{\partial~}{\partial t}(H\psi) +
\nabla\cdot(H\mathbf{u}\psi) = \psi_b (\mathcal{E} - \mathcal{D}) +
\mathcal{Q}_H\mathcal{Q}_{\psi},\\
\frac{\partial ~}{\partial t}(\rho H \mathbf{u}) + 
\nabla\cdot(\rho H\mathbf{u}\otimes\mathbf{u}) + 
\frac{g}{2\cos\theta}\nabla_s\left(\rho H^2\cos^2\theta \right) = 
-\rho H \nabla_s b - \frac{\rho\mathbf{u}}{|\mathbf{u}|}\mathcal{F} + \nabla\cdot(\nu \rho H \nabla \mathbf{u}), \\
\frac{\partial b}{\partial t} = \frac{\mathcal{D} - \mathcal{E}}{\cos\theta},
\label{eq:governing eqs 4}
\end{gather}
where $\psi_b$, $g$ and $\nu$ are user-defined modelling parameters,
$\theta(\mathbf{x},t)$ is the local slope angle between the bed normal and
gravity, and $\nabla_s = \nabla - \mathbf{s}(\mathbf{s}\cdot\nabla)$, 
with $\mathbf{s} \equiv \cos(\theta)\nabla b$.

The technical details of these equations, their derivation and our numerical
solution scheme are fully presented by @Langham:2023.  While most of the terms
are fixed by the underlying depth-averaged flow physics (and shall not be
discussed further), some parts of the right-hand sides are user-settable.  The
terms $\mathcal{F}$, $\mathcal{E}$ and $\mathcal{D}$ denote the basal friction,
erosion rate and deposition rate respectively. These are modelling closures,
assumed to be functions of the flow fields $H$, $\mathbf{u}$ and $\psi$. In each
case, the user may choose from different options, depending on the problem at
hand. For example, the friction $\mathcal{F}$ may be set either to a function
appropriate for turbulent fluids, to various models of purely granular flows, or
to a combined law that depends on the solids concentration.  This provides the
flexibility to simulate many different kinds of flow.  Furthermore, it is worth
noting that in many cases, the question of which closures most faithfully
capture the flow physics is an open problem that cannot easily be addressed
experimentally.  Using numerical simulations to investigate the effects of
different modelling choices is one way to approach this.

The remaining source terms $\mathcal{Q}_H$ and $\mathcal{Q}_{\psi}$ are
time-dependent functions that provide one way for a modeller to control fluxes
of material input into the simulation by specifying their form via time series
data and identifying regions over which these fluxes apply. The flux values may
be informed by geophysical measurements, expert judgement, or chosen according
to some other consideration. Simulations may also be initiated by constructing
an initial flow state which is evolved forward in time. `Kestrel` accepts
arbitrary initial conditions given in the same format as its result files, or
simple initial volumes of material (such as cubes and cylinders) can be
specified via an input file.

In deriving Eqs.\ \eqref{eq:governing eqs 1}--\eqref{eq:governing eqs 4}, some
physical assumptions are made that simplify our model compared to some of the
similar codes cited above, such as the neglect of interstitial pressure between
sediment particles and the assumption that the sediment is composed of particles
that are all roughly the same size.  Conversely, we include the effects of
morphodynamics, which not all codes support. The reasoning for these choices is
twofold: (1) to simplify the problem for modellers wishing to conduct
simulations, without compromising the most essential flow physics and (2) in
order to focus on making some important technical advances for morphodynamic
models, which are particular to our code and detailed by @Langham:2023. These
include improvements to standard numerical schemes, a careful treatment of the
effects of the basal geometry and a regularisation that ensures the model is
well-posed as an initial value problem.

Despite the relative simplicity of our approach, morphodynamic simulations
involve many free parameters, which are not all directly measurable for real
world scenarios. Proper calibration and an appreciation of the uncertainties
present at each stage of the modelling process is essential to obtain reliable
results from `Kestrel`.  Examples and guidance for getting started may be found
in our documentation.  Results from `Kestrel` simulations have thus far been
used in the following scientific publications:
[@Jenkins:2023;@Langham:2023]. `Kestrel` also forms the backend for the
`LaharFlow` volcanic hazard model [@Laharflow:2023].

# Acknowledgements

The development of this software spanned research grants from the the Newton
Fund (NE/S00274X/1), EPSRC Impact Acceleration Account (EP/X525674/1) and a NERC
Knowledge Exchange Fellowship (NE/R003890/1).  We thank our colleagues Andrew J.
Hogg, Luke T. Jenkins and Jeremy C. Phillips, who worked on the foundational
mathematical and geological ideas that underpin the code. We also thank users
who supported the development of `Kestrel`, especially Felipe Flores
(SERNAGEOMIN, Chile), Francisco Vasconez (IG-EPN, Ecuador), Selwyn Cabaluna
(PHIVOLCS, Phillipines) and other geological collaborators at SERNAGEOMIN,
IG-EPN, PHIVOLCS, INGEMMET (Perú), as well as all current and past `Kestrel`
users.

# References

