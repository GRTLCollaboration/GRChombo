---
title: 'GRChombo: An adaptable numerical relativity code for fundamental physics'
tags:
  - c++
  - MPI
  - Open MP
  - vector intrinsics
  - gravity
  - general relativity
  - numerical relativity
authors:
- name: Tomas Andrade
  orcid: 0000-0002-2472-114X
  affiliation: 1
- name: Llibert Areste Salo
  orcid: 0000-0002-3812-8523
  affiliation: 2
- name: Josu Aurrekoetxea
  orcid: 0000-0001-9584-5791
  affiliation: 3
- name: Jamie Bamber
  orcid: 0000-0001-7181-3365
  affiliation: 4
- name: Katy Clough
  orcid: 0000-0001-8841-1522
  affiliation: 4
- name: Robin Croft
  orcid: 0000-0002-1236-6566
  affiliation: 5
- name: Eloy de Jong
  orcid: 0000-0002-4505-0808
  affiliation: 3
- name: Amelia Drew
  orcid: 0000-0001-8252-602X
  affiliation: 5
- name: Alejandro Duran
  #orcid: tbc
  affiliation: 6
- name: Pedro G. Ferreira
  orcid: 0000-0002-3021-2851
  affiliation: 4
- name: Pau Figueras
  orcid: 0000-0001-6438-315X
  affiliation: 2
- name: Hal Finkel
  orcid: 0000-0002-7551-7122
  affiliation: 7
- name: Tiago Fran\c{c}a
  orcid: 0000-0002-1718-151X
  affiliation: 2
- name: Bo-Xuan Ge
  orcid: 0000-0003-0738-3473
  affiliation: 3    
- name: Chenxia Gu
  orcid: 0000-0001-9537-6139
  affiliation: 2
- name: Thomas Helfer
  orcid: 0000-0001-6880-1005
  affiliation: 8
- name: Juha Jäykkä
  orcid: 0000-0002-5929-3931
  affiliation: 5
- name: Cristian Joana
  orcid: 0000-0003-4642-3028
  affiliation: 9
- name: Markus Kunesch
  orcid: 0000-0003-3818-7897
  affiliation: 5
- name: Kacper Kornet
  orcid: 0000-0001-8376-8231
  affiliation: 5
- name: Eugene A. Lim
  orcid: 0000-0002-6227-9540
  affiliation: 3
- name: Francesco Muia
  orcid: 0000-0002-4872-6172
  affiliation: 5
- name: Zainab Nazari
  orcid: 0000-0002-2955-7262
  affiliation: 10, 11
- name: Miren Radia
  orcid: 0000-0001-8861-2025
  affiliation: 5
- name: Justin Ripley
  orcid: 0000-0001-7192-0021
  affiliation: 5
- name: Paul Shellard
  #orcid: tbc
  affiliation: 5
- name: Ulrich Sperhake
  orcid: 0000-0002-3134-7088
  affiliation: 5
- name: Dina Traykova
  orcid: 0000-0002-3451-0987
  affiliation: 4
- name: Saran Tunyasuvunakool
  orcid: 0000-0002-1620-6797
  affiliation: 5
- name: Zipeng Wang
  orcid: 0000-0002-4745-8209
  affiliation: 8
- name: James Widdicombe
  orcid: 0000-0003-2269-4544
  affiliation: 3
- name: Kaze Wong
  orcid: 0000-0001-8432-7788
  affiliation: 8
affiliations:
- name: Departament de Fısica Quantica i Astrofisica, Institut de Ciencies del Cosmos, Universitat de Barcelona, Marti i Franques 1, 08028 Barcelona, Spain
  index: 1
- name: School of Mathematical Sciences, Queen Mary University of London, Mile End Road, London E1 4NS, United Kingdom
  index: 2
- name: Theoretical Particle Physics and Cosmology, King's College London, Strand, London, WC2R 2LS, United Kingdom
  index: 3  
- name: Astrophysics, Oxford University, Denys Wilkinson Building, Keble Road, Oxford OX1 3RH, United Kingdom
  index: 4
- name: Department of Applied Mathematics and Theoretical Physics (DAMTP), University of Cambridge, Centre for Mathematical Sciences, Wilberforce Road, Cambridge CB3 0WA, United Kingdom
  index: 5
- name: Intel Iberia, Torre Picasso Plaza Pablo Ruiz Picasso 1 Madrid, 28020 Spain
  index: 6
- name: Argonne National Laboratory (ANL), 9700 S. Cass Avenue, Argonne, IL 60439-4815, United States
  index: 7
- name: Henry A. Rowland Department of Physics & Astronomy, Johns Hopkins University, 3701 San Martin Drive, Baltimore, Maryland (MD) 21218, United States
  index: 8
- name: Cosmology, Universe and Relativity at Louvain (CURL), Institut de Recherche en Mathematique et Physique, University of Louvain, 2 Chemin du Cyclotron, 1348 Louvain-la-Neuve, Belgium
  index: 9
- name: Department of Physics, Bogazici University, 34342 Bebek, 80820 Istanbul, Turkey
  index: 10
- name: HECAP Section, Abdus Salam International Centre for Theoretical Physics (ICTP), 34151, Trieste, Italy
  index: 11
date: 1 June 2021
bibliography: paper.bib

---

# Summary

The 2015 detection of gravitational waves (GWs) from a binary black hole merger [@Abbott:2016blz] was a breakthrough moment for science. More detections have since been made by the Advanced LIGO/Virgo network [@TheVirgo:2014hva;@TheLIGOScientific:2014jea;@Aasi:2013wya] and future ground and space based detectors [@Somiya:2011np;@Saleem:2021iwi;@Audley:2017drz;@Luo:2015ght;@Hu:2017mde] will further expand our reach. 

Strong gravity regimes are described by the *Einstein Field Equation* (EFE) of General Relativity [@Einstein:1916vd]. 
\begin{equation}
   R_{\mu \nu} - \frac{1}{2} R g_{\mu \nu} = 8 \pi G T_{\mu\nu} ~.
\end{equation}
Analytic solutions to the EFE only exist where there is a high degree of symmetry; in general the equations must be solved numerically. The need for observational predictions has thus led to the development of *numerical relativity* (NR), methods for numerically solving the above expression, typically utilising high performance computing (HPC) resources.
Expanding out the tensorial notation, the EFE is a set of second order partial differential equations for the metric tensor field $g_{\mu\nu}$, which describes the curvature of spacetime in the presence of matter with stress-energy $T_{\mu\nu}$, i.e.,
\begin{equation}
   \partial_t \partial_t g_{\mu\nu} \sim \partial_x \partial_x g_{\mu\nu} + \partial_y \partial_y g_{\mu\nu} + \partial_z \partial_z g_{\mu\nu} + {\rm non ~ linear ~ cross ~ terms} 
          + 8 \pi G T_{\mu\nu}
\end{equation}
where the indices $\mu, \nu$ run over the spacetime indices -- in 4 dimensions, $t, x, y, z$. Given that $g_{\mu\nu}$ is symmetric in its indices, this gives a set of ten coupled non-linear partial differential equations, sourced by the stress-energy of any matter present in the spacetime.

One common approach to NR is to specify an initial spatial distribution for the metric and matter fields (subject to certain constraints), and then solve a time evolution for all metric and matter quantities, thus populating their values thoughout the four dimensional spacetime. The canonical example of this is the simulation of two black holes in orbit around each other, which permits extraction of the gravitational wave signal produced during the merger. Such numerical results have been instrumental in discovering signals in the noisy LIGO/VIRGO detector data, as well as confirming the predictions of GR to a high precision in the strong field regime.

GRChombo is an open-source code for performing such NR time evolutions, built on top of the publicly available Chombo software [@Adams:2015kgr] for the solution of PDEs. Whilst GRChombo uses standard techniques in NR, it focusses on applications in theoretical physics where adaptability, both in terms of grid structure, and in terms of code modification, are key drivers. 

# Key features of GRChombo

Since its initial announcement in 2015 [@Clough:2015sqa], the GRChombo code has become a fully mature, open source NR resource.

The key features of GRChombo are as follows:

- BSSN/CCZ4 formalism with moving punctures: GRChombo evolves the Einstein equation in the BSSN [@Nakamura:1987zz;@Shibata:1995we;@Baumgarte:1998te] or CCZ4 [@Gundlach:2005eh;@Alic:2011gg] formalism with conformal factor $\chi = det(\gamma_{ij})^{-1/3}$. Singularities of black holes are managed using the moving puncture gauge conditions [@Campanelli:2005dd;@Baker:2005vv], and Kreiss-Oliger dissipation is used to control errors, both from truncation and the interpolation associated with regridding.

- Boundary Conditions: The code implements periodic, Sommerfeld (radiative), extrapolating and reflective boundary conditions.

- Initial Conditions: The current examples provide analytic or semi-analytic initial data for black hole binaries, Kerr black holes and scalar matter. The code also incorporates a standalone version of the TwoPunctures code [@Ansorg:2004ds] for accurate binary BH data of arbitrary spins, masses and momenta.

- Diagnostics: GRChombo has routines for finding black hole horizons, calculating spacetime masses, angular momenta, densities, fluxes and extracting gravitational waves. 

- C++ class structure: GRChombo is written in the C++ language, and makes heavy use of object oriented programming (OOP) and templating.

- Parallelism: GRChombo uses hybrid OpenMP/MPI  parallelism with explicit vectorisation of the evolution equations via intrinsics, and is AVX-512 compliant. Our code strong scales efficiently to several thousand CPU-cores for a typical BH binary problem, and further for larger problem sizes. 

- Adaptive Mesh Refinement: The underlying Chombo code provides Berger-Oliger style [@Berger:1984zza] AMR with block-structured  Berger-Rigoutsos grid generation [@BergerRigoutsos]. The tagging of refinement regions is fully flexible and can be based on truncation error or other user-defined measures.

The code continues to be actively developed with a number of ongoing projects to add new features.

# Statement of Need

Several 3+1D NR codes using the moving puncture formulation already exist and are under active development. The Einstein Toolkit (http://einsteintoolkit.org/), with its related Cactus (http://cactuscode.org) [@Loffler:2011ay;@Schnetter:2003rb], and Kranc (http://kranccode.org) [@Husa:2004ip] infrastructure used by LEAN [@Sperhake:2006cy;@Zilhao:2010sr] and Canuda (https://bitbucket.org/canuda) [@Witek:2018dmd]. Other notable but non public codes include \texttt{BAM} [@Marronetti:2007ya;@Brugmann:2008zz], AMSS-NCKU [@Galaviz:2010mx], PAMR/AMRD and HAD [@Neilsen:2007ua;@East:2011aa]. Codes such as SPeC [@Pfeiffer:2002wt] and bamps [@Hilditch:2015aba] implement the generalised harmonic formulation of the Einstein equations using a pseudospectral method, and discontinuous Galerkin methods are used in SpECTRE (https://spectre-code.org) [@deppe_nils_2021_4734670;@Kidder:2016hev] (see also [@Cao:2018vhw]). NRPy (http://astro.phys.wvu.edu/bhathome) [@Ruchlin:2017com] is a code aimed for use on non HPC systems, which generate C code from Python, and uses adapted coordinate systems to minimise computational costs. CosmoGRaPH (https://cwru-pat.github.io/cosmograph) [@Mertens:2015ttp] and GRAMSES [@Barrera-Hinojosa:2019mzo] are among several NR codes targeted at cosmological applications (see [@Adamek:2020jmr] for a comparison) and which also employ particle methods. Simflowny (https://bitbucket.org/iac3/simflowny/wiki/Home) [@Palenzuela:2018sly], like CosmoGRaPH, is based on the SAMRAI infrastructure, and has targeted fluid and MHD applications. GRAthena++ [@Daszuta:2021ecf] makes use of oct-tree AMR to maximise scaling.

While GRChombo is not the only open source NR code, its unique features (detailed above) have made it one of the premier codes for numerical relativity, especially in the study of fundamental physics beyond standard binary mergers. In particular, GRChombo's highly flexible adaptive mesh refinement scheme allows for complicated "many-boxes-in-many-boxes" topology , enabling users to simulate non-trivial systems, such as ring-like configurations [@Figueras:2015hkb,@Helfer:2018qgv] and inhomogeneous cosmological spacetimes [@Joana:2020rxm;@Aurrekoetxea:2019fhr;@Clough:2017efm;@Clough:2016ymm]. Nevertheless, with its efficient scalability and AMR capabilities, it can also play a leading role in the continuing efforts to simulate ``standard'' binary mergers to the required sensitivities required for the upcoming LISA space mission [@Radia:2021hjs]. Finally, GRChombo's object-oriented and template-based code can be rapidly modified for non-standard problems such as higher dimensional spacetimes [@Figueras:2015hkb;@Figueras:2017zwa;@Bantilan:2019bvf;@Andrade:2020dgc], modified gravity systems [@Figueras:2020dzx] and additional fundamental fields [@Nazari:2020fmk;@Muia:2019coe;@Widdicombe:2019woy;@Clough:2018exo;@Dietrich:2018bvi;@Helfer:2018vtq;@Helfer:2016ljl;@Bamber:2020bpu;@Clough:2019jpm;@Alexandre:2018crg].

# Key research projects using GRChombo

The wide range of fundamental physics problems for which the code has been used so far includes:

- the simulation of inhomogeneous pre-inflationary spacetimes, bubble collisions and preheating in early universe cosmology [@Joana:2020rxm;@Aurrekoetxea:2019fhr;@Clough:2017efm;@Clough:2016ymm].


# Acknowledgements

The GRChombo collaboration gratefully acknowledges support to its members by the ERC, UKRI/STFC, Intel, The Royal Society, PRACE and DiRAC. 

In particular, JB, KC, PF and DT acknowledge support from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (grant agreement No 693024).
AD has been supported by an EPSRC iCASE Studentship in partnership with Intel (EP/N509620/1, Voucher 16000206) and is currently supported by a Junior Research Fellowship at Homerton College, Cambridge. F.M. is funded by a UKRI/EPSRC Stephen Hawking fellowship, grant reference EP/T017279/1. This work has been partially supported by STFC consolidated grant ST/P000681/1.

We thank the developers of the Chombo code for their assistance and guidance on using their code, and the Intel Parallel Computing Centre at the University of Cambridge for their support of our code development. We acknowledge the support of the Intel Visualization team, led by Jim Jeffers, notably the collaboration on in-situ visualization with Carson Brownlee.

GRChombo users have benefited from the provision of HPC resources from:

  * DiRAC (Distributed Research utilising Advanced Computing) resources under the projects ACSP218, ACSP191, ACTP183 and ACTP186. Systems used include: 

    - Cambridge Service for Data Driven Discovery (CSD3), part of which is operated by the University of Cambridge Research Computing on behalf of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The DiRAC component of CSD3 was funded by BEIS capital funding via STFC capital grants ST/P002307/1 and ST/R002452/1 and STFC operations grant ST/R00689X/1. DiRAC is part of the National e-Infrastructure.

    - DiRAC Data Intensive service at Leicester, operated by the University of Leicester IT Services, which forms part of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The equipment was funded by BEIS capital funding via STFC capital grants ST/K000373/1 and ST/R002363/1 and STFC DiRAC Operations grant ST/R001014/1. DiRAC is part of the National e-Infrastructure.

    -  DiRAC at Durham facility managed by the Institute for Computational Cosmology on behalf of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The equipment was funded by BEIS capital funding via STFC capital grants ST/P002293/1 and ST/R002371/1, Durham University and STFC operations grant ST/R000832/1. DiRAC is part of the National e-Infrastructure.

    - DIRAC Shared Memory Processing system at the University of Cambridge, operated by the COSMOS Project at the Department of Applied Mathematics and Theoretical Physics on behalf of the STFC DiRAC HPC Facility (www.dirac.ac.uk). This equipment was funded by BIS National E-infrastructure capital grant ST/J005673/1, STFC capital grant ST/H008586/1, and STFC DiRAC Operations grant ST/K00333X/1. DiRAC is part of the National e-Infrastructure.

    - DiRAC Complexity system, operated by the University of Leicester IT Services, which forms part of the STFC DiRAC HPC Facility (www.dirac.ac.uk ). This equipment is funded by BIS National E-Infrastructure capital grant ST/K000373/1 and STFC DiRAC Operations grant ST/K0003259/1. DiRAC is part of the National e-Infrastructure.

  * PRACE (Partnership for Advanced Computing in Europe) resources under grant numbers 2018194669, 2020225359. Systems used include:

    - SuperMUCNG, Leibniz Supercomputing Center (LRZ), Germany

    - JUWELS, Juelich Supercomputing Centre (JSC), Germany

    - Cartesius (SURF), Netherlands

    - Marenostrum (BSC), Spain
    
  * the Texas Advanced Computing Center (TACC) at The University of Texas at Austin HPC and visualisation resources URL: http://www.tacc.utexas.edu under project number PHY20043.

  * Marconi (CINECA), Italy

  * the Glamdring cluster, Astrophysics, Oxford, UK

  * The Fawcett cluster, Faculty of Mathematics, Cambridge, UK

  * the Argo cluster at ICTP, Trieste, Italy

  * the Apocrita cluster at QMUL, UK
  
  * The Athena cluster at HPC Midlands Plus, UK

  * Consortium des Équipements de Calcul Intensif (CÉCI), funded by the Fonds de la Recherche Scientifique de Belgique (F.R.S.-FNRS), Belgium

# References

