---
title: 'GRChombo: A numerical relativity code with AMR'
tags:
  - c++
  - gravity
  - numerical relativity
authors:
  - name: Katy Clough
    orcid: 0000-0000-0000-xxxx
    affiliation: 1
  - name: Pau Figueras
    orcid: 0000-0000-0000-xxxx
    affiliation: 2
  - name: Eugene A. Lim
    orcid: 0000-0000-0000-xxxx
    affiliation: 3

affiliations:
 - name: Oxford University
   index: 1
 - name: Queen Mary University of London
   index: 2
 - name: King's College London
   index: 3
date: 1 June 2020
bibliography: paper.bib

# Summary

Gravity is described by the Einstein Field Equation. This tells things how to move. It is non trivial to solve analytically and so requires numerical methods for time domain evolutions... blah blah.

GRChombo is an open-source code for numerical general relativity simulations. It is developed and maintained by a collaboration of numerical relativists with a wide range of research interests, from early universe cosmology to astrophysics and mathematical general relativity, and has been used in many papers since its first release in 2015 (add refs like this [@Pearson:2017]).

GRChombo is written entirely in C++14, using hybrid MPI/OpenMP parallelism and vector intrinsics to achieve good performance on the latest architectures. Furthermore, it makes use of the Chombo library for adaptive mesh refinement to allow automatic increasing and decreasing of the grid resolution in regions of arbitrary shape and topology.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Fenced code blocks are rendered with syntax highlighting:
```python
for n in range(10):
    yield f(n)
```	

# Acknowledgements

We acknowledge contributions from xx and yy. Dirac, PRACE, Royal Society, ERC...

# References