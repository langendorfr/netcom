# netcom
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/netcom)](https://cran.r-project.org/package=netcom)
[![Build Status](https://travis-ci.org/langendorfr/netcom.svg?branch=master)](https://travis-ci.org/langendorfr/netcom)

`netcom` is an R package to align two networks using the dynamics of diffusion originating from their respective nodes.

Consider network alignment as trying to compare two hypothetical cities of houses connected by roads. The approach implemented here is to pairwise compare each house with those in the other city by creating a house-specific signature. This is accomplished by quantifying the predictability of the location of a person at various times after they left their house, assuming they move randomly. This predictability across all houses captures much of the way each city is organized and functions. This package can be used to align networks using this conceptual rationale, with nodes as houses, edges as roads, and random diffusion representing people leaving their houses and walking around the city to other houses. 

The mechanics of this, which are conceptually akin to flow algorithms and Laplacian dynamics, can be analytically expressed as a Markov chain raised to successive powers which are the durations of diffusion. The alignment algorithm implemented here uses a normalized entropy to compare the predictability of diffusion emanating from each node in each network at each time step. A distance matrix is created with the numerically integrated differences between every pair of nodes' entropy-over-time curves where rows are nodes in one network and columns are nodes in the other network. The Hungarian algorithm is then used to find the optimal way to pair each node in each network with at most one node in the other network. 


## Installation
1. Install the release version of `netcom` from CRAN. 
	```R
	install.packages("netcom").
	```

2. Install the development version of `netcom` from GitHub using `devtools`.
	```R
	install.packages("devtools")
	devtools::install_github("langendorfr/netcom")
	```
