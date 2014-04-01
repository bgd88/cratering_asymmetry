cratering_asymmetry
===================

Introduction
======= 
 
:Authors: Brent Delbridge, Ben Black, Ian Rose, Nick Swanson-Hysell 

:Web site: https://github.com/bgd88/cratering_asymetry

:Copyright: This document has been placed in the public domain.

:License:  Released under the Academic Free License.

Purpose
=======
The purpose of this Monte Carlo Model is to update the work of 
Zahnle et al., 2001 using modern observations and insight in order to 
explore differential cratering of synchronously rotating satelites.
Namely we hope to understand why the crater populations of the 
leading versus trailing hemisphere's are not correctly predicted
in several planetary satellites. This asymmetry is caused by the 
velocity of the satelite being much greater than the velocity of
the impactor, resulting in 15-50 times greater createring rate at
the apex (center of the leading hemisphere) than the anapex
(center of the trailing hempisphere). Our aim is understand what 
is happening when this is not the case. 

Code
====
To build our Monte Carlo Model we will utilize a Python Module
called PyMC that implements a Bayesian statistical models and
fitting algorithms. 
:Web site: https://github.com/pymc-devs/pymc
:Documentation: http://bit.ly/pymc_docs

*****
Model
*****


Random Parameters to be Estimated:
==================================

The following 6 Parameters will be estimated in our MC simulation.

First there are three parameters required to describe the orbit of 
an ecliptic comet making a close encounter with a giant planet.

	1) Encounter Velocity 

	2) Inclination

	3) Periapse Distance

There are two random numbers that will determine where on the hemi-
sphere the impact will occur.

	4) Azimuth of the impact site

	5) Incidence Angle of Impact

The last parameter is the mass of the impactor, which is necessary 
to predict an appropriate crater size.

	6) Impactor Mass

Calculated Quantities:
======================

From the first three quantites which determine comet's orbit we
can calculate:
	
	C1) The Impact Probability
	C2) The Impact Velocity
	C3) Apparent Radiant of Satelite

From which we may then calculate the crater diameter from the 
randomly sampled Mass and Incidence Angle and the calculated 
Impact Velocity
	C4) Crater Diameter


