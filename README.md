DG-OSPREY2
=====
Discontinuous Galerkin Off-gas SeParation and REcoverY model <Version 2.0>

This version was created soley for the purpose of conforming to new MOOSE standards. No other major changes. 

This is a MOOSE application designed to simulate mass and energy transport of gases through a packed-bed column reactor. It uses DG kernels to ensure conservation of mass and energy for problems that are highly advectively dominant, but does not currently employ any form of slope limiters to correct the over-shoot under-shoot problems that are notorious in finite element methods. As a consequence, the solutions tend to be slightly oscillatory around the wave fronts as the various concentration profiles move through the bed. 

Material and energy balances are fully coupled, along with solid phase adsorption compositions. Full coupling ensures that the mass and energy balances are not violated by the simulations. Parameters for fluid flow are determined by the molecular diffusivities of the species in the gas phase, mechanical mixing introduced by the particle sizes, fluid velocities, and bed diameter, and velocity is determined by the flow rate, bed porosity, and bed dimensions. 

Thermal conductivities in the bed are estimated through semi-emperical relationships between thermal properties of all gas species as well as properties of the adsorbents. Other thermal properties, such as bed-wall heat transfer coefficients must be given as input. Separate temperatures for wall and/or ambient air can be specified in the input files. 

Currently, this model employs the following adsorption expressions for mass transfer:

(1) Generalized Statistical Thermodynamic Adsorption (GSTA) isotherm
(2) Coupled GSTA isotherm with Linear Driving Force (LDF) kinetics
(3) Langmuir adsorption isotherm
(4) Coupled Langmuir adsorption isotherm with LDF kinetics
(5) Extended Langmuir model for multiple gas species
(6) Coupled Extended Langmuir model with LDF kinetics

As this model continues to be developed, more and more adsorption expressions will be added. 

This is application is still under active development. Use with caution.

Project is authored by Dr. Austin Ladshaw (aladshaw3@gatech.edu). Please contact Dr. Ladshaw for any questions regarding the software.

To cite this application, you may use the following citation format:

Ladshaw, A.P., "DG-OSPREY: A mass and energy transport modeling platform for adsorption in fixed-bed columns using Discontinuous Galerkin finite element methods," <Version (#)>, https://github.com/aladshaw3/dgosprey, Access (Month) (Day), (Year). 

"Fork Dgosprey" to create a new MOOSE-based application.

For more information see: [http://mooseframework.org/create-an-app/](http://mooseframework.org/create-an-app/)
