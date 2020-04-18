# RHCW18_Plate_Model

Finite difference code to compute thermal evolution of oceanic lithosphere and associated seafloor subsidence and heat flow. 
Scheme originally published in Richards et al. (JGR, 2018).

---------------------
Compilation and input
---------------------

To compile, run:

make

To clean directory of output files and executable before re-running code, run:

make clean

To execute code, run:

./cool [Tp] [zp] [zr]

where Tp is potential temperature in °C, zp is plate thickness in km and zr is zero-age ridge depth in m.

------------
Output files
------------

depth-[Tp]-[zp]-[zr].dat: Depth as a function of age. Column 1- age (Ma); column 2 - depth below sea-level (m).

hf-[Tp]-[zp].dat: Heat flow as a function of age. Column 1 - age (Ma); column 2 - heat flow (W m^-2).

Tgrid-[Tp]-[zp].dat: Temperature as a function of age and depth. Column 1 - age (Ma); column 2 - depth below top of plate 
					 (km); column 3 - depth below sea-level (km); column 4 - temperature (°C)
