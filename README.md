# RHCW18_Plate_Model

Fortran 90 finite difference code to compute thermal evolution of oceanic lithosphere and associated seafloor subsidence and heat flow. 
Scheme originally published in Richards et al. (JGR, 2018) and updated in Richards et al. (PEPI, 2020).

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

depth-[Tp]-[zp]-[zr].dat: Depth as a function of age:   Column 1 - Age (Myr)
                                                        Column 2 - Depth below sea-level (m)

hf-[Tp]-[zp].dat: Heat flow as a function of age:       Column 1 - Age (Myr)
                                                        Column 2 - Surface heat flow (W m^-2)

Tgrid-[Tp]-[zp].dat: Temperature structure:             Column 1 - Age (Myr)
                                                        Column 2 - Depth below top of plate (km)
                                                        Column 3 - Depth below sea-level (km)
                                                        Column 4 - Temperature (°C)

------------
Citations:
------------

Richards, F.D., M.J. Hoggard, L.R. Cowton & N.J. White (2018). Reassessing the thermal structure of oceanic lithosphere with revised global inventories of basement depths and heat flow measurements, Journal of Geophysical Research: Solid Earth, Vol. 123, 9136-9161.

Richards, F.D., M.J. Hoggard, A. Crosby, S. Ghelichkhan & N.J. White (2020). Structure and dynamics of the oceanic lithosphere-asthenosphere system, Physics of the Earth and Planetary Interiors, under review.
