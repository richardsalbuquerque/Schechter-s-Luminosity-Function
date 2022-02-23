# Schechter-s-Luminosity-Function

First we need to select the galaxies we want to work on.

Among the selection criteria we have:

• apparent magnitudes in the g filter less than 18

• redshifts between 0.001 and 0.02

The ADQL code used was:

select top 10000 p.ra, p.dec, p.u, p.g, p.r, p.i, p.z, s.z as redshift

from galaxy p, specobj s

where p.objid=s.bestobjid and p.g < 18 and s.z BETWEEN 0.001 AND 0.02

The objective of this project is to fit the Schechter function to the number density curve of galaxies.
