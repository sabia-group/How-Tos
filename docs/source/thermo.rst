############
Fluctuations
############


During a molecular dynamics simulation one generally focuses on three main ensembles, namely 
the :math:`NVT`, :math:`NPT` and :math:`NVE` ensembles.
 
In the case of :math:`NVT` one introduces a thermostat to **equilibrate** the temperaturexi, whereas in :math:`NPT` a barostat is also
used. Various methods for introducing both a thermostat as well as a barostat have been introduced.  

It is important to point out, however, that the target temperature and/or pressure are only constant **on average**, and throughout a simulation
it will oscilate.  


.. _Pressure fluctuations:

Pressure fluctuations
#####################

The pressure fluctuaions can be calculated from

:math:`\sigma_P^2 = -k_B T \left(\frac{\partial P}{\partial V}\right)_S = \frac{k_B T}{\kappa_S V}`

where :math:`\kappa_S` is the isentropic compressibility, which can be obtained by

:math:`\kappa_S = \frac{1}{\rho c^2}` 

where :math:`\rho` is the density, and :math:`c` is the speed of sound in the material.

From this we obtain that

:math:`\sigma_P = \sqrt{\frac{\rho c^2 k_B T}{V}}`


Estimates
*********

Let's take a box of water with volume :math:`V = 20  ~  nm^3` (512 water molecules) at :math:`T = 300 K` and  density :math:`\rho = 1 ~ g/cm^2`. Under these 
conditions the speed of sound is :math:`1500 ~ m/s`

Thus we expect :math:`\sigma_P \sim 25 ~ MPa = 250 ~ bar`
`




