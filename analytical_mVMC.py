# zlib license:

# Copyright (c) 2019 Maxime Charlebois

# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.

# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:

# 1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
# 2. Altered source versions must be plainly marked as such, and must not be
#   misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.

import numpy as np
import IO as io
from qState import *

def flavor(i,nSites):
  site = i%nSites
  spin = i//nSites
  return '%1d %1d '%(site,spin)   

def pretty(float1):
  if abs(float1)<1e-8:
    return '  . '#%(float1)   
  return '% 4.3e '%(float1)   

########################
##### main program #####
########################

nSites = 2

#f_[1][1] = 1.632588035
#f_[1][0] = 3.927111666
#f_[0][1] = 3.927111666
#f_[0][0] = 4.0

f_=io.Read_zqp_opt_dat(nSites,'zqp_opt.dat')
f_=io.Symmetrize_f_ij(nSites,f_)


#initialization of phi_pair:

phi_0 = qState(nSites,1.0)  #vacuum
phi_pair = qState(nSites,0.0) #zero (recipient)

for ii in range(0,nSites):
  for jj in range(0,nSites):
    phi_pair += f_[ii][jj] * (c_(nSites+ii) * (c_(jj) * phi_0))
print(phi_pair)

##############################
print('\nexample of use:', end='')
print('\n< c_0 a_1  c_1 a_0 > =', end='')

psi_2 = 1.0*(c_(0)*(a_(1)*(c_(1)*(a_(0)*(phi_pair)))))
val_2 = psi_2.scalarProd(phi_pair)/phi_pair.squareNorm()
    
print( (' % 4.3e '%val_2.real), end='\n')

exit()


