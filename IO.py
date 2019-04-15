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
import os

def ReadFile(fileName):
  fileExist = os.path.isfile(fileName)
  if not fileExist:
    print('\nerror: file '+fileName+' does not exist.')   
    print('terminated.')
    exit()
   
  f = open(fileName, 'r')
  file1 = f.readlines()
  f.close
  
  file2 = []
  for lines in file1:
    file2.append(lines.strip())
  return file2


def Read_zqp_opt_dat(nSites,fileName):
  f_ = np.zeros((nSites,nSites), dtype=complex)

  print('\nreading file '+fileName)
  data = ReadFile(fileName)
  data_list = data[0].split()
  print(len(data_list))

  for ii in range(nSites):
    for jj in range(nSites):
      data_ii_re = (data[0].split())[len(data_list)-3-3*(ii+jj*nSites)]
      data_ii_im = (data[0].split())[len(data_list)-2-3*(ii+jj*nSites)]
      f_[ii][jj] += float(data_ii_re) + 1.j * float(data_ii_im)
      print(ii,jj,f_[ii][jj])
      
  return f_
      
  #exit()

def Symmetrize_f_ij(nSites, f_):
  ftmp_ = np.zeros((nSites,nSites), dtype=complex)
  print('\nSymmetrizing f_ij from ./zqp_opt.dat')
  for ii in range(nSites):
    for jj in range(nSites):
      diff = (jj-ii)%nSites
      ftmp_[0][diff] += (1./nSites)*f_[ii][jj]
  for ii in range(nSites):
    for jj in range(nSites):
      diff = (jj-ii)%nSites
      ftmp_[ii][jj] = ftmp_[0][diff]
  return f_    