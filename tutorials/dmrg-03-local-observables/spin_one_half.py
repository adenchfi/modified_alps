# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Jan Gukelberger <gukelberger@phys.ethz.ch> 
# 
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************

import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot

#prepare the input parameters
parms = []
for sz in [0,1,2]:
    parms.append( { 
        'LATTICE'                               : "open chain lattice", 
        'MODEL'                                 : "spin",
        'CONSERVED_QUANTUMNUMBERS'              : 'N,Sz',
        'Sz_total'                              : sz,
        'J'                                     : 1,
        'SWEEPS'                                : 4,
        'NUMBER_EIGENVALUES'                    : 1,
        'L'                                     : 32,
        'MAXSTATES'                             : 40,
        'MEASURE_LOCAL[Local magnetization]'    : 'Sz'
       } )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one_half',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one_half'))

# extract local magnetization data
curves = []
for run in data:
    for s in run:
        if s.props['observable'] == 'Local magnetization':
            sz = s.props['Sz_total']
            s.props['label'] = '$S_z = ' + str(sz) + '$'
            s.y = s.y.flatten()
            curves.append(s)
        
# Plot local magnetization vs. site
plt.figure()
pyalps.plot.plot(curves)
plt.legend()
plt.title('Magnetization of antiferromagnetic Heisenberg chain (S=1/2)')
plt.ylabel('local magnetization')
plt.xlabel('site')


plt.show()
