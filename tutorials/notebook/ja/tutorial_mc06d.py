# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
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
import matplotlib.pyplot as plt
import pyalps.plot
import copy

#prepare the input parameters
parms=[]

for (l,c) in [(4,500),( 6,1000)]:
    parms.append(
       { 
          'LATTICE'        : "simple cubic lattice", 
          'MODEL'          : "spin",
          'local_S'        : 0.5,
          'L'              : l,
          'J'              : 1 ,
          'CUTOFF'         : c
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('mc06d',parms)
pyalps.runApplication('qwl',input_file)

#run the evaluation and load all the plots
results = pyalps.evaluateQWL(pyalps.getResultFiles(prefix='mc06d'),DELTA_T=0.05, T_MIN=0.5, T_MAX=1.5)

#extract just the staggered structure factor S(Q) and rescale it by L^{-2+\eta}
data = []
for s in pyalps.flatten(results):
  if s.props['ylabel']=='Staggered Structure Factor per Site':
    print 'yes'
    d = copy.deepcopy(s) # make a deep copy to not change the original
    l = s.props['L']
    d.props['label']='L='+str(l)
    d.y = d.y * pow(float(l),-1.97)
    data.append(d)

#make plot
plt.figure()
plt.title("Scaling plot for cubic lattice Heisenberg antiferromagnet")
pyalps.plot.plot(data)
plt.legend()
plt.xlabel('Temperature $T/J$')
plt.ylabel('$S(\pi,\pi,\pi) L^{-2+\eta}$')
plt.show()
