from __future__ import print_function
#/*****************************************************************************
#*
#* ALPS Project: Algorithms and Libraries for Physics Simulations
#*
#* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
#*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
#*                            Maximilian Poprawe <poprawem@ethz.ch>
#*
#* This software is part of the ALPS libraries, published under the ALPS
#* Library License; you can use, redistribute it and/or modify it under
#* the terms of the license, either version 1 or (at your option) any later
#* version.
#*
#* You should have received a copy of the ALPS Library License along with
#* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
#* available from http://alps.comp-phys.org/.
#*
#* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
#* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
#* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
#* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#* DEALINGS IN THE SOFTWARE.
#*
#*****************************************************************************/

import pyalps

# This is an example of how to easily calculate the variance of data stored in a hdf5 file.


filename = "testfile.h5"

# create the correct MCData object to load the data.
obs = pyalps.alea.MCScalarData()

# load the variable E saved in the file testfile.h5 into the mcdata object.
obs.load(filename, "simulation/results/" + pyalps.hdf5_name_encode("E"))

# calculate the mean
variance = pyalps.alea.variance(obs)

# print the result
print("The variance of E is: " + str(variance))
  
# write the result back to the file
ar = pyalps.hdf5.archive(filename, 1)
ar.write("simulation/results/" + pyalps.hdf5_name_encode("E") + "/variance/value", variance)


