/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Maximilian Poprawe <poprawem@ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#include <alps/alea.h>
#include <alps/alea/mcanalyze.hpp>
#include <alps/utility/encode.hpp>

#include <alps/hdf5.hpp>

#include <iostream>
#include <string>


// This is an example of how to estimate the error of the mean of data stored in a hdf5 file

int main() {
  // this using-statement makes life easier when entering the named parameters
  using namespace alps::alea;
 
  const std::string filename = "testfile.h5";

  // create mcdata object with the correct template parameter.
  mcdata<double> obs;

  // load the variable m saved in the file testfile.h5 into the mcdata object.
  obs.load(filename, "simulation/results/" + alps::hdf5_name_encode("m"));

  // calculate the autocorrelation until it has reached 20% of its start value
  mctimeseries<double> auto_corr = alps::alea::autocorrelation(obs, _limit = 0.2);
  
  // fit the autocorrelation exponentially between the values where it is at 80% and at 20% of the value at t = 1
  std::pair<double, double> fit = exponential_autocorrelation_time(auto_corr, _max=0.8, _min=0.2);

  // calculate the integrated autocorrelation time by summing up the autocorrelation and integrating the fit for the tail
  double int_autocorr_time = integrated_autocorrelation_time(cut_tail(auto_corr, _limit=0.2), fit);


  // calculate different error estimates:
    // assuming uncorrelated data:
    double error_uncorr = error(obs, uncorrelated);

    // accounting for autocorrelation using binning analysis:
    double error_binning = error(obs, binning);

    // accounting for autocorrelation using the integrated autocorrelation time:
    double error_corrtime = error_uncorr * std::sqrt(1. + 2. * int_autocorr_time);


  // write to std::cout
  std::cout << "The estimated integrated autocorrelation time is: " << int_autocorr_time << "\n";
  std::cout << "The different error estimates are:\n";
  std::cout << "uncorrelated:\t\t" << error_uncorr << "\n";
  std::cout << "with binning:\t\t" << error_binning << "\n";
  std::cout << "with correlation time:\t" << error_corrtime << "\n";


  // we can also write one of the errors back to the file
  alps::hdf5::archive ar(filename, "a");
  ar << alps::make_pvp("simulation/results/" + alps::hdf5_name_encode("m") + "/mean/error", error_corrtime);

  return 0;
}

