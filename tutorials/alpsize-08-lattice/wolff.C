/*****************************************************************************
*
* Tutorial: How to ALPSize your applications
*
* Copyright (C) 2005-2010 by Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person or organization
* obtaining a copy of the software and accompanying documentation covered by
* this license (the "Software") to use, reproduce, display, distribute,
* execute, and transmit the Software, and to prepare derivative works of the
* Software, and to permit third-parties to whom the Software is furnished to
* do so, all subject to the following:
*
* The copyright notices in the Software and this entire statement, including
* the above license grant, this restriction and the following disclaimer,
* must be included in all copies of the Software, in whole or in part, and
* all derivative works of the Software, unless such copies or derivative
* works are solely in the form of machine-executable object code generated by
* a source language processor.
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
#include <alps/parameter.h>
#include <alps/lattice.h>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <cmath>
#include <iostream>
#include <stack>
#include <vector>

int main() {

  alps::Parameters params(std::cin);
  const double T = params.value_or_default("T", 2.2);
  const int MCSTEP = params.value_or_default("SWEEPS", 1 << 15);
  const int MCTHRM = params.value_or_default("THERMALIZATION", MCSTEP >> 3);
  const unsigned int SEED = params.value_or_default("SEED", 93812);

  // setting up square lattice
  alps::graph_helper<> graph(params);
  const int N = graph.num_sites();

  // random number generator
  boost::mt19937 eng(SEED);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    random_01(eng, boost::uniform_real<>());

  // spin configuration
  std::vector<int> spin(N, 1);
  int sz = N;

  // stack for uninspected sites
  std::stack<int> stck;

  // connecting probability
  double pc = 1 - std::exp(-2./T);

  // measurement
  alps::ObservableSet measurements;
  measurements << alps::RealObservable("Magnetization");
  measurements << alps::RealObservable("Magnetization^2");
  measurements << alps::RealObservable("Magnetization^4");

  // timer
  boost::timer tm;

  for (int mcs = 0; mcs < MCSTEP + MCTHRM; ++mcs) {
    if (mcs == MCTHRM) measurements.reset(true);
    int s = static_cast<int>(random_01() * N);
    int so = spin[s];
    spin[s] = -so;
    stck.push(s);
    int cs = 0;
    while (!stck.empty()) {
      ++cs;
      int sc = stck.top();
      stck.pop();
      BOOST_FOREACH(alps::graph_helper<>::site_descriptor const& sn, graph.neighbors(sc)) {
        if (spin[sn] == so && random_01() < pc) {
          stck.push(sn);
          spin[sn] = -so;
        }
      }
    }
    sz -= 2 * so * cs;
    double dsz = sz / static_cast<double>(N);
    measurements["Magnetization"] << dsz;
    measurements["Magnetization^2"] << dsz * dsz;
    measurements["Magnetization^4"] << dsz * dsz * dsz * dsz;
  }

  // calculate Binder parameter
  alps::RealObsevaluator m2 = measurements["Magnetization^2"];
  alps::RealObsevaluator m4 = measurements["Magnetization^4"];
  alps::RealObsevaluator binder("Binder Ratio of Magnetization");
  binder = m2 * m2 / m4;
  measurements.addObservable(binder);

  // output results
  std::cout << measurements;
  std::cerr << "Elapsed time = " << tm.elapsed() << " sec\n";

  return 0;
}
