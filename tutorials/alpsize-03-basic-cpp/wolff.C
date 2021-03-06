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

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>

#define L 32
#define N (L*L)
#define T 2.2
#define MCSTEP (1 << 15)
#define MCTHRM (MCSTEP >> 3)
#define SEED 93812

double random_01() {
  return (double)(std::rand()) / RAND_MAX;
}

int main() {

  // setting up square lattice
  int nn[N][4];
  for (int y = 0; y < L; ++y)
    for (int x = 0; x < L; ++x) {
      nn[x+L*y][0] = ((x+L-1)%L) + L*y;
      nn[x+L*y][1] = ((x+1)%L) + L*y;
      nn[x+L*y][2] = x + L*((y+L-1)%L);
      nn[x+L*y][3] = x + L*((y+1)%L);
    }

  // random number generator
  srand(SEED);

  // spin configuration
  int spin[N];
  for (int s = 0; s < N; ++s) spin[s] = 1;
  int sz = N;

  // stack for uninspected sites
  int stck[N];
  int is = 0;

  // connecting probability
  double pc = 1 - std::exp(-2./T);

  // measurement
  double m = 0;
  double m2 = 0;
  double m4 = 0;

  // timer
  std::clock_t tm = std::clock();

  for (int mcs = 0; mcs < MCSTEP + MCTHRM; ++mcs) {
    int s = static_cast<int>(random_01() * N);
    int so = spin[s];
    spin[s] = -so;
    stck[0] = s;
    is = 1;
    int cs = 0;
    while (is) {
      ++cs;
      int sc = stck[--is];
      for (int k = 0; k < 4; ++k) {
        int sn = nn[sc][k];
        if (spin[sn] == so && random_01() < pc) {
          stck[is++] = sn;
          spin[sn] = -so;
        }
      }
    }
    sz -= 2 * so * cs;
    if (mcs >= MCTHRM) {
      double dsz = sz / static_cast<double>(N);
      m += dsz;
      m2 += dsz * dsz;
      m4 += dsz * dsz * dsz * dsz;
    }
  }

  // output results
  std::cout << "Magnetization = " << m / MCSTEP << std::endl;
  std::cout << "Magnetization^2 = " << m2 / MCSTEP << std::endl;
  std::cout << "Magnetization^4 = " << m4 / MCSTEP << std::endl;
  std::cout << "Binder Ratio of Magnetization = " << m2 * m2 / m4 / MCSTEP << std::endl;
  std::cerr << "Elapsed time = " << std::difftime(std::clock(), tm) / CLOCKS_PER_SEC << " sec\n";

  return 0;
}
