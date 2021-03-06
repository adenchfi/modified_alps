/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2007 - 2010 by Matthias Troyer <troyer@comp-phys.org>
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

#include "fleas.h"
#include <boost/random.hpp>
#include <valarray>

int main()
{
  const int N=50; // total number of fleas
  int M;          // number of samples
  std::cin >> M;
  
  typedef boost::mt19937 engine_type;
  typedef boost::uniform_int<> dist_type;
  
  engine_type engine(3242354);
  dist_type dist(0,1);
  boost::variate_generator<engine_type,dist_type> rng(engine,dist);
  
  std::valarray<unsigned long> histogram(N+1);
  
  for (int i=0;i<M;++i) {
    int n=0;
    for (int j=0;j<N;++j)
      if (rng()==0)
        ++n;
    histogram[n]++;
  }
  
  std::valarray<double> mean(N+1);
  std::valarray<double> error(N+1);
  for (int i=0; i<=N ; ++i) {
    mean[i] = static_cast<double>(histogram[i])/M;
    error[i] = std::sqrt((mean[i]-mean[i]*mean[i])/(M-1));
  }
  
  for (int i=0;i<=N;++i)
    std::cout << i << "\t" 
              << probability(N,i) << "\t"
              << mean[i] << "\t" 
              << error[i] << "\n";
  
  return 0;
}
