/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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

/* $Id: random.h 1747 2005-08-28 08:40:20Z troyer $ */

#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <cmath>

#ifndef ALPS_APPLICATIONS_WORM_RANDOM_H
#define ALPS_APPLICATIONS_WORM_RANDOM_H

template <class RNG, class INT>
inline INT make_uniform_int(RNG& rng, INT i)
{ 
  return INT(i*rng());
  //return boost::variate_generator< RNG&, boost::uniform_int<INT> >(rng,boost::uniform_int<INT>(0,i-1))();
}

//- Worm creation/annihilation --------------------------------------------

inline double worm_creation_probability(double lambda,double time)
{
#ifdef SIMPLE
  if(lambda!=0.)
    return 0.;
  else
    return time*time/2.;
#else
  if(lambda==0. || fabs(lambda*time)<1e-10)
    return time*time/2.;
  else
    return (exp(-lambda*time)-1.+lambda*time)/(lambda*lambda);
#endif
}   // worm_creation_probability

inline double worm_creation_probability_open(double lambda,double time)
{
#ifdef SIMPLE
if(lambda!=0.)
  return 0.;
else
  return time*time;
#else
  if(lambda==0. || fabs(lambda*time)<1e-10)
    return time*time;
  else
    return time/lambda*(1.-exp(-lambda*time));
#endif
//  return time*integrated_weight(lambda,time);
}   // worm_creation_probability_open

//- Random time from finite time interval ---------------------------------

class new_finite_exponential_random {
public:
  new_finite_exponential_random(double x, double l, double t)
   : random_x(x), lambda(l), time(t) {}

  double operator()() {
    if(fabs(lambda*time)<1e-8)
      return time*random_x;
    else if(lambda<0)
    {   
      if(lambda*time > std::log(std::numeric_limits<double>::min()))
      { 
        double factor=std::exp(lambda*time);
        return 1./lambda*std::log(factor+(1.-factor)*random_x);
      }   
      else
      {
        double t= 1./lambda*std::log(random_x);
        return (t < time ? t : time-std::numeric_limits<double>::epsilon());
       }
    }
    else 
    {
      if(-lambda*time > std::log(std::numeric_limits<double>::min()))
      {
        double factor=std::exp(-lambda*time);
        return time - 1./(-lambda)*std::log(factor+(1.-factor)*random_x);
        }
      else
      {
        double t = time -1./(-lambda)*std::log(random_x);
        return (t>0 ? t : std::numeric_limits<double>::epsilon());
      } 
    }     
  }

private:
  double random_x;
  double lambda;
  double time; 
};   // new_finite_exponential_random

template <class RNG>
class finite_exponential_random {
public:
  finite_exponential_random(RNG& r, double l, double t)
   : rng(r,boost::uniform_real<>()), lambda(l), time(t) {}

  double operator()() {
#ifdef SIMPLE
    return time*rng();
#else
  return new_finite_exponential_random(rng(), lambda, time)();
#endif
  }

private:
 boost::variate_generator<RNG&, boost::uniform_real<> > rng;
 double lambda;
 double time;
};   // finite_exponential_random 

#endif
