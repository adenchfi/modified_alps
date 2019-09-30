/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_TIME_H
#define LOOPER_TIME_H

#include <alps/config.h>
#include <alps/osiris.h>
#include <cmath>
#include <complex>

namespace looper {

inline std::complex<double> ctime(double t) {
  return std::exp(std::complex<double>(0, 2*M_PI*t));
}

inline std::complex<double> ctime2(double t) {
  return std::exp(std::complex<double>(0, 4*M_PI*t));
}

template<int STORE_CTIME = 1> class imaginary_time;

// STORE_CTIME = 0
// calculate complex time when needed
template<>
class imaginary_time<0> {
public:
  explicit imaginary_time(double t) : time_(t) {}
  operator double() const { return time_; }
  std::complex<double> ctime() const { return ::looper::ctime(time_); }
  std::complex<double> ctime2() const { return ::looper::ctime2(time_); }
private:
  double time_;
};

// STORE_CTIME = 1
// store (complex time) inside the object (default version)
template<>
class imaginary_time<1> {
public:
  explicit imaginary_time(double t) : time_(t), ctime_(::looper::ctime(t)) {}
  operator double() const { return time_; }
  std::complex<double> const& ctime() const { return ctime_; }
  std::complex<double> ctime2() const { return ctime_ * ctime_; }
private:
  double time_;
  std::complex<double> ctime_;
};

// STORE_CTIME = 2
// store both of (complex time) and (complex time)^2 inside the object
template<>
class imaginary_time<2> {
public:
  explicit imaginary_time(double t)
    : time_(t), ctime_(::looper::ctime(t)), ctime2_(ctime_ * ctime_) {}
  operator double() const { return time_; }
  std::complex<double> const& ctime() const { return ctime_; }
  std::complex<double> const& ctime2() const { return ctime2_; }
private:
  double time_;
  std::complex<double> ctime_, ctime2_;
};

inline std::complex<double> ctime(imaginary_time<0> const& t) {
  return t.ctime();
}

inline std::complex<double> ctime2(imaginary_time<0> const& t) {
  return t.ctime2();
}

inline std::complex<double> const& ctime(imaginary_time<1> const& t) {
  return t.ctime();
}

inline std::complex<double> ctime2(imaginary_time<1> const& t) {
  return t.ctime2();
}

inline std::complex<double> const& ctime(imaginary_time<2> const& t) {
  return t.ctime();
}

inline std::complex<double> const& ctime2(imaginary_time<2> const& t) {
  return t.ctime2();
}

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<int STORE_CTIME>
alps::ODump& operator<<(alps::ODump& dp, looper::imaginary_time<STORE_CTIME> const& t) {
  dp << static_cast<double>(t); return dp;
}

template<int STORE_CTIME>
alps::IDump& operator>>(alps::IDump& dp, looper::imaginary_time<STORE_CTIME>& t) {
  double t_in;
  dp >> t_in;
  t = looper::imaginary_time<STORE_CTIME>(t_in);
  return dp;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_TIME_H
