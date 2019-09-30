/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
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

/* $Id: set_negative_0.hpp 3932 2010-02-28 15:24:40Z troyer $ */

#ifndef ALPS_NUEMRIC_SET_NEGATIVE_0_HPP
#define ALPS_NUEMRIC_SET_NEGATIVE_0_HPP

#include <alps/type_traits/is_sequence.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>


namespace alps { namespace numeric {

template <class T>
inline typename boost::disable_if<is_sequence<T>,void>::type
set_negative_0(T& x)
{
  if (x<T()) 
    x=T();
}

template <class T>
inline void set_negative_0(std::complex<T>& x)
{ 
  if (std::real(x)<0. || std::imag(x)<0.) 
    x=0.;
}

template <class T>
inline typename boost::enable_if<is_sequence<T>,void>::type
set_negative_0(T& a) 
{
  for(std::size_t i=0; i!=a.size(); ++i)
    set_negative_0(a[i]);
}



} } // end namespace alps::numeric

#endif // ALPS_NUEMRIC_SET_NEGATIVE_0_HPP
