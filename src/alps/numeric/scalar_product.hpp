/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

/* $Id: scalar_product.hpp 7172 2013-09-24 15:00:50Z hehn $ */

#ifndef ALPS_NUMERIC_SCALAR_PRODUCT_HPP
#define ALPS_NUMERIC_SCALAR_PRODUCT_HPP

#include <alps/numeric/matrix/scalar_product.hpp>
#include <alps/functional.h>
#include <alps/type_traits/element_type.hpp>

#include <algorithm>
#include <functional>
#include <numeric>
#include <valarray>

namespace alps { namespace numeric {

// The generic implementation of the scalar_product moved to alps/numeric/matrix/scalar_product.hpp

/// \overload
template <class T>
inline T scalar_product(const std::valarray<T>& c1, const std::valarray<T>& c2) 
{
  return std::inner_product(data(c1),data(c1)+c1.size(),data(c2),T(), std::plus<T>(),conj_mult<T,T>());
}

} } // namespace alps::numeric

#endif // ALPS_NUMERIC_SCALAR_PRODUCT_HPP
