/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_TYPE_H
#define LOOPER_TYPE_H

#include <boost/mpl/bool.hpp>

namespace looper {

//
// QMC types
//

struct classical {};
struct path_integral {};
struct sse {};

//
// meta functions
//

template<typename QMC>
struct is_path_integral
{ typedef boost::mpl::false_ type; };

template<>
struct is_path_integral<path_integral>
{ typedef boost::mpl::true_ type; };

template<typename QMC>
struct is_sse
{ typedef boost::mpl::false_ type; };

template<>
struct is_sse<sse>
{ typedef boost::mpl::true_ type; };

} // end namepspace looper

#endif // LOOPER_TYPE_H
