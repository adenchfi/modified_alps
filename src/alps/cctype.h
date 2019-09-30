/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

/* $Id: cctype.h 1744 2005-08-26 07:55:30Z wistaria $ */

/// \file cctype.h
/// \brief A safe version of the standard cctype header
///
///  Some cctype headers do not undefine harmful macros, so undefine
///  them here.

#ifndef ALPS_CCTYPE_H
#define ALPS_CCTYPE_H

#include <cctype>

#ifdef isspace 
# undef isspace
#endif
#ifdef isprint
# undef isprint
#endif
#ifdef iscntrl
# undef iscntrl
#endif
#ifdef isupper
# undef isupper
#endif
#ifdef islower
# undef islower
#endif
#ifdef isalpha
# undef isalpha
#endif
#ifdef isdigit
# undef isdigit
#endif
#ifdef ispunct
# undef ispunct
#endif
#ifdef isxdigit
# undef isxdigit
#endif
#ifdef isalnum
# undef isalnum
#endif
#ifdef isgraph
# undef isgraph
#endif
#ifdef toupper
# undef toupper
#endif
#ifdef tolower
# undef tolower
#endif

#endif // ALPS_CCTYPE_H
