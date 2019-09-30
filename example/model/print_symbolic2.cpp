/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>
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

/* $Id: print_symbolic2.cpp 4026 2010-03-24 18:44:22Z troyer $ */

#include <alps/expression.h>

#include <alps/model/hamiltonian_matrix.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  using namespace boost::numeric::ublas ;

  alps::Parameters parms;
  std::cin >> parms;

  alps::hamiltonian_matrix<compressed_matrix<alps::Expression,row_major> > ham(parms);

  std::cout << ham.matrix() << "\n";

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
}

