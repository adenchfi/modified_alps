/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Axel Grzesik <axel@th.physik.uni-bonn.de>
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

/* $Id: example7.C 6149 2012-05-14 05:59:33Z wistaria $ */

#include <alps/model.h>
#include <fstream>
#include <iostream>
#include <string>

void write_two_site_basis(const std::string& name, const alps::ModelLibrary& lib, 
                          const alps::Parameters& p=alps::Parameters())
{
  alps::SiteBasisDescriptor<short> sitebasis=lib.get_site_basis(name),
                                   sitebasis2=lib.get_site_basis(name);
  sitebasis.set_parameters(p);
  std::cout << "States of basis " << name << "=" << alps::site_basis<short>(sitebasis);
  for(int i=0;i<sitebasis.size();++i)
    sitebasis[i]+=sitebasis2[i];
  std::cout << "States of basis " << name << "=" << alps::site_basis<short>(sitebasis);
}

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  
  std::ifstream in("../../lib/xml/models.xml");
  alps::ModelLibrary lib(in);
  alps::Parameters p;
  p["NMax"]=2;
  write_two_site_basis("spinful boson",lib,p);

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  //alps::comm_exit(true);
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
