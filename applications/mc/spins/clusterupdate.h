/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
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

/* $Id: clusterupdate.h 1956 2006-02-05 12:28:43Z troyer $ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_CLUSTERUPDATE_H_
#define ALPS_APPLICATIONS_MC_SPIN_CLUSTERUPDATE_H_

#include "connect.h"
#include "faststack.h"

#include <boost/random/uniform_int.hpp>
#include <vector>
#include <algorithm>

template <class Graph, class MomentMap, class RNG, class CouplingMap, class SpinFactorMap>
update_info_type cluster_update(const Graph& graph, MomentMap& moment,
                  double beta, RNG& rng, const CouplingMap& coupling, const SpinFactorMap& spinfactor_)
{
  // some typedefs for clarity
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<MomentMap>::value_type moment_type;
  typedef typename moment_type::update_type update_type;
  typedef typename CouplingMap::value_type MAT;

  // count cluster size
  update_info_type update_info;
  update_info.clustersize = 0;
  update_info.m2_ = 0;

  // data structures used for building the cluster
  fast_stack<vertex_descriptor> new_sites(boost::num_vertices(graph));
  std::vector<int> on_cluster(boost::num_vertices(graph),false);
  Connector<moment_type,RNG> connect(beta,rng);

  // step i)
  vertex_descriptor vertex=*(boost::vertices(graph).first + int(boost::num_vertices(graph)*rng()));
  new_sites.push(vertex);
  on_cluster[vertex]=true;
  
  // step ii)
  update_type update=moment[vertex].random_update(rng);

  // step iii)
  while (!new_sites.empty()) {
    // get current site
    vertex=new_sites.top();
    new_sites.pop();
    ++update_info.clustersize;

    // step iii.a)
    typename boost::graph_traits<Graph>::out_edge_iterator edge,end;
    for (boost::tie(edge,end)=boost::out_edges(vertex,graph);edge!=end; ++edge) {
      vertex_descriptor neighbor = boost::target(*edge,graph);
      double fact = -spinfactor_[boost::source(*edge,graph)]*spinfactor_[boost::target(*edge,graph)]; 
      MAT J = coupling[*edge];
      J = J * fact;
      if(!on_cluster[neighbor]&&connect(moment[vertex],moment[neighbor],update,J)) {
        new_sites.push(neighbor);
        on_cluster[neighbor]=true;
      } 
    }
    // step iii.b)
    update_info.m2_ += moment[vertex].project(update)*spinfactor_[vertex];
    moment[vertex].update(update);
  }
  return update_info;
}

#endif
