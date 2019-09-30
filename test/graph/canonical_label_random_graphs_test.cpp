/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "generate_random_graph.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <alps/lattice.h>
#include <alps/graph/canonical_properties.hpp>
#include <boost/random/mersenne_twister.hpp>

template <typename RNG, typename Graph>
bool canonical_label_test(RNG& rng, Graph const& g, unsigned int iterations = 100)
{
    using alps::graph::canonical_properties;
    typedef typename alps::graph::canonical_properties_type<Graph>::type canonical_properties_type;

    canonical_properties_type const gp(canonical_properties(g));
    bool ok = true;
    for(unsigned int i = 0; i < iterations; ++i)
    {
       Graph const ng = random_isomorphic_graph(rng,g);
       canonical_properties_type const ngp(canonical_properties(ng));
       ok = ok && (get<alps::graph::label>(gp) == get<alps::graph::label>(ngp));
       if(!ok)
           std::cout << "ERROR! graph label mismatch: "
               << get<alps::graph::label>(gp)
               << " != " << get<alps::graph::label>(ngp)
               << std::endl;
    }
    return ok;
}

int main()
{
    boost::random::mt19937 rng(23);

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;

    bool ok = true;
    for(unsigned int i=0; i < 50; ++i)
    {
        graph_type g = generate_random_simple_graph<graph_type>(rng, 20, 0.6);
        std::cout << "#v = " << num_vertices(g) << " #e = " << num_edges(g) << "\t";
        std::cout << get<alps::graph::label>(alps::graph::canonical_properties(g)) << std::endl;
        ok = ok && canonical_label_test(rng,g);
    }
    return ok ? 0 : -1;
}
