/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2014        by Andreas Hehn <hehn@phys.ethz.ch>                   *
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
#ifndef ALPS_GRAPH_DETAIL_ASSERT_HELPERS_HPP
#define ALPS_GRAPH_DETAIL_ASSERT_HELPERS_HPP

#include <algorithm>
#include <vector>
#include <boost/algorithm/minmax.hpp>
#include <boost/algorithm/cxx11/is_sorted.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/graph/graph_traits.hpp>

namespace alps {
namespace graph {
namespace detail {
namespace assert_helpers {

struct tuple_0th_equals_1st
{
    template <typename T1, typename T2>
    bool operator()(boost::tuple<T1,T2> const& t) { return get<0>(t) == get<1>(t); }
};

template <typename Graph>
bool edge_list_matches_graph(std::vector<typename boost::graph_traits<Graph>::edge_descriptor> edge_list, Graph const& g)
{
    using std::sort;
    using std::equal;
    typename boost::graph_traits<Graph>::edge_iterator it,end;
    boost::tie(it,end) = edges(g);
    std::vector<typename boost::graph_traits<Graph>::edge_descriptor> gel(it,end);
    sort(edge_list.begin(),edge_list.end());
    sort(gel.begin(),gel.end());
    return equal(edge_list.begin(),edge_list.end(),gel.begin());
}

template <typename Graph>
bool color_partitions_are_complete(typename color_partition<Graph>::type const& color_partition, Graph const& g)
{
    // Check if all edge colors occuring in the graph are also in the color_partition
    typename boost::graph_traits<Graph>::edge_iterator it, end;
    bool are_complete = true;
    for (boost::tie(it, end) = edges(g); it != end; ++it)
        are_complete = are_complete && color_partition.find(get(alps::edge_type_t(),g)[*it]) != color_partition.end();
    return are_complete;
}

//
// Checks if graph g is a simple graph.
// A simple graph is a graph where
//   - no edge is connected to the same vertex at both ends and
//   - no two vertices are directly connected by more than one edge.
//
template <typename Graph>
bool is_simple_graph(Graph const& g)
{
    using std::sort;
    using std::adjacent_find;
    using std::find_if;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typename boost::graph_traits<Graph>::edge_iterator it, end;
    std::vector<boost::tuple<vertex_descriptor,vertex_descriptor> > edge_list;
    edge_list.reserve(num_edges(g));
    for (boost::tie(it, end) = edges(g); it != end; ++it)
        edge_list.push_back(boost::minmax(source(*it,g), target(*it,g)));
    sort(edge_list.begin(),edge_list.end());
    return (find_if(edge_list.begin(),edge_list.end(), tuple_0th_equals_1st()) == edge_list.end() ) // no loop edge
         && (adjacent_find(edge_list.begin(),edge_list.end()) == edge_list.end()); // no more than one edge between two vertices
}

/**
  * Checks the Graph g has a Vertex labeled v
  * \parm g The graph to be checked
  * \parm v The vertex descriptor
  */
template <typename Graph>
bool graph_has_vertex(Graph const& g, typename boost::graph_traits<Graph>::vertex_descriptor v)
{
    typename boost::graph_traits<Graph>::vertex_iterator vb,ve;
    boost::tie(vb,ve) = vertices(g);
    return (std::find(vb,ve,v) != ve);
}

template <typename Graph>
bool graph_has_vertices(Graph const& g, std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const& v)
{
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    bool ok = true;
    for(typename std::vector<vertex_descriptor>::const_iterator it = v.begin(); it != v.end(); ++it)
        ok = ok && graph_has_vertex(g, *it);
    return ok;
}

/**
  * Checks if a partition has a valid structure
  * \parm p The partition to be checked
  * \parm g the graph the partition belongs to
  */
template <typename Graph>
bool partition_has_valid_structure(typename alps::graph::partition_type<Graph>::type const& p, Graph const& g)
{
    using std::sort;
    using std::equal;
    bool ok = true;
    // Check if each part is sorted
    for(typename alps::graph::partition_type<Graph>::type::const_iterator it = p.begin(), end = p.end(); it != end; ++it)
        ok = ok && boost::algorithm::is_sorted(it->begin(), it->end());
    typename boost::graph_traits<Graph>::vertex_iterator vit,vend;
    boost::tie(vit,vend) = vertices(g);

    // Check if all vertices are present
    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> vg(vit,vend);
    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> vp;
    for(typename alps::graph::partition_type<Graph>::type::const_iterator it = p.begin(), end = p.end(); it != end; ++it)
        vp.insert(vp.end(), it->begin(), it->end());
    sort(vg.begin(), vg.end());
    sort(vp.begin(), vp.end());
    ok = ok && (vg.size() == vp.size());
    ok = ok && equal(vg.begin(), vg.end(), vp.begin());
    return ok;
}

/**
  * Checks if a partition is a discrete (i.e. has only trivial parts)
  * \parm p The partition to be checked
  * \parm g the graph the partition belongs to
  */
template <typename Graph>
bool is_discrete_partition(typename alps::graph::partition_type<Graph>::type const& p, Graph const& g)
{
    bool ok = partition_has_valid_structure(p,g);
    for(typename alps::graph::partition_type<Graph>::type::const_iterator it = p.begin(), end = p.end(); it != end; ++it)
        ok = ok && (it->size() == 1);
    return ok;
}

template <typename Graph>
bool orbit_of_is_valid(std::vector<std::size_t> const& orbit_of, typename alps::graph::partition_type<Graph>::type const& orbit, Graph const& g)
{
    typename boost::graph_traits<Graph>::vertex_iterator vit, vend;
    if( orbit_of.size() != num_vertices(g) )
        return false;
    for(boost::tie(vit,vend) = vertices(g); vit != vend; ++vit)
    {
        if(*vit >= orbit_of.size())
            return false;
        if( std::find(orbit[orbit_of[*vit]].begin(), orbit[orbit_of[*vit]].end(), *vit) == orbit[orbit_of[*vit]].end() )
            return false;
    }
    return true;
}

} // end namespace assert_helpers
} // end namespace detail
} // end namespace graph
} // end namespace alps


#endif // ALPS_GRAPH_DETAIL_ASSERT_HELPERS_HPP
