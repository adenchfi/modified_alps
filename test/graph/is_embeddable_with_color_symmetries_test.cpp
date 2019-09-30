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

#include <alps/graph/is_embeddable.hpp>
#include <alps/graph/utils.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>


typedef boost::property<alps::edge_type_t,unsigned int> edge_props;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_props> graph_type;
typedef boost::graph_traits<graph_type>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<graph_type>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<graph_type>::edge_iterator edge_iterator;
typedef boost::property_map<graph_type,alps::edge_type_t>::type edge_color_map_type;
typedef alps::graph::graph_label<graph_type>::type label_type;
typedef alps::graph::canonical_properties_type<graph_type>::type canonical_properties_type;

void is_embeddable_with_color_symmetries_test()
{
    std::cout << "is_embeddable_with_color_symmetries_test()" << std::endl;
    using alps::graph::canonical_properties;
    using alps::graph::is_embeddable;

    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;

    // g: O---O+++O+++O---O
    graph_type g(5);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, g).first;
        edge_color[e] = 1;
        e = add_edge(2, 3, g).first;
        edge_color[e] = 1;
        e = add_edge(3, 4, g).first;
        edge_color[e] = 0;
    }
    canonical_properties_type gp = canonical_properties(g,color_symmetry);

    // h: O---O---O
    graph_type h(3);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),h);
        edge_descriptor e;
        e = add_edge(0, 1, h).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, h).first;
        edge_color[e] = 0;
    }
    canonical_properties_type hp = canonical_properties(h,color_symmetry);

    // i: O---O---O+++0
    graph_type i(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),i);
        edge_descriptor e;
        e = add_edge(0, 1, i).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, i).first;
        edge_color[e] = 0;
        e = add_edge(2, 3, i).first;
        edge_color[e] = 1;
    }
    canonical_properties_type ip = canonical_properties(i,color_symmetry);

    // j: O---O---O---O
    graph_type j(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),j);
        edge_descriptor e;
        e = add_edge(0, 1, j).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, j).first;
        edge_color[e] = 0;
        e = add_edge(2, 3, j).first;
        edge_color[e] = 0;
    }
    canonical_properties_type jp = canonical_properties(j,color_symmetry);

    std::cout << std::boolalpha << is_embeddable(h,g,get<alps::graph::partition>(hp),color_symmetry) << std::endl;
    std::cout << std::boolalpha << is_embeddable(i,g,get<alps::graph::partition>(ip),color_symmetry) << std::endl;
    std::cout << std::boolalpha << is_embeddable(j,g,get<alps::graph::partition>(jp),color_symmetry) << std::endl;
}

void is_embeddable_with_color_symmetries_test2()
{
    std::cout << "is_embeddable_with_color_symmetries_test2()" << std::endl;
    using alps::graph::canonical_properties;
    using alps::graph::is_embeddable;

    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;
    color_symmetry[2] = 0;

    // g:
    //     2    +++ c0
    //    / .   --- c1
    //   O+++1  ... c2
    //    . /
    //     3
    graph_type g(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, g).first;
        edge_color[e] = 1;
        e = add_edge(0, 3, g).first;
        edge_color[e] = 2;
        e = add_edge(1, 2, g).first;
        edge_color[e] = 2;
        e = add_edge(1, 3, g).first;
        edge_color[e] = 1;
    }
    canonical_properties_type gp = canonical_properties(g,color_symmetry);

    // h:
    //   0+++1   +++ c0
    //   |   |   --- c1
    //   2...3   ... c2
    graph_type h(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),h);
        edge_descriptor e;
        e = add_edge(0, 1, h).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, h).first;
        edge_color[e] = 1;
        e = add_edge(1, 3, h).first;
        edge_color[e] = 1;
        e = add_edge(2, 3, h).first;
        edge_color[e] = 2;
    }
    canonical_properties_type hp = canonical_properties(h,color_symmetry);
    std::cout << std::boolalpha << is_embeddable(h,g,get<alps::graph::partition>(hp),color_symmetry) << std::endl;
}

void is_embeddable_with_color_symmetries_test3()
{
    std::cout << "is_embeddable_with_color_symmetries_test3()" << std::endl;
    using alps::graph::canonical_properties;
    using alps::graph::is_embeddable;

    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;

    std::ifstream in("../../lib/xml/lattices.xml");
    alps::Parameters parm;
    parm["LATTICE"] = "anisotropic triangular lattice";
    parm["L"]       = 2*5+1;
    alps::graph_helper<> alps_lattice(in,parm);
    alps::coordinate_graph_type g = alps_lattice.graph();
    std::map<unsigned int, unsigned int> edge_map;
    edge_map[0] = 0;
    edge_map[1] = 0;
    edge_map[2] = 1;
    alps::graph::remap_edge_types(g, edge_map);

    // h:
    //     3
    //     |    +++ c0
    //     0    --- c1
    //    / +
    //   2+++1    
    graph_type h(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),h);
        edge_descriptor e;
        e = add_edge(0, 1, h).first;
        edge_color[e] = 1;
        e = add_edge(0, 2, h).first;
        edge_color[e] = 0;
        e = add_edge(0, 3, h).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, h).first;
        edge_color[e] = 1;
    }
    canonical_properties_type hp = canonical_properties(h,color_symmetry);
    std::vector<vertex_descriptor> pin(1, (2*5*2*5+1)/2);
    std::cout << std::boolalpha << is_embeddable(h, g, pin, get<alps::graph::partition>(hp)) << std::endl;
    std::cout << std::boolalpha << is_embeddable(h, g, pin, get<alps::graph::partition>(hp), color_symmetry) << std::endl;
}

void is_embeddable_with_color_symmetries_test4()
{
    std::cout << "is_embeddable_with_color_symmetries_test4()" << std::endl;
    using alps::graph::canonical_properties;
    using alps::graph::is_embeddable;

    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;

    // g:
    //     2        +++ c0
    //    / \       --- c1
    //   O+++1+++4 
    //    \ /
    //     3
    graph_type g(5);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, g).first;
        edge_color[e] = 1;
        e = add_edge(0, 3, g).first;
        edge_color[e] = 1;
        e = add_edge(1, 2, g).first;
        edge_color[e] = 1;
        e = add_edge(1, 3, g).first;
        edge_color[e] = 1;
        e = add_edge(1, 4, g).first;
        edge_color[e] = 0;
    }
    canonical_properties_type gp = canonical_properties(g,color_symmetry);

    // h:
    //     3
    //     |    +++ c0
    //     0    --- c1
    //    / +
    //   2+++1    
    graph_type h(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),h);
        edge_descriptor e;
        e = add_edge(0, 1, h).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, h).first;
        edge_color[e] = 1;
        e = add_edge(0, 3, h).first;
        edge_color[e] = 1;
        e = add_edge(1, 2, h).first;
        edge_color[e] = 0;
    }
    canonical_properties_type hp = canonical_properties(h,color_symmetry);
    std::vector<vertex_descriptor> pin(1, 1);
    std::cout << std::boolalpha << is_embeddable(h,g,pin,get<alps::graph::partition>(hp)) << std::endl;
    std::cout << std::boolalpha << is_embeddable(h,g,pin,get<alps::graph::partition>(hp),color_symmetry) << std::endl;
}

int main()
{
    is_embeddable_with_color_symmetries_test();
    is_embeddable_with_color_symmetries_test2();
    is_embeddable_with_color_symmetries_test3();
    is_embeddable_with_color_symmetries_test4();
    return 0;
}
