/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
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

#include <alps/graph/lattice_constant.hpp>

#include <boost/progress.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>

int main() {
    using boost::get;
    using alps::graph::canonical_properties;

    typedef unsigned int lc_type;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;

    alps::Parameters parm;
    unsigned int side_length = 40;

    std::ifstream in("../../lib/xml/lattices.xml");
    parm["LATTICE"] = "square lattice";
    parm["L"] = side_length;
    alps::graph_helper<> lattice(in,parm);

    std::vector<std::pair<graph_type,lc_type> > g;


    //
    //  0---1
    //  |   |
    //  2---3
    //
    g.push_back(std::make_pair(graph_type(), 1));
    add_edge(0, 1, g.back().first);
    add_edge(1, 3, g.back().first);
    add_edge(3, 2, g.back().first);
    add_edge(2, 0, g.back().first);

    //
    //  1---0---2
    //
    g.push_back(std::make_pair(graph_type(),6));
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);

    //
    //  3---1---0---2
    //
    g.push_back(std::make_pair(graph_type(),18));
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(1, 3,g.back().first);

    //
    //  1---0---2
    //      |
    //      3
    //
    g.push_back(std::make_pair(graph_type(),4));
    add_edge(0, 1, g.back().first);
    add_edge(0, 2, g.back().first);
    add_edge(0, 3, g.back().first);


    //
    //     3
    //     |
    // 1---0---2
    //     |
    //     4
    //
    g.push_back(std::make_pair(graph_type(),1));
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(0, 3,g.back().first);
    add_edge(0, 4,g.back().first);

    //
    //   2   5
    //   |   |
    //   0---1
    //   |   |
    //   3   4
    //
    g.push_back(std::make_pair(graph_type(),18));
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(0, 3,g.back().first);
    add_edge(1, 4,g.back().first);
    add_edge(1, 5,g.back().first);

    //
    //           8
    //           |
    //           4
    //           |
    //   6---2---0---1---5
    //           |
    //           3
    //           |
    //           7
    //
    g.push_back(std::make_pair(graph_type(),47));
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(0, 3,g.back().first);
    add_edge(0, 4,g.back().first);
    add_edge(1, 5,g.back().first);
    add_edge(2, 6,g.back().first);
    add_edge(3, 7,g.back().first);
    add_edge(4, 8,g.back().first);

    //
    //         3
    //         |
    //   6--2--0--1--4--7
    //            |
    //            5
    //
    g.push_back(std::make_pair(graph_type(),476));
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(0, 3,g.back().first);
    add_edge(1, 4,g.back().first);
    add_edge(1, 5,g.back().first);
    add_edge(2, 6,g.back().first);
    add_edge(4, 7,g.back().first);

    // graph no. 32340
    // 14    15    1    64272
    g.push_back(std::make_pair(graph_type(15),64272));
    add_edge(13, 14,g.back().first);
    add_edge(12, 14,g.back().first);
    add_edge( 9, 14,g.back().first);
    add_edge( 7, 14,g.back().first);
    add_edge( 8, 13,g.back().first);
    add_edge( 1, 13,g.back().first);
    add_edge( 6, 12,g.back().first);
    add_edge( 0, 12,g.back().first);
    add_edge(10,  9,g.back().first);
    add_edge( 2,  7,g.back().first);
    add_edge(11,  8,g.back().first);
    add_edge( 3,  6,g.back().first);
    add_edge( 4, 10,g.back().first);
    add_edge( 5, 11,g.back().first);

    //
    //  0---1---...---18---19
    //
//    g.push_back(std::make_pair(graph_type(10), 8134));
//    g.push_back(std::make_pair(graph_type(11), 22050));
//    g.push_back(std::make_pair(graph_type(12), 60146));
//    g.push_back(std::make_pair(graph_type(13), 162466));
//    g.push_back(std::make_pair(graph_type(14), 440750));
//    g.push_back(std::make_pair(graph_type(15), 1187222));
//    g.push_back(std::make_pair(graph_type(16), 3208298)); // 0.13 GB
    g.push_back(std::make_pair(graph_type(17), 8622666)); // 0.26 GB
//    g.push_back(std::make_pair(graph_type(18), 23233338)); // 0.64 GB 
//    g.push_back(std::make_pair(graph_type(19), 62329366)); // 2.13
//    g.push_back(std::make_pair(graph_type(20), 0)); // 
    add_edge( 0,  1, g.back().first);
    add_edge( 1,  2, g.back().first);
    add_edge( 2,  3, g.back().first);
    add_edge( 3,  4, g.back().first);
    add_edge( 4,  5, g.back().first);
    add_edge( 5,  6, g.back().first);
    add_edge( 6,  7, g.back().first);
    add_edge( 7,  8, g.back().first);
    add_edge( 8,  9, g.back().first);
    add_edge( 9, 10, g.back().first);
    add_edge(10, 11, g.back().first);
    add_edge(11, 12, g.back().first);
    add_edge(12, 13, g.back().first);
    add_edge(13, 14, g.back().first);
    add_edge(14, 15, g.back().first);
    add_edge(15, 16, g.back().first);
//    add_edge(16, 17, g.back().first);
//    add_edge(17, 18, g.back().first);
//    add_edge(18, 19, g.back().first);

    int success = 0;
    {
        boost::progress_timer timer;
        for(std::vector<std::pair<graph_type,lc_type> >::iterator it= g.begin(); it != g.end(); ++it)
        {
            lc_type lc = alps::graph::lattice_constant(
                  it->first
                , lattice.graph()
                , lattice.lattice()
                , alps::cell(std::vector<int>(2,side_length/2),lattice.lattice())
            );
            if (lc != it->second) {
                std::cerr<<"ERROR: lattice constant does not match!"<<std::endl;
                std::cerr<<"Calculated: "<<lc<<"\tReference: "<<it->second<<std::endl<<std::endl;
                success = -1;
            } else
                std::cerr<<"SUCCESS: "<<it->second<<std::endl;
        }
    }
    return success;
}
