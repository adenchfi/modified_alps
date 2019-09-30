/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2015 by Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_GRAPH_DETAIL_LATTICE_CONSTANT_IMPL_HPP
#define ALPS_GRAPH_DETAIL_LATTICE_CONSTANT_IMPL_HPP

#include <alps/ngs/stacktrace.hpp>

#include <alps/lattice.h>
#include <alps/numeric/vector_functions.hpp>
#include <alps/graph/canonical_properties.hpp>
#include <alps/numeric/matrix.hpp>
#include <alps/graph/detail/shared_queue.hpp>
#include <alps/graph/detail/helper_functions.hpp>
#include <alps/graph/detail/assert_helpers.hpp>

#include <boost/array.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/integer_traits.hpp>
#include <boost/algorithm/minmax.hpp>

#include <vector>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace alps {
    namespace graph {

        namespace detail {

            inline unsigned int num_bits_required_for(std::size_t number)
            {
                std::size_t bits = 0;
                while ((0x01u << ++bits) < number);
                return bits;
            }

            template <typename Graph>
            unsigned int num_label_bits(Graph const& g)
            {
                return num_bits_required_for(num_vertices(g));
            }

            template <typename VertexDataType>
            class embeddings_set
            {
                // This is an (unordered) set of embedding entries which allows insertion only.
                // Since all entries are of the same size (known only at runtime), we save memory
                // by storing the data of all entries in vertices_data_ and edges_data_
                // and only put the index of the entry into the actual unordered_set.
              public:
                struct embeddings_set_entry_equal
                {
                    embeddings_set_entry_equal(embeddings_set const& parent)
                    : parent_(parent)
                    {
                    }
                    bool operator()(std::size_t i, std::size_t j) const
                    {
                        assert( (i+1) * parent_.edges_entry_size_    <= distance(parent_.edges_data_.begin(), parent_.edges_data_.end()) );
                        assert( (j+1) * parent_.edges_entry_size_    <= distance(parent_.edges_data_.begin(), parent_.edges_data_.end()) );
                        assert( (i+1) * parent_.vertices_entry_size_ <= distance(parent_.vertices_data_.begin(), parent_.vertices_data_.end()) );
                        assert( (j+1) * parent_.vertices_entry_size_ <= distance(parent_.vertices_data_.begin(), parent_.vertices_data_.end()) );
                        return std::equal(
                                  parent_.edges_data_.begin() +    i  * parent_.edges_entry_size_
                                , parent_.edges_data_.begin() + (i+1) * parent_.edges_entry_size_
                                , parent_.edges_data_.begin() +    j  * parent_.edges_entry_size_
                            ) && std::equal(
                                  parent_.vertices_data_.begin() +    i  * parent_.vertices_entry_size_
                                , parent_.vertices_data_.begin() + (i+1) * parent_.vertices_entry_size_
                                , parent_.vertices_data_.begin() +    j  * parent_.vertices_entry_size_
                            );
                    }
                  private:
                    embeddings_set const& parent_;
                };

                struct embeddings_set_entry_hash
                {
                    embeddings_set_entry_hash(embeddings_set const& parent)
                    : parent_(parent)
                    {
                    }
                    std::size_t operator()(std::size_t i) const
                    {
                        assert( (i+1) * parent_.edges_entry_size_     <= distance(parent_.edges_data_.begin(), parent_.edges_data_.end()) );
                        assert( (i+1) * parent_.vertices_entry_size_  <= distance(parent_.vertices_data_.begin(), parent_.vertices_data_.end()) );
                        using boost::hash_combine;
                        std::size_t hash = 0;
                        typename std::vector<VertexDataType>::const_iterator vit (parent_.vertices_data_.begin() + i*parent_.vertices_entry_size_);
                        typename std::vector<VertexDataType>::const_iterator vend(vit + parent_.vertices_entry_size_);
                        while( vit != vend)
                        {
                            hash_combine(hash, *vit);
                            ++vit;
                        }
                        std::vector<boost::uint64_t>::const_iterator eit (parent_.edges_data_.begin() + i*parent_.edges_entry_size_);
                        std::vector<boost::uint64_t>::const_iterator eend(eit + parent_.edges_entry_size_);
                        while( eit != eend)
                        {
                            hash_combine(hash, *eit);
                            ++eit;
                        }
                        return hash;
                    }
                  private:
                    embeddings_set const& parent_;
                };

                explicit embeddings_set(std::size_t num_vertices)
                    : vertices_entry_size_(num_vertices)
                    , edges_entry_size_(0)
                    , set_(0, embeddings_set_entry_hash(*this), embeddings_set_entry_equal(*this))
                {
                    std::size_t num_bits = num_vertices * (num_vertices + 1) / 2;
                    edges_entry_size_ = ((num_bits >> 6) + ((num_bits & 0x3F) == 0 ? 0 : 1));
                }

                template <typename SubGraph>
                bool insert(std::vector< boost::tuple<boost::uint16_t, VertexDataType, unsigned int> > const& embedding_data, SubGraph const& S)
                {
                    // Append data for new entry and remove it again if it is already known
                    std::size_t const vsize = vertices_data_.size();
                    std::size_t const esize = edges_data_.size();
                    vertices_data_.resize(vsize + vertices_entry_size_);
                    edges_data_.resize(   esize + edges_entry_size_);
                    create_entry(&vertices_data_[vsize], &edges_data_[esize], embedding_data, S);
                    std::size_t const num_entries = set_.size();
                    assert( vertices_data_.size() / vertices_entry_size_ == num_entries + 1);
                    assert( edges_data_.size() / edges_entry_size_       == num_entries + 1);
                    if(!set_.insert(num_entries).second)
                    {
                        vertices_data_.resize(vsize);
                        edges_data_.resize(esize);
                        return false;
                    }
                    else
                        return true;
                }

                std::size_t size() const
                {
                    return set_.size();
                }

              private:
                template <typename SubGraph>
                void create_entry(VertexDataType * const vertices_entry, boost::uint64_t * const edges_entry, std::vector< boost::tuple<boost::uint16_t, VertexDataType, unsigned int> > const& embedding_data, SubGraph const& S)
                {
                    assert(num_vertices(S) == vertices_entry_size_);
                    // std::vector<std::size_t> index_of_subgraph_vertex_in_embedding_generic -> use buffer instead of reallocating each time
                    index_of_subgraph_vertex_in_embedding_generic_buffer_.clear();
                    create_vertex_part(vertices_entry, index_of_subgraph_vertex_in_embedding_generic_buffer_, embedding_data);
                    create_edge_part(edges_entry, index_of_subgraph_vertex_in_embedding_generic_buffer_, S);
                }

                static void create_vertex_part(VertexDataType * const vertices_entry, std::vector<std::size_t> & index_of_subgraph_vertex, std::vector< boost::tuple<boost::uint16_t, VertexDataType, unsigned int> > const& embedding_data)
                {
                    BOOST_STATIC_ASSERT(( boost::is_same<VertexDataType,boost::uint16_t>::value || boost::is_same<VertexDataType,boost::uint32_t>::value || boost::is_same<VertexDataType, boost::uint64_t>::value));
                    index_of_subgraph_vertex.clear();
                    index_of_subgraph_vertex.resize(embedding_data.size());
                    enum { orbit_number = 0 , compressed_embedding_data = 1 , subgraph_vertex = 2 };
                    for( typename std::vector< boost::tuple<boost::uint16_t, VertexDataType, unsigned int> >::const_iterator it = embedding_data.begin() ; it != embedding_data.end(); ++it)
                    {
                        std::size_t const pos = it - embedding_data.begin();
                        assert( get<subgraph_vertex>(*it) < index_of_subgraph_vertex.size() );
                        index_of_subgraph_vertex[ get<subgraph_vertex>(*it) ] = pos;
                        *(vertices_entry+pos)                                 = get<compressed_embedding_data>(*it);
                    }
                }

                template <typename SubGraph>
                static void create_edge_part(boost::uint64_t * const edges_entry, std::vector<std::size_t> const& index_of_subgraph_vertex, SubGraph const& S)
                {
                    typename boost::graph_traits<SubGraph>::edge_iterator s_ei, s_ee;
                    for( boost::tie(s_ei, s_ee) = edges(S); s_ei != s_ee; ++s_ei) {
                        std::size_t v1, v2;
                        boost::tie(v1,v2) = boost::minmax(index_of_subgraph_vertex[source(*s_ei, S)], index_of_subgraph_vertex[target(*s_ei, S)]);
                        std::size_t index = v1 * num_vertices(S) - (v1 - 1) * v1 / 2 + v2 - v1;
                        *(edges_entry + (index >> 6)) |= 0x01 << (index & 0x3F);
                    }
                }

                std::size_t vertices_entry_size_;
                std::size_t edges_entry_size_;
                std::vector<VertexDataType>  vertices_data_;
                std::vector<boost::uint64_t> edges_data_;
                boost::unordered_set<std::size_t, embeddings_set_entry_hash, embeddings_set_entry_equal> set_;
                std::vector<std::size_t> index_of_subgraph_vertex_in_embedding_generic_buffer_;
            };


            // Returns a table with the distances of the lattice graph vertices to the (periodic) boundary of the lattice in units of unit cell shifts:
            // distance of vertex v along dimension d: distance[d][vertexid]
            template <typename Graph, typename Lattice> std::vector<std::vector<boost::uint_t<8>::fast> > build_distance_table(
                  Graph const & graph
                , Lattice const & lattice
            ) {
                typedef typename alps::lattice_traits<Lattice>::cell_iterator cell_iterator;
                typedef typename alps::lattice_traits<Lattice>::offset_type offset_type;
                typedef typename alps::lattice_traits<Lattice>::size_type cell_index_type;

                unsigned const invalid = num_vertices(graph);
                std::vector<std::vector<boost::uint_t<8>::fast> > distance_to_boarder(dimension(lattice), std::vector<boost::uint_t<8>::fast>(num_vertices(graph),0));
                unsigned vtcs_per_ucell = num_vertices(alps::graph::graph(unit_cell(lattice)));
                for(std::size_t d = 0; d < dimension(lattice); ++d) {
                    std::vector<unsigned> translated_vertex_along_d(num_vertices(graph), invalid);
                    for(std::pair<cell_iterator,cell_iterator> c = cells(lattice); c.first != c.second; ++c.first) {
                        const cell_index_type cellidx = index(*c.first,lattice);
                        offset_type ofst = offset(*c.first,lattice);
                        offset_type move(dimension(lattice));
                        move[d] = -1;
                        std::pair<bool,bool> on_lattice_pbc_crossing = shift(ofst,move,lattice);
                        // first = shifted ofst is on lattice (i.e. ofst is valid)
                        // second = (converted to bool) true if boundary was crossed, false otherwise
                        if(on_lattice_pbc_crossing.first && !on_lattice_pbc_crossing.second) {
                            // if no boundary was crossed,
                            // save the mapping 'vertex' -> 'translated vertex' in translated_vertex_along_d
                            const cell_index_type neighboridx = index(cell(ofst, lattice), lattice);
                            for(unsigned v = 0; v < vtcs_per_ucell; ++v)
                                translated_vertex_along_d[cellidx * vtcs_per_ucell + v] = neighboridx * vtcs_per_ucell + v;
                        }
                    }
                    for (std::vector<unsigned>::const_iterator it = translated_vertex_along_d.begin(); it != translated_vertex_along_d.end(); ++it) {
                        if (*it != invalid)
                        {
                            unsigned int v = *it;
                            while ((v = translated_vertex_along_d[v]) != invalid)
                            {
                                assert( distance_to_boarder[d][*it] < std::numeric_limits<boost::uint_t<8>::fast>::max() && "Overflow in distance_to_boarder!");
                                ++distance_to_boarder[d][*it];
                            }
                        }
                    }
                }
                return distance_to_boarder;
            }


            template <typename SubGraph>
            class count_translational_invariant_embeddings
            {
                typedef SubGraph subgraph_type;
              public:

                template <typename Graph, typename Lattice>
                count_translational_invariant_embeddings(subgraph_type const& S, Graph const& G, Lattice const& L, typename partition_type<subgraph_type>::type const & subgraph_orbit)
                : orbit_of_( build_vertex_to_partition_index_map(subgraph_orbit,S) )
                , distance_to_boarder_( build_distance_table(G,L) )
                , matches_( num_vertices(S))
                , unit_cell_size_( num_vertices(alps::graph::graph(unit_cell(L))) )
                {
                    assert(( get<alps::graph::partition>(canonical_properties(S)) == subgraph_orbit ));
                    assert(( assert_helpers::partition_has_valid_structure(subgraph_orbit, S) ));
                    if( !( (0x01ull << (distance_to_boarder_.size() * num_label_bits(S) + num_bits_required_for(unit_cell_size_))) < static_cast<unsigned long long>(boost::integer_traits<boost::uint16_t>::const_max) ) )
                        throw std::runtime_error("Subgraph, lattice dimension or unit cell is too large. uint16_t is not large enough to store full embedding data.");
                }

                template <typename Graph>
                void operator()(
                      subgraph_type const& S
                    , Graph const& G
                    , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const & pinning
                ) {
                    assert(( assert_helpers::orbit_of_is_valid(orbit_of_, get<alps::graph::partition>(canonical_properties(S)), S) ));

                    enum { orbit_number = 0 , compressed_embedding_data = 1 , subgraph_vertex = 2 };
                    embedding_data_buffer_.clear();

                    // How many bits do we need to store the maximal distance in a subgraph
                    // This is the same number of bits as we need for a unique vertex label of the subgraph
                    unsigned int const bits_per_dim = num_label_bits(S);

                    // Consider translational invariance and get a unique embedding representative
                    // embedding_compressed contains the full embedding information for a subgraph vertex (i.e. index within unit cell, "canonical" position in all dimensions)
                    for (typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::const_iterator it = pinning.begin(); it != pinning.end(); ++it) {
                        unsigned int    const subgraph_vertex      = it - pinning.begin();
                        boost::uint16_t const embedding_compressed = calc_compressed_transinv_embedding(bits_per_dim, *it, pinning, G);
                        embedding_data_buffer_.push_back( boost::make_tuple( static_cast<boost::uint16_t>(orbit_of_[subgraph_vertex]), embedding_compressed, subgraph_vertex) );
                    }

                    sort(embedding_data_buffer_.begin(), embedding_data_buffer_.end());

                    matches_.insert(embedding_data_buffer_,S);
                }

                std::size_t get_count() const
                {
                    return matches_.size();
                }
              private:

                /// Gets the (minimal) distance to boarder in each dimension of an embedding/pinning
                template <typename Graph>
                void get_distance_to_boarder(std::vector<boost::uint_t<8>::fast> & distance, std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const& pinning, Graph const& G) const
                {
                    distance.resize(distance_to_boarder_.size(), std::numeric_limits<boost::uint_t<8>::fast>::max());
                    for(std::size_t d = 0; d < distance_to_boarder_.size(); ++d)
                    {
                        boost::uint_t<8>::fast dist = std::numeric_limits<boost::uint_t<8>::fast>::max();
                        for (typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::const_iterator it = pinning.begin(); it != pinning.end(); ++it)
                            dist = (std::min)(dist, distance_to_boarder_[d][*it]);
                        distance[d] = dist;
                    }
                }

                template <typename Graph>
                boost::uint16_t calc_compressed_transinv_embedding(
                      unsigned int const bits_per_dim
                    , typename boost::graph_traits<Graph>::vertex_descriptor v_id
                    , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const& pinning
                    , Graph const& G
                ) const {
                    assert((0x01u << (distance_to_boarder_.size() * bits_per_dim + num_bits_required_for(unit_cell_size_))) < boost::integer_traits<boost::uint16_t>::const_max && "boost::uint16_t is not large enough to store full embedding data.");
                    get_distance_to_boarder(distance_buffer_, pinning, G);
                    // data =  bitfield looking like: unit_cell_vtx_idx|d[0]|d[1]|...
                    boost::uint16_t data = v_id % unit_cell_size_;
                    for(std::size_t d = 0; d < distance_to_boarder_.size(); ++d)
                    {
                        data <<= bits_per_dim;
                        assert( distance_to_boarder_[d][v_id] >= distance_buffer_[d] );
                        boost::uint_t<8>::fast const dist = distance_to_boarder_[d][v_id] - distance_buffer_[d];
                        assert( dist < (0x01u << bits_per_dim) );
                        data += dist;
                    }
                    return data;
                }

                std::vector<std::size_t>                                        orbit_of_;
                std::vector<std::vector<boost::uint_t<8>::fast> > const         distance_to_boarder_;
                std::vector< boost::tuple<boost::uint16_t, boost::uint16_t, unsigned int> > embedding_data_buffer_;
                embeddings_set<boost::uint16_t>                                 matches_;
                mutable std::vector<boost::uint_t<8>::fast>                     distance_buffer_;
                std::size_t const                                               unit_cell_size_;
            };

            template <typename SubGraph>
            class count_how_subgraphs_are_embedded
            {
                typedef SubGraph subgraph_type;
              public:

                template <typename Graph, typename Lattice>
                count_how_subgraphs_are_embedded(
                      subgraph_type const& S
                    , Graph const& G
                    , Lattice const& L
                    , typename partition_type<subgraph_type>::type const & subgraph_orbit
                    , alps::numeric::matrix<unsigned int> & orbit_mapping_matrix
                    , typename boost::graph_traits<subgraph_type>::vertex_descriptor breaking_vertex
                )
                : orbit_of_( build_vertex_to_partition_index_map(subgraph_orbit,S) )
                , distance_to_boarder_( build_distance_table(G,L) )
                , matches_(num_vertices(S))
                , unit_cell_size_( num_vertices(alps::graph::graph(unit_cell(L))) )
                , orbit_mapping_matrix_(orbit_mapping_matrix)
                , breaking_vertex_(breaking_vertex)
                {
                    assert(( get<alps::graph::partition>(canonical_properties(S,breaking_vertex_)) == subgraph_orbit ));
                    assert(( assert_helpers::partition_has_valid_structure(subgraph_orbit, S) ));
                    if( !( (0x01ull << (distance_to_boarder_.size() * num_label_bits(G) + num_bits_required_for(unit_cell_size_))) < static_cast<unsigned long long>(boost::integer_traits<boost::uint64_t>::const_max) ) )
                        throw std::runtime_error("Lattice graph, lattice dimension or unit cell is too large. uint64_t is not large enough to store full embedding data.");
                }

                template <typename Graph>
                void operator()(
                      subgraph_type const& S
                    , Graph const& G
                    , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const & pinning
                ) {
                    assert(( assert_helpers::orbit_of_is_valid(orbit_of_, get<alps::graph::partition>(canonical_properties(S,breaking_vertex_)), S) ));

                    enum { orbit_number = 0 , compressed_embedding_data = 1 , subgraph_vertex = 2 };
                    embedding_data_buffer_.clear();

                    // How many bits do we need to store the maximal distance to the boundary of the lattice (in hops)
                    // This is the same number of bits as we need for a unique vertex label of the graph
                    unsigned int const bits_per_dim = num_label_bits(G);

                    // Build lattice_pinning data.
                    // embedding_compressed contains the full embedding information for a subgraph vertex (i.e. index within unit cell, position in all dimensions)
                    for (typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::const_iterator it = pinning.begin(); it != pinning.end(); ++it) {
                        unsigned int    const subgraph_vertex      = it - pinning.begin();
                        boost::uint64_t const embedding_compressed = calc_compressed_embedding(bits_per_dim, *it, pinning, G);
                        embedding_data_buffer_.push_back( boost::make_tuple( static_cast<boost::uint16_t>(orbit_of_[subgraph_vertex]), embedding_compressed, subgraph_vertex) );
                    }

                    sort(embedding_data_buffer_.begin(), embedding_data_buffer_.end());

                    bool const inserted = matches_.insert(embedding_data_buffer_, S);
                    if (inserted)
                        add_pinning_to_orbit_mapping_matrix(pinning, G);
                }

              private:

                template <typename Graph>
                void add_pinning_to_orbit_mapping_matrix(std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const& pinning, Graph const& ) const
                {
                    for (typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::const_iterator it = pinning.begin(); it != pinning.end(); ++it)
                        ++orbit_mapping_matrix_(*it, orbit_of_[it - pinning.begin()]);
                }

                template <typename Graph>
                boost::uint64_t calc_compressed_embedding(
                      unsigned int const bits_per_dim
                    , typename boost::graph_traits<Graph>::vertex_descriptor v_id
                    , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const& pinning
                    , Graph const& G
                ) const {
                    assert((0x01u << (distance_to_boarder_.size() * bits_per_dim + num_bits_required_for(unit_cell_size_))) < boost::integer_traits<boost::uint64_t>::const_max && "boost::uint64_t is not large enough to store full embedding data.");
                    // data =  bitfield looking like: unit_cell_vtx_idx|d[0]|d[1]|...
                    boost::uint64_t data = v_id % unit_cell_size_;
                    for(std::size_t d = 0; d < distance_to_boarder_.size(); ++d)
                    {
                        data <<= bits_per_dim;
                        boost::uint_t<8>::fast const dist = distance_to_boarder_[d][v_id];
                        assert( dist < (0x01u << bits_per_dim) );
                        data += dist;
                    }
                    return data;
                }


                std::vector<std::size_t>                                        orbit_of_;
                std::vector<std::vector<boost::uint_t<8>::fast> > const         distance_to_boarder_;
                std::vector< boost::tuple<boost::uint16_t, boost::uint64_t, unsigned int> > embedding_data_buffer_; // TODO: this calls for bad alignment...
                embeddings_set<boost::uint64_t>                                 matches_;
                std::size_t const                                               unit_cell_size_;
                alps::numeric::matrix<unsigned int> &                           orbit_mapping_matrix_;
                typename boost::graph_traits<subgraph_type>::vertex_descriptor const  breaking_vertex_;
            };

            struct embedding_found {};

            struct throw_on_embedding_found
            {
                template <typename SubGraph, typename Graph>
                void operator()(SubGraph const&, Graph const&, std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const&)
                {
                    throw embedding_found();
                }
            };

            template <typename Subgraph>
            struct edge_equal_simple
            {
              private:
                template <typename Graph>
                bool impl(
                  typename boost::graph_traits<Subgraph>::edge_descriptor const & s_e
                , typename boost::graph_traits<Graph>::edge_descriptor const & g_e
                , Subgraph const & S
                , Graph const & G
                , boost::mpl::true_
                ) const {
                    return get(alps::edge_type_t(), S)[s_e] == get(alps::edge_type_t(), G)[g_e];
                }

                template <typename Graph>
                bool impl(
                  typename boost::graph_traits<Subgraph>::edge_descriptor const & s_e
                , typename boost::graph_traits<Graph>::edge_descriptor const & g_e
                , Subgraph const & S
                , Graph const & G
                , boost::mpl::false_
                ) const {
                    return true;
                }

              public:
                template <typename Graph>
                bool operator()(
                  typename boost::graph_traits<Subgraph>::edge_descriptor const & s_e
                , typename boost::graph_traits<Graph>::edge_descriptor const & g_e
                , Subgraph const & S
                , Graph const & G
                ) const {
                    return impl(s_e,g_e,S,G,boost::mpl::bool_<has_property<alps::edge_type_t,Subgraph>::edge_property>());
                }

            };

            template <typename Subgraph>
            struct edge_equal_mapped_colors
            {
              public:
                edge_equal_mapped_colors(std::vector<alps::type_type> const& map)
                : map_(map)
                {
                }

                template <typename Graph>
                bool operator()(
                  typename boost::graph_traits<Subgraph>::edge_descriptor const & s_e
                , typename boost::graph_traits<Graph>::edge_descriptor const & g_e
                , Subgraph const & S
                , Graph const & G
                ) const {
                    assert( get(alps::edge_type_t(), S, s_e) < map_.size() );
                    return map_[get(alps::edge_type_t(), S, s_e)] == get(alps::edge_type_t(), G, g_e);
                }
              private:
                std::vector<alps::type_type> const& map_;
            };

            template <typename Subgraph>
            struct vertex_equal_simple
            {
              private:
                template <typename Graph>
                bool impl(
                  typename boost::graph_traits<Subgraph>::vertex_descriptor const & s_v
                , typename boost::graph_traits<Graph>::vertex_descriptor const & g_v
                , Subgraph const & S
                , Graph const & G
                , boost::mpl::true_
                ) const {
                    return get(alps::vertex_type_t(), S)[s_v] == get(alps::vertex_type_t(), G)[g_v];
                }

                template <typename Graph>
                bool impl(
                  typename boost::graph_traits<Subgraph>::vertex_descriptor const & s_v
                , typename boost::graph_traits<Graph>::vertex_descriptor const & g_v
                , Subgraph const & S
                , Graph const & G
                , boost::mpl::false_
                ) const {
                    return true;
                }

              public:
                template <typename Graph>
                bool operator()(
                  typename boost::graph_traits<Subgraph>::vertex_descriptor const & s_v
                , typename boost::graph_traits<Graph>::vertex_descriptor const & g_v
                , Subgraph const & S
                , Graph const & G
                ) const {
                    return impl(s_v,g_v,S,G,boost::mpl::bool_<has_property<alps::vertex_type_t,Subgraph>::vertex_property>());
                }
            };

            // TODO: make an object out of walker
            template<typename Subgraph, typename Graph, typename VertexEqual, typename EdgeEqual, typename EmbeddingFoundPolicy> void lattice_constant_walker(
                  typename boost::graph_traits<Subgraph>::vertex_descriptor const & s
                , typename boost::graph_traits<Graph>::vertex_descriptor const & g
                , Subgraph const & S
                , Graph const & G
                , shared_queue_view<
                      typename boost::graph_traits<Subgraph>::vertex_descriptor
                    , typename boost::graph_traits<Graph>::vertex_descriptor
                  > const& queue
                , boost::dynamic_bitset<> & visited
                , unsigned int visited_cnt
                , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> & pinning
                , VertexEqual & vertex_equal
                , EdgeEqual & edge_equal
                , EmbeddingFoundPolicy& register_embedding
            ) {
                typedef typename boost::graph_traits<Subgraph>::vertex_descriptor subgraph_vertex_descriptor;

                // Check if the vertex mapping s->g is valid by checking if...
                // ... the degrees are compatible
                if (out_degree(s, S) > out_degree(g, G))
                    return;
                // ... the vertex types are equal (with respect to symmetries maybe)
                if (!vertex_equal(s, g, S, G))
                    return;
                // ... the existing edges from s are compatible with those of g.
                typename boost::graph_traits<Subgraph>::adjacency_iterator s_ai, s_ae;
                for (boost::tie(s_ai, s_ae) = adjacent_vertices(s, S); s_ai != s_ae; ++s_ai)
                    if (pinning[*s_ai] != num_vertices(G)) {
                        typename boost::graph_traits<Graph>::edge_descriptor e;
                        bool is_e;
                        boost::tie(e, is_e) = edge(g, pinning[*s_ai], G);
                        if (!is_e || !edge_equal( edge(s, *s_ai, S).first , e, S, G) )
                            return;
                    }

                // s->g seems legit => pin s->g.
                pinning[s] = g;
                visited.set(g);
                ++visited_cnt;
                assert(visited_cnt == visited.count()); // visited.count is called frequently and is rather slow -> cached
                // If not all vertices are mapped yet
                if (visited_cnt < num_vertices(S)) {
                    // queue mapping adjecent vertices of s
                    shared_queue_view<
                          typename boost::graph_traits<Subgraph>::vertex_descriptor
                        , typename boost::graph_traits<Graph>::vertex_descriptor
                    > local_queue(queue);

                    typename boost::graph_traits<Graph>::adjacency_iterator g_ai, g_ae;
                    for (boost::tie(s_ai, s_ae) = adjacent_vertices(s, S); s_ai != s_ae; ++s_ai)
                        if ( !local_queue.was_queued(*s_ai) ) {
                            local_queue.push_back(std::make_pair(*s_ai, g));
                        }
                    // take the first entry of the queue
                    // and check if entry.s can be mapped to adjacent vertices of entry.g
                    subgraph_vertex_descriptor t = local_queue.front().first;
                    boost::tie(g_ai, g_ae) = adjacent_vertices(local_queue.front().second, G);
                    local_queue.pop_front();
                    for (; g_ai != g_ae; ++g_ai)
                        if (!visited[*g_ai])
                            detail::lattice_constant_walker(
                                  t
                                , *g_ai
                                , S
                                , G
                                , local_queue
                                , visited
                                , visited_cnt
                                , pinning
                                , vertex_equal
                                , edge_equal
                                , register_embedding
                            );
                } else
                    register_embedding(S,G,pinning);
                pinning[s] = num_vertices(G);
                visited[g] = false;
            }

            // Input: Subgraph, Graph, vertices of G contained in mapping of S on G
            // Output: lattice_constant of S in G containing v
            template<typename Subgraph, typename Graph, typename VertexEqual, typename EdgeEqual, typename EmbeddingFoundPolicy> void lattice_constant_impl(
                  Subgraph const & S
                , Graph const & G
                , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const & V
                , typename partition_type<Subgraph>::type const & subgraph_orbit
                , VertexEqual & vertex_equal
                , EdgeEqual & edge_equal
                , EmbeddingFoundPolicy& embedding_found_policy
            ) {
                // Assume the vertex desciptor is an unsigned integer type (since we want to use it as an index for a vector)
                BOOST_STATIC_ASSERT((boost::is_unsigned<typename alps::graph_traits<Subgraph>::vertex_descriptor>::value));
                assert(num_vertices(S) > 0);
                // if larger, extend the space
                assert(num_vertices(S) < 21);
                assert(num_edges(S) < 21);

                BOOST_STATIC_ASSERT((boost::is_unsigned<typename alps::graph_traits<Graph>::vertex_descriptor>::value));
                assert(num_vertices(G) > 0);

                // make sure, that a distance in one direction fits in a boost::uint8_t
                assert(std::size_t(num_vertices(G)) < 256 * 256);

                boost::dynamic_bitset<> visited; // (num_vertices(G));
                shared_queue_data<
                      typename boost::graph_traits<Subgraph>::vertex_descriptor
                    , typename boost::graph_traits<Graph>::vertex_descriptor
                > queue(num_vertices(S));
                std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> pinning;    // (num_vertices(S), num_vertices(G));
                for (typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::const_iterator it = V.begin(); it != V.end(); ++it)
                    for (typename partition_type<Subgraph>::type::const_iterator jt = subgraph_orbit.begin(); jt != subgraph_orbit.end(); ++jt)
                        if (out_degree(jt->front(), S) <= out_degree(*it, G)) {
                            
                            // Reset all data
                            queue.reset();
                            visited.clear();
                            visited.resize(num_vertices(G));
                            pinning.clear();
                            pinning.resize(num_vertices(S), num_vertices(G));

                            queue.mark_queued(jt->front());
                            shared_queue_view<
                                  typename boost::graph_traits<Subgraph>::vertex_descriptor
                                , typename boost::graph_traits<Graph>::vertex_descriptor
                            > queue_view(queue);
                            lattice_constant_walker(
                                  jt->front()
                                , *it
                                , S
                                , G
                                , queue_view
                                , visited
                                , 0
                                , pinning
                                , vertex_equal
                                , edge_equal
                                , embedding_found_policy
                            );
                            break;
                        }
            }

            // Input: Subgraph, Graph, vertices of G contained in mapping of S on G
            // Output: lattice_constant of S in G containing v
            template<typename Subgraph, typename Graph, typename VertexEqual, typename EdgeEqual, typename EmbeddingFoundPolicy> void lattice_constant_impl_geometric(
                  Subgraph const & S
                , Graph const & G
                , typename boost::graph_traits<Graph>::vertex_descriptor const & lattice_pin
                , typename partition_type<Subgraph>::type const & subgraph_orbit
                , VertexEqual & vertex_equal
                , EdgeEqual & edge_equal
                , EmbeddingFoundPolicy& embedding_found_policy
                , typename boost::graph_traits<Subgraph>::vertex_descriptor breaking_vertex
            ) {
                // Assume the vertex desciptor is an unsigned integer type (since we want to use it as an index for a vector)
                BOOST_STATIC_ASSERT((boost::is_unsigned<typename alps::graph_traits<Subgraph>::vertex_descriptor>::value));
                assert(num_vertices(S) > 0);
                // if larger, extend the space
                assert(num_vertices(S) < 21);
                assert(num_edges(S) < 21);

                BOOST_STATIC_ASSERT((boost::is_unsigned<typename alps::graph_traits<Graph>::vertex_descriptor>::value));
                assert(num_vertices(G) > 0);

                // make sure, that a distance in one direction fits in a boost::uint8_t
                assert(std::size_t(num_vertices(G)) < 256 * 256);

                boost::dynamic_bitset<> visited; // (num_vertices(G));
                shared_queue_data<
                      typename boost::graph_traits<Subgraph>::vertex_descriptor
                    , typename boost::graph_traits<Graph>::vertex_descriptor
                > queue(num_vertices(S));
                std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> pinning;    // (num_vertices(S), num_vertices(G));
                if (out_degree(breaking_vertex, S) <= out_degree(lattice_pin, G)) {
                    // Reset all data
                    queue.reset();
                    visited.clear();
                    visited.resize(num_vertices(G));
                    pinning.clear();
                    pinning.resize(num_vertices(S), num_vertices(G));

                    queue.mark_queued(breaking_vertex);
                    shared_queue_view<
                          typename boost::graph_traits<Subgraph>::vertex_descriptor
                        , typename boost::graph_traits<Graph>::vertex_descriptor
                    > queue_view(queue);
                    lattice_constant_walker(
                          breaking_vertex
                        , lattice_pin
                        , S
                        , G
                        , queue_view
                        , visited
                        , 0
                        , pinning
                        , vertex_equal
                        , edge_equal
                        , embedding_found_policy
                    );
                }
            }

        } // end namespace detail
    } // end namespace graph
} // end namespace alps

#endif // ALPS_GRAPH_DETAIL_LATTICE_CONSTANT_IMPL_HPP
