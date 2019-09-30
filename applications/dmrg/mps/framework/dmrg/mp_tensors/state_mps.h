/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2012-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MAQUIS_DMRG_STATE_MPS_H
#define MAQUIS_DMRG_STATE_MPS_H

#include "dmrg/mp_tensors/mps.h"

#include <boost/tuple/tuple.hpp>

template <class Matrix, class SymmGroup, class Charge>
MPS<Matrix, SymmGroup> state_mps(std::vector<boost::tuple<Charge, std::size_t> > const & state,
                                 std::vector<Index<SymmGroup> > const& phys_dims, std::vector<int> const& site_type)
{
    typedef typename SymmGroup::charge charge;
    typedef boost::tuple<charge, size_t> local_state;
    
    MPS<Matrix, SymmGroup> mps(state.size());
    
    Index<SymmGroup> curr_i;
    curr_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
    size_t curr_b = 0;
    for (int i=0; i<state.size(); ++i)
    {
        charge newc = SymmGroup::fuse(curr_i[0].first, boost::get<0>(state[i]));
        size_t news = 1;
        Index<SymmGroup> new_i;
        new_i.insert(std::make_pair(newc, news));
        ProductBasis<SymmGroup> left(phys_dims[site_type[i]], curr_i);
        mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], curr_i, new_i, false, 0);
        size_t b_in = left(boost::get<0>(state[i]), curr_i[0].first) + boost::get<1>(state[i]) * curr_i[0].second + curr_b;
        size_t b_out = 0;
        
        mps[i].make_left_paired();
        block_matrix<Matrix, SymmGroup> & block = mps[i].data();
        Matrix & m = block(SymmGroup::fuse(curr_i[0].first, boost::get<0>(state[i])), new_i[0].first);
        m(b_in, b_out) = 1.;
        
        curr_i = new_i;
        curr_b = b_out;
    }
    return mps;
}


#endif
