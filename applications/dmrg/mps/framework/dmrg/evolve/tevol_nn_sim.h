/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_DMRG_TEVOL_NN_SIM_H
#define APP_DMRG_TEVOL_NN_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>

#include "dmrg/evolve/te_utils.hpp"
#include "dmrg/evolve/trotter_decomposer.h"
#include "dmrg/utils/results_collector.h"


// ******     HELPER OBJECTS     ******
template<class Matrix, class SymmGroup>
struct trotter_gate {
    typedef std::vector<long> idx_t;
    typedef std::vector<block_matrix<Matrix, SymmGroup> > vgates_t;
    
    std::size_t pfirst;
	idx_t idx;
	vgates_t vgates;
    
    trotter_gate(std::size_t L) : idx(L, -1) { }
    
    void add_term(std::size_t p, block_matrix<Matrix, SymmGroup> const & block)
    {
        assert(idx[p] == -1);
        vgates.push_back(block);
        idx[p] = vgates.size()-1;
    }
    
    void clear()
    {
        vgates.clear();
        idx = idx_t(idx.size(), -1);
        pfirst = 0;
    }
};



// ******   SIMULATION CLASS   ******
template <class Matrix, class SymmGroup>
class nearest_neighbors_evolver {
    typedef trotter_decomposer::steps_iterator steps_iterator;
public:
    nearest_neighbors_evolver(DmrgParameters * parms_, MPS<Matrix, SymmGroup> * mps_,
                              Lattice const& lattice_, Model<Matrix, SymmGroup> const& model_,
                              int init_sweep=0)
    : parms(parms_)
    , mps(mps_)
    , lattice(lattice_) // shallow copy
    , model(model_) // shallow copy
    , L(lattice.size())
    , decomp(2,  // NN time evolution has only two types of operators
             (*parms)["te_order"],
             (*parms)["te_optim"])
    , block_terms(hamil_to_blocks(lattice, model))
    {
        maquis::cout << "Using nearest-neighbors time evolution." << std::endl;
        maquis::cout << "Using " << decomp.description() << std::endl;
        
        /// compute the time evolution gates
        prepare_te_terms(init_sweep);
    }
    

    void prepare_te_terms(unsigned sweep)
    {
        double dt = (*parms)["dt"];
        typename Matrix::value_type I;
        if (sweep < (*parms)["nsweeps_img"])
            I = maquis::traits::real_identity<typename Matrix::value_type>::value;
        else
            I = maquis::traits::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I * dt;
        
        Uterms.resize(decomp.size(), trotter_gate<Matrix, SymmGroup>(L));
        for (size_t i=0; i<decomp.size(); ++i) {
            Uterms[i].clear();
            Uterms[i].pfirst = decomp.trotter_term(i).first;
            for (size_t p=decomp.trotter_term(i).first; p<L-1; p+=2){
                int type1 = lattice.get_prop<int>("type", p);
                int type2 = lattice.get_prop<int>("type", p+1);
                if ((*parms)["expm_method"] == "heev")
                    Uterms[i].add_term(p, op_exp_hermitian(model.phys_dim(type1)*model.phys_dim(type2), block_terms[p], decomp.trotter_term(i).second*alpha));
                else
                    Uterms[i].add_term(p, op_exp(model.phys_dim(type1)*model.phys_dim(type2), block_terms[p], decomp.trotter_term(i).second*alpha));
            }
        }
    }
    
    void operator()(unsigned sweep, unsigned nsteps)
    {
        iteration_results_.clear();
        
        for (steps_iterator it=decomp.steps_begin(nsteps); it != steps_iterator(); ++it)
            evolve_time_step(*it);
    }
    
    results_collector const& iteration_results() const
    {
        return iteration_results_;
    }
    
private:

    void evolve_time_step(std::vector<std::size_t> const & gates_i)
    {
        assert(gates_i.size() > 0);
        
        for (size_t i=0; i<gates_i.size(); ++i) {
            if (mps->canonization(true) < mps->length()/2)
                evolve_l2r(gates_i[i]);
            else
                evolve_r2l(gates_i[i]);
        }
    }
    
    void evolve_l2r(std::size_t gate_index)
    {
        std::size_t L = mps->length();
        std::vector<block_matrix<Matrix, SymmGroup> > const & ops = Uterms[gate_index].vgates;
        std::vector<long> const & idx = Uterms[gate_index].idx;
        int pfirst = Uterms[gate_index].pfirst;
        std::size_t Mmax=(*parms)["max_bond_dimension"];
        double cutoff=(*parms)["truncation_final"];
        MPS<Matrix, SymmGroup> const& constmps = *mps;
        
        assert(L == idx.size());

        if (mps->canonization() != pfirst + 1)
            mps->canonize(pfirst + 1);
        
        // TODO: remove pfirst!
        for (std::size_t p = pfirst; p <= L-1; p += 2)
        {
            if (idx[p] != -1)
            {
                constmps[p].make_left_paired();
                constmps[p+1].make_right_paired();
                
                block_matrix<Matrix, SymmGroup> v0, v1;
                gemm(constmps[p].data(), constmps[p+1].data(), v0); // outer product of two sites
                
                v1 = contraction::multiply_with_twosite(v0, ops[idx[p]],
                                                        constmps[p].row_dim(), constmps[p+1].col_dim(),
                                                        constmps[p].site_dim());
                truncation_results trunc = compression::replace_two_sites_l2r(*mps, Mmax, cutoff, v1, p);
                iteration_results_["BondDimension"]     << trunc.bond_dimension;
                iteration_results_["TruncatedWeight"]   << trunc.truncated_weight;
                iteration_results_["TruncatedFraction"] << trunc.truncated_fraction;
                iteration_results_["SmallestEV"]        << trunc.smallest_ev;
            }
            mps->move_normalization_l2r(p+1, p+3, DefaultSolver());
        }
        mps->canonization(true);
        assert(mps->canonization() == L-1);
        // maquis::cout << "Norm loss " << i << ": " << trace(t) << " " << -log(trace(t)) << std::endl;
    }
    
    void evolve_r2l(std::size_t gate_index)
    {
        std::size_t L = mps->length();
        std::vector<block_matrix<Matrix, SymmGroup> > const & ops = Uterms[gate_index].vgates;
        std::vector<long> const & idx = Uterms[gate_index].idx;
        int pfirst = Uterms[gate_index].pfirst;
        std::size_t Mmax=(*parms)["max_bond_dimension"];
        double cutoff=(*parms)["truncation_final"];
        MPS<Matrix, SymmGroup> const& constmps = *mps;

        assert(L == idx.size());
        
        int startpos = std::min(L-2-(L-pfirst)%2, L-2);
        if (mps->canonization() != startpos)
            mps->canonize(startpos);
        
        for (int p = std::min(L-2-(L-pfirst)%2, L-2); p >= pfirst; p -= 2)
        {
            if (idx[p] != -1)
            {
                constmps[p].make_left_paired();
                constmps[p+1].make_right_paired();
                
                block_matrix<Matrix, SymmGroup> v0, v1;
                gemm(constmps[p].data(), constmps[p+1].data(), v0); // outer product of two sites
                
                v1 = contraction::multiply_with_twosite(v0, ops[idx[p]],
                                                        constmps[p].row_dim(), constmps[p+1].col_dim(),
                                                        constmps[p].site_dim());
                truncation_results trunc = compression::replace_two_sites_r2l(*mps, Mmax, cutoff, v1, p);
                iteration_results_["BondDimension"]     << trunc.bond_dimension;
                iteration_results_["TruncatedWeight"]   << trunc.truncated_weight;
                iteration_results_["TruncatedFraction"] << trunc.truncated_fraction;
                iteration_results_["SmallestEV"]        << trunc.smallest_ev;
            }
            mps->move_normalization_r2l(p, std::max(static_cast<long>(p)-2,0L), DefaultSolver());
        }
        
        
        mps->canonization(true);
        assert(mps->canonization() == 0);
        // maquis::cout << "Norm loss " << i << ": " << trace(t) << " " << -log(trace(t)) << std::endl;
    }

private:
    DmrgParameters * parms;
    MPS<Matrix, SymmGroup> * mps;
    Lattice lattice;
    Model<Matrix, SymmGroup> model;
    size_t L;
    
    results_collector iteration_results_;
    
    trotter_decomposer decomp;
    std::vector<block_matrix<Matrix, SymmGroup> > block_terms;
    std::vector<trotter_gate<Matrix, SymmGroup> > Uterms;
};

#endif
