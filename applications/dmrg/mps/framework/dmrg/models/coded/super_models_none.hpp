/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *               2012      by Jan Gukelberger <gukelberger@phys.ethz.ch>
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

#ifndef SUPER_MODELS_NONE_HPP
#define SUPER_MODELS_NONE_HPP

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"


/*
 * Time evolution for Bose-Hubbard density operators.
 *
 * The density operator rho is represented as a superstate whose time evolution
 * is generated by a superoperator according to the Lindblad master equation
 *      d rho/dt = (-i ad H + lind L) rho .
 * In this framework the time evolution of the density operator is computed
 * just like that of a pure state
 *      rho(t) = exp{ -i (ad H + i lind L) t } rho(0) .
 *
 * Operators are represented in the occupation number basis
 *      O_{ij} = <i|O|j>
 * with basis states {|i>}, i=0,...,N-1 for each site.
 * The corresponding superstate is obtained by the mapping
 *      O_{ij} -> O_n, n=i*N+j ,
 * i.e. the basis for operator states is given by the matrix elements
 *      |i><j|
 *
 * Superoperators represent the map
 *      rho -> d rho/dt .
 * We mostly need the adjoint Hamiltonian superoperator
 *      ad H: rho -> [H,rho]
 * and the dissipative superoperator
 *      lind L: rho -> 2 L rho L^\dag - L^\dag L rho - rho L^\dag L .
 * The superoperator representation of any linear map is conveniently computed
 * by applying the map to all N^2 operator basis states, filling the
 * corresponding columns of the superoperator.
 */


/// Take a Hamilton operator H and construct its adjoint Hamiltonian
/// superoperator, corresponding to the map
///   ad H: rho -> [H,rho]
template<class Matrix>
Matrix adjoint_hamiltonian(const Matrix& h)
{
    const size_t N = num_rows(h);
    const size_t N2 = N*N;
    assert( num_cols(h) == N );

    Matrix adH(N2,N2);
    for( size_t i = 0; i < N; ++i )
    {
        for( size_t j = 0; j < N; ++j )
        {
            const size_t n = i + j*N;
            Matrix phi (N,N);   phi(i,j) = 1.;      // n'th operator basis state phi
            Matrix hphi(N,N);   gemm(h,phi,hphi);   // H.phi
            Matrix phih(N,N);   gemm(phi,h,phih);   // phi.H
            Matrix comm = hphi - phih;              // [H,phi]
            
            // operator state adH.phi -> n'th column of adH
            for( size_t ii = 0; ii < N; ++ii )
                for( size_t jj = 0; jj < N; ++jj )
                    adH(ii+jj*N, n) = comm(ii,jj);
        }
    }
    
    return adH;
}

/// Take a Lindblad operator L and construct its contribution to the Liouville 
/// superoperator, i.e. the map
///   lind L: rho -> 2 L rho L^\dag - L^\dag L rho - rho L^\dag L
template<class Matrix>
Matrix super_lindblad(const Matrix& l)
{
    const size_t N = num_rows(l);
    const size_t N2 = N*N;
    assert( num_cols(l) == N );

    Matrix ldag = transpose(conj(l));       // L^\dag
    Matrix ldagl(N,N); gemm(ldag,l,ldagl);  // L^\dag L
    Matrix lindL(N2,N2);

    for( size_t i = 0; i < N; ++i )
    {
        for( size_t j = 0; j < N; ++j )
        {
            const size_t n = i + j*N;
            Matrix phi (N,N);   phi(i,j) = 1.;  // n'th operator basis state phi

            // phi -> lindL.phi = 2*l.phi.ldag - ldag.l.phi - phi.ldag.l
            Matrix lphi    (N,N); gemm(l,phi,lphi);
            Matrix lphildag(N,N); gemm(lphi,ldag,lphildag);
            Matrix ldaglphi(N,N); gemm(ldag,lphi,ldaglphi);
            Matrix phildagl(N,N); gemm(phi,ldagl,phildagl);
            Matrix ll = 2*lphildag - ldaglphi - phildagl;
            
            // operator state lindL.phi -> n'th column of lindL
            for( size_t ii = 0; ii < N; ++ii )
                for( size_t jj = 0; jj < N; ++jj )
                    lindL(ii+jj*N, n) = ll(ii,jj);
        }
    }

    return lindL;
}

/// Take a Lindblad operator L and construct the map corresponding to left
/// multiplication with L
///   rho -> L rho
template<class Matrix>
Matrix super_left(const Matrix& l)
{
    const size_t N = num_rows(l);
    const size_t N2 = N*N;
    assert( num_cols(l) == N );

    Matrix superL(N2,N2);

    for( size_t i = 0; i < N; ++i )
    {
        for( size_t j = 0; j < N; ++j )
        {
            const size_t n = i + j*N;
            Matrix phi(N,N);   phi(i,j) = 1.;  // n'th operator basis state phi

            // phi -> L.phi
            Matrix lphi(N,N); gemm(l,phi,lphi);
            for( size_t ii = 0; ii < N; ++ii )
                for( size_t jj = 0; jj < N; ++jj )
                    superL(ii+jj*N, n) = lphi(ii,jj);
        }
    }

    return superL;
}

/// Take a Lindblad operator L and construct the map corresponding to right
/// multiplication with L
///   rho -> rho L
template<class Matrix>
Matrix super_right(const Matrix& l)
{
    const size_t N = num_rows(l);
    const size_t N2 = N*N;
    assert( num_cols(l) == N );

    Matrix superL(N2,N2);

    for( size_t i = 0; i < N; ++i )
    {
        for( size_t j = 0; j < N; ++j )
        {
            const size_t n = i + j*N;
            Matrix phi (N,N);   phi(i,j) = 1.;  // n'th operator basis state phi

            // phi -> phi.L
            Matrix phil(N,N); gemm(phi,l,phil);
            for( size_t ii = 0; ii < N; ++ii )
                for( size_t jj = 0; jj < N; ++jj )
                    superL(ii+jj*N, n) = phil(ii,jj);
        }
    }

    return superL;
}


/// Fuse indices n[i] into one p = \sum_i n[i] d^i
template<class T, class A>
T fuse(const A& ind, T d)
{
    T fused = 0;
    T stride = 1;
    for( T i = 0; i < ind.size(); ++i )
    {
        fused += ind[i] * stride;
        stride *= d;
    }
    return fused;
}

///
template<class T, class A>
void unfuse(T fused, T d, A& ind)
{
    for( T i = 0; i < ind.size(); ++i )
    {
        ind[i] = fused % d;
        fused /= d;
    }
    assert( fused == 0 );
}



///
template<class Matrix>
Matrix reshape_bond2site(const Matrix& a)
{
    typedef typename Matrix::size_type size_type;
    size_type d4 = num_rows(a);
    size_type d = sqrt(sqrt(double(d4)));
    assert( d4 == num_cols(a) );
    assert( d4 == d*d*d*d );
    
    Matrix b(d4,d4);
    boost::array<size_type,4> ii, jj, kk, ll;
    for( size_type i = 0; i < d4; ++i )
    {
        unfuse(i,d,ii);
        for( size_type j = 0; j < d4; ++j )
        {
            unfuse(j,d,jj);
            kk[0] = ii[0]; kk[1] = ii[2]; kk[2] = jj[0]; kk[3] = jj[2];
            ll[0] = ii[1]; ll[1] = ii[3]; ll[2] = jj[1]; ll[3] = jj[3];
            b(fuse(kk,d),fuse(ll,d)) = a(i,j);
        }
    }
    return b;
}

///
template<class Op,class Matrix>
std::vector< std::pair<Op,Op> > decompose_bond_super(const Matrix& bondop, const Index<TrivialGroup>& phys)
{
    Matrix rbond = reshape_bond2site(bondop);
    
    Matrix U, V;
    typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type S, Ssqrt;
    svd(rbond, U, V, S);
    Ssqrt = sqrt(S);
    
    TrivialGroup::charge C = TrivialGroup::IdentityCharge;
    Op left, right;
    {
        Matrix tmp;
        gemm(U, Ssqrt, tmp);
        left.insert_block(tmp, C, C);
    }
    {
        Matrix tmp;
        gemm(Ssqrt, V, tmp);
        right.insert_block(tmp, C, C);
    }
    
    std::vector<Op> leftops  = reshape_right_to_list(phys, left);
    std::vector<Op> rightops = reshape_left_to_list (phys, right);
    assert(leftops.size() == rightops.size());
    
    // discard terms with no weight
    std::vector< std::pair<Op,Op> > terms;
    for( unsigned i = 0; i < num_rows(S) && std::abs(S(i,i)) > 1e-10; ++i )
    {
        leftops[i].transpose_inplace();
        rightops[i].transpose_inplace();
        terms.push_back(std::pair<Op,Op>( leftops[i], rightops[i] ));
    }
    return terms;
}


/* ****************** BOSE-HUBBARD */
template<class Matrix>
class SuperBoseHubbardNone : public model_impl<Matrix, TrivialGroup>
{
    typedef model_impl<Matrix, TrivialGroup> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurement_term_desc_type measurement_term_desc_type;
    
    typedef typename base::size_t size_t;
    typedef typename Matrix::value_type value_type;
    
public:
    // Dissipation needs complex types, that's why we forward to do_init with a tag
    SuperBoseHubbardNone(const Lattice& lat, BaseParameters & model_)
    : model(model_)
    , lattice(lat)
    , tag_handler(new table_type())
    {
        do_init(lat,model,typename Matrix::value_type());
    }
    
    void do_init(const Lattice& lat, BaseParameters & model_, double)
    {
        throw std::runtime_error("need complex value type");
    }
    
    void do_init(const Lattice& lat, BaseParameters & model_, std::complex<double>)
    {
        // retrieve model parameters
        int Nmax      = model["Nmax"];
        double t      = model["t"];
        double U      = model["U"];
        double mu     = model["mu"];
        double omega  = model["omega"];
        double V      = model["V"];
        double Lambda = model["Lambda"];
        double Delta  = model["Delta"];
        double Gamma1a = model["Gamma1a"];
        double Gamma1b = model["Gamma1b"];
        double Gamma2  = model["Gamma2"];
        double nbar    = model["nbar"];
        
        TrivialGroup::charge C = TrivialGroup::IdentityCharge;
        const size_t N = Nmax+1;
        const size_t N2 = N*N;
        const std::complex<double> I(0,1);
        
        phys_psi.insert(std::make_pair(C, N));
        
        phys.insert(std::make_pair(C, N2));
        ident.insert_block(Matrix::identity_matrix(N2), C, C);
        
        // construct basic on-site operators
        mcount.resize(N,N); minteraction.resize(N,N); mcreate.resize(N,N); mdestroy.resize(N,N);
        for( int n = 1; n <= Nmax; ++n )
        {
            mcount(n,n) = n;
            if ((n*n-n) != 0)
                minteraction(n,n) = n*n-n;
            
            mcreate(n,n-1) = std::sqrt(value_type(n));   // input n-1, output n
            mdestroy(n-1,n) = std::sqrt(value_type(n));  // input n,   output n-1
        }
        Matrix mcreate2 (N,N);   gemm(mcreate ,mcreate ,mcreate2 );
        Matrix mdestroy2(N,N);   gemm(mdestroy,mdestroy,mdestroy2);
        
        // construct on-site superoperators
        Matrix screate      = adjoint_hamiltonian(mcreate);
        Matrix sdestroy     = adjoint_hamiltonian(mdestroy);
        Matrix sdrive       = adjoint_hamiltonian(mcreate+mdestroy);
        Matrix spump        = adjoint_hamiltonian(mcreate2+mdestroy2);
        Matrix scount       = adjoint_hamiltonian(mcount);
        Matrix sinteraction = adjoint_hamiltonian(minteraction);
        
        Matrix ldestroy  = super_lindblad(mdestroy );
        Matrix lcreate   = super_lindblad(mcreate  );
        Matrix ldestroy2 = super_lindblad(mdestroy2);
        
        // cast superoperators to op_t
        create     .insert_block(transpose(screate     ), C,C);
        destroy    .insert_block(transpose(sdestroy    ), C,C);
        drive      .insert_block(transpose(sdrive      ), C,C);
        pump       .insert_block(transpose(spump       ), C,C);
        count      .insert_block(transpose(scount      ), C,C);
        interaction.insert_block(transpose(sinteraction), C,C);

        lindDestroy .insert_block(transpose(ldestroy     ), C,C);
        lindCreate  .insert_block(transpose(lcreate      ), C,C);
        lindDestroy2.insert_block(transpose(ldestroy2    ), C,C);
        
        leftDestroy.insert_block(transpose(super_left(mdestroy)), C,C);
        rightCreate.insert_block(transpose(super_right(mcreate)), C,C);
        
        std::vector< std::pair<op_t,op_t> > hopops = decompose_bond_super<op_t>(adjoint_hamiltonian(kron(mcreate, mdestroy)),phys);
        std::vector< std::pair<op_t,op_t> > Vops = decompose_bond_super<op_t>(adjoint_hamiltonian(kron(mcount, mcount)),phys);
        
        
#define REGISTER(op, kind) op ## _tag = tag_handler->register_op(op, kind);
        REGISTER(ident,       tag_detail::bosonic)
        REGISTER(create,      tag_detail::bosonic)
        REGISTER(destroy,     tag_detail::bosonic)
        REGISTER(drive,       tag_detail::bosonic)
        REGISTER(pump,        tag_detail::bosonic)
        REGISTER(count,       tag_detail::bosonic)
        REGISTER(interaction, tag_detail::bosonic)
        REGISTER(lindDestroy, tag_detail::bosonic)
        REGISTER(lindCreate,  tag_detail::bosonic)
        REGISTER(lindDestroy2,tag_detail::bosonic)
        REGISTER(leftDestroy, tag_detail::bosonic)
        REGISTER(rightCreate, tag_detail::bosonic)
        
        std::vector< std::pair<tag_type,tag_type> > hopops_tag(hopops.size());
        for (size_t i=0; i<hopops.size(); ++i)
            hopops_tag[i] = std::make_pair( tag_handler->register_op(hopops[i].first,  tag_detail::bosonic),
                                            tag_handler->register_op(hopops[i].second, tag_detail::bosonic) );
        std::vector< std::pair<tag_type,tag_type> > Vops_tag(Vops.size());
        for (size_t i=0; i<Vops.size(); ++i)
            Vops_tag[i] = std::make_pair( tag_handler->register_op(Vops[i].first,  tag_detail::bosonic),
                                          tag_handler->register_op(Vops[i].second, tag_detail::bosonic) );
        
        
#undef REGISTER

        // insert superoperators for each site
        for( int p=0; p < lat.size(); ++p ) 
        {
            // interaction H_U = U/2 n_i (n_i - 1)
            if( U != 0 )
            {
                term_descriptor term;
                term.coeff = 0.5*U;
                term.push_back( boost::make_tuple(p, interaction_tag) );
                this->terms_.push_back(term);
            }
            
            // site-dependent chemical potential H_mu = -mu n_i
            // harmonic trap H_omega = (omega x_i)^2/2 n_i
            double x = p - 0.5*lat.size();
            double mup = -mu + 0.5*omega*omega*x*x;
            if( mup != 0 )
            {
                term_descriptor term;
                term.coeff = mup;
                term.push_back( boost::make_tuple(p, count_tag) );
                this->terms_.push_back(term);
            }
            
            // drive H_Lambda = Lambda (b_i^\dag + b_i)
            if( Lambda != 0 )
            {
                term_descriptor term;
                term.coeff = Lambda;
                term.push_back( boost::make_tuple(p, drive_tag) );
                this->terms_.push_back(term);
            }

            // pump H_Delta = Delta (b_i^\dag^2 + b_i^2)
            if( Delta != 0 )
            {
                term_descriptor term;
                term.coeff = Delta;
                term.push_back( boost::make_tuple(p, pump_tag) );
                this->terms_.push_back(term);
            }

            // one-boson dissipation L_{1a} = Gamma_{1a} (1+\bar{n}) lind b_i
            //   = Gamma_{1a} (1+\bar{n}) (2 b_i rho b_i^\dag - b_i^\dag b_i rho - rho b_i^\dag b_i)
            if( Gamma1a != 0 )
            {
                term_descriptor term;
                term.coeff = I*Gamma1a*(1+nbar);
                term.push_back( boost::make_tuple(p, lindDestroy_tag) );
                this->terms_.push_back(term);
            }
            
            // one-boson dissipation at finite temperature (thermal population \bar{n})
            // L_{1a} = Gamma_{1a} \bar{n} lind b_i^\dag
            //   = Gamma_{1a} \bar{n} (2 b_i^\dag rho b_i - b_i b_i^\dag rho - rho b_i b_i^\dag)
            if( Gamma1a*nbar != 0 )
            {
                term_descriptor term;
                term.coeff = I*Gamma1a*nbar;
                term.push_back( boost::make_tuple(p, lindCreate_tag) );
                this->terms_.push_back(term);
            }
            
            // two-boson dissipation L_2 = Gamma_2/2 lind b_i^2
            //   = Gamma_2/2 (2 b_i^2 rho b_i^\dag^2 - b_i^\dag^2 b_i^2 rho - rho b_i^\dag^2 b_i^2)
            if( Gamma2 != 0 )
            {
                term_descriptor term;
                term.coeff = I*0.5*Gamma2;
                term.push_back( boost::make_tuple(p, lindDestroy2_tag) );
                this->terms_.push_back(term);
            }
            
            // bond terms
            std::vector<int> neighs = lat.forward(p);
            for( int n = 0; n < neighs.size(); ++n ) 
            {
                // hopping H_J = -J (b_i^\dag b_{i+1} + b_i b_{i+1}^\dag)
                for( unsigned i = 0; i < hopops_tag.size(); ++i )
                {
                    {
                        term_descriptor term;
                        term.coeff = -t;
                        term.push_back( boost::make_tuple(p,         hopops_tag[i].first) );
                        term.push_back( boost::make_tuple(neighs[n], hopops_tag[i].second) );
                        this->terms_.push_back(term);
                    }
                    {
                        term_descriptor term;
                        term.coeff = -t;
                        term.push_back( boost::make_tuple(p,         hopops_tag[i].second) );
                        term.push_back( boost::make_tuple(neighs[n], hopops_tag[i].first) );
                        this->terms_.push_back(term);
                    }
                }
                
                // nearest-neighbor interaction H_V = V n_i n_{i+1}
                if( V != 0 )
                {
                    for( unsigned i = 0; i < Vops_tag.size(); ++i )
                    {
                        {
                            term_descriptor term;
                            term.coeff = V;
                            term.push_back( boost::make_tuple(p,        Vops_tag[i].first) );
                            term.push_back( boost::make_tuple(neighs[n],Vops_tag[i].second) );
                            this->terms_.push_back(term);
                        }
                        {
                            term_descriptor term;
                            term.coeff = V;
                            term.push_back( boost::make_tuple(p,        Vops_tag[i].second) );
                            term.push_back( boost::make_tuple(neighs[n],Vops_tag[i].first) );
                            this->terms_.push_back(term);
                        }
                    }
                }
            
                // one-boson dissipation L_{1b} = -Gamma_{1b}/2 (
                //          2 b_{i+1} rho b_i^\dag - b_i^\dag b_{i+1} rho - rho b_i^\dag b_{i+1}
                //        + 2 b_i rho b_{i+1}^\dag - b_i b_{i+1}^\dag rho - rho b_i b_{i+1}^\dag )
                //  = Gamma_{1b}/2 (
                //          (ad b^\dag)_i (b rho)_{i+1} + (b rho)_i (ad b^\dag)_{i+1}
                //        - (ad b)_i (rho b^\dag)_{i+1} - (rho b^\dag)_i (ad b)_{i+1} )
                if( Gamma1b != 0 )
                {
                    term_descriptor term;
                    term.coeff = I*Gamma1b/2.;
                    term.push_back( boost::make_tuple(p,         create_tag) );
                    term.push_back( boost::make_tuple(neighs[n], leftDestroy_tag) );
                    this->terms_.push_back(term);
                }
                if( Gamma1b != 0 )
                {
                    term_descriptor term;
                    term.coeff = I*Gamma1b/2.;
                    term.push_back( boost::make_tuple(p,         leftDestroy_tag) );
                    term.push_back( boost::make_tuple(neighs[n], create_tag) );
                    this->terms_.push_back(term);
                }
                if( Gamma1b != 0 )
                {
                    term_descriptor term;
                    term.coeff = -I*Gamma1b/2.;
                    term.push_back( boost::make_tuple(p,         destroy_tag) );
                    term.push_back( boost::make_tuple(neighs[n], rightCreate_tag) );
                    this->terms_.push_back(term);
                }
                if( Gamma1b != 0 )
                {
                    term_descriptor term;
                    term.coeff = -I*Gamma1b/2.;
                    term.push_back( boost::make_tuple(p,         rightCreate_tag) );
                    term.push_back( boost::make_tuple(neighs[n], destroy_tag) );
                    this->terms_.push_back(term);
                }
            }
        }
        
    }
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }
    
    Index<TrivialGroup> const& phys_dim(size_t type) const
    {
        return phys;
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ident_tag;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return identity_matrix_tag(type);
    }
    typename TrivialGroup::charge total_quantum_numbers(BaseParameters & parms) const
    {
        return TrivialGroup::IdentityCharge;
    }
    
//    measurements_type measurements() const
//    {
//        TrivialGroup::charge C = TrivialGroup::IdentityCharge;
//        
//        op_t ident_psi = identity_matrix<Matrix>(phys_psi);
//        op_t count_psi, create_psi, destroy_psi;
//        count_psi.insert_block(mcount, C, C);
//        create_psi.insert_block(transpose(mcreate), C, C);
//        destroy_psi.insert_block(transpose(mdestroy), C, C);
//        
//        typedef std::vector<block_matrix<Matrix, TrivialGroup> > op_vec;
//        typedef std::vector<std::pair<op_vec, bool> > bond_element;
//        
//        measurements_type meas;
//
//        if (model["MEASURE[Density]"]) {
//            meas.push_back( new measurements::average<Matrix, TrivialGroup>("Density", lattice,
//                                                                            op_vec(1,ident_psi), op_vec(1,ident_psi),
//                                                                            op_vec(1,count_psi)) );
//            meas[meas.size()-1].set_super_meas(phys_psi);
//        }
//        
//        if (model["MEASURE[Local density]"]) {
//            meas.push_back( new measurements::local<Matrix, TrivialGroup>("Local density", lattice,
//                                                                          op_vec(1,ident_psi), op_vec(1,ident_psi),
//                                                                          op_vec(1,count_psi)) );
//            meas[meas.size()-1].set_super_meas(phys_psi);
//        }
//        
//        if (model["MEASURE[Local density^2]"]) {
//            op_t count2_psi;
//            gemm(count_psi, count_psi, count2_psi);
//            meas.push_back( new measurements::local<Matrix, TrivialGroup>("Local density^2", lattice,
//                                                                          op_vec(1,ident_psi), op_vec(1,ident_psi),
//                                                                          op_vec(1,count2_psi)) );
//            meas[meas.size()-1].set_super_meas(phys_psi);
//        }
//        
//        if (model["MEASURE[Onebody density matrix]"]) {
//            bond_element ops;
//            ops.push_back( std::make_pair(op_vec(1,create_psi), false) );
//            ops.push_back( std::make_pair(op_vec(1,destroy_psi), false) );
//            meas.push_back( new measurements::correlations<Matrix, TrivialGroup>("Onebody density matrix", lattice,
//                                                                                 op_vec(1,ident_psi), op_vec(1,ident_psi),
//                                                                                 ops, true, false) );
//            meas[meas.size()-1].set_super_meas(phys_psi);
//        }
//        
//        if (model["MEASURE[Density correlation]"]) {
//            bond_element ops;
//            ops.push_back( std::make_pair(op_vec(1,count_psi), false) );
//            ops.push_back( std::make_pair(op_vec(1,count_psi), false) );
//            meas.push_back( new measurements::correlations<Matrix, TrivialGroup>("Density correlation", lattice,
//                                                                                 op_vec(1,ident_psi), op_vec(1,ident_psi),
//                                                                                 ops, true, false) );
//            meas[meas.size()-1].set_super_meas(phys_psi);
//        }
//        
//        return meas;
//    }
    
    bool has_operator(std::string const & name, size_t type) const
    {
        return false;
    }
    
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        throw std::runtime_error("Operator not defined for this model.");
        return 0;
    }

    std::vector<measurement_term_desc_type> unpack_measurement_terms(std::string const & name) const
    {
        throw std::runtime_error("Function unpack_measurement_terms not implemented for this model class.");
        return std::vector<measurement_term_desc_type>();
    }

    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
private:
    BaseParameters & model;
    Lattice lattice;
    
    op_t ident;
    Matrix mcount, minteraction, mcreate, mdestroy;
    op_t create, destroy, drive, pump, count, interaction;
    op_t lindDestroy, lindCreate, lindDestroy2, leftDestroy, rightCreate;
    Index<TrivialGroup> phys, phys_psi;

    boost::shared_ptr<TagHandler<Matrix, TrivialGroup> > tag_handler;
    tag_type ident_tag;
    tag_type create_tag, destroy_tag, drive_tag, pump_tag, count_tag, interaction_tag;
    tag_type lindDestroy_tag, lindCreate_tag, lindDestroy2_tag, leftDestroy_tag, rightCreate_tag;
};



#endif
