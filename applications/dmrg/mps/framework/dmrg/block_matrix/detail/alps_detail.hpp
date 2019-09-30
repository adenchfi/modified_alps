/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Timothee Ewart <timothee.ewart@gmail.com>
 *                            Alexandr Kosenkov <alex.kosenkov@gmail.com>
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

#ifndef MAQUIS_BLOCK_MATRIX_DEATAIL_ALPS_MATRIX_DETAIL_HPP
#define MAQUIS_BLOCK_MATRIX_DEATAIL_ALPS_MATRIX_DETAIL_HPP

template<class T, class SymmGroup>
class block_matrix;

namespace maquis { namespace dmrg { namespace detail {
        
    template<class InputIterator, class OutputIterator, class T>
    void iterator_axpy(InputIterator in1, InputIterator in2,
                       OutputIterator out1, T val)
    {
        std::transform(in1, in2, out1, out1, boost::lambda::_1*val+boost::lambda::_2);
    }
    
    inline void iterator_axpy(double const * in1, double const * in2,
                              double * out1, double val)
    {
        fortran_int_t one = 1, diff = in2-in1;
#if defined(BIND_FORTRAN_LOWERCASE) || defined(__xlC__)
        daxpy(&diff, &val, in1, &one, out1, &one);
#else
        daxpy_(&diff, &val, in1, &one, out1, &one);
#endif
    }
    
    inline void iterator_axpy(std::complex<double> const * in1, std::complex<double> const * in2,
                              std::complex<double> * out1, double val)
    {
        throw std::runtime_error("Not implemented.");
    }

    //template <typename T, class A>
    //void op_kron(alps::numeric::matrix<T,A>& out, const alps::numeric::matrix<T,A>& in, const alps::numeric::matrix<T,A>& alfa,
    template <class Matrix1, class Matrix2>
    void op_kron(Matrix2& out, const Matrix1& in, const Matrix1& alfa,
                 size_t out_y_offset, size_t out_x_offset, 
                 size_t ldim1, size_t ldim2, 
                 size_t rdim1, size_t rdim2)
    {
            for(int l1 = 0; l1 < ldim1; ++l1)
            for(int r1 = 0; r1 < rdim1; ++r1)
            for(int l2 = 0; l2 < ldim2; ++l2)
            for(int r2 = 0; r2 < rdim2; ++r2)
                out(out_y_offset + l1*ldim2 + l2, out_x_offset + r1*rdim2 + r2) = 
                in(l2, r2)*alfa(l1, r1);
    }
    
    template <typename T, class A>
    void reshape_l2b(alps::numeric::matrix<T,A>& out, const alps::numeric::matrix<T,A>& in,
                     size_t in_left_offset, size_t in_phys_offset, 
                     size_t out_left_offset, size_t out_right_offset,
                     size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
        for(size_t ss2 = 0; ss2 < sdim2; ++ss2)
        {
            size_t ss_out = in_phys_offset + ss1*sdim2 + ss2;
            for(size_t rr = 0; rr < rdim; ++rr)
                for(size_t ll = 0; ll < ldim; ++ll)
                    out(out_left_offset + ss1*ldim + ll, out_right_offset + ss2*rdim + rr) =
                    in(in_left_offset + ss_out*ldim + ll, rr);
        }
    }

    template <typename T, class A>
    void reshape_b2l(alps::numeric::matrix<T,A>& out, const alps::numeric::matrix<T,A>& in,
                     size_t in_left_offset, size_t in_right_offset, 
                     size_t out_left_offset, size_t out_phys_offset,
                     size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
        for(size_t ss2 = 0; ss2 < sdim2; ++ss2)
        {
            size_t ss_out = out_phys_offset + ss1*sdim2 + ss2;
            for(size_t rr = 0; rr < rdim; ++rr)
                for(size_t ll = 0; ll < ldim; ++ll)
                    out(out_left_offset + ss_out*ldim + ll, rr) = 
                    in(in_left_offset + ss1*ldim+ll, in_right_offset + ss2*rdim+rr);
        }
    }
    
    template <typename T, class A1, class A2>
    void reshape_r2l(alps::numeric::matrix<T,A1>& left, const alps::numeric::matrix<T,A2>& right,
                     size_t left_offset, size_t right_offset, 
                     size_t sdim, size_t ldim, size_t rdim)
    {
        for (size_t ss = 0; ss < sdim; ++ss)
            for (size_t rr = 0; rr < rdim; ++rr)
                for(size_t ll = 0; ll < ldim; ++ll)
                    left(left_offset + ss*ldim+ll, rr) = 
                    right(ll, right_offset + ss*rdim+rr);
        // memcpy(&left(left_offset + ss*ldim, rr),
        //        &right(0, right_offset + ss*rdim+rr),
        //        sizeof(T) * ldim);
    }
    
    template <typename T, class A1, class A2>
    void reshape_l2r(const alps::numeric::matrix<T,A1>& left, alps::numeric::matrix<T,A2>& right,
                     size_t left_offset, size_t right_offset, 
                     size_t sdim, size_t ldim, size_t rdim)
    {
        for (size_t ss = 0; ss < sdim; ++ss)
            for (size_t rr = 0; rr < rdim; ++rr)
                for (size_t ll = 0; ll < ldim; ++ll)
                    right(ll, right_offset + ss*rdim+rr) = left(left_offset + ss*ldim+ll, rr);
    }
    
    template <typename T1, class A1,
              typename T2, class A2,
              typename T3, class A3>
    void lb_tensor_mpo(alps::numeric::matrix<T1,A1>& out, const alps::numeric::matrix<T2,A2>& in, const alps::numeric::matrix<T3,A3>& alfa,
                       size_t out_offset, size_t in_offset, 
                       size_t sdim1, size_t sdim2, size_t ldim, size_t rdim, T2 alfa_scale)
    {
        for(size_t rr = 0; rr < rdim; ++rr) {
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1) {
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2) {
                    T3 alfa_t = alfa(ss1, ss2) * alfa_scale;
                    iterator_axpy(&in(0, in_offset + ss1*rdim + rr),
                                  &in(0, in_offset + ss1*rdim + rr) + ldim, // bugbug
                                  &out(out_offset + ss2*ldim, rr),
                                  alfa_t);
                }
            }
        }
    }
    
    template <typename T1, class A1,
              typename T2, class A2,
              typename T3, class A3>
    void rb_tensor_mpo(alps::numeric::matrix<T1,A1>& out, const alps::numeric::matrix<T2,A2>& in, const alps::numeric::matrix<T3,A3>& alfa,
                       size_t out_offset, size_t in_offset, 
                       size_t sdim1, size_t sdim2, size_t ldim, size_t rdim, T2 alfa_scale)
    {
        for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
            for(size_t ss2 = 0; ss2 < sdim2; ++ss2) {
                T3 alfa_t = alfa(ss1, ss2) * alfa_scale;
                for(size_t rr = 0; rr < rdim; ++rr)
                    for(size_t ll = 0; ll < ldim; ++ll) {
                        out(ll, out_offset + ss2*rdim+rr) += in(ll + in_offset + ss1*ldim, rr) * alfa_t;
                    }
            }
    }
    
    template <typename T1, class A1,
              typename T2, class A2,
              typename T3, class A3>
    void mwt(alps::numeric::matrix<T1,A1>& out, const alps::numeric::matrix<T2,A2>& in, const alps::numeric::matrix<T3,A3>& alfa,
             size_t out_y_offset,  size_t out_x_offset,
             size_t in_y_offset,   size_t in_x_offset,
             size_t alfa_y_offset, size_t alfa_x_offset,
             size_t ldim,   size_t rdim,
             size_t lpdim,  size_t rpdim,
             size_t ilpdim, size_t irpdim)
    {
        for(size_t ll = 0; ll < ldim; ++ll)
            for(size_t rr = 0; rr < rdim; ++rr)
                for(size_t lp = 0; lp < lpdim; ++lp)
                    for(size_t rp = 0; rp < rpdim; ++rp)
                        for(size_t ilp = 0; ilp < ilpdim; ++ilp)
                            for(size_t irp = 0; irp < irpdim; ++irp)
                                out(out_y_offset + lp*ldim + ll, out_x_offset + rp*rdim + rr) +=  
                                in(in_y_offset + ilp*ldim + ll, in_x_offset + irp*rdim + rr) *
                                alfa(alfa_y_offset + ilp*irpdim + irp, alfa_x_offset + lp*rpdim + rp);
    }
    
    template<class T, class SymmGroup>
    std::vector<double> bond_renyi_entropies(const block_matrix<alps::numeric::diagonal_matrix<T>, SymmGroup>& set){
        std::vector<double> ret;
        for(std::size_t k = 0; k < set.n_blocks(); ++k){
            for (typename alps::numeric::diagonal_matrix<T>::const_diagonal_iterator it = diagonal(set[k]).first;
                 it != diagonal(set[k]).second; ++it)
            {
                double a = std::abs(*it);
                if (a > 1e-10)
                    ret.push_back(a*a);
            }
        }
        return ret;
    }
    
    template <typename T, class A>
    void left_right_boundary_init(alps::numeric::matrix<T,A> & M){
        //            memset((void*)&M(0,0),1,num_rows(M)*num_cols(M)*sizeof(T));
        for_each(elements(M).first,elements(M).second, boost::lambda::_1 = 1); // boost::lambda ^^' because iterable matrix concept 
    }
    
} } } // namespace maquis::dmrg::detail

#endif // MAQUIS_BLOCK_MATRIX_DEATAIL_ALPS_MATRIX_DETAIL_HPP
