/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef BLOCK_MATRIX_H
#define BLOCK_MATRIX_H

#include <sstream>
#include <algorithm>
#include <numeric>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"

#include "utils/timings.h"
#include "utils/traits.hpp"
#include "dmrg/utils/storage.h"

#include <boost/ptr_container/ptr_vector.hpp>

template<class Matrix, class SymmGroup>
class block_matrix
{
    friend class block_matrix<typename storage::constrained<Matrix>::type, SymmGroup>;
private:
    typedef typename SymmGroup::charge charge;
public:
    typedef typename Matrix::size_type size_type;
    typedef typename Matrix::value_type value_type;
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef typename maquis::traits::real_type<Matrix>::type real_type;
   
    block_matrix();

    block_matrix(Index<SymmGroup> rows,
                 Index<SymmGroup> cols);
    
    block_matrix(block_matrix const&);

    template <class OtherMatrix>
    block_matrix(block_matrix<OtherMatrix,SymmGroup> const&);

    block_matrix& operator=(block_matrix rhs);
    template<class OtherMatrix>
    block_matrix& operator=(const block_matrix<OtherMatrix, SymmGroup>& rhs);

    Index<SymmGroup> const & left_basis() const;
    Index<SymmGroup> const & right_basis() const;

    void shift_basis(charge diff);

//  Remove by Tim 06/08/2012, presently not used in any DMRG/TE code
//  block_matrix(charge rc, charge cc, Matrix& m);

    std::string description() const;
    std::size_t num_elements() const;
    
    Matrix &             operator[](size_type c);
    Matrix const &       operator[](size_type c) const;
    value_type &         operator()(std::pair<charge, size_type> const & r,
                                    std::pair<charge, size_type> const & c);
    value_type const &   operator()(std::pair<charge, size_type> const & r,
                                    std::pair<charge, size_type> const & c) const;
    block_matrix &       operator+=(block_matrix const & rhs);
    block_matrix &       operator-=(block_matrix const & rhs);
    block_matrix const & operator*=(const scalar_type& v);
    block_matrix const & operator/=(const scalar_type& v);

    size_type n_blocks() const;
    size_type find_block(charge r, charge c) const;
    bool has_block(charge r, charge c) const;
    bool has_block(std::pair<charge, size_type> const & r,
                   std::pair<charge, size_type> const & c) const;
    
    size_type insert_block(Matrix const &, charge, charge);
    size_type insert_block(Matrix *, charge, charge);
    void remove_block(charge r, charge c);
    void remove_block(std::size_t which);

// Remove by Tim 06/08/2012, presently not used in any DMRG/TE code
//  void remove_rows_from_block(size_type block, size_type r, size_type k = 1);
//  void remove_cols_from_block(size_type block, size_type r, size_type k = 1);
    
    scalar_type trace() const;
    real_type norm() const;
    void transpose_inplace();
    void conjugate_inplace();
    void adjoint_inplace();
    void clear();
    template<class Generator>
    void generate(Generator g);

    void match_and_add_block(Matrix const &, charge, charge);
    
    void reserve(charge, charge, std::size_t, std::size_t);
    inline void reserve_pos(charge, charge, std::size_t, std::size_t);
    void allocate_blocks();
    
    void resize_block(charge r, charge c,
                      size_type new_r, size_type new_c,
                      bool pretend = false);
    
    friend void swap(block_matrix & x, block_matrix & y)
    {
        swap(x.data_, y.data_);
        swap(x.rows_, y.rows_);
        swap(x.cols_, y.cols_);
    }
    
    Matrix const & operator()(charge r, charge c) const
    {
        assert( has_block(r, c) );
        assert( rows_.position(r) == cols_.position(c) );
        return data_[rows_.position(r)];
    }
    
    Matrix & operator()(charge r, charge c)
    {
        assert( has_block(r, c) );
        assert( rows_.position(r) == cols_.position(c) );
        return data_[rows_.position(r)];
    }
    
    template<class Archive> void load(Archive & ar);
    template<class Archive> void save(Archive & ar) const;
    
    template <class Archive>
    inline void serialize(Archive & ar, const unsigned int version);
    
    bool reasonable() const;
    
private:
    Index<SymmGroup> rows_, cols_;
    boost::ptr_vector<Matrix> data_;
};    

#include "dmrg/block_matrix/block_matrix.hpp"
template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> operator*(const typename block_matrix<Matrix,SymmGroup>::scalar_type& v,
                                          block_matrix<Matrix, SymmGroup> bm)
{
    bm *= v;
    return bm;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> operator*(block_matrix<Matrix, SymmGroup> bm,
                                          const typename block_matrix<Matrix,SymmGroup>::scalar_type& v)
{
    bm *= v;
    return bm;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> operator+(block_matrix<Matrix,SymmGroup> b1, block_matrix<Matrix, SymmGroup> const& b2)
{
    b1 += b2;
    return b1;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> operator-(block_matrix<Matrix,SymmGroup> b1, block_matrix<Matrix, SymmGroup> const& b2)
{
    b1 -= b2;
    return b1;
}


template<class Matrix, class SymmGroup>
std::size_t size_of(block_matrix<Matrix, SymmGroup> const & m)
{
    size_t r = 0;
    for (size_t i = 0; i < m.n_blocks(); ++i)
        r += size_of(m[i]);
    return r;
}


#endif
