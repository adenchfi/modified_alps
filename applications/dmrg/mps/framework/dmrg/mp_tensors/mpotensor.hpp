/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Sebastian Keller <sebkelle@phys.ethz.ch>
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

#include "dmrg/mp_tensors/reshapes.h"

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup>::MPOTensor(index_type ld,
                                        index_type rd,
                                        prempo_t const & tags,
                                        op_table_ptr tbl_)
: left_i(ld)
, right_i(rd)
, col_tags(ld, rd)
, operator_table(tbl_)
{
    using namespace boost::tuples;
    typedef boost::tuple<index_type, index_type, tag_type, value_type> prempo_descriptor;
    typedef std::vector<prempo_descriptor> converted_prempo_t;

    row_index.resize(ld);

    if (tags.size() > 0 && operator_table.get() != NULL) {
        converted_prempo_t tmp_tags;
        
        // copy (due to const &) and convert to index_type
        for (typename prempo_t::const_iterator it = tags.begin(); it != tags.end(); ++it) {
            index_type row_i = (left_i == 1) ? 0 : index_type(get<0>(*it));
            index_type col_i = (right_i == 1) ? 0 : index_type(get<1>(*it));
            tmp_tags.push_back( prempo_descriptor(row_i, col_i, get<2>(*it), get<3>(*it)) );
        }


        std::sort(tmp_tags.begin(), tmp_tags.end(), MPOTensor_detail::col_cmp<prempo_descriptor>());

        for (typename converted_prempo_t::const_iterator it = tmp_tags.begin(); it != tmp_tags.end(); ++it) {
            col_tags(get<0>(*it), get<1>(*it)) = std::make_pair(get<2>(*it), get<3>(*it));
            row_index[get<0>(*it)].insert(get<1>(*it));
        }
    }
    else {
        // Initialize a private operator table
        operator_table = op_table_ptr(new OPTable<Matrix, SymmGroup>());
    }
}

/*
template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & MPOTensor<Matrix, SymmGroup>::operator()(index_type left_index,
                                                                         index_type right_index) const
{
    throw std::runtime_error("operator() doesn't work for MPOTensors anymore!\n");
    assert( left_index < left_i );
    assert( right_index < right_i );
    return (*operator_table)[col_tags(left_index, right_index).first];
}


template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & MPOTensor<Matrix, SymmGroup>::operator()(index_type left_index,
                                                                         index_type right_index)
{
    throw std::runtime_error("operator() doesn't work for MPOTensors anymore!\n");
    assert( left_index < left_i );
    assert( right_index < right_i );
    typename CSCMatrix::value_type const & p = col_tags(left_index, right_index);
    return (*operator_table)[p.first];
}
*/

template<class Matrix, class SymmGroup>
bool MPOTensor<Matrix, SymmGroup>::has(index_type left_index,
                                       index_type right_index) const
{
    assert(left_index < left_i && right_index < right_i);
    return col_tags.find_element(left_index, right_index) != NULL;
}

// warning: this method allows to (indirectly) change the op in the table, all tags pointing to it will
//          get a modified matrix!
//          better design needed
template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::set(index_type li, index_type ri, op_t const & op, value_type scale_){
    if (this->has(li, ri)) {
        col_tags.find_element(li, ri)->second = scale_;
        (*operator_table)[col_tags.find_element(li, ri)->first] = op;
    }
    else {
        tag_type new_tag = operator_table->register_op(op);
        col_tags(li, ri) = internal_value_type(new_tag, scale_);
        row_index[li].insert(ri);
    }
}

template<class Matrix, class SymmGroup>
MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup>
MPOTensor<Matrix, SymmGroup>::at(index_type left_index, index_type right_index) const {
    assert(this->has(left_index, right_index));
    typename CSCMatrix::value_type const & p = col_tags(left_index, right_index);
    return MPOTensor_detail::make_const_term_descriptor((*operator_table)[p.first], p.second);
}

// warning: this method allows to (indirectly) change the op in the table, all tags pointing to it will
//          get a modified matrix!
//          better design needed
template<class Matrix, class SymmGroup>
MPOTensor_detail::term_descriptor<Matrix, SymmGroup>
MPOTensor<Matrix, SymmGroup>::at(index_type left_index, index_type right_index) {
    if (!this->has(left_index, right_index))
        this->set(left_index, right_index, op_t(), 1.);
    typename CSCMatrix::value_type & p = col_tags(left_index, right_index).ref();
    return MPOTensor_detail::term_descriptor<Matrix, SymmGroup>((*operator_table)[p.first], p.second);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::row_proxy MPOTensor<Matrix, SymmGroup>::row(index_type row_i) const
{  
    return row_proxy(row_index[row_i].begin(), row_index[row_i].end());
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::col_proxy MPOTensor<Matrix, SymmGroup>::column(index_type col_i) const
{  
    return col_proxy(col_tags, col_i);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::tag_type
MPOTensor<Matrix, SymmGroup>::tag_number(index_type left_index, index_type right_index) const {
    return col_tags(left_index, right_index).first;
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::multiply_by_scalar(const scalar_type& v)
{
    for (typename CSCMatrix::iterator2 it2 = col_tags.begin2(); it2 != col_tags.end2(); ++it2)
        for (typename CSCMatrix::iterator1 it1 = it2.begin(); it1 != it2.end(); ++it1)
            it1->second *= v; 
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::divide_by_scalar(const scalar_type& v)
{
    for (typename CSCMatrix::iterator2 it2 = col_tags.begin2(); it2 != col_tags.end2(); ++it2)
        for (typename CSCMatrix::iterator1 it1 = it2.begin(); it1 != it2.end(); ++it1)
            it1->second /= v; 
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::op_table_ptr MPOTensor<Matrix, SymmGroup>::get_operator_table() const
{
    return operator_table;
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::index_type MPOTensor<Matrix, SymmGroup>::row_dim() const
{
    return left_i;
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::index_type MPOTensor<Matrix, SymmGroup>::col_dim() const
{
    return right_i;
}
