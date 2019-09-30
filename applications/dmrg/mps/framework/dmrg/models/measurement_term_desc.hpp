/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MEASUREMENT_TERM_DESC_HPP
#define MEASUREMENT_TERM_DESC_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>

template <class T>
struct MeasurementTermDescriptor {
    typedef std::vector<std::string> op_prod_type;
    T coeff;
    std::vector<op_prod_type> op_names;
    
    MeasurementTermDescriptor() { }
    MeasurementTermDescriptor(T const& coeff_, std::vector<op_prod_type> const& op_names_)
    : coeff(coeff_)
    , op_names(op_names_)
    { }
};

template <class T>
std::ostream & operator<<(std::ostream & os, MeasurementTermDescriptor<T> const& meas_term_desc)
{
    os << "coeff=" << meas_term_desc.coeff;
    os << ", op_names=[";
    for (unsigned i=0; i<meas_term_desc.op_names.size(); ++i) {
        os << "[";
        std::copy(meas_term_desc.op_names[i].begin(), meas_term_desc.op_names[i].end()-1, std::ostream_iterator<std::string const&>(os, ","));
        os << meas_term_desc.op_names[i].back();
        os << "]";
        if (i < meas_term_desc.op_names.size()-1)
            os << ", ";
    }
    os << "]";
    return os;
}

template <class T>
std::ostream & operator<<(std::ostream & os, std::vector<MeasurementTermDescriptor<T> > const& meas_terms_desc)
{
    os << "[\n  ";
    std::copy(meas_terms_desc.begin(), meas_terms_desc.end()-1, std::ostream_iterator<MeasurementTermDescriptor<T> const&>(os, ",\n  "));
    os << meas_terms_desc.back();
    os << "\n]\n";
    return os;
}

#endif