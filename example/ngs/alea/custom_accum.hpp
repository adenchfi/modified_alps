/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
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


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how one can turn a std::vector into an accumulator           *
 * by impementing all functions, that the concept requires                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef EXAMPLE_NGS_CUSTOM_ACCUM_HEADER
#define EXAMPLE_NGS_CUSTOM_ACCUM_HEADER

#include <alps/ngs.hpp>

#include <iostream>
#include <vector>
#include <boost/cstdint.hpp> //for the count-function

#include <numeric> //for accumulate

// = = = = = = = = = = C U S T O M   A C C U M = = = = = = = = = = = = =

struct myVector: public std::vector<int>
{
};

/* it's better to typedef a derived struct because 
 * typedef std::vector<int> custom_accum 
 * creates namespace problems
 */
typedef myVector custom_accum;


namespace alps
{
    namespace accumulator
    {
        // with this specialisation the detail::accumulator_wrapper knows, that custom_accum has a mean function
        template<> struct has_feature<tag::mean, custom_accum> {
            enum{value = true};
        };


        // all traits (has_X) that are not specialised are set to false 
        // unless custom_accum has a memberfunction that is called X() const

        // the value_type trait must be specialised
        template<> struct value_type<custom_accum> {
            typedef int type;
        };
        template<> struct weight_value_type<custom_accum> {
            typedef void type;
        };

        
    }//end namespace accumulator
}//end namespace alps

// this is the free mean function for the custom_accum (return value via mean trait)
alps::accumulator::mean_type<alps::accumulator::value_type<custom_accum>::type>::type mean(const custom_accum & arg)
{
    typedef alps::accumulator::mean_type<alps::accumulator::value_type<custom_accum>::type>::type mean_type;
    
    return mean_type(std::accumulate(arg.begin(), arg.end(), 0))/arg.size();
}

// the count op returns the count of the accum
boost::uint64_t count(custom_accum const & arg)
{
    return arg.size();
}

void reset(custom_accum & arg)
{
    return arg.clear();
}

// the two stream ops must also be provided. This one takes in a value of type value_type
void add_value(custom_accum & acc, int const & a)
{
    acc.push_back(a);
}

// normal ostream operator
std::ostream& operator<<(std::ostream& os,  custom_accum const & arg)
{
    os << "custom_accum(";
    os << "count = ";
    os << count(arg);
    os << " mean = ";
    os << mean(arg);
    os << ")";
    return os;
}

#endif //EXAMPLE_NGS_CUSTOM_ACCUM_HEADER
