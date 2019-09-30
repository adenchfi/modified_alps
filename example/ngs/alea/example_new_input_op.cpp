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

#include <iostream>
#include <alps/ngs.hpp>

using namespace std;
using namespace alps::accumulator;

template<typename T> 
struct has_method //checks if T has a method mean
{
    template <typename U, void (U::*)(int) > struct helper {};
    template<typename U> static char check(helper<U, &U::operator()>*);
    template<typename U> static double check(...);
    enum { value = (sizeof(char) == sizeof(check<T>(0))) };
};

int main()
{
    
    //accumulator with all the features so far
    typedef accumulator<  int
                        , features<  tag::fixed_size_binning
                                , tag::max_num_binning
                                , tag::log_binning
                                , tag::autocorrelation
                                >
                        , double
                        > accum;

    accum acc;
    detail::accumulator_wrapper w(acc);
    
    
    acc << 1;
    acc << 2;
    acc << 3;
    cout << acc << endl;
    
    acc(4);
    acc(5);
    acc(10, 100);
    acc(10, Weight = 100);
    acc(1000);
    
    w(4);
    cout << "vt :" << typeid(value_type<accum>::type).name() << endl;
    cout << "wvt:" << typeid(weight_value_type<accum>::type).name() << endl;
    cout << "test()op: " << has_method<accum>::value << endl;
    cout << "Type acc: " << typeid(accum).name() << endl;
    w(5, 1.5);
    w(10, Weight = double(100.));
    
    std::cout << mean(acc) << std::endl;
    std::cout << mean(w.get<int>()) << std::endl;
    
    cout << sizeof(weight_type<accum>::type) << endl;
    
}
