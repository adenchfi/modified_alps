/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
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

#include "ising.hpp"

#include <boost/lambda/lambda.hpp>

ising_sim::ising_sim(parameters_type const & params)
    : parameters(params)
    , random(boost::mt19937((parameters["SEED"] | 42)), boost::uniform_real<>())
    , length(parameters["L"])
    , sweeps(0)
    , thermalization_sweeps(int(parameters["THERMALIZATION"]))
    , total_sweeps(int(parameters["SWEEPS"]))
    , beta(1. / double(parameters["T"]))
    , spins(length)
{
    for(int i = 0; i < length; ++i)
        spins[i] = (random() < 0.5 ? 1 : -1);
    measurements
        << alps::accumulator::RealObservable("Energy")
        << alps::accumulator::RealObservable("Magnetization")
        << alps::accumulator::RealObservable("Magnetization^2")
        << alps::accumulator::RealObservable("Magnetization^4")
        << alps::accumulator::RealVectorObservable("Correlations")
    ;
}

void ising_sim::update() {
    for (int j = 0; j < length; ++j) {
        using std::exp;
        int i = int(double(length) * random());
        int right = ( i + 1 < length ? i + 1 : 0 );
        int left = ( i - 1 < 0 ? length - 1 : i - 1 );
        double p = exp( 2. * beta * spins[i] * ( spins[right] + spins[left] ));
        if ( p >= 1. || random() < p )
            spins[i] = -spins[i];
    }
}

void ising_sim::measure() {
    sweeps++;
    if (sweeps > thermalization_sweeps) {
        double tmag = 0;
        double ten = 0;
        double sign = 1;
        std::vector<double> corr(length);
        for (int i = 0; i < length; ++i) {
            tmag += spins[i];
            sign *= spins[i];
            ten += -spins[i] * spins[ i + 1 < length ? i + 1 : 0 ];
            for (int d = 0; d < length; ++d)
                corr[d] += spins[i] * spins[( i + d ) % length ];
        }
        std::transform(corr.begin(), corr.end(), corr.begin(), boost::lambda::_1 / double(length));
        ten /= length;
        tmag /= length;
        measurements["Energy"] << ten;
        measurements["Magnetization"] << tmag;
        measurements["Magnetization^2"] << tmag * tmag;
        measurements["Magnetization^4"] << tmag * tmag * tmag * tmag;
        measurements["Correlations"] << corr;
    }
}

double ising_sim::fraction_completed() const {
    return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}

bool ising_sim::run(boost::function<bool ()> const & stop_callback) {
    bool stopped = false;
    do {
        update();
        measure();
    } while(!(stopped = stop_callback()) && fraction_completed() < 1.);
    return !stopped;
}

// implement a nice keys(m) function
ising_sim::result_names_type ising_sim::result_names() const {
    result_names_type names;
    for(accumulators_type::const_iterator it = measurements.begin(); it != measurements.end(); ++it)
        names.push_back(it->first);
    return names;
}

ising_sim::result_names_type ising_sim::unsaved_result_names() const {
    return result_names_type(); 
}

ising_sim::results_type ising_sim::collect_results() const {
    return collect_results(result_names());
}

ising_sim::results_type ising_sim::collect_results(result_names_type const & names) const {
    results_type partial_results;
    for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it)
        partial_results[*it] = measurements[*it].result());
    return partial_results;
}

void ising_sim::save(boost::filesystem::path const & filename) const {
    alps::hdf5::archive ar(filename, "w");
    ar["/"] << *this;
}

void ising_sim::load(boost::filesystem::path const & filename) {
    alps::hdf5::archive ar(filename);
    ar["/"] >> *this;
}

void ising_sim::save(alps::hdf5::archive & ar) const {
    std::string context = ar.get_context();

    ar["/parameters"] << parameters;

    ar.set_context("/simulation/realizations/0/clones/0");
    ar["measurements"] << measurements;

    ar.set_context("checkpoint");
    ar["sweeps"] << sweeps;
    ar["spins"] << spins;

    {
        std::ostringstream os;
        os << random.engine();
        ar["engine"] << os.str();
    }

    ar.set_context(context);
}

void ising_sim::load(alps::hdf5::archive & ar) {
    std::string context = ar.get_context();

    ar["/parameters"] >> parameters;
    length = int(parameters["L"]);
    thermalization_sweeps = int(parameters["THERMALIZATION"]);
    total_sweeps = int(parameters["SWEEPS"]);
    beta = 1. / double(parameters["T"]);

    ar.set_context("/simulation/realizations/0/clones/0");
    ar["measurements"] >> measurements;

    ar.set_context("checkpoint");
    ar["sweeps"] >> sweeps;
    ar["spins"] >> spins;

    {
        std::string state;
        ar["engine"] >> state;
        std::istringstream is(state);
        is >> random.engine();
    }

    ar.set_context(context);
}
