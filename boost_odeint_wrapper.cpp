//      Copyright Julio Daniel Machado Silva 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "boost_odeint_wrapper.hpp"

#include <boost/math/constants/constants.hpp>

#include <array>
#include <iostream>

/////////////////////////////////////////////////////////////////////////////////////////
// User definitions: type of state vectors and ode system to be integrated
/////////////////////////////////////////////////////////////////////////////////////////

using StateType = std::array<double, 2>;

// Got this class from boost::numeric::odeint examples
class HarmonicOscillator {
    double m_gam;

    public:
        HarmonicOscillator(double gam = 0.0) : m_gam(gam) {}

        StateType operator()(const StateType& x)
        {
            StateType dxdt;

            dxdt[0] = x[1];
            dxdt[1] = -x[0] - m_gam * x[1];

            return dxdt;
        }
};

/////////////////////////////////////////////////////////////////////////////////////////
// Main
/////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    //
    // Define variables:
    //

    // ode system
    HarmonicOscillator oscillator;

    // numerical integrator
    odeintWrapper<HarmonicOscillator, StateType> integrator(oscillator);

    StateType x0 = {1.0, 0.0};  // Initial state
    std::cout << "\nInitial state is (" << x0[0] << ", " << x0[1] << "):\n\n";

    double t0 = 0.0;                                   // Initial time
    double t1 = boost::math::constants::pi<double>();  // Final time
    double dt = 0.1;                                   // Time step

    //
    // Integrate the flux `HarmonicOscillator` with initial data `x0` from `t0` to `t1`
    // with suggested time-steps of size `dt`.
    //
    // Computed states and times are stored on the respective variables.
    //

    // C++17 style
    auto [states, times] = integrator.integrate(x0, t0, t1, dt);
    auto y = integrator.map(x0, t0, t1, dt);

    // C++14 style
    /* using odeIntegrator = odeintWrapper<HarmonicOscillator, StateType>; */
    /* typename odeIntegrator::ObservedStates states;  // std::vector of StateType */
    /* typename odeIntegrator::ObservedTimes times;    // std::vector of double */
    /* std::tie(states, times) = integrator.integrate(x0, t0, t1, dt); */

    //
    // Output states and times:
    //

    for (size_t i = 0; i < states.size(); i++) {
        std::cout << std::fixed << "times[i] = " << times[i]
            << ",\t states[i][0] = " << states[i][0]
            << ",\t states[i][1] = " << states[i][1] << '\n';
    }

    //
    // Check solution:
    //

    bool success = states.size() > 0;
    if (!success) {
        std::cerr << "Failure: output is empty.\n" << std::flush;

        return 1;
    }

    // The solution is x(t) = (cos(t), sin(t)), thus:
    // x(pi) = (-1.0, 0.0).
    double tol = 1E-4;
    std::size_t last = states.size() > 0 ? states.size() - 1 : 0;
    success &= (std::fabs(states[last][0] + 1.0) < tol);
    success &= (std::fabs(states[last][1] - 0.0) < tol);

    if (!success) {
        std::cerr << std::scientific;
        std::cerr << "\nFailure:";
        std::cerr << "\nExpected solution: x[0] = " << -1.0 << ", x[1] = " << 0.0;
        std::cerr << "\nComputed solution: x[0] = " << states[last][0]
            << ", x[1] = " << states[last][1];
        std::cerr << "\n" << std::flush;

        return 2;
    }

    // final state and mapped state `y` must also be identical
    success &= (std::fabs(y[0] + 1.0) < tol);
    success &= (std::fabs(y[1] - 0.0) < tol);

    if (!success) {
        std::cerr << std::scientific;
        std::cerr << "\nFailure:";
        std::cerr << "\nExpected solution: y[0] = " << -1.0 << ", y[1] = " << 0.0;
        std::cerr << "\nComputed solution: y[0] = " << y[0]
            << ", y[1] = " << y[1];
        std::cerr << "\n" << std::flush;

        return 3;
    }

    std::cout << "\nSuccess.\n";

    return 0;
}
