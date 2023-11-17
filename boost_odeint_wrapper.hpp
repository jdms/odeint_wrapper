//      Copyright Julio Daniel Machado Silva 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#ifndef SIMPLE_BOOST_ODEINT_WRAPPER
#define SIMPLE_BOOST_ODEINT_WRAPPER

#include <functional>
#include <tuple>
#include <vector>

#include <boost/numeric/odeint.hpp>

/// This class wraps the default integrator of boost::numeric::odeint, currently
/// the Dormand-Prince scheme RK5(4) to integrate a given flux
///
/// It expects a `State` type for the physical state and a `System` type for the
/// vector-field, such that, for a `State x` and a `System F`, d/dt x(t) = F(x(t))
///
/// Reference:
/// Dormand, J. R.; Prince, P. J. (1980), "A family of embedded Runge-Kutta formulae",
/// Journal of Computational and Applied Mathematics", 6 (1): 19–26
///
template <typename System, typename State>
class odeintWrapper {
    public:
        /// Type alias to hold a vector of states
        using ObservedStates = std::vector<State>;

        /// Type alias to hold a vector of times
        using ObservedTimes = std::vector<double>;

        /// Initialize the flux to be integrated with `system`
        odeintWrapper(System& system) : flux_(AutonomousFlux(system)) {}

        /// Initialize the flux to be integrated with `system`
        odeintWrapper(System&& system) : flux_(AutonomousFlux(system)) {}

        /// Default dtor
        ~odeintWrapper() = default;

        /// Map initial state `x` at `time0` to a new state `y` at `time1` with the flux;
        /// returns: (return value) mapped state `y`; and (forward reference) number of steps `num_steps`.
        template<typename T = std::size_t>
        State map(const State& x /**< Initial state/Mapped state */,
                double time0 = 0.0 /**< Initial time */,
                double time1 = 1.0 /**< Final time */,
                double dt = 1.0 /**< Initial discretization (scheme is adaptive) */,
                T&& num_steps = 0) const
        {
            State y{x}; // mapped state

            AutonomousFlux f = const_cast<odeintWrapper*>(this)->flux_;
            if (time0 > time1) {
                num_steps = boost::numeric::odeint::integrate(-f, y, time1, time0, dt);
            }
            else {
                num_steps = boost::numeric::odeint::integrate(f, y, time0, time1, dt);
            }

            return y;
        }

        /// Integrate given flux from `state0` at `time0` to `time1` with `dt` increments;
        /// returns: tuple `std::tie(states, times)` computed during integration
        std::tuple<ObservedStates, ObservedTimes> integrate(
                State state0 /**< Initial state */,
                double time0 /**< Initial time */,
                double time1 /**< Final time */,
                double dt = 1.0 /**< Initial discretization (scheme is adaptive) */) const
        {
            ObservedStates states;
            ObservedTimes times;

            AutonomousFlux f = const_cast<odeintWrapper*>(this)->flux_;
            if (time0 > time1) {
                boost::numeric::odeint::integrate(
                        -f, state0, time1, time0, dt, Observer(states, times));
            }
            else {
                boost::numeric::odeint::integrate(
                        f, state0, time0, time1, dt, Observer(states, times));
            }

            return std::tie(states, times);
        }

    private:
        /// Helper class to hold a reference to vector-field, satisfies odeint's usage
        class AutonomousFlux {
            public:
                /// Initialize the flux to be integrated with `system`
                AutonomousFlux(System& system) : system_(system) {}

                /// Initialize the flux to be integrated with `system`
                AutonomousFlux(System&& system) : system_(system) {}

                /// Default dtor
                ~AutonomousFlux() = default;

                /// Default copy ctor
                AutonomousFlux(AutonomousFlux&) = default;

                /// Default copy assignment
                AutonomousFlux& operator=(AutonomousFlux&) = default;

                /// Default move ctor
                AutonomousFlux(AutonomousFlux&&) = default;

                /// Default move assignment
                AutonomousFlux& operator=(AutonomousFlux&&) = default;

                /// Evaluate vector field, `odeint` expects this signature for its integrators
                void operator()(const State& p /**< Current state */,
                        State& Fp /**< Derivative F(p) */,
                        double /**< Time (not used) */) const
                {
                    Fp = system_(p);

                    if (invert_flux_) {
                        for (auto& Fp_i : Fp) {
                            Fp_i = -Fp_i;
                        }
                    }
                }

                /// Invert direction of vector field
                AutonomousFlux operator-() const
                {
                    AutonomousFlux minus_flux = *const_cast<AutonomousFlux*>(this);
                    minus_flux.invertFlux();

                    return minus_flux;
                }

            private:
                std::reference_wrapper<System> system_;  ///< wraps the given ode system
                bool invert_flux_ = false;               ///< operator() evaluates -system

                /// Invert direction of vector field
                void invertFlux() { invert_flux_ = !invert_flux_; }

        } flux_;  ///< Flux to be integrated */

        /// Helper class used to store states and times computed during integration
        class Observer {
            public:
                /// Ctor stores references to states and times computed during integration
                Observer(ObservedStates& states, ObservedTimes& times)
                    : states_(states), times_(times){};

                /// Default dtor
                ~Observer() = default;

                /// Store new pair (state, time)
                void operator()(const State& s /**< state to be stored */,
                        double t /**< time to be stored */)
                {
                    states_.push_back(s);
                    times_.push_back(t);
                }

            private:
                ObservedStates& states_;  ///< Vector of observed states
                ObservedTimes& times_;    ///< Vector of observed times
        };
};

#endif /* ifndef SIMPLE_BOOST_ODEINT_WRAPPER */
