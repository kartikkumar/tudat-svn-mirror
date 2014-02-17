/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120203    B. Tong Minh      Copied RungeKutta4Stepsize unit test.
 *      120207    K. Kumar          Adapted to use modified benchmark functions in Tudat Core.
 *      120213    K. Kumar          Modified getCurrentInterval() to getIndependentVariable();
 *                                  transferred to Boost unit test framework.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>

#include <limits>

#include <boost/test/unit_test.hpp>

#include <TudatCore/Mathematics/NumericalIntegrators/UnitTests/benchmarkFunctions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/NumericalIntegrators/bulirschStoerVariableStepsizeIntegrator.h"
#include "Tudat/Mathematics/NumericalIntegrators/bulirschStoerRationalFunctionSequences.h"

namespace tudat
{
namespace unit_tests
{

const double TUDAT_MACHINE_PRECISION = std::numeric_limits< double >::epsilon( );

BOOST_AUTO_TEST_SUITE( test_bulirsch_stoer_variable_step_size_integrator )

using mathematics::numerical_integrators::BulirschStoerVariableStepSizeIntegratorXd;
using mathematics::numerical_integrators::BulirschStoerRationalFunctionSequences;

//! Test the validity of the BurlischStoerVariableStepSize integrator.
/*!
 * Tests the validity of the RungeKuttaVariableStepSize integrator.
 * \param sequence The rational function sequence used for the Bulirsch-Stoer integrator being
 *          tested.
 * \param stateDerivativeFunction Function pointer to the state derivative function.
 * \param intervalStart The start of the integration interval.
 * \param intervalEnd The end of the integration interval.
 * \param stepSize The step size to take.
 * \param initialState The initial state.
 * \param expectedState Expected final state.
 * \param tolerance Tolerance when comparing.
 * \return True if actual final state equals the expected final state within the specified
 *          tolerance.
 */
bool testValidityOfBulirschStoerVariableStepSizeIntegrator(
        const BulirschStoerRationalFunctionSequences sequence,
        const BulirschStoerVariableStepSizeIntegratorXd::StateDerivativeFunction&
        stateDerivativeFunction,
        const double intervalStart, const double intervalEnd, const double stepSize,
        const Eigen::VectorXd& initialState, const Eigen::VectorXd expectedState,
        const double tolerance )
{
    // Create forward BulirschStoerVariableStepSize integrator.
    {
        BulirschStoerVariableStepSizeIntegratorXd integrator( sequence, stateDerivativeFunction,
                                                              intervalStart, initialState,
                                                              1.0e-15, tolerance / 10.0 );

        Eigen::VectorXd finalState = integrator.integrateTo( intervalEnd, stepSize );

        // Compute differences between computed and expected interval end and generate
        // cerr statement if test fails.
        if ( std::fabs( integrator.getCurrentIndependentVariable( ) - intervalEnd ) / intervalEnd >
             10.0 * TUDAT_MACHINE_PRECISION )
        {
            return false;
        }

        // Compute differences between computed and expected results and generate
        // cerr statement if test fails.
        if ( !expectedState.isApprox( finalState, tolerance ) )
        {
            return false;
        }
    }

    // Try the same again, but in two steps.
    {
        BulirschStoerVariableStepSizeIntegratorXd integrator( sequence, stateDerivativeFunction,
                                                              intervalStart, initialState,
                                                              1.0e-15, tolerance / 10.0 );

        const double intermediateIndependentVariable
                = intervalStart + ( intervalEnd - intervalStart ) / 2.0;

        const Eigen::VectorXd intermediateState = integrator.integrateTo(
                    intermediateIndependentVariable, stepSize );

        // Compute differences between computed and expected interval end.
        if ( std::fabs( integrator.getCurrentIndependentVariable( )
                       - intermediateIndependentVariable ) /
             intermediateIndependentVariable > TUDAT_MACHINE_PRECISION )
        {
            return false;
        }

        // Integrate to the end
        Eigen::VectorXd finalState = integrator.integrateTo( intervalEnd, stepSize );

        // Compute differences between computed and expected interval end.
        if ( std::fabs( integrator.getCurrentIndependentVariable( ) - intervalEnd ) / intervalEnd >
             10.0 * TUDAT_MACHINE_PRECISION )
        {
            return false;
        }

        // Compute differences between computed and expected results and generate
        // cerr statement if test fails.
        if ( !expectedState.isApprox( finalState, tolerance ) )
        {
            return false;
        }

        integrator.performIntegrationStep( stepSize );
        if ( !integrator.rollbackToPreviousState( ) )
        {
            return false;
        }

        // No need to check machine precision, because this interval is stored exact.
        if ( std::fabs( integrator.getCurrentIndependentVariable( ) - intervalEnd ) / intervalEnd >
             10.0 * TUDAT_MACHINE_PRECISION )
        {
            return false;
        }

        // This result should be exactly the same.
        if ( integrator.getCurrentState( ) != finalState )
        {
            return false;
        }

        if ( integrator.rollbackToPreviousState( ) )
        {
            return false;
        }
    }

    return true;
}

//! Test different types of states and state derivatives.
/*!
 * Tests if different types of states and state derivatives work. If something
 * is broken, then a compile time error will be generated.
 * \return Unconditionally true
 */
bool testDifferentStateAndStateDerivativeTypes( )
{
    using tudat::unit_tests::computeZeroStateDerivative;
    tudat::mathematics::numerical_integrators::BulirschStoerVariableStepSizeIntegrator
            < double, Eigen::Vector3d, Eigen::VectorXd > integrator(
                BulirschStoerRationalFunctionSequences( ),  &computeZeroStateDerivative,
                0.0, Eigen::Vector3d::Zero( ) );

    integrator.integrateTo( 0.0, 0.1 );

    // No need to test anything, this is just to check compile time errors.
    return true;
}

//! Test the validity of the BulirschStoerVariableStepSize integrator.
/*!
 * Tests the validity of the BulirschStoerVariableStepSize integrator.
 * \param sequence Bulirsch-Stoer, variable-step size rational function sequence.
 * \return True if actual final state equals the expected final state within the specified
 *          tolerance.
 */
bool testValidityOfBulirschStoerVariableStepSizeIntegrator(
        const BulirschStoerRationalFunctionSequences& sequence )
{
    // Test result initialised to false.
    bool testBulirschStoerVariableStepSizeIsOk = true;

    using namespace tudat::unit_tests;
    std::map< BenchmarkFunctions, BenchmarkFunction >& benchmarkFunctions =
            getBenchmarkFunctions( );

    // Case 1: test with x_dot = 0, which results in x_f = x_0.
    {
        testBulirschStoerVariableStepSizeIsOk
                &= testValidityOfBulirschStoerVariableStepSizeIntegrator(
                    sequence,
                    benchmarkFunctions[ Zero ].pointerToStateDerivativeFunction_,
                    benchmarkFunctions[ Zero ].intervalStart_,
                    benchmarkFunctions[ Zero ].intervalEnd_, 0.2,
                    benchmarkFunctions[ Zero ].initialState_,
                    benchmarkFunctions[ Zero ].finalState_,
                    TUDAT_MACHINE_PRECISION );
    }

    // Case 2: test with x_dot = 1, which results in x_f = x_0 + t_f.
    {
        testBulirschStoerVariableStepSizeIsOk
                &= testValidityOfBulirschStoerVariableStepSizeIntegrator(
                    sequence,
                    benchmarkFunctions[ Constant ].pointerToStateDerivativeFunction_,
                    benchmarkFunctions[ Constant ].intervalStart_,
                    benchmarkFunctions[ Constant ].intervalEnd_, 0.2,
                    benchmarkFunctions[ Constant ].initialState_,
                    benchmarkFunctions[ Constant ].finalState_, 1.0e-14 );
    }

    // Case 3: test with x_dot = x, which results in x_f = x0 * exp( t_f )
    {
        testBulirschStoerVariableStepSizeIsOk
                &= testValidityOfBulirschStoerVariableStepSizeIntegrator(
                    sequence,
                    benchmarkFunctions[ Exponential ].pointerToStateDerivativeFunction_,
                    benchmarkFunctions[ Exponential ].intervalStart_,
                    benchmarkFunctions[ Exponential ].intervalEnd_, 1.0,
                    benchmarkFunctions[ Exponential ].initialState_,
                    benchmarkFunctions[ Exponential ].finalState_, 1.0e-5 );
    }

    // Case 4: test with an example from Burden and Faires.
    {
        testBulirschStoerVariableStepSizeIsOk
                &= testValidityOfBulirschStoerVariableStepSizeIntegrator(
                    sequence,
                    benchmarkFunctions[ BurdenAndFaires ]
                    .pointerToStateDerivativeFunction_,
                    benchmarkFunctions[ BurdenAndFaires ].intervalStart_,
                    benchmarkFunctions[ BurdenAndFaires ].intervalEnd_, 0.25,
                    benchmarkFunctions[ BurdenAndFaires ].initialState_,
                    benchmarkFunctions[ BurdenAndFaires ].finalState_, 1.0e-11 );
    }

    return testBulirschStoerVariableStepSizeIsOk;
}

BOOST_AUTO_TEST_CASE( testBulirschStoerVariableStepSizeIntegrator )
{
    // Case 1: test if difference in type between state and state derivative works.
    BOOST_CHECK( testDifferentStateAndStateDerivativeTypes( ) );

    // Case 2: test if Bulirsch-Stoer sequence works.
    BOOST_CHECK( testValidityOfBulirschStoerVariableStepSizeIntegrator(
                     BulirschStoerRationalFunctionSequences::get(
                         BulirschStoerRationalFunctionSequences::bulirschStoer ) ) );

//    // Case 2: test if Deufelhard sequence works.
//    BOOST_CHECK( testValidityOfBulirschStoerVariableStepSizeIntegrator(
//                     BulirschStoerRationalFunctionSequences::get(
//                         BulirschStoerRationalFunctionSequences::deufelhard ) ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
