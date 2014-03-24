/*    Copyright (c) 2010-2014, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *                D. Dirkx          First version of code created.
 *      140324    K. Kumar          Cleaned up code.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_LAGRANGE_INTERPOLATOR_H
#define TUDAT_LAGRANGE_INTERPOLATOR_H

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

// #include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"

namespace tudat
{
namespace interpolators
{

//! Get value one.
template< typename ScalarType >
ScalarType getOne( )
{
  return static_cast< ScalarType >( 1.0 );
}

//! Class to perform Lagrange polynomial interpolation.
/*!
 * Class to perform Lagrange polynomial interpolation from a set of independent and dependent 
 * values, as well as the order of the interpolation. Note that this class is optimized for many 
 * function calls to interpolate, since the denominators for the interpolations are pre-computed 
 * for all interpolation intervals.
 */
template< typename IndependentVariableType, typename DependentVariableType, 
          typename ScalarType = IndependentVariableType >
class LagrangeInterpolator 
    : public OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >
{
public:

    //! Using statements to prevent having to put 'this' everywhere in the code.
    using OneDimensionalInterpolator< 
        IndependentVariableType, DependentVariableType >::dependentValues_;
    using OneDimensionalInterpolator< 
        IndependentVariableType, DependentVariableType >::independentValues_;
    using OneDimensionalInterpolator< 
        IndependentVariableType, DependentVariableType >::lookUpScheme_;

    //! Constructor from vectors of independent/dependent data.
    /*!
     * This constructor initializes the interpolator from two vectors containing the independent
     * variables and dependent variables. A look-up scheme can be provided to
     * override the given default.
     * \param independentVariables Vector of values of independent variables that are used.
     * \param dependentVariables Vector of values of dependent variables that are used.
     * \param numberOfStages Number of data points that are used to calculate the interpolating
     *          polynomial (must be even).
     * \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *          to find the nearest lower data point in the independent variables when requesting
     *          interpolation.
     */
    LagrangeInterpolator( const std::vector< IndependentVariableType >& independentVariables,
                          const std::vector< DependentVariableType >& dependentVariables,
                          const int numberOfStages,
                          const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm )
        : numberOfStages_( numberOfStages )
    {
        independentValues_ = independentVariables;
        dependentValues_ = dependentVariables;          
        numberOfIndependentValues_ = static_cast< int >( independentValues_.size( ) );

        // Verify that the initialization variables are not empty.
        if ( numberOfIndependentValues_ == 0 || dependentValues_.size( ) == 0 )
        {
            throw std::runtime_error(
                "ERROR: The vectors used in the Lagrange interpolator initialization are empty!" );
        }

        // Check consistency of input data.
        if ( static_cast< int >( dependentValues_.size( ) ) != numberOfIndependentValues_ )
        {
            throw std::runtime_error( 
                "ERROR: independent and dependent variables not of same size in Lagrange \
                    interpolator!" );
        }

        // Define zero entry for dependent variable.
        zeroEntry_ = dependentVariables[ 0 ] - dependentVariables[ 0 ];

        // Create lookup scheme from independent variable values.
        this->makeLookupScheme( selectedLookupScheme );

        // Calculate denominators for each interval, to prevent recalculations during each 
        // interpolation call.
        initializeDenominators( );
        initializeBoundaryInterpolators( );

        independentVariableDifferenceCache.resize( 2 * offsetEntries_ + 2 );
    }

    //! Constructor from map of independent/dependent data.
    /*!
     * This constructor initializes the interpolator from a map containing independent variables
     * as key and dependent variables as value. A look-up scheme can be provided to override the
     * given default.
     * \param dataMap Map containing independent variables as key and dependent variables as value.
     * \param numberOfStages Number of data points that are used to calculate the interpolating
     *          polynomial (must be even).
     * \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *          to find the nearest lower data point in the independent variables when requesting
     *          interpolation.
     */
    LagrangeInterpolator( 
        const std::map< IndependentVariableType, DependentVariableType >& dataMap,
        const int numberOfStages,
        const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm )
        : numberOfStages_( numberOfStages )
    {
        numberOfIndependentValues_ = dataMap.size( );

        // Verify that the initialization variables are not empty.
        if ( dataMap.size( ) == 0 )
        {
            throw std::runtime_error(
                "ERROR: The vectors used in the Lagrange interpolator initialization are empty!" );
        }

        // Resize data vectors of independent/dependent values.
        independentValues_.resize( dataMap.size( ) );
        dependentValues_.resize( dataMap.size( ) );

        // Fill data vectors with data from map.
        int counter = 0;
        for ( typename std::map< IndependentVariableType, 
                                 DependentVariableType >::const_iterator mapIterator 
                                    = dataMap.begin( );
             mapIterator != dataMap.end( ); mapIterator++ )
        {
            independentValues_[ counter ] = mapIterator->first;
            dependentValues_[ counter ] = mapIterator->second;
            counter++;
        }

        // Define zero entry for dependent variable.
        zeroEntry_ = dependentValues_[ 0 ] - dependentValues_[ 0 ];

        // Create lookup scheme from independent variable data points.
        this->makeLookupScheme( selectedLookupScheme );

        // Calculate denominators for each interval, to prevent recalculations during each 
        // interpolation call.
        initializeDenominators( );
        initializeBoundaryInterpolators( );

        independentVariableDifferenceCache.resize( 2 * offsetEntries_ + 2 );
    }

    // Using statement to prevent compiler warning.
    using Interpolator< IndependentVariableType, DependentVariableType, 1 >::interpolate;

    //! Interpolate dependent variable value at given independent variable value.
    /*!
     * Interpolates dependent variable value at given independent variable value.
     * The polynomial centered on the requested interval is used for the interpolation. Note that 
     * the number of data points, not the values of the independent variables are used for 
     * determining the center (in case of a non-equispaced grid). If the required interpolating 
     * polynimial goes beyond the independent variable bondaries, a cubic spline with natural 
     * boundary conditions is used.
     * \param targetIndependentVariableValue Value of independent variable at which interpolation 
     * is to take place.
     * \return Interpolated value of dependent variable.
     */
    DependentVariableType interpolate( 
        const IndependentVariableType targetIndependentVariableValue )
    {
        using std::pow;

        // Determine the lower entry in the table corresponding to the target independent variable
        // value.
        DependentVariableType interpolatedValue = zeroEntry_;

        // Find interpolation interval.
        int lowerEntry = lookUpScheme_->findNearestLowerNeighbour( 
            targetIndependentVariableValue );

        // Check if requested interval is inside region in which centered lagrange interpolation 
        // can be used.
        if ( lowerEntry < offsetEntries_ )
        {
            std::cerr << "Below, time diff " << 
                         targetIndependentVariableValue - independentValues_[ 0 ] << std::endl;

            // Use start cubic spline interpolator.
            interpolatedValue = beginInterpolator_->interpolate( targetIndependentVariableValue );
        }

        else if ( lowerEntry >= numberOfIndependentValues_ - offsetEntries_ - 1 )
        {
            // Use end cubic spline interpolator.
            std::cout << "Above, time diff " 
                      << targetIndependentVariableValue 
                         - independentValues_[ independentValues_.size( ) - 1 ] << std::endl;

            interpolatedValue = endInterpolator_->interpolate( targetIndependentVariableValue );
        }

        else
        {
            // Initialize repeated numerator to 1.
            ScalarType repeatedNumerator = getOne< ScalarType >( );

            // Check if requested independent variable is equal to data point.
            if ( independentValues_[ lowerEntry ] == targetIndependentVariableValue )
            {
                interpolatedValue = dependentValues_[ lowerEntry ];
            }

            else if ( independentValues_[ lowerEntry + 1 ] == targetIndependentVariableValue )
            {
                interpolatedValue = dependentValues_[ lowerEntry + 1 ];
            }

            else if ( independentValues_[ lowerEntry - 1 ] == targetIndependentVariableValue )
            {
                interpolatedValue = dependentValues_[ lowerEntry - 1 ];
            }

            else
            {
                // Set up repeated numerator and cache of independent variable values from which
                // interpolant is created.
                int counter = 0;
                for( int i = lowerEntry - offsetEntries_; 
                     i <= lowerEntry + offsetEntries_ + 1; 
                     i++ )
                {
                    independentVariableDifferenceCache[ counter ] 
                        = static_cast< ScalarType >( 
                            targetIndependentVariableValue - independentValues_[ i ] );
                    repeatedNumerator *= independentVariableDifferenceCache[ counter ];
                    counter++;
                }

                counter = 0;

                // Evaluate interpolating polynomial at requested data point.
                for( int i = lowerEntry - offsetEntries_; 
                     i <= lowerEntry + offsetEntries_ + 1; 
                     i++ )
                {
                    interpolatedValue += dependentValues_[ i ]  
                        * ( repeatedNumerator / 
                            ( independentVariableDifferenceCache[ counter ] 
                              * denominators[ lowerEntry ][ i - lowerEntry + offsetEntries_ ] ) );
                    counter++;
                }
            }
        }

        return interpolatedValue;
    }

protected:

private:

    //! Pre-compute the denominators of the interpolants at each interval.
    /*!
     * Pre-computes at initialization the denominators of the interpolants at each interval, 
     * i.e. each interval between two subsequent independent variable values.
     */
    void initializeDenominators( )
    {
        // Check validity of requested number of stages".
        if ( numberOfStages_% 2 != 0 )
        {
            throw std::runtime_error(
             "ERROR: Lagrange interplator currently only implemented for even number of stages!" );
        }

        if ( numberOfStages_ < 2 )
        {
            throw std::runtime_error( 
                "ERROR: Lagrange interplator number of stages must be greater than 2!" );
        }

        // Determine offset from boundary of interpolation interval where interpolant is valid.
        offsetEntries_ = ( numberOfStages_ / 2 ) - 1;

        // Iterate over all intervals and calculate denominators.
        int currentIterationStart;
        denominators.resize( numberOfIndependentValues_ );
        for ( int i = offsetEntries_; i <= numberOfIndependentValues_ - offsetEntries_; i++ )
        {
            // Determine start index in independent variables for current polynomial.
            currentIterationStart = i - offsetEntries_;

            denominators[ i ].resize( 2 * offsetEntries_ + 2 );

            // Calculate all denominators for single interval.
            for ( int j = 0; j <= 2 * offsetEntries_ + 1; j++ )
            {
                denominators[ i ][ j ] = getOne< ScalarType >( );
                for ( int k = 0; k <= 2 * offsetEntries_ + 1; k++ )
                {
                    if ( k != j )
                    {
                        denominators[ i ][ j ] 
                            *= static_cast< ScalarType >( 
                                independentValues_[ j + currentIterationStart ] 
                                - independentValues_[ k + currentIterationStart ] );
                    }
                }
            }
        }
    }

    //! Create interpolators used at the boundaries of the interpolation domain.
    /*!
     * At initialization, creates the interpolators used at the boundaries of the interpolation 
     * domain. Cubic spline interpolators with natural boundary conditions are used.
     */
    void initializeBoundaryInterpolators( )
    {
        // Set input maps for interpolators.
        std::map< IndependentVariableType, DependentVariableType > startMap;
        for ( int i = 0; i <= offsetEntries_; i++ )
        {
            startMap[ independentValues_[ i ] ] = dependentValues_[ i ];
        }

        std::map< IndependentVariableType, DependentVariableType > endMap;
        for ( int i = numberOfIndependentValues_ - offsetEntries_ - 1; 
              i < numberOfIndependentValues_; 
              i++ )
        {
            endMap[ independentValues_[ i ] ] = dependentValues_[ i ];
        }

        // Create interpolators.
        beginInterpolator_ 
            = boost::make_shared< CubicSplineInterpolator< 
                IndependentVariableType, DependentVariableType > >( startMap );

        endInterpolator_ 
            = boost::make_shared< CubicSplineInterpolator< 
                IndependentVariableType, DependentVariableType > >( endMap );
    }

    std::vector< std::vector< ScalarType > > denominators;

    //! Zero entry for dependent variables.
    /*!
     * Zero entry for dependent variables, i.e. algebraic identity element for addition of 
     * dependent variables.
     */
    DependentVariableType zeroEntry_;

    //! Number of stages of interpolator.
    /*!
     * Number of stages of interpolator, i.e. order of interpolating polynomial/number of 
     * intervals used to create interpolant.
     */
    int numberOfStages_;

    int offsetEntries_;

    std::vector< ScalarType > independentVariableDifferenceCache;

    boost::shared_ptr< OneDimensionalInterpolator< 
        IndependentVariableType, DependentVariableType > > beginInterpolator_;

    boost::shared_ptr< OneDimensionalInterpolator< 
        IndependentVariableType, DependentVariableType > > endInterpolator_;

    int numberOfIndependentValues_;
};

typedef LagrangeInterpolator< double, double > LagrangeInterpolatorDouble;

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_LAGRANGE_INTERPOLATOR_H
