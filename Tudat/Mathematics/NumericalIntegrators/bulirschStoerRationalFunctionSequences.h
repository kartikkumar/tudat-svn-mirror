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
 *      120316    K. Kumar          File created
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#ifndef TUDAT_BULIRSCH_STOER_RATIONAL_FUNCTION_SEQUENCES_H
#define TUDAT_BULIRSCH_STOER_RATIONAL_FUNCTION_SEQUENCES_H

#include <vector>

namespace tudat
{
namespace mathematics
{
namespace numerical_integrators
{

//! Struct that defines the rational function sequence for a Bulirsch-Stoer integrator.
/*!
 * Struct that defines the rational function sequence for a Bulirsch-Stoer integrator.
 */
struct BulirschStoerRationalFunctionSequences
{

    //! Default constructor.
    /*!
     * Default constructor that initializes without setting the rational function sequence.
     */
    BulirschStoerRationalFunctionSequences( ) { }

    //! Constructor.
    /*!
     * Constructor that sets the rational function sequence.
     * \param sequence Rational function sequence.
     */
    BulirschStoerRationalFunctionSequences( const std::vector< unsigned int > sequence ) :
        rationalFunctionSequence( sequence ) { }

    //! Rational function sequence.
    /*!
     * Rational function sequence.
     */
    std::vector< unsigned int > rationalFunctionSequence;

    //! Enum of predefined rational function sequences.
    enum RationalFunctionSequences
    {
        bulirschStoer,
        deufelhard
    };

    //! Get rational function sequence.
    /*!
     * Returns requested rational function sequence.
     * \param rationalFunctionSequence The rational function sequence.
     * \return The requested rational function sequence.
     */
    static const BulirschStoerRationalFunctionSequences& get(
            RationalFunctionSequences sequence, const unsigned int lengthOfSequence = 100 );
};

} // namespace integrators
} // namespace mathematics
} // namespace tudat

#endif // TUDAT_BULIRSCH_STOER_RATIONAL_FUNCTION_SEQUENCES_H
