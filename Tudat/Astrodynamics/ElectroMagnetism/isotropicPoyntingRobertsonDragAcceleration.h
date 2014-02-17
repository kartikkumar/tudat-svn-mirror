/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      131016    S. Hirsh          File Created.
 *
 *    References
 *      Eigen. Structures having Eigen members,
 *          http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html, last accessed: 5th
 *          March, 2013.
 *
 *    Notes
 *
 */

#ifndef TUDAT_ISOTROPIC_POYNTING_ROBERTSON_DRAG_ACCELERATION_H
#define TUDAT_ISOTROPIC_POYNTING_ROBERTSON_DRAG_ACCELERATION_H 

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureForce.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute isotropic Poynting-Robertson drag acceleration.
/*!
 * Computes Poynting-Robertson drag acceleration on a particle. Assumes the particle radiates 
 * radiates isotropically in the particle's local frame. 
 * \param velocityOfAcceleratedBody Velocity vector of the body which is accelerated by the 
 *          Poynting-Robertson drag. This vector must be with respect to the inertial 
 *          (barycentric) frame.                                                              [m/s]                                                                                              
 * \param vectorFromSource Vector pointing from source to target. N.B: this must be a unit
 *          vector! To compute the unit vector based on a given position vector, you can
 *          use the .normalize() or .normalized() member functions of an Eigen::Vector3d
 *          object.                                                                             [m]
 * \param radiationPressureForce Force due to radiation pressure from source. The function 
 *          computeCannonBallRadiationPressureForce() can be used to compute this value for
 *          a cannon-ball model. Note: Input is the magnitude of the force, not the vector.
 *          To compute the unit vector based on a given position vector, you can use the 
 *          .norm() member function of an Eigen::Vector3d object.                               [N]
 * \param mass Mass of accelerated body.                                                       [kg]                          
 * \return Acceleration due to Poynting-Robertson drag.                                     [m/s^2]
 * \sa computeCannonBallRadiationPressureForce().
 */
Eigen::Vector3d computeIsotropicPoyntingRobertsonDragAcceleration(
        const Eigen::Vector3d& velocityOfAcceleratedBody,
        const Eigen::Vector3d& vectorFromSource,
        const double radiationPressureForce,
        const double mass );

//! Isotropic Poynting-Robertson drag acceleration model class.
/*!
 * Class that can be used to compute the isotropic Poynting-Robertson drag. The particle subject to
 * the Poynting-Robertson drag is modelled as a cannon-ball.
 */
class IsotropicPoyntingRobertsonDrag : public basic_astrodynamics::AccelerationModel3d
{
private:

    //! Typedef for double-returning function.
    typedef boost::function< double( ) > DoubleReturningFunction;

    //! Typedef for Eigen::Vector3d-returning function.
    typedef boost::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor taking function pointers for all variables.
    /*!
     * Constructor taking function pointers for all variables.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing
     *          Poynting-Robertson drag acceleration.
     * \param acceleratedBodyVelocityFunction Function returning velocity of body undergoing
     *          Poynting-Robertson drag acceleration.     
     * \param radiationPressureForceFunction Function returning current radiation pressure force.
     * \param massFunction Function returning current mass of body undergoing acceleration.
     */
    IsotropicPoyntingRobertsonDrag(
            Vector3dReturningFunction sourcePositionFunction,
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            Vector3dReturningFunction acceleratedBodyVelocityFunction,            
            DoubleReturningFunction radiationPressureForceFunction,
            DoubleReturningFunction massFunction )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          acceleratedBodyVelocityFunction_( acceleratedBodyVelocityFunction ),          
          radiationPressureForceFunction_( radiationPressureForceFunction ),
          massFunction_( massFunction )
    {
        this->updateMembers( );
    }

    //! Constructor taking functions pointers and constant values for parameters.
    /*!
     * Constructor taking function pointers for position and velocity vectors and
     * constant values for all other parameters.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing
     *          Poynting-Robertson drag acceleration.
     * \param acceleratedBodyVelocityFunction Function returning velocity of body undergoing
     *          Poynting-Robertson drag acceleration.     
     * \param radiationPressureForce Constant radiation pressure.                               [N]
     * \param mass Constant mass of body undergoing acceleration.                              [kg]
     */
    IsotropicPoyntingRobertsonDrag(
            Vector3dReturningFunction sourcePositionFunction,
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            Vector3dReturningFunction acceleratedBodyVelocityFunction, 
            const double radiationPressureForce,
            const double mass )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          acceleratedBodyVelocityFunction_( acceleratedBodyVelocityFunction ),  
          radiationPressureForceFunction_( 
            boost::lambda::constant ( radiationPressureForce ) ),
          massFunction_( boost::lambda::constant( mass ) )
    {
        this->updateMembers( );
    }

    //! Get isotropic Poynting-Robertson drag acceleration.
    /*!
     * Returns the isotropic Poynting-Robertson drag acceleration. No arguments are passed to this 
     * function. Instead, all data required for computation is to be obtained from pointers to functions/
     * classes/structs, etc., which are to be set in a derived class and evaluated by the
     * updateMembers() function below. This function is essentially a wrapper for the free
     * function that computes the Poynting-Roberston drag acceleration.
     * \return Isotropic Poynting-Robertson drag acceleration.
     * \sa computeIsotropicPoyntingRobertsonDragAcceleration().
     */
    Eigen::Vector3d getAcceleration( );

    //! Update member variables used by the isotropic Poynting-Robertson drag acceleration model.
    /*!
     * Updates member variables used by the acceleration model. This function evaluates all
     * dependent variables to the 'current' values of these parameters. Only these current values,
     * not the function-pointers are then used by the getAcceleration() function.
     */
    void updateMembers( );

private:

    //! Function pointer returning position of source.
    /*!
     * Function pointer returning position of source (3D vector).
     */
    const Vector3dReturningFunction sourcePositionFunction_;

    //! Function pointer returning position of accelerated body.
    /*!
     * Function pointer returning position of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyPositionFunction_;

    //! Function pointer returning velocity of accelerated body.
    /*!
     * Function pointer returning velocity of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyVelocityFunction_;    

    //! Function pointer returning radiation pressure force.
    /*!
     * Function pointer returning radiation pressure force [N].
     */
    const DoubleReturningFunction radiationPressureForceFunction_;

    //! Function pointer returning mass of accelerated body.
    /*!
     * Function pointer returning mass of accelerated body [kg].
     */
    const DoubleReturningFunction massFunction_;

    //! Current vector from source to accelerated body.
    /*!
     * Current vector from source to accelerated body (3D vector).
     */
    Eigen::Vector3d currentVectorFromSource_;

    //! Current vector of velocity of accelerated body.
    /*!
     * Current vector of velocity of accelerated body (3D vector).
     */
    Eigen::Vector3d currentVelocityOfAcceleratedBody_;

    //! Current radiation pressure force.
    /*!
     *  Current radiation pressure force[N].
     */
    double currentRadiationPressureForce_;

    //! Current mass of accelerated body.
    /*!
     * Current mass of accelerated body [kg].
     */
    double currentMass_;
};

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_ISOTROPIC_POYNTING_ROBERTSON_DRAG_ACCELERATION_H 
