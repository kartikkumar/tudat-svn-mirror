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
 *		131016	  S. Hirsh			File Created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_ISOTROPIC_POYNTING_ROBERTSON_DRAG_FORCE_H
#define TUDAT_ISOTROPIC_POYNTING_ROBERTSON_DRAG_FORCE_H

#include <Eigen/Core>

namespace tudat
{
namespace electro_magnetism
{

//! Compute isotropic Poynting-Robertson drag force.
/*!
 * Computes Poynting-Robertson drag force on a particle. Assumes the particle radiates 
 * radiates isotropically with respect to the particle's local frame. 
 * \param velocityOfAcceleratedBody Velocity vector of the body which is accelerated by the 
 *			Poynting-Robertson drag. This vector must be with respect to the inertial 
 *          (barycentric) frame.                                                              [m/s]															
 * \param vectorFromSource Vector pointing from source to target. N.B: this must be a unit
 *          vector! To compute the unit vector based on a given position vector, you can
 *          use the .normalize() or .normalized() member functions of an Eigen::Vector3d
 *          object.                                                                             [m]
 * \param radiationPressureForce Force due to radiation pressure from source. The function 
 *			computeCannonBallRadiationPressureForce() can be used to compute this value for
 *			a cannon-ball model. Note: Input is the magnitude of the force, not the vector.
 *			To compute the unit vector based on a given position vector, you can use the 
 *			.norm() member function of an Eigen::Vector3d object.                               [N]		
 * \return Force due to Poynting-Robertson drag.                                                [N]
 * \sa computeCannonBallRadiationPressureForce().
 */
Eigen::Vector3d computeIsotropicPoyntingRobertsonDragForce(
        const Eigen::Vector3d& velocityOfAcceleratedBody,
        const Eigen::Vector3d& vectorFromSource,
        const double radiationPressureForce );

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_ISOTROPIC_POYNTING_ROBERTSON_DRAG_FORCE_H
