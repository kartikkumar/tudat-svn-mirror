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
 *	     Klacka, J. Poynting-Robertson Effect II. Perturbation Equations, Earth, Moon, and Planets,
 *          1992. 
 *
 *    Notes
 *
 */

#include "TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/isotropicPoyntingRobertsonDragForce.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute isotropic Poynting-Robertson drag force. 
Eigen::Vector3d computeIsotropicPoyntingRobertsonDragForce(
        const Eigen::Vector3d& velocityOfAcceleratedBody,
        const Eigen::Vector3d& vectorFromSource,
        const double radiationPressureForce )
{
    // Compute velocity in the radial direction [m/s].
    const double radialVelocity = velocityOfAcceleratedBody.dot( vectorFromSource );

    // Return force due to Poynting-Robertson drag [N].
    return radiationPressureForce 
        * ( ( 1.0 - radialVelocity / tudat::physical_constants::SPEED_OF_LIGHT ) 
            * vectorFromSource - velocityOfAcceleratedBody 
            / tudat::physical_constants::SPEED_OF_LIGHT );
}

} // namespace electro_magnetism
} // namespace tudat
