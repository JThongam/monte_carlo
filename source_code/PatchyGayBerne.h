/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+vmmc@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _PATCHYDISC_H
#define _PATCHYDISC_H

#include "Model.h"

/*! \file PatchyGayBerne.h
*/

//! Class defining the Patchy-Disc potential.
class PatchyGayBerne : public Model
{
public:
    //! Constructor.
    /*! \param box_
            A reference to the simulation box object.

        \param particles_
            A reference to the particle list.

        \param cells_
            A reference to the cell list object.

        \param maxInteractions_
            The maximum number of interactions per particle (number of patches).

        \param interactionRange_
            The potential cut-off distance (patch diameter).

        \param aa_
            The length of the major axis of the ellipse(oid) i.e., 2a.

        \param epsilonP_
            The multiplicative factor for the patch interaction strength.
     */
    PatchyGayBerne(Box&, std::vector<Particle>&, CellList&, unsigned int, double, double, double, double, std::vector<double>&);

    //! Calculate the pair energy between two particles.
    /*! \param particle1
            The index of the first particle.

        \param position1
            The position vector of the first particle.

        \param orientation1
            The orientation vector of the first particle.

        \param particle2
            The index of the second particle.

        \param position2
            The position vector of the second particle.

        \param orientation2
            The orientation vector of the second particle.

        \return
            The pair energy between particles 1 and 2.
     */
    double computePairEnergy(unsigned int, const double*, const double*, unsigned int, const double*, const double*);

private:
    double aa;                      // the long axis, i.e., 2a; short axis is taken to be 1.0
    double chi;                     // the eccentricity of the ellipse(oid)
    double epsilonP;                // interaction strength of the patchy regions
    double alpha;                   // the parameter alpha that determines the patch interaction range
    std::vector<double> cosTheta;   //!< Lookup table for cosine rotation matrix components.
    std::vector<double> sinTheta;   //!< Lookup table for sine rotation matrix components.
};

#endif  /* _PATCHYDISC_H */
