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

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "Box.h"
#include "Particle.h"
#include "CellList.h"
#include "PatchyGayBerne.h"

PatchyGayBerne::PatchyGayBerne(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionRange_, double aa_,  double epsilonP_, double patchAngle_, std::vector<double>& patchTheta_) :
    Model(box_, particles_, cells_, maxInteractions_, interactionRange_)
{
    aa = aa_;
    epsilonP = epsilonP_;
    chi = (aa*aa - 1.0) / (aa*aa + 1.0);
    alpha = patchAngle_;

    unsigned int dimn = patchTheta_.size();

    // Resize rotation matrix lookup tables.
    cosTheta.resize(dimn);
    sinTheta.resize(dimn);
    // Populate lookup tables.
    for (unsigned int i=0;i<dimn;i++)
    {
        double theta = atan2(aa* sin(patchTheta_[i]), cos(patchTheta_[i]));
        cosTheta[i] = 0.5 * aa * cos(theta);
        sinTheta[i] = 0.5 * sin(theta);
    }
}

const double mu = 1.0;
const double nu = 3.0;
const double Kp = 1.0; // is 1.0 isotropic? doesn't look like from simulations
const double Kp_rmu = pow(Kp, 1.0 / mu);
const double chi_p = (Kp_rmu - 1.0) / (Kp_rmu + 1.0);


double PatchyGayBerne::computePairEnergy(unsigned int particle1, const double* position1,
    const double* orientation1, unsigned int particle2, const double* position2, const double* orientation2)
{
    int dim = box.dimension;
    // Separation vector.
    std::vector<double> sep(dim);

    // Calculate particle separation.
    for(int idim = 0; idim < dim; ++idim)
        sep[idim] = position1[idim] - position2[idim];

    // Enforce minimum image.
    box.minimumImage(sep);

    double ssep = 0.0;
    // calculate the scalar separation
    for(int idim = 0; idim < dim; ++idim)
        ssep += sep[idim] * sep[idim];
    ssep = sqrt(ssep);

    // normalize r12
    for(int idim = 0; idim < dim; ++idim)
        sep[idim] /= ssep;

    double u1dotu2 = 0.0; double sepdotu1 = 0.0; double sepdotu2 = 0.0;
    for(int idim = 0; idim < dim; ++idim)
    {
        u1dotu2 += orientation1[idim] * orientation2[idim];
        sepdotu1 += sep[idim] * orientation1[idim];
        sepdotu2 += sep[idim] * orientation2[idim];
    }

    double Xu1dotu2 = chi * u1dotu2;

    double epsilon2 = 1.0 / sqrt(1.0 - Xu1dotu2 * Xu1dotu2);

    double Xpu1dotu2 = chi_p * u1dotu2;
    double epsilon1f1 = (sepdotu1 + sepdotu2) * (sepdotu1 + sepdotu2) / (1.0 + Xpu1dotu2);
    double epsilon1f2 = (sepdotu1 - sepdotu2) * (sepdotu1 - sepdotu2) / (1.0 - Xpu1dotu2);

    double epsilon1 = 1.0 - 0.5 * chi_p * (epsilon1f1 + epsilon1f2);

    # ifdef REPULSIVE
    double epsilon = 1.0;
    # else
    // the strength parameter
    double epsilon = pow(epsilon1,mu) * pow(epsilon2,nu);
    # endif
    
    // the two terms in the expression for sigma - the range parameter
    double sigmaf1 = (sepdotu1 + sepdotu2) * (sepdotu1 + sepdotu2) / (1.0 + Xu1dotu2);
    double sigmaf2 = (sepdotu1 - sepdotu2) * (sepdotu1 - sepdotu2) / (1.0 - Xu1dotu2);

    double sigma = 1.0 / sqrt(1.0 - 0.5 * chi * (sigmaf1 + sigmaf2));

    double s0bysep = 1.0 / (ssep - sigma + 1.0);

    double s0bysep6 = s0bysep * s0bysep * s0bysep;
    s0bysep6 *= s0bysep6;

    // non-patchy GB energy
    double energy1 = 4.0 * epsilon * (s0bysep6 * s0bysep6 - s0bysep6);
    
    double p1 = (particle1 % 2) ? 1.0 : -1.0;
    double p2 = (particle2 % 2) ? 1.0 : -1.0;
    double energy2=0.0;
    
    unsigned int dimn = cosTheta.size();

    // Test interactions between all patch pairs.
    for (unsigned int i=0;i<dimn;i++)
    {
        // Compute position of patch i on first disc.
        std::vector<double> patchdir1(2);
        patchdir1[0] = orientation1[0]*cosTheta[i] - orientation1[1]*sinTheta[i]*p1;
        patchdir1[1] = p1*orientation1[0]*sinTheta[i] + orientation1[1]*cosTheta[i];

        std::vector<double> patchpos1(2);
        patchpos1[0] = position1[0] + patchdir1[0];
        patchpos1[1] = position1[1] + patchdir1[1];

        int j = i; // i.e., only compute the interactions between patches of the same kind
        // well, technically I could have used i instead of j, but whatever...
        {
            // Compute position of patch j on second disc.
            std::vector<double> patchdir2(2);
            patchdir2[0] = orientation2[0]*cosTheta[j] - orientation2[1]*sinTheta[j]*p2;
            patchdir2[1] = p2*orientation2[0]*sinTheta[j] + orientation2[1]*cosTheta[j];

            std::vector<double> patchpos2(2);
            patchpos2[0] = position2[0] + patchdir2[0];
            patchpos2[1] = position2[1] + patchdir2[1];

            std::vector<double> patchsep(2);
            patchsep[0] = patchpos2[0] - patchpos1[0];
            patchsep[1] = patchpos2[1] - patchpos1[1];

            double patchsep2 = patchsep[0] * patchsep[0] + patchsep[1] * patchsep[1];

            energy2 -= epsilonP * exp(-0.5 * patchsep2 / alpha / alpha); // i.e., alpha = 0.1
        }
    }
    if((particle1 + particle2) % 2 == 0){
                energy2 *= 0.0;
    }
    double energy=energy1+energy2;

    return energy;
}

