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

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "Initialise.h"
#include "MersenneTwister.h"

Initialise::Initialise()
{
}

void Initialise::random(std::vector<Particle>& particles, CellList& cells, Box& box, MersenneTwister& rng, double aa_)
{
    // Copy box dimensions.
    boxSize = box.boxSize;

    // Calculate chi.
    chi = (aa_*aa_ - 1.0) / (aa_*aa_ + 1.0);

    unsigned int dim = box.dimension;

    for (unsigned i=0;i<particles.size();i++)
    {
        // Current number of attempted particle insertions.
        unsigned int nTrials = 0;

        // Whether particle overlaps.
        bool isOverlap = true;

        // Temporary vector.
        std::vector<double> vec(dim);

        // Set particle index.
        particles[i].index = i;

        // Keep trying to insert particle until there is no overlap.
        while (isOverlap)
        {
            nTrials++;

            // Generate a random position.
            for (unsigned int j=0;j<dim;j++)
                vec[j] = rng()*box.boxSize[j];

            particles[i].position = vec;

            // Generate a random orientation.
            for (unsigned int j=0;j<dim;j++)
                vec[j] = rng.normal();

            // Calculate vector norm.
            double norm = 0;
            for (unsigned int j=0;j<dim;j++)
                norm += vec[j]*vec[j];
            norm = sqrt(norm);

            // Convert orientation to a unit vector.
            for (unsigned int j=0;j<dim;j++)
                vec[j] /= norm;

            particles[i].orientation = vec;

            // Calculate the particle's cell index.
            particles[i].cell = cells.getCell(particles[i]);

            // See if there is any overlap between particles.
            isOverlap = checkOverlap(particles[i], particles, cells, box);

            // Check trial limit isn't exceeded.
            if (nTrials == MAX_TRIALS)
            {
                std::cout << "[ERROR] Initialise: Maximum number of trial insertions reached.\n";
                exit(EXIT_FAILURE);
            }
        }

        // Update cell list.
        cells.initCell(particles[i].cell, particles[i]);
    }
}

bool Initialise::checkOverlap(Particle& particle, std::vector<Particle>& particles, CellList& cells, Box& box)
{
    unsigned int cell, neighbour;
    int dim = box.dimension;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        cell = cells[particle.cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure particles are different.
            if (neighbour != particle.index)
            {
                // Particle separtion vector.
                std::vector<double> sep(box.dimension);

                // Compute separation.
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] = particle.position[k] - particles[neighbour].position[k];

                // Compute minimum image.
                box.minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<box.dimension;k++)
                    normSqd += sep[k]*sep[k];

                double norm = sqrt(normSqd);
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] /= norm;

                double u1dotu2 = 0.0; double sepdotu1 = 0.0; double sepdotu2 = 0.0;
                for(int idim = 0; idim < dim; ++idim)
                {
                    u1dotu2 += particle.orientation[idim] * particles[neighbour].orientation[idim];
                    sepdotu1 += sep[idim] * particle.orientation[idim];
                    sepdotu2 += sep[idim] * particles[neighbour].orientation[idim];
                }
                double Xu1dotu2 = chi * u1dotu2;

                // the two terms in the expression for sigma - the range parameter
                double sigmaf1 = (sepdotu1 + sepdotu2) * (sepdotu1 + sepdotu2) / (1.0 + Xu1dotu2);
                double sigmaf2 = (sepdotu1 - sepdotu2) * (sepdotu1 - sepdotu2) / (1.0 - Xu1dotu2);

                double sigmaSqd = 1.0 / (1.0 - 0.5 * chi * (sigmaf1 + sigmaf2));

                // Overlap if normSqd is less than particle diameter (box is scaled in diameter units).
                if (normSqd < sigmaSqd) return true;
            }
        }
    }

    // If we get this far, no overlaps.
    return false;
}
