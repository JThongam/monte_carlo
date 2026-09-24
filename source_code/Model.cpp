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
#include <limits>

#include "Box.h"
#include "CellList.h"
#include "Model.h"
#include "Particle.h"

double INF = std::numeric_limits<double>::infinity();

Model::Model(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionRange_) :

    box(box_),
    particles(particles_),
    cells(cells_),
    maxInteractions(maxInteractions_),
    interactionRange(interactionRange_)
{
    // Work out squared cut-off distance.
    squaredCutOffDistance = interactionRange * interactionRange;
}

double Model::computeEnergy(unsigned int particle, const double* position, const double* orientation)
{
    // N.B. This method is somewhat redundant since the same functionality
    // could be achieved by using a combination of the computeInteractions
    // and model specific computePairEnergy methods.

    // Energy counter.
    double energy = 0;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        // Cell index.
        unsigned int cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            // Index of neighbouring particle.
            unsigned int neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                // Calculate model specific pair energy.
                energy += computePairEnergy(particle, position, orientation,
                          neighbour, &particles[neighbour].position[0],
                          &particles[neighbour].orientation[0]);

                // Early exit test for hard core overlaps and large finite energy repulsions.
                if (energy > 1e6) return INF;
            }
        }
    }

    return energy;
}

double Model::computePairEnergy(unsigned int particle1, const double* position1, const double* orientation1,
    unsigned int particle2, const double* position2, const double* orientation2)
{
    std::cout << "[ERROR] Model: Virtual function Model::computePairEnergy() must be defined.\n";
    exit(EXIT_FAILURE);
}

void Model::volumeMoveUpdates(std::vector<double> boxArgument)
{
    // Update/Restore the box size.
    box.updateBox(boxArgument);

    // Build/Restore the cell list.
    cells.reset();
    cells.setDimension(box.dimension);
    cells.initialise(box.boxSize, interactionRange);
    cells.initCellList(particles);
}

double Model::getEnergy()
{
    double energy = 0;

    for (unsigned int i=0;i<particles.size();i++)
        energy += computeEnergy(i, &particles[i].position[0], &particles[i].orientation[0]);

    return energy/2;
}
