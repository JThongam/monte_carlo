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
#include <fstream>
#include <iostream>

#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "InputOutput.h"
#include "PatchyGayBerne.h"

InputOutput::InputOutput() {}

void InputOutput::loadConfiguration(std::string fileName, Box& box,
    std::vector<Particle>& particles, CellList& cells, bool isIsotropic)
{
    std::ifstream dataFile;

    // Reset cell list.
    cells.reset();

    // Attempt to read data file.
    dataFile.open(fileName.c_str());

    // Check that the file is valid.
    if (dataFile.good())
    {
        std:: string firstline;

        // Skip the first line.
        std::getline(dataFile, firstline);
        
        for (unsigned int i=0;i<particles.size();i++)
        {
            // Set particle index.
            particles[i].index = i;

            // Resize position and orientation vectors.
            particles[i].position.resize(box.dimension);
            particles[i].orientation.resize(box.dimension);

            // Load position.
            for (unsigned int j=0;j<box.dimension;j++)
                dataFile >> particles[i].position[j];

            if (!isIsotropic)
            {
                // Load orientation.
                for (unsigned int j=0;j<box.dimension;j++)
                    dataFile >> particles[i].orientation[j];
            }
            else
            {
                // Assign dummy orientation.
                for (unsigned int j=0;j<box.dimension;j++)
                    particles[i].orientation[j] = 1.0/sqrt(box.dimension);
            }

            // Enforce periodic boundary conditions.
            box.periodicBoundaries(particles[i].position);

            // Calculate the particle's cell index.
            particles[i].cell = cells.getCell(particles[i]);

            // Update cell list.
            cells.initCell(particles[i].cell, particles[i]);
        }
    }
    else
    {
        std::cout << "[ERROR] InputOutput: Invalid restart file!\n";
        exit(EXIT_FAILURE);
    }

    // Close file stream.
    dataFile.close();
}


////////////////////////////////////////////////////////
void InputOutput::saveConfiguration(std::string fileName, Box& box, std::vector<Particle>& particles, 
    std::vector<double>& patchTheta, double aa, CellList& cells, unsigned int maxInteractions,
    double interactionEnergy, double interactionRange, double epsilonP, double patchAngle)
{
    double cut_off = -5*interactionEnergy;

    // Create file pointer.
    FILE *pFile = fopen(fileName.c_str(), "a");
    
    PatchyGayBerne patchyGB(box, particles, cells, maxInteractions, interactionRange, aa, epsilonP, patchAngle, patchTheta);
    
    int numAtoms = particles.size()*3;
      
    fprintf(pFile, "ITEM: TIMESTEP\n");
    fprintf(pFile, "0\n");
    fprintf(pFile, "ITEM: NUMBER OF ATOMS\n");
    fprintf(pFile, "%d\n", numAtoms);
    fprintf(pFile, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(pFile, "0.0000000000000000e+00 %5.4f \n", box.boxSize[0]);
    fprintf(pFile, "0.0000000000000000e+00 %5.4f \n", box.boxSize[1]);
    fprintf(pFile, "0.0000000000000000e+00 1\n");
    fprintf(pFile, "ITEM: ATOMS id type x y z c_q[1] c_q[2] c_q[3] c_q[4] c_sh[1] c_sh[2] c_sh[3] \n");
    
    std::vector<double> num, hand, pos, ow, ori, oz;
    num.resize(1); hand.resize(1); pos.resize(2); ow.resize(1); ori.resize(2); oz.resize(1);
    double erg;
    
    

    
    unsigned int dimn = patchTheta.size();
    for (unsigned int i=0;i<particles.size();i++){
        double handedness = 1.0 + 3*(i % 2);
        
        
        
        ///////////////////////////////////////////////////////  
        double min=0.0;
        for (unsigned int j=0;j<particles.size();j++){
            double ri = sqrt(pow(particles[i].position[0],2.0)+pow(particles[i].position[1],2.0));
            double rj = sqrt(pow(particles[j].position[0],2.0)+pow(particles[j].position[1],2.0));
            if(i != j){
                if(abs(ri-rj)<=3){
                
                    double position1[2]={particles[i].position[0],particles[i].position[1]};
                    double orientation1[2]={particles[i].orientation[0],particles[i].orientation[1]};
                    
                    double position2[2]={particles[j].position[0],particles[j].position[1]};
                    double orientation2[2]={particles[j].orientation[0],particles[j].orientation[1]};
                    
                    erg = interactionEnergy*(patchyGB.computePairEnergy(i,position1,orientation1, j, position2, orientation2));
                    if(erg < min){
                        min = erg;
                    }
                }
            }
        }
        if(min<cut_off){
            handedness = 7;
        }
        /////////////////////////////////////////////////////
        
        
        for(int j=0;j<2;++j)
        {
            pos[j] = particles[i].position[j];
            ori[j] = particles[i].orientation[j];
        }
        
        // Write particle position.
        fprintf(pFile, "%1.0f %1.0f %5.4f %5.4f %5.4f", 3*i+1.0, handedness, pos[0], pos[1], 0.0);
        if (box.dimension == 3) fprintf(pFile, " %5.4f", pos[2]);
        // Write particle orientation.
        fprintf(pFile, " %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f \n", sqrt((1+ori[0])/2), 0.0, 0.0, ori[1]/(sqrt(2*(1+ori[0]))), 3.0, 1.0, 1.0);
        if (box.dimension == 3) fprintf(pFile, " %5.4f", ori[2]);

        for(unsigned int i0 = 0; i0 < dimn; ++i0)
        {
            double theta = atan2(aa* sin(patchTheta[i0]), cos(patchTheta[i0]));
            double raw0 = 0.5 * aa * cos(theta);
            double raw1 = 0.5 * sin(theta);
            double p1 = (i % 2) ? 1.0 : -1.0;
            double xx = pos[0] + raw0 * ori[0] - raw1 * ori[1] * p1;
            double yy = pos[1] + raw0 * ori[1] + raw1 * ori[0] * p1;
            fprintf(pFile, "%d %d %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f \n", i0+2+3*i, i0+2+3*(i % 2), xx, yy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        
        // Terminate line.
        //fprintf(pFile, "\n");
    }
    // Close file pointer.
    fclose(pFile);
}

void InputOutput::restartfile(std::string fileName, Box& box, std::vector<Particle>& particles){
    // Create file pointer.
    FILE *pFile = fopen(fileName.c_str(), "w");
        
    // Print the box size
    fprintf(pFile, "%.5f %.5f", box.boxSize[0], box.boxSize[1]);
    if (box.dimension == 3) fprintf(pFile, "%.5f", box.boxSize[2]);
    fprintf(pFile, "\n");
    
    for (unsigned int i=0;i<particles.size();i++)
    {
        // Write particle position.
        fprintf(pFile, "%5.4f %5.4f", particles[i].position[0], particles[i].position[1]);
        if (box.dimension == 3) fprintf(pFile, " %5.4f", particles[i].position[2]);

        // Write particle orientation.
        fprintf(pFile, " %.16f %.16f", particles[i].orientation[0], particles[i].orientation[1]);
        if (box.dimension == 3) fprintf(pFile, " %.16f", particles[i].orientation[2]);

        // Terminate line.
        fprintf(pFile, "\n");
        
    }
    
    fclose(pFile);
}

////////////////////////////////////////////////////////
void InputOutput::save_PairEnergy(std::string fileName, Box& box, std::vector<Particle>& particles, CellList& cells, unsigned int maxInteractions, double interactionEnergy, double interactionRange, double aa,  double epsilonP, double patchAngle, std::vector<double>& patchTheta){
    
    double cut_off = -5*interactionEnergy;

    // Create file pointer.
    FILE *pFile = fopen(fileName.c_str(), "a");
    
    PatchyGayBerne patchyGB(box, particles, cells, maxInteractions, interactionRange, aa, epsilonP, patchAngle, patchTheta);
    
    //double number_of_lines = 0.0;
    //std::string line;
    
    double erg;
    double cout=0.0;
    for (unsigned int i=0;i<particles.size();i++){
        double min=0.0;
        for (unsigned int j=0;j<particles.size();j++){
            double ri = sqrt(pow(particles[i].position[0],2.0)+pow(particles[i].position[1],2.0));
            double rj = sqrt(pow(particles[j].position[0],2.0)+pow(particles[j].position[1],2.0));
            if(i != j){
                if(abs(ri-rj)<=3){
                
                    double position1[2]={particles[i].position[0],particles[i].position[1]};
                    double orientation1[2]={particles[i].orientation[0],particles[i].orientation[1]};
                    
                    double position2[2]={particles[j].position[0],particles[j].position[1]};
                    double orientation2[2]={particles[j].orientation[0],particles[j].orientation[1]};
                    
                    erg = interactionEnergy*(patchyGB.computePairEnergy(i,position1,orientation1, j, position2, orientation2));
                    if(erg < min){
                        min = erg;
                    }
                }
            }
        }
        //fprintf(pFile, "%5.4f ", min) ;
        //fprintf(pFile, "\n");
        
        if(min<cut_off){
            cout += 1.0;
        }
    }
    //while (std::getline(fileName.c_str(), line))
        //++number_of_lines;
    fprintf(pFile, "%1.0f", cout) ;
    fprintf(pFile, "\n");
    //}
    
    // Close file pointer.
    fclose(pFile);
}


void InputOutput::appendXyzTrajectory(std::string fileName, unsigned int step, Box& box, std::vector<Particle>& particles, 
    std::vector<double>& patchTheta, double aa, CellList& cells, unsigned int maxInteractions,double interactionEnergy, 
    double interactionRange, double epsilonP, double patchAngle, bool clearFile)
{
    double cut_off = -5*interactionEnergy;
    
    FILE* pFile;

    // Wipe existing trajectory file.
    if (clearFile)
    {
        pFile = fopen(fileName.c_str(), "w");
        fclose(pFile);
    }

    PatchyGayBerne patchyGB(box, particles, cells, maxInteractions, interactionRange, aa, epsilonP, patchAngle, patchTheta);
    
    int numAtoms = particles.size()*3;

    pFile = fopen(fileName.c_str(), "a");
    fprintf(pFile, "ITEM: TIMESTEP\n");
    fprintf(pFile, "%d\n", step);
    fprintf(pFile, "ITEM: NUMBER OF ATOMS\n");
    fprintf(pFile, "%d\n", numAtoms);
    fprintf(pFile, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(pFile, "0.0000000000000000e+00 %5.4f \n", box.boxSize[0]);
    fprintf(pFile, "0.0000000000000000e+00 %5.4f \n", box.boxSize[1]);
    fprintf(pFile, "0.0000000000000000e+00 1\n");
    fprintf(pFile, "ITEM: ATOMS id type x y z c_q[1] c_q[2] c_q[3] c_q[4] c_sh[1] c_sh[2] c_sh[3] \n");

    std::vector<double> num, hand, pos, ow, ori, oz;
    num.resize(1); hand.resize(1); pos.resize(2); ow.resize(1); ori.resize(2); oz.resize(1);
    double erg;
    
    

    
    unsigned int dimn = patchTheta.size();
    for (unsigned int i=0;i<particles.size();i++){
        double handedness = 1.0 + 3*(i % 2);
        
        
        
        ///////////////////////////////////////////////////////  
        double min=0.0;
        for (unsigned int j=0;j<particles.size();j++){
            double ri = sqrt(pow(particles[i].position[0],2.0)+pow(particles[i].position[1],2.0));
            double rj = sqrt(pow(particles[j].position[0],2.0)+pow(particles[j].position[1],2.0));
            if(i != j){
                if(abs(ri-rj)<=3){
                
                    double position1[2]={particles[i].position[0],particles[i].position[1]};
                    double orientation1[2]={particles[i].orientation[0],particles[i].orientation[1]};
                    
                    double position2[2]={particles[j].position[0],particles[j].position[1]};
                    double orientation2[2]={particles[j].orientation[0],particles[j].orientation[1]};
                    
                    erg = interactionEnergy*(patchyGB.computePairEnergy(i,position1,orientation1, j, position2, orientation2));
                    if(erg < min){
                        min = erg;
                    }
                }
            }
        }
        if(min<cut_off){
            handedness = 7;
        }
        /////////////////////////////////////////////////////
        
        
        for(int j=0;j<2;++j)
        {
            pos[j] = particles[i].position[j];
            ori[j] = particles[i].orientation[j];
        }
        
        // Write particle position.
        fprintf(pFile, "%1.0f %1.0f %5.4f %5.4f %5.4f", 3*i+1.0, handedness, pos[0], pos[1], 0.0);
        if (box.dimension == 3) fprintf(pFile, " %5.4f", pos[2]);
        // Write particle orientation.
        fprintf(pFile, " %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f \n", sqrt((1+ori[0])/2), 0.0, 0.0, ori[1]/(sqrt(2*(1+ori[0]))), 3.0, 1.0, 1.0);
        if (box.dimension == 3) fprintf(pFile, " %5.4f", ori[2]);

        for(unsigned int i0 = 0; i0 < dimn; ++i0)
        {
            double theta = atan2(aa* sin(patchTheta[i0]), cos(patchTheta[i0]));
            double raw0 = 0.5 * aa * cos(theta);
            double raw1 = 0.5 * sin(theta);
            double p1 = (i % 2) ? 1.0 : -1.0;
            double xx = pos[0] + raw0 * ori[0] - raw1 * ori[1] * p1;
            double yy = pos[1] + raw0 * ori[1] + raw1 * ori[0] * p1;
            fprintf(pFile, "%d %d %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f \n", i0+2+3*i, i0+2+3*(i % 2), xx, yy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    
    }

    fclose(pFile);
}