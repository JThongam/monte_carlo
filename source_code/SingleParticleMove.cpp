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

#include "Box.h"
#include "CellList.h"
#include "Model.h"
#include "Particle.h"
#include "SingleParticleMove.h"

SingleParticleMove::SingleParticleMove(
    Model* model_,
    double interactionEnergy_,
    double pressure_,
    double maxTrialTranslation_,
    double maxTrialRotation_,
    double maxTrialVolumeMove_,
    double probTranslate_,
    std::string ensemble_) :
    

    model(model_),
    interactionEnergy(interactionEnergy_),
    pressure(pressure_),
    maxTrialTranslation(maxTrialTranslation_),
    maxTrialRotation(maxTrialRotation_),
    maxTrialVolumeMove(maxTrialVolumeMove_),
    probTranslate(probTranslate_),
    ensemble(ensemble_)
{
    // Check dimensionality.
    if (model->box.dimension == 3) is3D = true;
    else is3D = false;

    // Reset the statistics.
    reset();

    // Allocate memory.
    moveParams.trialVector.resize(model->box.dimension);

    // Check that orientation is a unit vector.
    for (unsigned int i = 0; i < model->particles.size(); i++)
    {
        if (std::abs(1.0 - computeNorm(model->particles[i].orientation)) > 1e-6)
        {
            std::cout << "[ERROR] VMMC: Particle orientations must be unit vectors!\n";
            exit(EXIT_FAILURE);
        }
    }
    
}

void SingleParticleMove::step(const int nSteps)
{
    if (ensemble == "nvt") {
        for (int i=0;i<nSteps;i++) step();
    }
    else{    
        // Perform a volume move for every (nParticles+1) steps.
        for (int i=0;i<nSteps;i++)
        {
            unsigned int volumeMoveIndex = rng() * (nSteps+1);
            if (volumeMoveIndex <= model->particles.size()) step();
            else volumeMove();
            
        }
    }
}

void SingleParticleMove::operator ++ (const int)
{
    step();
}

void SingleParticleMove::operator += (const int nSteps)
{
    step(nSteps);
}

void SingleParticleMove::step()
{
    // Propose a trial move.
    proposeMove();

    // Check whether move was accepted.
    if (accept())
    {
        // Increment number of rotations.
        nRotationAccepts += moveParams.isRotation;

        if (!moveParams.isRotation)
        {
            // Translation is accepted.
            nTranslationAccepts++;

            unsigned int oldCell = moveParams.preMoveParticle.cell;
            unsigned int newCell = model->particles[moveParams.seed].cell;

            // update cell list
            if (oldCell != newCell)
            {
                model->particles[moveParams.seed].cell = oldCell;
                model->cells.updateCell(newCell, model->particles[moveParams.seed], model->particles);
            }
        }
    }
    else
    {
        // Revert particle to pre-move state.
        model->particles[moveParams.seed] = moveParams.preMoveParticle;
    }
}

unsigned long long SingleParticleMove::getTranslationAttempts() const
{
    return nTranslationAttempts;
}

unsigned long long SingleParticleMove::getTranslationAccepts() const
{
    return nTranslationAccepts;
}

unsigned long long SingleParticleMove::getRotationAttempts() const
{
    return nRotationAttempts;
}

unsigned long long SingleParticleMove::getRotationAccepts() const
{
    return nRotationAccepts;
}

unsigned long long SingleParticleMove::getVolumeMoveAttempts() const
{
    return nVolumeMoveAttempts;
}

unsigned long long SingleParticleMove::getVolumeMoveAccepts() const
{
    return nVolumeMoveAccepts;
}

double SingleParticleMove::getMaxTranslation() const
{
    return maxTrialTranslation;
}

double SingleParticleMove::getMaxRotation() const
{
    return maxTrialRotation;
}

double SingleParticleMove::getMaxVolumeMove() const
{
    return maxTrialVolumeMove;
}

void SingleParticleMove::reset()
{
    nTranslationAttempts = nTranslationAccepts = 0;
    nRotationAttempts = nRotationAccepts = 0;
    nVolumeMoveAttempts = nVolumeMoveAccepts = 0;
    nTAttempt = nTAccept = 0;
    nRAttempt = nRAccept = 0;
    nVAttempt = nVAccept = 0;
}

void SingleParticleMove::proposeMove()
{
    // Choose a seed particle.
    moveParams.seed = rng.integer(0, model->particles.size()-1);

    // Choose a random point on the surface of the unit sphere/circle.
    for (unsigned int i=0;i<model->box.dimension;i++)
        moveParams.trialVector[i] = rng.normal();

    // Normalise the trial vector.
    double norm = computeNorm(moveParams.trialVector);
    for (unsigned int i=0;i<model->box.dimension;i++)
        moveParams.trialVector[i] /= norm;

    // Choose the move type.
    if (rng() < probTranslate)
    {
        // Translation.
        nTranslationAttempts++;
        moveParams.isRotation = false;

        // Scale step-size to uniformly sample unit sphere/circle.
        if (is3D) moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/3.0);
        else moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/2.0);
    }
    else
    {
        // Rotation.
        nRotationAttempts++;
        moveParams.isRotation = true;
        moveParams.stepSize = maxTrialRotation*(2.0*rng()-1.0);
    }

    // Calculate pre-move energy.
    double initialEnergy = model->computeEnergy(moveParams.seed,
        &model->particles[moveParams.seed].position[0],
        &model->particles[moveParams.seed].orientation[0]);

    // Store initial coordinates/orientation.
    moveParams.preMoveParticle = model->particles[moveParams.seed];

    // Execute the move.
    if (!moveParams.isRotation) // Translation.
    {
        for (unsigned int i=0;i<model->box.dimension;i++)
            model->particles[moveParams.seed].position[i] += moveParams.stepSize*moveParams.trialVector[i];

        // Apply periodic boundary conditions.
        model->box.periodicBoundaries(model->particles[moveParams.seed].position);

        // Work out new cell index.
        model->particles[moveParams.seed].cell = model->cells.getCell(model->particles[moveParams.seed]);
    }
    else                        // Rotation.
    {
        std::vector<double> vec(model->box.dimension);

        // Calculate orientation rotation vector.
        if (is3D) rotate3D(model->particles[moveParams.seed].orientation, moveParams.trialVector, vec, moveParams.stepSize);
        else rotate2D(model->particles[moveParams.seed].orientation, vec, moveParams.stepSize);

        // Update orientation.
        for (unsigned int i=0;i<model->box.dimension;i++)
            model->particles[moveParams.seed].orientation[i] += vec[i];
    }

    // Calculate post-move energy.
    double finalEnergy = model->computeEnergy(moveParams.seed,
        &model->particles[moveParams.seed].position[0],
        &model->particles[moveParams.seed].orientation[0]);

    energyChange = finalEnergy - initialEnergy;
}

bool SingleParticleMove::accept()
{
    double exponentGuard = 75.0;
    if ((energyChange == INF) or (energyChange > exponentGuard)) return false;
    if ((energyChange == 0) or (rng() < exp(-interactionEnergy*energyChange))) return true;
    else return false;
}

void SingleParticleMove::volumeMove()
{
    nVolumeMoveAttempts++;

    // Old Energy
    double oldEnergy = model->getEnergy();

    // Old Box Size.
    std::vector<double> oldBoxSize;
    if (is3D) oldBoxSize = {model->box.boxSize[0], model->box.boxSize[1], model->box.boxSize[2]};
    else oldBoxSize = {model->box.boxSize[0], model->box.boxSize[1]};

    // Old Volume
    double oldVolume = oldBoxSize[0] * oldBoxSize[1];
    if (is3D) oldVolume *= oldBoxSize[2];

    // Random walk in lnv
    double lnvn = log(oldVolume) + (rng() - 0.5) * maxTrialVolumeMove;

    // New volume
    double newVolume = exp(lnvn);

    // New Box Size
    std::vector<double> newBoxSize;
    if(ensemble == "isonpt"){
        if (is3D) newBoxSize = {pow(newVolume, 1.0/3.0), pow(newVolume, 1.0/3.0), pow(newVolume, 1.0/3.0)};
        else newBoxSize = {pow(newVolume, 1.0/2.0), pow(newVolume, 1.0/2.0)};
    }
    else if (ensemble == "anisonpt"){
        // Choose Which Axis to Move
        unsigned int axis;
        if (is3D) axis = static_cast<int>(rng() * 3);
        else axis = static_cast<int>(rng() * 2);

        double newSide = newVolume;
        for (unsigned int i = 0; i < model->box.dimension; i++)
            if (i != axis) newSide /= oldBoxSize[i];

        if (is3D){
            if (axis == 0) newBoxSize = {newSide, oldBoxSize[1], oldBoxSize[2]};
            else if (axis == 1) newBoxSize = {oldBoxSize[0], newSide, oldBoxSize[2]};
            else if (axis == 2) newBoxSize = {oldBoxSize[0], oldBoxSize[1], newSide};}
        else{
            if (axis == 0) newBoxSize = {newSide, oldBoxSize[1]};
            else if (axis == 1) newBoxSize = {oldBoxSize[0], newSide};}
    }
    

    // Perform volume move
    executeMove(newBoxSize, oldBoxSize);

    // New energy
    double newEnergy = model->getEnergy();

    // Appropriate weight function.
    double argument = -interactionEnergy * (
        (newEnergy - oldEnergy) + pressure * (newVolume - oldVolume) - (model->particles.size() + 1) * log(newVolume/oldVolume) / interactionEnergy
    );

    double exponentGuard = 75.0;
    if ((argument > exponentGuard) or (rng() >= exp(argument))) executeMove(oldBoxSize, newBoxSize); // Reverse the volume move
    else nVolumeMoveAccepts++; // Move accepted
}

void SingleParticleMove::executeMove(std::vector<double> newBox, std::vector<double> oldBox)
{
    for (unsigned int i = 0; i < model->particles.size(); i++)
        for (unsigned int j = 0; j < model->box.dimension; j++) model->particles[i].position[j] *= newBox[j]/oldBox[j];

    // Update/Restore the Box Size.
    model->volumeMoveUpdates(newBox);
}

void SingleParticleMove::rotate3D(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& v3, double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    double v1Dotv2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

    v3[0] = ((v1[0] - v2[0]*v1Dotv2))*(c - 1) + (v2[2]*v1[1] - v2[1]*v1[2])*s;
    v3[1] = ((v1[1] - v2[1]*v1Dotv2))*(c - 1) + (v2[0]*v1[2] - v2[2]*v1[0])*s;
    v3[2] = ((v1[2] - v2[2]*v1Dotv2))*(c - 1) + (v2[1]*v1[0] - v2[0]*v1[1])*s;
}

void SingleParticleMove::rotate2D(std::vector<double>& v1, std::vector<double>& v2, double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    v2[0] = (v1[0]*c - v1[1]*s) - v1[0];
    v2[1] = (v1[0]*s + v1[1]*c) - v1[1];
}


double SingleParticleMove::computeNorm(std::vector<double>& vec)
{
    double normSquared = 0;

    for (unsigned int i=0;i<vec.size();i++)
        normSquared += vec[i]*vec[i];

    return sqrt(normSquared);
}

void SingleParticleMove::adjust()
{
    if (nTranslationAttempts == 0)
    {
        nTAttempt = nTranslationAttempts;
        nTAccept = nTranslationAccepts;
        nRAttempt = nRotationAttempts;
        nRAccept = nRotationAccepts;
        nVAttempt = nVolumeMoveAttempts;
        nVAccept = nVolumeMoveAccepts;
    }
    else
    {
        // Calculate the current translation acceptance rate.
        double fracT = static_cast<double>(nTranslationAccepts - nTAccept) / static_cast<double>(nTranslationAttempts - nTAttempt);

        // Calculate the current rotation acceptance rate.
        double fracR = static_cast<double>(nRotationAccepts - nRAccept) / static_cast<double>(nRotationAttempts - nRAttempt);

        // Calculate the current volume move acceptance rate.
        double fracV = static_cast<double>(nVolumeMoveAccepts - nVAccept) / static_cast<double>(nVolumeMoveAttempts - nVAttempt);

        // Store the old maximum trial translation value.
        double oldMaxTrialTranslation = maxTrialTranslation;

        // Store the old maximum trial rotation value.
        double oldMaxTrialRotation = maxTrialRotation;

        // Store the old maximum trial volume move value.
        double oldMaxTrialVolumeMove = maxTrialVolumeMove;

        // Calculate how much trial translation value should be to achieve 40% (0.4) acceptance rate.
        maxTrialTranslation = maxTrialTranslation * std::abs(fracT / 0.4);

        // Calculate how much trial rotation value should be to achieve 40% (0.4) acceptance rate.
        maxTrialRotation = maxTrialRotation * std::abs(fracR / 0.4);

        // Calculate how much trial volume move value should be to achieve 40% (0.4) acceptance rate.
        maxTrialVolumeMove = maxTrialVolumeMove * std::abs(fracV / 0.4);

        // Limit the change to avoid drastic adjustments.
        if (maxTrialTranslation/oldMaxTrialTranslation > 1.5) maxTrialTranslation = oldMaxTrialTranslation * 1.5;
        if (maxTrialTranslation/oldMaxTrialTranslation < 0.5) maxTrialTranslation = oldMaxTrialTranslation * 0.5;

        double smallestSide = model->box.boxSize[0];
        if (smallestSide < model->box.boxSize[1]) smallestSide = model->box.boxSize[1];
        if(is3D) if (smallestSide < model->box.boxSize[2]) smallestSide = model->box.boxSize[2];
        
        if (maxTrialTranslation > smallestSide/2) maxTrialTranslation = smallestSide/2;

        // Limit the change to avoid drastic adjustments.
        if (maxTrialRotation/oldMaxTrialRotation > 1.5) maxTrialRotation = oldMaxTrialRotation * 1.5;
        if (maxTrialRotation/oldMaxTrialRotation < 0.5) maxTrialRotation = oldMaxTrialRotation * 0.5;
        if (maxTrialRotation > M_PI) maxTrialRotation = M_PI;        // Max 180 degrees
        if (maxTrialRotation < 1e-6) maxTrialRotation = 1e-6;        // Minimum sensible rotation

        // Limit the change to avoid drastic adjustments.
        if (maxTrialVolumeMove/oldMaxTrialVolumeMove > 1.5) maxTrialVolumeMove = oldMaxTrialVolumeMove * 1.5;
        if (maxTrialVolumeMove/oldMaxTrialVolumeMove < 0.5) maxTrialVolumeMove = oldMaxTrialVolumeMove * 0.5;

        nTAttempt = nTranslationAttempts;
        nTAccept = nTranslationAccepts;

        nRAttempt = nRotationAttempts;
        nRAccept = nRotationAccepts;

        nVAttempt = nVolumeMoveAttempts;
        nVAccept = nVolumeMoveAccepts;
    }
}