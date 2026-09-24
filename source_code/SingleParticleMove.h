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

#ifndef _SINGLEPARTICLEMOVE_H
#define _SINGLEPARTICLEMOVE_H

#include <vector>

#include "MersenneTwister.h"

/*! \file SingleParticleMove.h
    \brief A class for executing single particle Move Monte Carlo moves (translations and rotations).
*/

// FORWARD DECLARATIONS

class  Model;
struct Particle;

//! Container for storing move parameters.
struct MoveParams
{
    unsigned int seed;                          //!< Index of the seed particle.
    bool isRotation;                            //!< Whether the move is a rotation.
    double stepSize;                            //!< The magnitude of the trial move.
    std::vector<double> trialVector;            //!< Vector for trial move.
    Particle preMoveParticle;                   //!< Particle state before the trial move.
};

class SingleParticleMove
{
public:
    //! Constructor.
    /*! \param model_
            A pointer to the model object.

        \param interactionEnergy_
            The inverse temperature of the system.

        \param pressure_
            The pressure of the simulation.

        \param maxTrialTranslation_
            The maximum trial translation (in units of the reference particle diameter).

        \param maxTrialRotation_
            The maximum trial rotation.

        \param maxTrialVolumeMove_
            The maximum trial volume move.

        \param probTranslate_
            The probability of performing a translation move (versus a rotation).
        
        \param ensemble_
            Whether to use NVT, anisoNPT, or isoNPT ensemble.
     */
    SingleParticleMove(Model*, double, double, double, double, double, double, std::string);

    //! Overloaded ++ operator. Perform a single step.
    void operator ++ (const int);

    //! Overloaded += operator. Perform "n" steps.
    void operator += (const int);

    //! Perform a single trial move.
    void step();

    //! Perform a specified number of trial moves.
    /*! \param nSteps
            The number of attempted trial moves.
     */
    void step(const int);

    //! Get the number of attempted translation moves.
    /*! \return
            The number of attempted translation moves.
     */
    unsigned long long getTranslationAttempts() const;

    //! Get the number of accepted translation moves.
    /*! \return
            The number of accepted translation moves.
     */
    unsigned long long getTranslationAccepts() const;

    //! Get the number of attempted rotation moves.
    /*! \return
            The number of attempted rotation moves.
     */
    unsigned long long getRotationAttempts() const;

    //! Get the number of accepted rotation moves.
    /*! \return
            The number of accepted rotation moves.
     */
    unsigned long long getRotationAccepts() const;

    //! Get the number of attempted volume moves.
    /*! \return
            The number of attempted volume moves.
     */
    unsigned long long getVolumeMoveAttempts() const;

    //! Get the number of accepted volume moves.
    /*! \return
            The number of accepted volume moves.
     */
    unsigned long long getVolumeMoveAccepts() const;

    //! Get the maximum translation value.
    /*! \return
            The maximum translation value.
    */
    double getMaxTranslation() const;

    //! Get the maximum translation value.
    /*! \return
            The maximum rotation value.
    */
    double getMaxRotation() const;

    //! Get the maximum translation value.
    /*! \return
            The maximum volume move value.
    */
    double getMaxVolumeMove() const;

    //! Reset statistics.
    void reset();

    // Adjust the maximum translation/rotation/volume value.
    void adjust();


    MersenneTwister rng;                        //!< Random number generator.

private:
    MoveParams moveParams;                      //!< Parameters for the trial move.
    Model* model;                               //!< A pointer to the model object.
    unsigned long long nTranslationAttempts;    //!< Number of attempted translation moves.
    unsigned long long nTranslationAccepts;     //!< Number of accepted translations moves.
    unsigned long long nRotationAttempts;       //!< Number of attempted rotation moves.
    unsigned long long nRotationAccepts;        //!< Number of accepted rotation moves.
    unsigned long long nVolumeMoveAttempts;    //!< Number of attempted volume moves.
    unsigned long long nVolumeMoveAccepts;     //!< Number of accepted volume moves.

    unsigned long long nTAttempt;
    unsigned long long nTAccept;
    unsigned long long nRAttempt;
    unsigned long long nRAccept;
    unsigned long long nVAttempt;
    unsigned long long nVAccept;

    double interactionEnergy;                   //!< The inverse energy of the systewm.
    double pressure;                            //!< The pressure of the system. 
    double maxTrialTranslation;                 //!< The maximum trial translation (in units of the reference diameter).
    double maxTrialRotation;                    //!< The maximum trial rotation.
    double maxTrialVolumeMove;                  //!< The maximum trial volume move.
    double probTranslate;                       //!< The relative probability of translational moves (vs rotations).
    std::string ensemble;                       //!< Whether to use NVT, anisoNPT, or isoNPT ensemble.
    bool is3D;                                  //!< Whether the simulation is three-dimensional.
    double energyChange;                        //!< Energy change resulting from trial move.

    //! Propose a trial particle translation/rotation.
    void proposeMove();

    //! Determine whether move is accepted.
    bool accept();

    // Prepare for a volume move.
    void volumeMove();

    //! Execute the volume move.
    /*! \param newBox
            The new simulation boxSize.
            
        \param oldBox
            The old simulation boxSize.
     */
    void executeMove(std::vector<double>, std::vector<double>);

    //! Calculate an unbiased rotation vector in 3D (Beard & Schlick, BJ 85 2973 (2003)).
    /*! \param v1
            The vector about which to rotate (either the position or orientation).

        \param v2
            Rotation unit vector.

        \param v3
            The rotation vector.

        \param angle
            Trial rotation angle.
     */
    void rotate3D(std::vector<double>&, std::vector<double>&, std::vector<double>&, double);

    //! Calculate a simple in plane rotatation vector.
    /*! \param v1
            The vector to rotate (either the position or orientation).

        \param v2
            The rotation vector.

        \param angle
            Trial rotation angle.
     */
    void rotate2D(std::vector<double>&, std::vector<double>&, double);

    //! Compute the norm of a vector.
    /*! \param vec
            A reference to the vector

        \return
            The norm of the vector.
     */
    double computeNorm(std::vector<double>&);
};

#endif  /* _SINGLEPARTICLEMOVE_H */
