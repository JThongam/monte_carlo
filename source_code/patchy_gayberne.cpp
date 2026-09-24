/*
  Copyright (c) 2022 Arunkumar Bupathy <arunbupathy@gmail.com>

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

#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>

#include "Box.h"
#include "CellList.h"
#include "Initialise.h"
#include "InputOutput.h"
#include "Particle.h"
#include "PatchyGayBerne.h"
#include "SingleParticleMove.h"

struct SimParams
{
    unsigned int productRun                         = 400000;     // Run length of the simulation.
    unsigned int equilibRun                         = 10000;     // Equilibration length.
    unsigned int dimension                          = 2;        // Dimension of the simulation.
    unsigned int nParticles                         = 300;      // Number of particles.
    std::string ensemble                            = "isonpt"; // "nvt", "isonpt", "anisonpt".
    double pressure                                 = 1.0;      // Pressure of the system.
    double temperature                              = 1.0;      // Temperature of the system.
    double density                                  = 0.2;      // Initial Density of the system
    bool useRestart                                 = false;    // Whether the simulation is starting from a restart file.
    std::string rname                               = "Location of Restart File"; // Location of restart
};

[[noreturn]] void print_usage_and_exit(const char* prog, const std::string& msg = "") {
    if (!msg.empty()) std::cerr << "error: " << msg << "\n\n";
    std::cerr
        << "Usage: " << prog << " [options]\n"
        << "Options (defaults in [ ]):\n"
        << "  -e, --ensemble <nvt|isonpt|anisonpt>   [" << "isonpt" << "]\n"
        << "  -P, --pressure <double>                [" << 1.0 << "]\n"
        << "  -T, --temperature <double>             [" << 1.0 << "]  (interactionEnergy = 1/T)\n"
        << "  -d, --density <double>                 [" << 0.2 << "]\n"
        << "  -N, --nParticles <uint>                [" << 300 << "]\n"
        << "  -D, --dimension <2|3>                  [" << 2 << "]\n"
        << "  -r, --productRun <uint>                [" << 1000000 << "]\n"
        << "  -q, --equilibRun <uint>                [" << 0 << "]\n"
        << "  -R, --useRestart <bool>                [" << false << "]\n"
        << "  -L, --rname <std::string>              [ Location of Restart File ]\n"
        << "  -h, --help\n";
    std::exit(msg.empty() ? 0 : 2);
}

static bool is_flag(const std::string& s, const char* longf, const char* shortf) {
    return s == longf || s == shortf;
}

template <class T>
T read_value(int& i, int argc, char** argv, const char* name) {
    if (i + 1 >= argc) print_usage_and_exit(argv[0], std::string("missing value for ") + name);
    std::istringstream iss(argv[++i]);
    T v{};
    if (!(iss >> v)) print_usage_and_exit(argv[0], std::string("invalid value for ") + name + ": " + argv[i]);
    return v;
}

void parse_args(int argc, char** argv, SimParams& p) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (is_flag(a, "--help", "-h")) print_usage_and_exit(argv[0]);

        else if (is_flag(a, "--ensemble", "-e"))        p.ensemble         = read_value<std::string>(i, argc, argv, "--ensemble");
        else if (is_flag(a, "--pressure", "-P"))        p.pressure         = read_value<double>(i, argc, argv, "--pressure");
        else if (is_flag(a, "--temperature", "-T"))     p.temperature      = read_value<double>(i, argc, argv, "--temperature");
        else if (is_flag(a, "--density", "-d"))         p.density          = read_value<double>(i, argc, argv, "--density");
        else if (is_flag(a, "--nParticles", "-N"))      p.nParticles       = read_value<unsigned int>(i, argc, argv, "--nParticles");
        else if (is_flag(a, "--dimension", "-D"))       p.dimension        = read_value<unsigned int>(i, argc, argv, "--dimension");
        else if (is_flag(a, "--productRun", "-r"))      p.productRun       = read_value<unsigned int>(i, argc, argv, "--productRun");
        else if (is_flag(a, "--equilibRun", "-q"))      p.equilibRun       = read_value<unsigned int>(i, argc, argv, "--equilibRun");
        else if (is_flag(a, "--useRestart", "-R"))      p.useRestart       = read_value<bool>(i, argc, argv, "--useRestart");
        else if (is_flag(a, "--rname", "-L"))           p.rname            = read_value<std::string>(i, argc, argv, "--rname");
        else print_usage_and_exit(argv[0], "unknown option: " + a);
    }

    if (p.dimension != 2 && p.dimension != 3) print_usage_and_exit(argv[0], "dimension must be 2 or 3");
    if (p.temperature <= 0.0) print_usage_and_exit(argv[0], "temperature must be > 0");
}

static void print_summary(const SimParams& p) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Ensemble           : " << p.ensemble << "\n"
              << "Pressure           : " << p.pressure << "\n"
              << "Temperature (kBT)  : " << p.temperature << "\n"
              << "Density            : " << p.density << "\n"
              << "nParticles         : " << p.nParticles << "\n"
              << "Dimension          : " << p.dimension << "\n"
              << "productRun         : " << p.productRun << "\n"
              << "equilibRun         : " << p.equilibRun << "\n"
              << "Use Restart?       : " << p.useRestart << "\n"
              << "Restart Location   : " << p.rname << "\n";
}


int main(int argc, char** argv)
{
    SimParams params;
    parse_args(argc, argv, params);
    print_summary(params);

    // Simulation parameters.
    const unsigned int productRun               = params.productRun;
    const unsigned int equilibRun               = params.equilibRun;
    const unsigned int dimension                = params.dimension;
    const unsigned int nParticles               = params.nParticles;
    const double pressure                       = params.pressure;
    const double density                        = params.density;
    const std::string ensemble                  = params.ensemble;
    const bool useRestart                       = params.useRestart;
    std::string rname                           = params.rname;

    unsigned int fileInterval = 1;
    if (params.productRun > 1000) fileInterval  = params.productRun/1000; 
    double interactionEnergy                    = 1 / params.temperature; // = 1/temperature    
    double interactionRange                     = 4.0;
    unsigned int maxInteractions                = 100;

    // MC move sizes
    double maxT                                 = 0.05;
    double maxR                                 = 0.06;
    double maxV                                 = 0.15;
    double probT                                = 0.5; // translation probability

    double baseLengthX;                         // base length of simulation box along x axis.
    double baseLengthY;                         // base length of simulation box along y axis.
    double baseLengthZ;                         // base length of simulation box along z axis.


    double aa                                   = 3.0; // the long axis, i.e., 2a; short axis is taken to be 1.0
    double epsilonP                             = 15; // interaction strength of the patchy regions relative to base energy
    double patchAngle                           = 0.1; // the parameter alpha's value which, in a way decides, the patch angle

    std::vector<double> patchTheta;
    patchTheta.emplace_back(5 * M_PI / 180.0);  // 5 degrees from the pole
    patchTheta.emplace_back(10 * M_PI / 180.0); // 10 degrees from the pole

    // Open output files for erg ans statistics.
    std::string fileName = "dump.txt";          // File to store system snapshots.
    std::ofstream ergFile("energy.txt");        // File to store crystal energy, lattice energy etc.
    std::ofstream outputFile("systemStat.txt"); // File to store the system statistics.

    ///////////////////////////////////////////////////////////
    // Data structures.
    std::vector<Particle> particles(nParticles);// particle container
    CellList cells;                             // cell list

    // Resize particle container.
    particles.resize(nParticles);

    // Path of the restart file.
    if (useRestart){
        // Read box dimensions from box.txt
        std::ifstream restart_File(rname);
        if (!restart_File.is_open()) {
            std::cerr << "Error: Failed to open " << rname << std::endl;
            return 1;
        }

        std::string line;

        // Read the first line.
        if (std::getline(restart_File, line)) {
            std::istringstream iss(line);
            if (dimension == 2){
               double readBaseLengthX, readBaseLengthY;
                if (!(iss >> readBaseLengthX >> readBaseLengthY)) {
                    std::cerr << "Error: Invalid format in " << rname << std::endl;
                    return 1;}
                if (ensemble == "isonpt" && readBaseLengthX != readBaseLengthY){
                    std::cerr << "Error: Square Box required in isotropic NPT ensemble." << std::endl;
                    return 1;}
                
                baseLengthX = readBaseLengthX;
                baseLengthY = readBaseLengthY;
            }else{
                double readBaseLengthX, readBaseLengthY, readBaseLengthZ;
                if (!(iss >> readBaseLengthX >> readBaseLengthY >> readBaseLengthZ)) {
                    std::cerr << "Error: Invalid format in " << rname << std::endl;
                    return 1;}
                if (ensemble == "isonpt" && readBaseLengthX != readBaseLengthY != readBaseLengthZ){
                    std::cerr << "Error: Cubic Box required in isotropic NPT ensemble." << std::endl;
                    return 1;}

                baseLengthX = readBaseLengthX;
                baseLengthY = readBaseLengthY;
                baseLengthZ = readBaseLengthZ;
            }
        }
        else {
            std::cerr << "Error: File is empty or cannot read first line in " << rname << std::endl;
            return 1;
        }
    }
    else
    {
        // Work out base length of simulation box (particle diameter is one).
        if (dimension == 2) {
            baseLengthX = std::pow((nParticles*M_PI*aa)/(4.0*1.0*density), 1.0/2.0);
            baseLengthY = baseLengthX;}
        else {
            baseLengthX = std::pow((nParticles*M_PI*aa)/(6.0*3.0*density), 1.0/3.0);
            baseLengthY = baseLengthX;
            baseLengthZ = baseLengthX;
        }
    }

    std::vector<double> boxSize;
    boxSize.push_back(baseLengthX);
    boxSize.push_back(baseLengthY);
    outputFile << "Base LengthX: " << baseLengthX << std::endl;
    outputFile << "Base LengthY: " << baseLengthY << std::endl;
    if(dimension == 3) {
        boxSize.push_back(baseLengthZ);
        outputFile << "Base LengthZ: " << baseLengthZ << std::endl;
    }

    // Initialise simulation box object.
    Box box(boxSize);

    // Initialise input/output class,
    InputOutput io;

    // Initialise cell list.
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // Initialise the patchy GB model.
    PatchyGayBerne patchyGB(box, particles, cells, maxInteractions, interactionRange, aa, epsilonP, patchAngle, patchTheta);

    // Initialise random number generator.
    MersenneTwister rng;
    outputFile << "Seed:\t" << rng.getSeed() << "\n\n";

    // Initialise particle initialisation object.
    Initialise initialise;

    if (useRestart) io.loadConfiguration(rname, box, particles, cells, false);      // Read a restart file.
    else initialise.random(particles, cells, box, rng, aa);                         // Generate a random particle configuration.
    
    // Initialise SingleParticleMove object.
    SingleParticleMove singleParticleMove(&patchyGB, interactionEnergy, pressure, maxT, maxR, maxV, probT, ensemble);

    // Store the initial snapshot.
    io.appendXyzTrajectory(fileName, 0, box, particles, patchTheta, aa, cells, maxInteractions, interactionEnergy, interactionRange, epsilonP, patchAngle, true);

    // Initial crystal energy
    double erg = patchyGB.getEnergy();  
    double den = (nParticles*M_PI*3) / (4 * box.boxSize[0] * box.boxSize[1]); 
    outputFile << std::fixed << std::setprecision(7);
    outputFile << "\nTotal energy of initial configuration  : " << erg/nParticles << std::endl;
    outputFile << "\nInitial density of the configuration  : " << den << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Execute the equilibration process.
    outputFile << "\nSystem Equilibration\n####################" << std::endl;

    // Variable declaration
    unsigned long long attemptedTranslation; unsigned long long acceptedTranslation; 
    unsigned long long attemptedRotation; unsigned long long acceptedRotation;
    unsigned long long attemptedVolumeMove; unsigned long long acceptedVolumeMove;
    double maxTranslation; double maxRotation; double maxVolumeMove;
    double translationAcceptanceRate; double rotationAcceptanceRate; double volumeMoveAcceptanceRate;

    // Call the adjust function to initialise it.
    singleParticleMove.adjust();

    for (unsigned int i = 0; i < equilibRun; i++)
    {
        unsigned int step = i + 1;

        // Increment simulation by 1 Monte Carlo Sweeps.
        singleParticleMove += (nParticles+1);

        if (step % (equilibRun / 100) == 0)
        {
            outputFile << "\nAdjust at " << step << " of " << equilibRun << " equilibration run.\n";
            singleParticleMove.adjust();
            attemptedTranslation = singleParticleMove.getTranslationAttempts();
            acceptedTranslation = singleParticleMove.getTranslationAccepts();
            attemptedRotation = singleParticleMove.getRotationAttempts();
            acceptedRotation = singleParticleMove.getRotationAccepts();
            maxTranslation = singleParticleMove.getMaxTranslation();
            maxRotation = singleParticleMove.getMaxRotation();
            translationAcceptanceRate = static_cast<double>(acceptedTranslation) / attemptedTranslation;
            rotationAcceptanceRate = static_cast<double>(acceptedRotation) / attemptedRotation; 

            if(ensemble == "isonpt" || ensemble == "anisonpt"){
                attemptedVolumeMove = singleParticleMove.getVolumeMoveAttempts();
                acceptedVolumeMove = singleParticleMove.getVolumeMoveAccepts();
                maxVolumeMove = singleParticleMove.getMaxVolumeMove();
                volumeMoveAcceptanceRate = static_cast<double>(acceptedVolumeMove) / attemptedVolumeMove;
            }
            outputFile << "Attempted Translation : " << attemptedTranslation << "\n" << "Accepted Translation  : " << acceptedTranslation << "\n";
            outputFile << "Attempted Rotation : " << attemptedRotation << "\n" << "Accepted Rotation  : " << acceptedRotation << "\n";
            if(ensemble == "isonpt" || ensemble == "anisonpt")
                outputFile << "Attempted Volume Move : " << attemptedVolumeMove << "\n" << "Accepted VolumeMove  : " << acceptedVolumeMove << "\n";

            outputFile << std::fixed << std::setprecision(6);
            outputFile << "\nMaximum Translation : " << maxTranslation << "\n";
            outputFile << "Maximum Rotation : " << maxRotation << "\n";
            if(ensemble == "isonpt" || ensemble == "anisonpt")
                outputFile << "Maximum Volume Move : " << maxVolumeMove << "\n";

            outputFile << std::fixed << std::setprecision(2);
            outputFile << "\nTranslation Acceptance Rate : " << translationAcceptanceRate*100 << "\n";
            outputFile << "Rotation Acceptance Rate : " << rotationAcceptanceRate*100 << "\n";
            if(ensemble == "isonpt" || ensemble == "anisonpt")
                outputFile << "Volume Move Acceptance Rate : " << volumeMoveAcceptanceRate*100 << "\n";
        }
    }

    // Energy at the end of equilibration run.
    erg = patchyGB.getEnergy();
    den = (nParticles*M_PI*3) / (4 * box.boxSize[0] * box.boxSize[1]); 
    // Print the energy at the end of equilibration.
    outputFile << std::fixed << std::setprecision(16);
    outputFile << "\nTotal crystal energy at end of equilibration : " << erg/nParticles << std::endl;
    outputFile << "\nInitial density of the configuration  : " << den << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Execute the production process.
    outputFile << "\nSystem Production\n#################" << std::endl;

    // Reset the statistics.
    singleParticleMove.reset();
    
    double nSamp = 0;
    double avgErg = 0;
    double avgDen = 0;
    
    // This is for cosmetic purpose, so that the step and energy columns are arranged neatly.
    int numDigits = std::log10(productRun) + 1;
    unsigned int avgPoint = 2*productRun/3;

    // Execute the simulation.
    for (unsigned int j=0;j<productRun;j++) {

        // Increment simulation by 1 Monte Carlo Sweeps.
        singleParticleMove += (nParticles+1);

        unsigned int step = j+1;
        
        // Report.
        if (step % fileInterval == 0){
            // Record the snapshots.
            io.appendXyzTrajectory(fileName, step, box, particles, patchTheta, aa, cells, maxInteractions, interactionEnergy, interactionRange, epsilonP, patchAngle, false);
            // Production system energy.
            erg = patchyGB.getEnergy();
            ergFile << std::fixed << std::setprecision(7) 
                << std::setw(numDigits) << std::setfill('0') << step << "\t" 
                << erg/nParticles << "\t";
            if(ensemble == "isonpt" || ensemble == "anisonpt"){
                den = (nParticles*M_PI*3) / (4 * box.boxSize[0] * box.boxSize[1]); 
                ergFile << den;
            }
            ergFile << std::endl;
            if( step > avgPoint) 
            {
                avgErg += erg/(nParticles);
                avgDen += den;
                nSamp++;
            }
        }    
        // Write out the restart file.
        if(step % productRun == 0){
            std::ostringstream fsufx;
            fsufx << (step);
            std::string fname = "restart_" + fsufx.str();
            io.restartfile(fname, box, particles);

        }  
    }

    attemptedTranslation = singleParticleMove.getTranslationAttempts();
    acceptedTranslation = singleParticleMove.getTranslationAccepts();
    attemptedRotation = singleParticleMove.getRotationAttempts();
    acceptedRotation = singleParticleMove.getRotationAccepts();
    maxTranslation = singleParticleMove.getMaxTranslation();
    maxRotation = singleParticleMove.getMaxRotation();
    translationAcceptanceRate = static_cast<double>(acceptedTranslation) / attemptedTranslation;
    rotationAcceptanceRate = static_cast<double>(acceptedRotation) / attemptedRotation;

    if(ensemble == "isonpt" || ensemble == "anisonpt"){
        attemptedVolumeMove = singleParticleMove.getVolumeMoveAttempts();
        acceptedVolumeMove = singleParticleMove.getVolumeMoveAccepts();
        maxVolumeMove = singleParticleMove.getMaxVolumeMove();
        volumeMoveAcceptanceRate = static_cast<double>(acceptedVolumeMove) / attemptedVolumeMove;
    }
    outputFile << "Attempted Translation : " << attemptedTranslation << "\n" << "Accepted Translation : " << acceptedTranslation << "\n";
    outputFile << "Attempted Rotation : " << attemptedRotation << "\n" << "Accepted Rotation : " << acceptedRotation << "\n";
    if(ensemble == "isonpt" || ensemble == "anisonpt")
        outputFile << "Attempted Volume Move : " << attemptedVolumeMove << "\n" << "Accepted Volume Move : " << acceptedVolumeMove << "\n";

    outputFile << std::fixed << std::setprecision(6);
    outputFile << "\nMaximum Translation : " << maxTranslation << "\n";
    outputFile << "Maximum Rotation : " << maxRotation << "\n";
    if(ensemble == "isonpt" || ensemble == "anisonpt")
        outputFile << "Maximum Volume Move : " << maxVolumeMove << "\n";
    
    outputFile << std::fixed << std::setprecision(6);
    outputFile << "\nTranslation Acceptance Rate : " << translationAcceptanceRate*100 << "\n";
    outputFile << "Rotation Acceptance Rate : " << rotationAcceptanceRate*100 << "\n";
    if(ensemble == "isonpt" || ensemble == "anisonpt")
        outputFile << "Volume Move Acceptance Rate : " << volumeMoveAcceptanceRate*100 << "\n";
    outputFile << "\nAverage Energy/particle : " << avgErg/nSamp << "\n";
    outputFile << "Average Density : " << avgDen/nSamp << "\n";
    outputFile << "\nComplete!\n";

    // Close the output files
    ergFile.close();
    outputFile.close();

    (void)productRun;(void)equilibRun; (void)dimension; (void)nParticles; 
    (void)pressure; (void)density; (void)ensemble; (void)useRestart; (void)rname;

    // We're done!
    return (EXIT_SUCCESS);
}
