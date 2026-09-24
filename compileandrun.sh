#!/bin/sh
#SBATCH -J isTm1.35
#SBATCH -N 1
#SBATCH --account=p70623
#SBATCH --qos=p70623_0512
#SBATCH --time=120:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jenis.thongam@tuwien.ac.at


# module purge
# module load intel/19.1.3.304

folder_name=temp
ensemble="anisonpt" # or anisonpt or isonpt.
SCRIPT_DIR="$(pwd)"

while read value; do
    mkdir "${folder_name}_${value}"

    cd source_code

    g++ -Wall -O2 -march=native -std=c++11 Box.cpp CellList.cpp InputOutput.cpp Model.cpp \
    Particle.cpp Initialise.cpp PatchyGayBerne.cpp SingleParticleMove.cpp patchy_gayberne.cpp -lm -o ${folder_name}_${value}.out
    
    mv ./${folder_name}_${value}.out ./../${folder_name}_${value}
    cd ./../${folder_name}_${value}
    time -p ./${folder_name}_${value}.out       \
    --ensemble ${ensemble}                      \
    --temperature "${value}"                    \
    --pressure 0                                \
    --useRestart 1                              \
    --rname "${SCRIPT_DIR}/restart_press_500"   \
    1>params 2>err 
    cd ..
    
    cp ./${folder_name}_${value}/restart_* ./restart_press_500

done < statePoints.txt
wait