#!/bin/bash

#MSUB -v id=2
#MSUB -N Ram_PaperCode_2
#MSUB -o /home/999999889/Ram_PaperCode/run_sim/Ram_PaperCode_2.out
#MSUB -e /home/999999889/Ram_PaperCode/run_sim/Ram_PaperCode_2.error

#MSUB -l walltime=60:00:00:00
#MSUB -m abe
#MSUB -M hkannan@ucsd.edu
#MSUB -l nodes=1:ppn=12
#MSUB -d /home/999999889/Ram_PaperCode/run_sim/
 
rm -rf /research/CNSM-Sun/Harish_crossfeeding/Ram_PaperCode/Ram_PaperCode_$id
mkdir -p /research/CNSM-Sun/Harish_crossfeeding/Ram_PaperCode/Ram_PaperCode_$id/Code

cp ./in_$id.txt /research/CNSM-Sun/Harish_crossfeeding/Ram_PaperCode/Ram_PaperCode_$id
cp ./submitJob.sh /research/CNSM-Sun/Harish_crossfeeding/Ram_PaperCode/Ram_PaperCode_$id
cp -r ../*.cpp /research/CNSM-Sun/Harish_crossfeeding/Ram_PaperCode/Ram_PaperCode_$id/Code
cp -r ../*.h /research/CNSM-Sun/Harish_crossfeeding/Ram_PaperCode/Ram_PaperCode_$id/Code

./Ram_PaperCodeA ./in_$id.txt 12 /research/CNSM-Sun/Harish_crossfeeding/Ram_PaperCode/Ram_PaperCode_$id
