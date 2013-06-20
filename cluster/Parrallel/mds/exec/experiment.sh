#!/bin/bash
#$ -j y
#$ -l h_rt=1:0:0
#$ -wd /home/Regimantas/errors 
#S -P WWW
#$ -pe orte 1
/opt/openmpi/bin/mpirun -np 1 /home/Regimantas/cpp/mds/mds.out /home/Regimantas/www/experimental_data/newexp350/data.txt 100 2 /home/Regimantas/www/experimental_data/newexp350/viz000/viz_
