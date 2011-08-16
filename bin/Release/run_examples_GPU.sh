#!/bin/sh

SAMPLES_DIR=../../data/3D/example_tests
#RUN="mpiexec.hydra -f myhosts -n 3 ./FluidSolver3D"
RUN=./FluidSolver3D

box_pipe_2D_data_txt=$SAMPLES_DIR/box_pipe/box_pipe_2D_data.txt
box_pipe_2D_config_txt=$SAMPLES_DIR/box_pipe/box_pipe_2D_config.txt
white_sea_config_txt=$SAMPLES_DIR/white_sea/white_sea_config.txt

ext=".tmp"

for file in $box_pipe_2D_data_txt $box_pipe_2D_config_txt $white_sea_config_txt
do
	awk '{ sub("\r$", ""); print }' $file > $file$ext
done

$RUN $box_pipe_2D_data_txt$ext box_pipe_example $box_pipe_2D_config_txt$ext align GPU
#transpose
$RUN $SAMPLES_DIR/white_sea/white_sea_data.nc white_sea_example $white_sea_config_txt$ext align GPU
#transpose  


