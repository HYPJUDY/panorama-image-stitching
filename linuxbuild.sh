#!/bin/bash
# put all vlfeat .o files in 'objs' folder
g++ -O3 -c stitching_main.cpp
g++ stitching_main.o objs/* -lm -fopenmp -lX11 -o panorama
