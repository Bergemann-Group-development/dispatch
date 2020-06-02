#!/bin/bash

rm -rf $1                               # remove any existing folder
mkdir $1                                # create folder
mkfifo $1/pipe1 $1/pipe2                # create pipes
x=$PWD/../mesa_template/star            # path to executable
i=$PWD/../mesa_template/inlist          # path to input file
cd $1                                   # change directory to the folder
$x $i > log.txt &                       # start MESA
