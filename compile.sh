#!/bin/bash

# last argument is always sourcecode's filename to compile
# preceeding arguments will be passed to compiler
# output executable has same name as source w/o suffix

narg=${#@}
source=${@:$narg:$narg}
exe=${source%.*}
clopts=${@:1:$narg-1}

if [ ! -e "$source" ];
then
    echo "Usage:  compile.sh [compiler clopts] c-file"
    exit
fi

#compiler=gcc
compiler=g++
#compiler=gfortran

opts=' -W -Wall -Wshadow -Wstrict-aliasing'
#opts=' -lm'

rootflags=$(root-config --cflags --libs )
#echo $rootflags
#echo $opts
#echo $clopts

$compiler \
    $opts \
    $clopts \
    $rootflags \
    -I$HOME/.root \
    -o $exe $source

