#!/bin/bash
# usage:
# ./compile.sh source.cpp -r - Compile and run source.cpp
# Dependencias eigen

outFile=(${1//./ })
outFile=${outFile[0]}
g++ $1 -I ./eigen/ -o $outFile -O3 -std=c++17
echo copile command:
echo g++ $1 -I ./eigen/ -o $outFile -O3 -std=c++17
chmod +x $outFile
echo executable file:
echo $outFile
while test $# -gt 0; do
    case "$1" in
        -r|--run)
            echo run command:
            echo ./$outFile
            ./$outFile
            shift
            ;;
        -rv)
            echo run command:
            echo ./$outFile --verbose
            ./$outFile --verbose
            shift
            ;;
        
        *)
            shift
            ;;
    esac
done
