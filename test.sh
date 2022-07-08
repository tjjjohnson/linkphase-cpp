#!/bin/bash


script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

set -euo pipefail

#echo "Compiling"
#g++ -o bin/linkphase src/linkphase.cpp -lboost_program_options

mkdir -p test/output/linkphase
cd test/output/linkphase

#cp ${script_dir}/resources/linkin.txt ./
#cp ${script_dir}/resources/linkphase_parameters.cfg ./
#${script_dir}/lic/unsorted.ped
#valgrind --tool=callgrind 
${script_dir}/bin/linkphase --pedigree-file ${script_dir}/lic/unsorted.ped.codes  \
  --genotype-file /data/small/thjoh0/dev/linkphase/test/output/GMK50k_chr29_10000000-11000000_old/gtFile.txt \
  --marker-file /data/small/thjoh0/dev/linkphase/test/output/GMK50k_chr29_10000000-11000000_old/markerFile.txt \
  --halfsib-phasing \
  --hmm-phasing \
  --templates 50 \
  --check-prephasing \
  --columns 

