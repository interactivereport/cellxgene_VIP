#!/bin/bash
# source ~/bin/VIP.sh to export the environment to the current shell

# module load Arrow/0.17.1-fosscuda-2020b

source /edgehpc/apps/gb/anaconda3/4.9.2/etc/profile.d/conda.sh
conda activate VIP
cd ~/tool/cellxgene_VIP

export PYTHONNOUSERSITE=1
unset R_LIBS_USER

#export LIBARROW_MINIMAL=false
