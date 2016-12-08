#!/bin/bash

#==========================================================
#
#         Configuration file for
#
# OCEAN MONITORING for NEMO v3.6 ORCA1 on 46 levels and LIM2
#
#            HPC: voima.fmi.fi 
#
#        L. Brodeau, 2015
#
#===========================================================

. ./configs/config_ORCA1_L46_v36_voima.sh

# LIM2 exceptions
export NN_ICEF="ice_pres" ; # name of ice fraction in "FILE_ICE_SUFFIX" file...
export NN_ICET="sit" ; # ice thickness but 'sit' is only in icemod file !!!
export JTITLE="NEMO v3.6 ${CONF} (L${NBL}) - LIM2 / ocean-only experiment"

