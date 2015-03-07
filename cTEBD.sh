#!/bin/bash
# Execution wrapper script
# Combine stdout and sterr into one output file:
#$ -j y
# Use "bash" shell:
#$ -S /bin/bash
# Change the output directory for stdout:
#$ -o sge_output
# Name the job:
#$ -N cTEBD
# Use current directory as working root:
#$ -cwd
# Set default memory request:
#$ -l h_vmem=3G
#$ -l medium
# Send mail: n=none, b=begin, e=end, a=abort, s=suspend
#$ -m a 
#$ -M tomi.johnson@physics.ox.ac.uk
 
/share/apps/MATLAB/run_mcc2012b.sh /home/johnson/WorkDoneOnIsingModel/TEBD/cTEBD $*

