#!/bin/bash
# Cluster loading wrapper script
 
echo Loading job
qsub cTEBD.sh $*
qstat -ne -u johnson
