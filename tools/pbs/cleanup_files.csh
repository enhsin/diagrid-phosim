#!/bin/csh
#
# Usage: cleanup_files.csh compute_node_id local_scratch_dir
#   -- Runs "rm -r" on directory "local_scratch_dir" on Athena 
# compute node "compute_node_id"
#   -- If no command line arguments are given, arguments are taken from 
# the following environment variables:
#         CLEAN_MASTER_NODE_ID    (sets compute_node_id)
#         CLEAN_LOCAL_SCRATCH_DIR (sets local_scratch_dir)
#
#   Example: "cleanupfiles.csh compute-5-1 /state/partition1/username/tmp.12345"
#
# The following lines are useful when this script is submitted to PBS.  
# Otherwise,they are ignored.

#PBS -j oe
#PBS -l walltime=1:00,nodes=1:ppn=1
#PBS -m abe

# Get arguments from either the command-line or environment variables
if ($#argv == 2) then
  set master_node_id = $1
  set local_scratch_dir = $2
  echo Cleanup script received these arguments from the command line:
else
  set master_node_id = $CLEAN_MASTER_NODE_ID
  set local_scratch_dir = $CLEAN_LOCAL_SCRATCH_DIR
  echo Cleanup script received these arguments from environment variables:
endif
 echo "   master_node_id:    $master_node_id"
 echo "   local_scratch_dir: $local_scratch_dir"

# Comment the following command when you are finished testing...
#ssh $master_node_id "ls -al $local_scratch_dir"
# ...and uncomment the following line for actual usage
ssh $master_node_id "rm -r $local_scratch_dir"
echo Cleanup script complete!
