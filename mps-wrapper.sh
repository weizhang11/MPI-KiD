#!/bin/bash

# select_gpu_device wrapper script

export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log

{
   echo quit | nvidia-cuda-mps-control
} || {
   echo "Failed to quit nvidia-mps daemon"
}

# Launch MPS daemon for local ID 0,1,2,3
if [ $SLURM_LOCALID -lt 4 ]; then
   nvidia-cuda-mps-control -d
fi

# Round-Robin way
export CUDA_VISIBLE_DEVICES=$(( SLURM_LOCALID % 4 ))

# Wait for MPS to start
sleep 5

# Run the command
"$@"

# Quit MPS control daemon before exiting
if [ $SLURM_LOCALID -lt 4 ]; then
   echo quit | nvidia-cuda-mps-control
fi
