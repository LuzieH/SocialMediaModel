we need to run julia on a node with

-t 30-0
-c 64
to have a time limit of 30 days and 64 cores

this is possible on the big partition, i.e.
--partition = big

to run this from a bash shell use

srun  --partition = big -t 30-0 -c 64 --pty bash


what is the maximum time limit I can ask for? -t 30-0?
probably


remember to set the number of blas threads accordingly
anything else?

probably writing an sbatch file incorporating all these tips is worthwile.

srun  --partition=big -t 30-0 -c 64 -C Gold6338 --pty bash