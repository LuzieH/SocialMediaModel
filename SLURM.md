we need to run julia on a node with

-t 1-0
-c 40
to have a time limit of 1 day and 40 cores

this is possible on the big partition, i.e.
-p big

to run this from a bash shell use

srun -p big -t 1-0 -c 40 --pty bash


what is the maximum time limit I can ask for? -t 30-0?
probably


remember to set the number of blas threads accordingly
anything else?

probably writing an sbatch file incorporating all these tips is worthwile.

