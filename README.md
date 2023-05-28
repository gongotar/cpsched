# cpsched
Checkpoint Scheduling for Shared Usage of Burst-Buffers in Supercomputers

This project is a fork from the [Scalable Checkpoint/Restart](https://github.com/LLNL/scr) library and enables it to coordinate 
checkpointing operations with a global scheduler. Checkpoints are only performed if the scheduler approves them, otherwise, 
the checkpoint requests are ignored and the application continues with the computation. As a result, checkpoint conflicts 
of different applications are optimised. This work is published as follows:

ICPP Workshops '18: Workshop Proceedings of the 47th International Conference on Parallel ProcessingAugust 2018Article No.: 44Pages 1–10https://doi.org/10.1145/3229710.3229755
