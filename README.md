# CellularAutomaton

This is a standalone track reconstruction algorithm for the LHCb experiment based on a Cellular Automaton approach.
The master repo contains everything needed, including 20 input files in the correct data structure.

You can use cmake and make to compile the repo; you'll need GCC 5.4+.

1. cmake CMakeLists.txt
2. make
3. run using: *./CellularAutomaton pathToFolderContaintingData #FilesToBeProcessed*

The thesis to the repo can be found [here](http://www.scriptiesonline.uba.uva.nl/scriptie/652350).


## PCM
We used Processor Counter Monitor (PCM) to read several of the hardware counters.
Its possible to add hardware counters for performance indicators such as FLOPs, Cache hit ratio, bytes read and written.
You'll need the PCM repo from [here](https://github.com/opcm/pcm).

Additionally, you'll need to uncomment 4 lines in the CMakeLists.txt file.

## Project structure
The src folder contains:
1. main.cpp
2. main_pcm.cpp (main file incase you want to use PCM)
3. runCellularAutomaton.cpp (the actual CA algorithm)
4. Tools.cpp (several functions for reading the input data by [Daniel Campora](https://gitlab.cern.ch/dcampora/search_by_triplet).

The input folder contains 20 simulated event files.

The include folder contains:
1. Common.h (defines things like the Event, doublet, track object)
2. computation.h (header for runCellularAutomaton)
3. Tools.h

The FLOPS folder contains a file detailing the umask of the used FLOP hardware counter.