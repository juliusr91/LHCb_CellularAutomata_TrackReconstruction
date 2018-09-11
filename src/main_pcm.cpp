//
// Created by julius on 1-3-18.
//
//You just have to get the PCM repo from here:
//https://github.com/opcm/pcm
//for all the measurements to work


#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>
#include "../include/Tools.h"
#include "../include/computation.h"
#include <time.h>
#include <chrono>

#include <cpucounters.h>

#define NB_EVENTS 4

int main(int argc, char *argv[]){


    if (argc < 3){
        std::cerr<<"need more arguments; folder containing ./bin files and number of files to process"
                 <<std::endl;
        return 0;
    }

    std::string folderName {std::string(argv[1])};
    int fileNumber {std::stoi(argv[2])};

    if (folderName.empty()){
        std::cerr << "Either didn't specify the folder correctly or there is no such folder" << std::endl;
        return -1;
    }

    std::vector<std::vector<unsigned char>> input {readFolder(folderName, fileNumber)};


    //to measure time properly
    int numberOfIterations = 20;
    std::vector<std::chrono::duration<double>> times;
    for (int i =0; i < numberOfIterations; i++){

        ulong FLOP_count = 0;
        std::vector<int> multiplier ({1,1,2,4,4,8,8,16});

        //pcm custom events that can be used to retrieve different hardware counters.
        //pcm actually has the ability to retrieve ht FLOPs with an internal function.
        //TODO: replace the custom events with the internal function
        PCM::CustomCoreEventDescription events[NB_EVENTS];
//        double values[NB_EVENTS+1];
        // FP_ARITH_INST_RETIRED.SCALAR_DOUBLE
        events[0].event_number = 0xC7;
        events[0].umask_value  = 0x01;
        // FP_ARITH_INST_RETIRED.SCALAR_SINGLE
        events[1].event_number = 0xC7;
        events[1].umask_value  = 0x02;
//        // FP_ARITH_INST_RETIRED.128B_PACKED_DOUBLE
//        events[2].event_number = 0xC7;
//        events[2].umask_value  = 0x04;
        // FP_ARITH_INST_RETIRED.128B_PACKED_SINGLE
        events[2].event_number = 0xC7;
        events[2].umask_value  = 0x08;
        // FP_ARITH_INST_RETIRED.256B_PACKED_DOUBLE
//        events[2].event_number = 0xC7;
//        events[2].umask_value  = 0x10;
        // FP_ARITH_INST_RETIRED.256B_PACKED_SINGLE
        events[3].event_number = 0xC7;
        events[3].umask_value  = 0x20;




        PCM * m = PCM::getInstance();
            m->resetPMU();

        PCM::ErrorCode returnResult = m->program(PCM::CUSTOM_CORE_EVENTS, events);
        if (returnResult != PCM::Success){
            std::cerr << "Intel's PCM couldn't start" << std::endl;
            std::cerr << "Error code: " << returnResult << std::endl;
            exit(1);
        }

        SystemCounterState before_sstate = getSystemCounterState();


        auto start = std::chrono::system_clock::now();
        runCellularAutomaton(input);
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = end-start;
        std::cout << diff.count() <<std::endl;
        times.emplace_back(diff);


        SystemCounterState after_sstate = getSystemCounterState();

        for ( int i=0; i < NB_EVENTS; i++ ) {
            uint64 value = getNumberOfCustomEvents(i, before_sstate, after_sstate);
            printf("Event %0d: 0x%04x0x%04x: %lld\n", i+1, events[i].event_number, events[i].umask_value, value);
            FLOP_count += multiplier[i] * value;
        }


        std::cout << "Instructions per clock:" << getIPC(before_sstate,after_sstate) << std::endl;
        std::cout << "Bytes read:" << getBytesReadFromMC(before_sstate,after_sstate) <<std::endl;
        std::cout << "Bytes write:" << getBytesWrittenToMC(before_sstate,after_sstate) <<std::endl;
        std::cout << "Energy consumed:" << getConsumedEnergy(before_sstate,after_sstate) <<std::endl;
        std::cout << "L2 Hits: " << getL2CacheHits(before_sstate,after_sstate) <<std::endl;
        std::cout << "L2 Misses: " << getL2CacheMisses(before_sstate,after_sstate) <<std::endl;
        std::cout << "L2 Hit Ratio: " << getL2CacheHitRatio(before_sstate,after_sstate) <<std::endl;
        std::cout << "L3 Hits: " << getL3CacheHits(before_sstate,after_sstate) <<std::endl;
        std::cout << "L3 Misses: " << getL3CacheMisses(before_sstate,after_sstate) <<std::endl;
        std::cout << "L3 Hit Ratio: " << getL3CacheHitRatio(before_sstate,after_sstate) <<std::endl;
        std::cout << "Instr retired: " << getInstructionsRetired(before_sstate,after_sstate) <<std::endl;
//        std::cout << "Instr retired: " << getL3CacheHitRatio(before_sstate,after_sstate) <<std::endl;

        std::cout << "GFLOPs:" << FLOP_count / 1e9 <<std::endl;
        std::cout << "GFLOPs/sec:" << FLOP_count / 1e9 /diff.count()  <<std::endl;

        outfile << "Final Optimisation - New Physics Parameters"  << "," << getBytesReadFromMC(before_sstate,after_sstate)      << "," <<getBytesWrittenToMC(before_sstate,after_sstate)<< "," << getL2CacheHitRatio(before_sstate,after_sstate) << " ," << getL3CacheHitRatio(before_sstate,after_sstate) << std::endl;

        m->cleanup();
    }


    std::map<std::string, double> results = calcTimeResults(times);
    std::cout << (fileNumber) / results["mean"] << " events / s, "
              << results["mean"] << " s (st dev " << results["standardDeviation"] << "), "
              << results["mean"]/(fileNumber) << " seconds per event)"
              << std::endl;
}
