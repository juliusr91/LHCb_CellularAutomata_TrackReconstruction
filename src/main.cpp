//
// Created by julius on 1-3-18.
//
#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <algorithm>
#include "../include/Tools.h"
#include "../include/computation.h"
#include <time.h>
#include <chrono>



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
    int numberOfIterations =5;

    std::vector<std::chrono::duration<double>> times;

    for (int i =0; i < numberOfIterations; i++){
        auto start = std::chrono::system_clock::now();
        runCellularAutomaton(input);
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = end-start;
        std::cout << diff.count() <<std::endl;
        times.emplace_back(diff);
    }


    std::map<std::string, double> results = calcTimeResults(times);
    std::cout << (fileNumber) / results["mean"] << " events / s, "
              << results["mean"] << " s (st dev " << results["standardDeviation"] << "), "
              << results["mean"]/(fileNumber) << " seconds per event)"
              << std::endl;
}
