//
// Created by julius on 4-3-18.
//
#include <vector>
#include <stdint.h>
#include <string>
#include <iostream>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <map>

std::vector<std::vector<uint8_t>> readFolder(
        const std::string & folderName,
        int fileNumber
);


void getFolderContents(
        std::vector<std::string>& folderContents,
        const std::string & folderName
);

void writeBinaryTrack();

std::map<std::string, double> calcTimeResults(std::vector<std::chrono::duration<double>> &times);