//
// All functions by Daniel Campora (linked below)
// minor fix on line 30 to sort the input files.
// uses same functions and event data structure as:
// https://gitlab.cern.ch/dcampora/search_by_triplet


#include "../include/Tools.h"
#include "../include/Common.h"


void getFolderContents(std::vector<std::string> & folderContents, const std::string & folderName){
    DIR *directory {opendir(folderName.c_str())};
    struct dirent *entry;

    if (directory != NULL){
        while ((entry = readdir(directory))!=NULL){
            std::string filename {entry->d_name};
            if (filename.find(".dat") != std::string::npos || filename.find(".bin") != std::string::npos){
                folderContents.push_back(filename);
            }
        }
        closedir(directory);
        if (folderContents.size() == 0) {
            std::cerr << "No binary files found in the folder." << std::endl;
            exit(-1);
        }else{
            //for Test files by daniel
//            std::sort(folderContents.begin(), folderContents.end(), sortFileNames());

            //for the new 10K files
            std::sort(folderContents.begin(), folderContents.end());
//            for (int i = 0; i < folderContents.size(); ++i) {
//                std::cout << folderContents[i] <<std::endl;
//            }
            std::cout << "Found " << folderContents.size() << " binary files." << std::endl;
        }
    } else {
        closedir(directory);
        std::cerr << "Folder couldn't be opened" << std::endl;
        exit(-1);
    }
}

void readFileIntoVector(std::vector<uint8_t>& output, const std::string & currentFile){
    std::ifstream infile(currentFile.c_str(), std::ifstream::binary | std::ios::ate);

//    infile.seekg(0, std::ios::end);
    auto end = infile.tellg();
    infile.seekg(0, std::ios::beg);
    auto dataSize = end - infile.tellg();

    output.resize(dataSize);
    infile.read((char*) &output[0], dataSize);
    infile.close();
}


std::vector<std::vector<uint8_t>> readFolder(const std::string & folderName, int fileNumber){
    std::vector<std::string> folderContents;
    getFolderContents(folderContents, folderName);

    std::cout << "Requested " << fileNumber << " files"<< std::endl;
    std::vector<std::vector<uint8_t >> input;

    int readFiles {0};

    for (int i=0; i <fileNumber; i++){
        std::string currentFile {folderContents[i % folderContents.size()]};

        std::vector<uint8_t > fileContents;
        readFileIntoVector(fileContents, folderName+ "/" + currentFile);

        auto eventInfo = EventInfo(fileContents);
        if (eventInfo.numberOfModules == NUMBER_OF_SENSORS){
            fileContents.resize(eventInfo.size);
            input.push_back(fileContents);
        }

        readFiles++;
        if ((readFiles % 100) == 0){
            std::cout << "." << std::flush;
        }
    }

    std::cout << std::endl << input.size() << " files read" << std::endl << std::endl;

    return  input;
}


std::map<std::string, double> calcTimeResults(std::vector<std::chrono::duration<double>> &times){
    std::map<std::string, double> results;
    double standardDeviation = 0.0, variance = 0.0, mean = 0.0;

    for (int i = 0; i < times.size(); i++){
        mean += times[i].count();
        variance += times[i].count() * times[i].count();
    }

    mean /= times.size();
    variance = (variance/times.size()) - (mean * mean);
    standardDeviation = std::sqrt(variance);
    results["variance"] = variance;
    results["mean"] = mean;
    results["standardDeviation"] = standardDeviation;

    return results;
};
