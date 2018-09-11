#also includes the make options for the intel PCM tool.
#You just have to get the PCM repo from here:
#https://github.com/opcm/pcm
#and adjust the folders below.

cmake_minimum_required(VERSION 3.9.6)

project(CellularAutomata)

set(CMAKE_CXX_STANDARD 11)

#TODO: use msse4 onwards
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -std=c++11 -Wall -fopenmp -lpthread -mavx -msse3")

#INCLUDE_DIRECTORIES(/home/julius/Projects/pcm)
#LINK_DIRECTORIES(/home/julius/Projects/pcm)

#set(PCM_STUFF /home/julius/Projects/pcm/ /home/julius/Projects/pcm/cpucounters.o /home/julius/Projects/pcm/pci.o /home/julius/Projects/pcm/msr.o /home/julius/Projects/pcm/client_bw.o /home/julius/Projects/pcm/cpucounters.o)

add_executable(CellularAutomata src/main.cpp src/Tools.cpp include/Common.h src/runCellularAutomaton.cpp include/computation.h)
#add_executable(CellularAutomataRoofline ${PCM_STUFF} src/main_pcm.cpp src/Tools.cpp include/Common.h src/runCellularAutomaton.cpp include/computation.h )
