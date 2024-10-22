PROJECT(SPECT)
cmake_minimum_required(VERSION 2.8)
# add the folder containing the dlib folder to the include path
#INCLUDE_DIRECTORIES(./lib)
#include_directories(./src)


ADD_EXECUTABLE(SPECT main.cpp bgeot_ftool.cpp acqphy.cpp image.cpp mlem.cpp)

set(CMAKE_CXX_FLAGS "-O3 -std=c++11")
find_package(OpenMP)
if(OPENMP_FOUND)
  message("OpenMP is found.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
set(OMP_NUM_THREADS 4)

