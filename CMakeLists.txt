cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(HDalitzEle LANGUAGES CXX)

# require C++17 in all targets(the same as the C++ version ROOT 6.24 used)
# show all of the warning information
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-O3 -Wall -Wextra -pedantic)

# find all of the shared libraries of ROOT in config file
# execute 'root-config --libs' in command line to check
# TMVA is also needed (not in the root-config)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(ROOT 6.24 CONFIG REQUIRED COMPONENTS TMVA XMLIO)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(fmt CONFIG REQUIRED) # txt format
find_package(yaml-cpp REQUIRED)
find_package(xgboost REQUIRED)

# create a shared libarary
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
add_library(HDalitzEle
    SHARED
        src/tdrstyle.cpp
        src/CMS_lumi_modification.cpp
        src/Utilities.cpp
        src/PUWeightCalculator.cpp
        src/EnCalibrater.cpp
        src/XGBReader.cpp
        src/GsfTracks.cpp
        src/WeightHandler.cpp
        src/MergedIDPred.cpp
        src/MergedEnCorrector.cpp
        src/ElectronSel.cpp
        src/PhotonSel.cpp
        src/PhoR9SSCorrection.cpp
        src/Categorizer.cpp
        src/Generator.cpp
)
target_include_directories(HDalitzEle PUBLIC include)
target_link_libraries(HDalitzEle
    PUBLIC
        ${ROOT_LIBRARIES}
        Boost::program_options
        xgboost::xgboost
        yaml-cpp
        fmt::fmt
        fmt::fmt-header-only
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
# script to perform the main HDalitz analysis
add_executable(runAnalysis src/runAnalysis.cpp)
target_link_libraries(runAnalysis
    PRIVATE
        HDalitzEle
)
# script to skim the signal tree
add_executable(skimTree src/skimTree.cpp)
target_link_libraries(skimTree
    PRIVATE
        HDalitzEle
)
# script to draw kinematics distributions
add_executable(drawKinematics src/drawKinematics.cpp)
target_link_libraries(drawKinematics
    PRIVATE
        HDalitzEle
)
# script to calculate the significance for each category
add_executable(runSignificance src/runSignificance.cpp)
target_link_libraries(runSignificance
    PRIVATE
        HDalitzEle
)



