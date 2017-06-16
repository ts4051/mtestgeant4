# mtestgeant4

This is a Geant4 model of the Fermilab MTest proton beam. The beam structure is described in mtestgeant4/src/MTestProtonSource.cc

An example Geant4 project is included for testing.

## Dependencies

GEANT4 (tested with V9.6.4)
ROOT (tested with V6)
cmake (V3.1 or higher)

## Building the example

1) Clone the git repository

2) Setup a standard GEANT4 and ROOT environemnt (ensure G4INSTALL and ROOTSYS variables are set)

3) Create a build directory

4) cd to build directory

5) cmake path/to/mtestgeant4

6) make

This will build the executable 'example' in the build directory.

## Running the example

cd path/to/build/dir
./example example.mac

An output ROOT file (MTestProtonSource.root) contains some histograms than be used to verify the operation.