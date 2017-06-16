//Particle source utilising geant4 particle gun for protons at the Fermilab MTest facility.
//Takes into account beam path, user requests to control room, etc to determine proton momentum, arrival times, etc
//Each event is one proton
//Only designed to generate a single spill (proton times generated in ascending order until end of spill). Run multiple jobs and combine results with arts for multiple spill statistics.
//Tom Stuttard, thomas.stuttard.13@fnal.gov

#ifndef MTestProtonSource_hh
#define MTestProtonSource_hh

// Get the Geant4 code
#include "G4Event.hh" 
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"

// ROOT includes
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"

//C++ includes
#include <memory>
#include <vector>

//Use "override" keyword if using modern enough compiler
#if __cplusplus >= 201103L
#define OVERRIDE override
#else
#define OVERRIDE
#endif


//Define namespace
namespace mtest
{

  class MTestProtonPGA; //Forward declare service for friendship purposes

  //Define MTest proton source class
  class MTestProtonSource : public G4VPrimaryGenerator
  {
  public:
    MTestProtonSource();
    ~MTestProtonSource();
    
    void Initialize();
    void GeneratePrimaryVertex(G4Event*) OVERRIDE ;

    friend class MTestProtonPGA; //Define art service wrapper as friend so that it can directly write to private member params
    
    void SetDebug(bool v) { debug_ = v; }

    void EnablePlotting();


    //Getters
    unsigned int GetNumBucketsPerBatch() const;
    unsigned int GetNumBucketsPerTurn() const;
    unsigned int GetNumBucketsPerSpill() const;


  private:

    //G4ParticleGun* protonGun_;
    std::unique_ptr<G4ParticleGun> protonGun_;

    bool debug_;

    //Core accelerator parameters
    G4double bosTime_; //Start of spill time
    G4double rfBucketFrequency_;
    G4double bunchWidth_; //Time width of proton bunch within RF bucket
    G4double bunchOffset_; //Offset from start of RF bucket period to start of proton bunch
    unsigned int numBucketsPerBatch_;
    unsigned int numBatchesPerTurn_;
    unsigned int numMainInjectorTurnsPerSpill_;

    //Beam parameters at MTest
    G4ThreeVector beamDirection_;
    G4ThreeVector aperturePosition_;
    G4double beamNorm_;
    G4double beamMeanX_;
    G4double beamMeanY_;
    G4double beamRMSX_;
    G4double beamRMSY_;
    G4double beamNormWide_;
    G4double beamRMSXwide_;
    G4double beamRMSYwide_;
    G4double beamMeanXwide_;
    G4double beamMeanYwide_;
    G4double deltaP_;
    std::unique_ptr<TF2> beamProfilePDF_; //This will be used to create a beam profile PDF to sample from

    //Parameters defined by SeaQuest operation or by user request
    unsigned int firstFilledBucketInBatch_;
    unsigned int numFilledBucketsInBatch_;
    unsigned int firstFilledBatchInTurn_;
    unsigned int numFilledBatchesInTurn_;
    unsigned int numProtonsRequestedPerSpill_;
    G4double multipleOccupancyProb_; //Fraction
    std::vector<G4double> multipleOccupancyWeights_; //Relative weights of numbers of protons in multiple occupancy (first elem is 2 protons, second is 3, etc)
    G4double protonEnergy_;

    //Derived parameters
    G4double protonInTurnProb_;
    G4double rfBucketPeriod_;
    G4double batchPeriod_;
    G4double turnPeriod_;

    //State members
    bool isInitialized_;
    unsigned int currentTurn_;
    unsigned int currentSpill_;
    std::vector<G4double> protonTimesToGen_;

    //Debug plotting
    bool plot_;
    std::unique_ptr<TFile> rootFile_;
    std::map<std::string,TH1*> histoMap_;



  }; // End MTestProtonSource class

} //namespace end

#endif //MTestProtonSource_hh

