#include "MTestProtonSource.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4UnitsTable.hh"

#include "TF1.h"

//#include <algorithm>
#include <numeric> //for std::accumulate

using CLHEP::RandGeneral;

// Constructor
mtest::MTestProtonSource::MTestProtonSource()

  : protonGun_()
  , debug_(false)

  //Beam structure
  , rfBucketFrequency_(53.1*megahertz)
  , bunchWidth_(1.5*ns) //1.5ns proton bunch width within 19ns RF bucket is typical
  , bunchOffset_(0.) //Kind of arbitrary since haven't defined whether start is peak or zero
  , numBucketsPerBatch_(84)
  , numBatchesPerTurn_(7)
  , numMainInjectorTurnsPerSpill_(360000)

  //MTest beam hall
  , beamDirection_(0,1.,0) //Unit vector pointing in direction of beam
  , aperturePosition_(0*mm,0*mm,0*mm) //Location for particles to start from

  //Beam profile
  , beamNorm_(1.)
  , beamMeanX_(0.*mm)
  , beamMeanY_(0.*mm)
  , beamRMSX_(2.2*mm)
  , beamRMSY_(3.5*mm)
  , beamNormWide_(0.00103)
  , beamMeanXwide_(beamMeanX_)
  , beamMeanYwide_(beamMeanY_)
  , beamRMSXwide_(24.1*mm)
  , beamRMSYwide_(25.2*mm)
  , beamProfilePDF_()

  //Beam energy/momentum
  , protonEnergy_(120.0*GeV)
  , deltaP_(0.02) //Fractional momentum spread, 2%

  //Accelerator operations
  , firstFilledBucketInBatch_(0)
  , numFilledBucketsInBatch_(82)
  , firstFilledBatchInTurn_(0)
  , numFilledBatchesInTurn_(6)
  , numProtonsRequestedPerSpill_(50000)

  //Multiple occupancy
  , multipleOccupancyProb_(0.35) //The probability of an RF bucket event at MTest being multiply occupied (e.g. prob. of multiple protons at same time)
  , multipleOccupancyWeights_({1.}) //Weight of each N num particles in multiple occupancy event, e.g. {3.,1.} means 2 protons 75% of the time, 3 protons the other 25%

  //Internal states
  , protonInTurnProb_(0.)
  , rfBucketPeriod_(0.)
  , batchPeriod_(0.)
  , turnPeriod_(0.)
  , isInitialized_(false)
  , currentTurn_(0)
  , currentSpill_(0)
  , protonTimesToGen_()

  //Plotting
  , plot_(false)
  , rootFile_()
  , histoMap_()
{

  //Instantiate gun
  protonGun_.reset( new G4ParticleGun() );

}

mtest::MTestProtonSource::~MTestProtonSource() 
{ 

  //
  // Write output ROOT plots file
  //

  if(plot_) {
    rootFile_->Write();
    rootFile_->Close();
    G4cout << "Dumped MTest proton source debugging plots to ROOT file : " << rootFile_->GetName() << G4endl;
  }


}


void mtest::MTestProtonSource::EnablePlotting() {

  //Check initialization completed successfully
  if(!isInitialized_) {
    //TODO throw
    std::cerr << "MTestProtonSource : Cannot enable plotting, class not initialized" << std::endl;
    return;
  }

  plot_ = true;


  //
  // Create ROOT file
  //

  rootFile_.reset( new TFile("MTestProtonSource.root","RECREATE") );


  //
  // Book histos
  //

  //These histos show data that can't be accessed from the art records, so i useful for debugging

  std::string name = "";

  //Energy spread
  name = "h_E"; 
  double energy5Sigma = 5. *  deltaP_ * protonEnergy_; //5 sigma range
  histoMap_[name] = new TH1F(name.c_str(), ";Proton E [GeV]", 
    100, (protonEnergy_-energy5Sigma)/GeV, (protonEnergy_+energy5Sigma)/GeV ); 


  //MI turn
  name = "h_MITurnTime"; 
  histoMap_[name] = new TH1F(name.c_str(), ";Time MI turn started for proton [ms]", 
    numMainInjectorTurnsPerSpill_, bosTime_/ms, (turnPeriod_*numMainInjectorTurnsPerSpill_)/ms );

  name = "h_MITurnNumber"; 
  histoMap_[name] = new TH1F(name.c_str(), ";MI turn counter for proton", 
    numMainInjectorTurnsPerSpill_/100, 0., numMainInjectorTurnsPerSpill_);

  name = "h_protonTimeInTurn"; 
  histoMap_[name] = new TH1F(name.c_str(), ";Time of proton relative to start of MI turn [#mus]", 
    this->GetNumBucketsPerTurn(), 0./microsecond, turnPeriod_/microsecond );


  //Booster batch
  name = "h_boosterBatch"; 
  histoMap_[name] = new TH1F(name.c_str(), ";Booster batch number (within MI turn) of proton", 
    numBatchesPerTurn_, -0.5, numBatchesPerTurn_-0.5);

  name = "h_protonTimeInBatch"; 
  histoMap_[name] = new TH1F(name.c_str(), ";Time of proton relative to start of booster batch [ns]", this->GetNumBucketsPerBatch(), 0/ns, batchPeriod_/ns );


  //RF buckets and bunches
  name = "h_bucket"; 
  histoMap_[name] = new TH1F(name.c_str(), ";RF bucket number (within booster batch) of proton", 
    numBucketsPerBatch_, -0.5, numBucketsPerBatch_-0.5);

  name = "h_protonTimeInBucket"; 
  histoMap_[name] = new TH1F(name.c_str(), ";Time of proton relative to start of RF bucket", 100, 0, rfBucketPeriod_);

  name = "h_protonTimeInBunch"; 
  histoMap_[name] = new TH1F(name.c_str(), ";Time of proton relative to start of proton bunch", 100, 0, bunchWidth_);

  name = "h_bucketOccupancy"; 
  histoMap_[name] = new TH1F(name.c_str(), ";RF bucket occupancy at MTest", 10, -0.5, 9.5);


  //Spill
  name = "h_spill"; 
  histoMap_[name] = new TH1F(name.c_str(), ";Spill number of proton", 
    100, -0.5, 99.5);

  name = "h_protonTimeInSpill"; 
  auto spillDuration = this->GetNumBucketsPerSpill() * rfBucketPeriod_;
  histoMap_[name] = new TH1F(name.c_str(), ";Time of proton relative to start of spill [ms]", 
    1000, bosTime_/ms, (bosTime_+spillDuration)/ms );


  //Beam profile PDF
  name = "h_beamProfilePDFX"; 
  histoMap_[name] = new TH1F(name.c_str(), "Input beam profile PDF in x (beam frame, horizontal transverse);x [mm]", 
    100, -100., 100.);

  name = "h_beamProfilePDFY"; 
  histoMap_[name] = new TH1F(name.c_str(), "Input beam profile PDF in y (beam frame, vertical transverse);y [mm]", 
    100, -100., 100.);

  name = "h_beamProfilePDFXY"; 
  histoMap_[name] = new TH2F(name.c_str(), "Input beam profile PDF in x,y (beam frame, transverse);y [mm];x [mm]", 
    100, -100., 100., 100, -100., 100.);
  histoMap_[name]->SetOption("COLZ");


  //Generated beam profile
  name = "h_beamProfileX"; 
  histoMap_[name] = new TH1F(name.c_str(), "Generated beam profile in x (beam frame, horizontal transverse);x [mm]", 
    100, -100., 100.);

  name = "h_beamProfileY"; 
  histoMap_[name] = new TH1F(name.c_str(), "Generated beam profile in y (beam frame, vertical transverse);y [mm]", 
    100, -100., 100.);

  name = "h_beamProfileXY"; 
  histoMap_[name] = new TH2F(name.c_str(), "Generated beam profile in x,y (beam frame, transverse);y [mm];x [mm]", 
    100, -100., 100., 100, -100., 100.);
  histoMap_[name]->SetOption("COLZ");



  //
  // Fill beam PDF histo
  //

  //Plot the input beam profile PDF (sample a random slection of values from it)
  if(plot_) {
    for( unsigned long i = 0 ; i < 1000000l ; i++ ) {
      double xRand, yRand;
      beamProfilePDF_->GetRandom2(xRand,yRand);
      histoMap_.find("h_beamProfilePDFX")->second->Fill( xRand/mm ); 
      histoMap_.find("h_beamProfilePDFY")->second->Fill( yRand/mm ); 
      static_cast<TH2F*>(histoMap_.find("h_beamProfilePDFXY")->second)->Fill( xRand/mm , yRand/mm ); 
    }
  }

}


void mtest::MTestProtonSource::Initialize() {


  //
  // Set derived params from the user specified ones
  //

  //This computation is done here rather than in GeneratePrimaryVertex to avoid repetition (e.g. for performance)

  //Determine probability of getting proton during a single MI turn. The probabolity is uniform throughout 
  //the spill and tuned by the accelerator team to give the requested number of protons during the spill)
  protonInTurnProb_ = (double)numProtonsRequestedPerSpill_ / (double)numMainInjectorTurnsPerSpill_;

  //Determine time period of one RF bucket
  rfBucketPeriod_ = 1. / rfBucketFrequency_ ;

  //Determine time period of one booster batch
  batchPeriod_ = rfBucketPeriod_ * numBucketsPerBatch_;

  //Determine time period of one MI turn
  turnPeriod_ = rfBucketPeriod_ * numBucketsPerBatch_ * numBatchesPerTurn_;


  //
  // Check user defined params
  //

  //Check did not request more protons that there are MI turns
  if( numProtonsRequestedPerSpill_ > numMainInjectorTurnsPerSpill_ ) {
    std::cerr << "MTestProtonSource : This code does not support requesting more protons per spill than there are turns\n";
    return;
  }

  //Check cannot get unphysical number of buckets
  if( (firstFilledBucketInBatch_+numFilledBucketsInBatch_) >= numBucketsPerBatch_ ) {
    std::cerr << "MTestProtonSource : Error in determining numbers of filled buckets per batch\n";
    return;
  }

  //Check cannot get unphysical number of batches
  if( (firstFilledBatchInTurn_+numFilledBatchesInTurn_) >= numBatchesPerTurn_ ) {
    std::cerr << "MTestProtonSource : Error in determining numbers of filled batches per turn\n";
    return;
  }

  //Check not requesting more protons than there are main innjector turns. This is not a use case for us, and code 
  //below isn't designed to handle it.
  if( numProtonsRequestedPerSpill_ > numMainInjectorTurnsPerSpill_ ) {
    std::cerr << "MTestProtonSource : Error: Requested more protons per spill than there are main injector turns\n";
    return;
  }

  //Check bunch width and offset are compatible with RF bucket width
  if( (bunchOffset_ < 0.) || ( bunchOffset_+bunchWidth_ > rfBucketPeriod_ ) ) {
    std::cerr << "MTestProtonSource : Error in bunch definition compared to bucket definition\n";
    return;
  }

  //Define transverse spatial profile PDF
  beamProfilePDF_.reset( new TF2("distribution","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4]) + [5]*TMath::Gaus(x,[6],[7])*TMath::Gaus(y,[8],[9])",-100,100,-100,100) );
  beamProfilePDF_->SetParameters(beamNorm_, beamMeanX_, beamRMSX_, beamMeanY_, beamRMSY_, beamNormWide_, beamMeanXwide_, beamRMSXwide_, beamMeanYwide_, beamRMSYwide_);
  beamProfilePDF_->SetNpx(100); //1000
  beamProfilePDF_->SetNpy(100); //1000


  //
  // Normalize multiple occupancy # proton weights to a probability
  //

  double sumWeights = std::accumulate(multipleOccupancyWeights_.begin(),multipleOccupancyWeights_.end(),(double)0.);//3rd arg gives type
  for( auto& weight : multipleOccupancyWeights_ ) weight /= sumWeights;


  //
  // Register initialization complete
  //

  isInitialized_ = true;

}


void mtest::MTestProtonSource::GeneratePrimaryVertex(G4Event* evt) {


  //Check initialization completed successfully
  if(!isInitialized_) {
    //TODO Move to better error mechanism, whilst maintaining stand-alone capabilities (e.g. no art)
    std::cerr << "MTestProtonSource : Cannot generate particle, class not initialized" << std::endl;
    return;
  }


  //
  // Create proton
  //

  //Create proton with MTest parameters (at aperture to beam hall)...

  //Get proton definition
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  protonGun_->SetParticleDefinition(particle);


  //
  // Set proton energy
  //

  //Use Gaussian spread about target energy to represent momentum spread of beam (dp/p = 2% from MTest webpage)
  double sigma = deltaP_*protonEnergy_;
  G4double generatedEnergy = G4RandGauss::shoot( protonEnergy_, sigma ); //TODO Truncate at 3sigma (common accerlator modelling practise)?
  protonGun_->SetParticleEnergy( generatedEnergy );
  if(plot_) histoMap_.find("h_E")->second->Fill( (generatedEnergy)/GeV ); 


  //
  // Set proton direction and position
  //

  //Generate proton transverse position relative to beam mean within beam local coords (where +z is beam direction, x is horizintal and y is vertical )
  double xRand, yRand;
  beamProfilePDF_->GetRandom2(xRand,yRand);
  G4ThreeVector protonPositionBeamFrame( xRand*mm , yRand*mm, 0.*mm );
 
  //Now transform to world frame. Currently have particle at z=0 (beam dir = +z), y,x = generated transverse position. 
  //(e.g. in beam frame). Need to rotate and translate this to  put the particle at the beam aperture in MTest, and 
  //pointing in the specified direction...

  //Start from beam frame
  G4ThreeVector protonPosWorldFrame(protonPositionBeamFrame);

  //Rotate to point in specified direction in parent (world) frame
  protonPosWorldFrame.rotateUz(beamDirection_.unit());

  //Translate proton position to 3D initial position in world (e.g. beam aperture)
  protonPosWorldFrame += aperturePosition_;

  //Write position and direction to particle
  protonGun_->SetParticlePosition(protonPosWorldFrame);
  protonGun_->SetParticleMomentumDirection(beamDirection_);


  //
  // Set proton arrival time
  //

  /* Beam path is as follows... (see http://www-ppd.fnal.gov/MTBF-w/beam/delivery.html)
    
     - Protons are accelerated in 53MHz RF buckets. The protons within one bucket are called a bunch (typically
       smaller width that period of one RF bucket, as want to move bunch within bucket depending on whether you
       are accelerating or storing particles).
     - Protons from linac enter the booster (circular), where one turn of particles consists of 84 RF buckets 
       (1.6us total width), and is called a boost batch. Not all RF buckets are filled with protons however 
       (depends on accelerator operations, 82 is common). The particles are accelerated to 8GeV.
     - Booster batches enter main injector (MI), which holds exactly 7 consecutive booster batches (11.2us total
       width, e.g. one MI turn). 
     - Protons are resonantly extracted from the MI into a beam line that serves SeaQuest and MTest (at least in
       2015). The number of booster batches this extraction is drawn from depends on SeaQuest operations, in 2015 
       (full SeaQuest running) this will be 6 of the 7 batches in the MI.
     - Exactly what is extracted from these batches depends on SeaQuest operations, but during full SeaQuest 
       running typically protons are extracted from ~approx. all filled RF buckets in the used batches (so 6 of 7),
       and then a small amount of these are filtered off for MTest (the rest goes to SeaQuest). 
     - The protons for MTest are then passed through a pin hole collimator before the MTest beam hall where most 
       are lost, resulting in an approx. 2mm Gaussian beam at MTest.
     - This extraction to MTest occurs during a slow spill, which is for ~4s every 60s. More precisely, the slow 
       spill is ~360,000 MI turns (each turn is ~11.2us). Between slows spills no particles arrive at MTest.
     - The number of protons that reach MTest can be controlled by the accelerator team. The MTest user can request 
       this number from the control room, with 50,000-150,000 per spill being typical values.
     - The number of protons per spill is less than the number of MI turns, so the user should not expect a proton 
       every turn. The acclerator team control the time distribution of these protons to MTest such that their arrival
       times are approximately unfiform throughout the spill. They also try to only deliver a single proton in a given
       turn to MTest. However the difficulties in extracting such a small amount of beam and the uncertainty involved 
       with using a collimator mean that there is an ~35% probability of a multiple occupancy event where more than 
       one (normally 2) protons arrive during a single RF bucket (probability increases at large intensities). 
     - Note that these double/multiple occupancy protons could be from any filled bucket in the turn (e.g. they are 
       not doubly extracted from the same bucket or batch necessarily or anything like that).
     - Note also that there is a learn-in time for the extraction process each time the accelerator restarts. This
       affects the numbers mentioned above until stable.
     - Due to the aforementioned issues with extracting the small amount of beam, there may be significant variation 
       in the number of particles arriving at MTest relative to the number requested by the user. Note that the number
       arriving at MTest is measured by a counter somewhere in MTest.
   */

  //Note that this function is called once per event, and only want one proton per event. Generation scheme is to 
  //increment MI turns until generate a proton in one, and then generates arrival times for all protons that turn 
  //(e.g. includes double occupancy). Each new event uses one of these times, and when all have been used a new set 
  //is generated for the next cycle where protons are generated.

  //If no remaining protons times from last event, generate new one(s)
  if( protonTimesToGen_.size() == 0 ) {

    //No buffered times, generate new one(s)...

    //Generate next MI turn a proton arrives at MTest during
    //Increment turn number, checking each time whetherto generate a particle, and if not repeat 
    while( G4RandFlat::shoot(0.,1.) > protonInTurnProb_ ) {
      currentTurn_++;
    }

    //Generate arrival RF bucket at MTest for this event proton(s) (multiple occupancy drawm from same bucket)
    //Uniformly randomly sampling from filled bunches and batches...

    //Get start time of current turn
    G4double startOfTurnTime = bosTime_ + (currentTurn_ * turnPeriod_);
      
    //Generate booster batch p+ is in and start time of this batch within turn
    unsigned int batch = firstFilledBatchInTurn_ + floor( G4RandFlat::shoot(0.,numFilledBatchesInTurn_) );
    double batchTimeInTurn = batch * batchPeriod_;

    //Generate bucket/bunch (a bunch is the protons within the RF bucket) p+ is in (within batch) 
    //and start times of these within batch
    unsigned int bucket = firstFilledBucketInBatch_ + floor( G4RandFlat::shoot(0.,numFilledBucketsInBatch_) );
    double bucketTimeInBatch = bucket * rfBucketPeriod_;
    double bunchTimeInBatch = bucketTimeInBatch + bunchOffset_;

    //Get start time of this bucket
    G4double bucketStartTime = startOfTurnTime + (batch*batchPeriod_) + (bucket*rfBucketPeriod_);

    //Generate whether or not this is a multiple occupancy RF bucket
    int numProtonsThisBucket = 1;
    if( G4RandFlat::shoot(0.,1.) < multipleOccupancyProb_ ) {

      //Generate number of protons in this multiple occupancy event
      double random = G4RandFlat::shoot(0.,1.);
      double runningWeightSum = 0.;
      for( auto weight : multipleOccupancyWeights_ ) { //Loop through probabilities for various # protons
        numProtonsThisBucket++;
        runningWeightSum += weight;
        if(random < runningWeightSum) break; //Use this # protons
      }

    }//if multiple occupancy

    if(plot_) histoMap_.find("h_bucketOccupancy")->second->Fill( numProtonsThisBucket ); 

    //Generate the times within the RF bucket for all protons (either single or multiple occupancy)
    for(int i(0) ; i<numProtonsThisBucket ; i++) {

      //Uniformly randomly sample a time from the bunch within the bucket to get this proton time (note 
      //that bunch width typically is smaller than total width of bucket)
      G4double protonTimeInBunch = G4RandFlat::shoot(0.,bunchWidth_);

      //Get the time of this proton within the spill
      G4double protonTime = startOfTurnTime + batchTimeInTurn + bunchTimeInBatch + protonTimeInBunch;

      //Add this generated time to vector so that it can be used to a proton event
      protonTimesToGen_.push_back(protonTime);

      //Fill some histos relating to this proton for debugging (this information isn't in tracking
      //action art record so cannot be directly accessed by analyzers later
      if(plot_) {
        histoMap_.find("h_MITurnTime")->second->Fill( startOfTurnTime/ms ); 
        histoMap_.find("h_MITurnNumber")->second->Fill( currentTurn_ ); 
        histoMap_.find("h_protonTimeInTurn")->second->Fill( (protonTime-startOfTurnTime)/microsecond ); 
        histoMap_.find("h_boosterBatch")->second->Fill( batch ); 
        histoMap_.find("h_protonTimeInBatch")->second->Fill( (bunchTimeInBatch+protonTimeInBunch)/ns ); 
        histoMap_.find("h_bucket")->second->Fill( bucket ); 
        histoMap_.find("h_protonTimeInBucket")->second->Fill( (bunchOffset_+protonTimeInBunch)/ns ); 
        histoMap_.find("h_protonTimeInBunch")->second->Fill( protonTimeInBunch/ns ); 
        histoMap_.find("h_spill")->second->Fill( currentSpill_ ); 
        histoMap_.find("h_protonTimeInSpill")->second->Fill( protonTime/ms ); 
        histoMap_.find("h_beamProfileX")->second->Fill( protonPositionBeamFrame.x()/mm ); 
        histoMap_.find("h_beamProfileY")->second->Fill( protonPositionBeamFrame.y()/mm ); 
        static_cast<TH2F*>(histoMap_.find("h_beamProfileXY")->second)->Fill( protonPositionBeamFrame.x()/mm , protonPositionBeamFrame.y()/mm ); 
      }

    }//num protons this turn loop

    //Sort the generated times into descending order(e.g. want back of vector to be earliest time)
    std::sort(protonTimesToGen_.begin(), protonTimesToGen_.end(), std::greater<G4double>());

    //Set current turn read for next event where proton time generation takes place
    currentTurn_++; 

    //Check if have reached end of spill
    if( currentTurn_ >= numMainInjectorTurnsPerSpill_ ) {
      currentTurn_ = 0;
      currentSpill_++;
      G4cout << "Starting spill " << currentSpill_ << G4endl;
    }

  }//end if size == 0

  //Now have at least proton time for this event (either because have leftover multiple occupancy 
  //protons from last event, or have generated a new turn this event). Use last one in vector (e.g. 
  //earliest time as vector sorted in descending order) as proton time for this event
  protonGun_->SetParticleTime(protonTimesToGen_.back());
  protonTimesToGen_.pop_back(); //Remove this time now have used it for a proton


  //TODO COmbine mutliple occupancy event sinto single event?



  //
  // Generate proton with this properties and fire it from particle gun
  //

  if (debug_) {
    G4cout
       << "Generated " << protonGun_->GetParticleDefinition()->GetParticleName() << " : \n"
       << "  Direction : " << protonGun_->GetParticleMomentumDirection() << "\n"
       << "  Energy    : " << G4BestUnit(protonGun_->GetParticleEnergy() ,"Energy") << "\n"
       << "  Time      : " << G4BestUnit(protonGun_->GetParticleTime(),"Time")
       << G4endl;
  }

  protonGun_->GeneratePrimaryVertex( evt );


}


unsigned int mtest::MTestProtonSource::GetNumBucketsPerBatch() const {
  return numBucketsPerBatch_;
}


unsigned int mtest::MTestProtonSource::GetNumBucketsPerTurn() const {
  return numBatchesPerTurn_ * this->GetNumBucketsPerBatch();
}


unsigned int mtest::MTestProtonSource::GetNumBucketsPerSpill() const {
  return numMainInjectorTurnsPerSpill_ * this->GetNumBucketsPerTurn();
}
