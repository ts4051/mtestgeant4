//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
/// \file ExamplePrimaryGeneratorAction.hh
/// \brief Definition of the ExamplePrimaryGeneratorAction class

#ifndef ExamplePrimaryGeneratorAction_h
#define ExamplePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
//#include "G4ParticleGun.hh"
#include "globals.hh"

#include <memory>

#include "MTestProtonSource.hh"

//class G4ParticleGun;
class G4Event;
class ExampleDetectorConstruction;

/// The primary generator action class with particle gum.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued 
/// in front of the phantom across 80% of the (X,Y) phantom size.

class ExamplePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExamplePrimaryGeneratorAction();    
    virtual ~ExamplePrimaryGeneratorAction();

    // static access method
    static const ExamplePrimaryGeneratorAction* Instance();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         
  
    // method to access particle gun
    //const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
    //const mtest::MTestProtonSource* getParticleSource() const { return fParticleSource; }
  
  private:
    static ExamplePrimaryGeneratorAction* fgInstance;
   
    std::unique_ptr<mtest::MTestProtonSource> fParticleSource; // pointer to the particle generator

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


