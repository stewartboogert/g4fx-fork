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
/// \file persistency/gdml/G01/load_gdml.cc
/// \brief Main program of the persistency/gdml/G01 example
//
//
//
//
// --------------------------------------------------------------
//      GEANT 4 - load_gdml
//
// --------------------------------------------------------------

#include <vector>

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"

#include "G01PrimaryGeneratorAction.hh"
#include "G01DetectorConstruction.hh"
#include "G01ActionInitialization.hh"

#include "FTFP_BERT.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4GDMLParser.hh"

#include "G4BooleanProcessorCGAL.hh"
#include "G4BooleanSolid.hh"

#include "G4TessellatedSolid.hh"
#include "G4TessellatedNew.hh"
#include "G4VtkSurfaceMeshLoader.hh"
// --------------------------------------------------------------

int main(int argc,char **argv)
{

    G4TessellatedSolid *tess_solid = new G4TessellatedSolid();
    G4TessellatedNew   *tess_solid2 = new G4TessellatedNew();


    G4VtkSurfaceMeshLoader *mesh_load = new G4VtkSurfaceMeshLoader();
    mesh_load->Load("teapot.stl");
    mesh_load->Fill(0, tess_solid);
    mesh_load->Fill(0, tess_solid2);

    tess_solid->SetSolidClosed(true);

    auto extent = tess_solid->GetExtent();

    auto mX = extent.GetXmin();
    auto mY = extent.GetYmin();
    auto mZ = extent.GetZmin();

    auto dX = extent.GetXmax() - extent.GetXmin();
    auto dY = extent.GetYmax() - extent.GetYmin();
    auto dZ = extent.GetZmax() - extent.GetZmin();

    dX *= 2;
    dY *= 2;
    dZ *= 2;

    // Get starting timepoint
    auto start = std::chrono::high_resolution_clock::now();

    auto nTest = 10000;
    for(int i=0;i<nTest;i++) {
        auto rX = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        auto rY = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        auto rZ = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        auto dX = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        auto dY = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        auto dZ = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        auto dNorm = sqrt(pow(dX,2) + pow(dY,2) + pow(dZ,2));
        dX /= dNorm;
        dY /= dNorm;
        dZ /= dNorm;

        auto x = rX * dX + mX;
        auto y = rY * dY + mY;
        auto z = rZ * dZ + mZ;
        tess_solid->DistanceToIn(G4ThreeVector(x,y,z));
        tess_solid->DistanceToOut(G4ThreeVector(x,y,z));
        tess_solid->DistanceToIn(G4ThreeVector(x,y,z), G4ThreeVector(dX,dY,dZ));
        tess_solid->Inside(G4ThreeVector(x,y,z));
        tess_solid->SafetyFromInside(G4ThreeVector(x,y,z),true);
        tess_solid->SafetyFromOutside(G4ThreeVector(x,y,z),true);
    }

    // Get starting timepoint
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    G4cout << "duration " << static_cast <float> (duration.count())/nTest << G4endl;

}
