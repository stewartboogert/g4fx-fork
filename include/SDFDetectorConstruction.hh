//
// Created by Stewart Boogert on 11/04/2023.
//

#ifndef SDFDETECTORCONSTRUCTION_HH
#define SDFDETECTORCONSTRUCTION_HH

#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SignedDistanceField.hh"
#include "G4Box.hh"
#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Orb.hh"

class SDFDetectorConstruction : public G4VUserDetectorConstruction
{
public:

    SDFDetectorConstruction() {}

    virtual G4VPhysicalVolume *Construct()
    {
        // Get nist material manager
        G4NistManager* nist = G4NistManager::Instance();

        G4double world_sizeXY = 20*m;
        G4double world_sizeZ  = 20*m;
        G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
        G4Material* sdf_mat = nist->FindOrBuildMaterial("G4_Al");

        auto solidWorld = new G4Box("World",
                                    0.5 * world_sizeXY,
                                    0.5 * world_sizeXY,
                                    0.5 * world_sizeZ);

        auto logicWorld = new G4LogicalVolume(solidWorld,
                                              world_mat,
                                              "World");

        auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                           G4ThreeVector(),
                                           logicWorld,
                                           "World",
                                           nullptr,
                                           false,
                                           0,
                                           false);

        auto solidSdf1  = new G4SphereSDF("solidSdf1",1*m);
        auto solidSdf2  = new G4BoxSDF("solidSdf2",3*m,1*m,1*m);
        auto solidSdf3  = new G4BoxRoundSDF("solidSdf3",3*m, 3*m,3*m, 0.5*m);
        auto solidSdf4  = new G4BoxFrameSDF("solidSdf4",0.5*m, 0.3*m,0.5*m,0.1*m);
        auto solidSdf5  = new G4TorusSDF("solidSdf5",1.5*m,0.25*m);
        auto solidSdf6  = new G4TorusCappedSDF("solidSdf6",5*m,0.5*m, cos(0.5), sin(0.5));
        auto solidSdf7  = new G4LinkSDF("solidSdf7",2.5*m,0.5*m,2.5*m);

        auto solidUnion      = new G4UnionSDF("union",solidSdf1,solidSdf6,nullptr,G4ThreeVector(-2*m,0,0));
        auto solidIntersection = new G4IntersectionSDF("intersection",solidSdf6,solidSdf2,nullptr,G4ThreeVector(-2*m,0,0));
        auto solidSubtraction = new G4IntersectionSDF("difference",solidSdf6,solidSdf2,nullptr,G4ThreeVector(-2*m,0,0));

        auto solidUnionSmooth = new G4UnionSmoothSDF("unionSmooth",solidSdf7,solidSdf2,nullptr,G4ThreeVector(0*m,0,0),1*m);
        auto solidIntersectionSmooth = new G4IntersectionSmoothSDF("intersectionSmooth",solidSdf7,solidSdf2,nullptr,G4ThreeVector(0*m,0,0),1*m);
        auto solidSubtractionSmooth = new G4SubtractionSmoothSDF("subtractionSmooth",solidSdf7,solidSdf2,nullptr,G4ThreeVector(0*m,0,0),1*m);

        auto logicSdf = new G4LogicalVolume(solidSubtractionSmooth,sdf_mat,"sdf");

        auto physSdf = new G4PVPlacement(nullptr,  // no rotation
                                         G4ThreeVector(0,0,0),
                                         logicSdf,
                                         "sdfPhys",
                                         logicWorld,
                                         false,
                                         0,
                                         false);

        return physWorld;
    }

private:
    G4VPhysicalVolume *fWorld;
};

#endif //SDFDETECTORCONSTRUCTION_HH
