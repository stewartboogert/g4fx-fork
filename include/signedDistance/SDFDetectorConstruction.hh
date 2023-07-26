//
// Created by Stewart Boogert on 11/04/2023.
//

#ifndef SDFDETECTORCONSTRUCTION_HH
#define SDFDETECTORCONSTRUCTION_HH

#include <math.h>

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

        auto solidSdf1  = new G4SphereSDF("solidSdf1",1.0*m);
        auto solidSdf2  = new G4BoxSDF("solidSdf2",1*m,1*m,1*m);
        auto solid2     = new G4Box("solid",1*m,1*m,1*m);
        auto solidSdf3  = new G4BoxRoundSDF("solidSdf3",3*m, 3*m,3*m, 0.5*m);
        auto solidSdf4  = new G4BoxFrameSDF("solidSdf4",0.5*m, 0.3*m,0.5*m,0.1*m);
        auto solidSdf5  = new G4TorusSDF("solidSdf5",1.5*m,0.25*m);
        auto solidSdf6  = new G4TorusCappedSDF("solidSdf6",5*m,0.5*m, cos(0.5), sin(0.5));
        auto solidSdf7  = new G4LinkSDF("solidSdf7",2.5*m,0.5*m,2.5*m);

        auto solidDisplaced = new G4DisplacedSDF("displaced",solidSdf1,new G4RotationMatrix(G4ThreeVector(0,1,0),0.0),G4ThreeVector(0*m,0*m,0.0*m));
        auto solidScaled = new G4ScaledSDF("scaled",solidSdf2,G4Scale3D(1,2,3));

        auto solidUnion      = new G4UnionSDF("union",solidSdf2,solidSdf1,nullptr,G4ThreeVector(1*m,0,0));
        auto solidIntersection = new G4IntersectionSDF("intersection",solidSdf2,solidSdf1,nullptr,G4ThreeVector(1*m,0,0));
        auto solidSubtraction = new G4SubtractionSDF("subtraction",solidSdf2,solidSdf1,nullptr,G4ThreeVector(1*m,0,0));

        auto multiUnion = new G4MultiUnionSDF("multiunion");
        for(auto i = 0; i<1;i++) {
            auto rx = (double)rand() / RAND_MAX;
            auto ry = (double)rand() / RAND_MAX;
            auto rz = sqrt(pow(rx,2)+pow(ry,2));
            auto dx = 1000*(double)rand() / RAND_MAX;
            auto dy = 1000*(double)rand() / RAND_MAX;
            auto dz = 1000*(double)rand() / RAND_MAX;
            auto a = M_PI*(double)rand() / RAND_MAX;

            rx = 0;
            ry = 0;
            rz = 0;
            dx = 0;
            dy = 0;
            dz = 0;
            a  = 0;

            multiUnion->AddNode(solidSdf2,G4Transform3D(G4RotationMatrix(G4ThreeVector(rx,ry,rz),a), G4ThreeVector(dx,dy,dz)));

        }

        auto solidUnionSmooth = new G4UnionSmoothSDF("unionSmooth",solidSdf2,solidSdf1,nullptr,G4ThreeVector(1*m,0,0),0.25*m);
        auto solidIntersectionSmooth = new G4IntersectionSmoothSDF("intersectionSmooth",solidSdf2,solidSdf1,nullptr,G4ThreeVector(1*m,0,0),0.25*m);
        auto solidSubtractionSmooth = new G4SubtractionSmoothSDF("subtractionSmooth",solidSdf2,solidSdf1,nullptr,G4ThreeVector(1*m,0,0),0.25*m);

        auto logicSdf = new G4LogicalVolume(solidDisplaced,sdf_mat,"sdf");

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
