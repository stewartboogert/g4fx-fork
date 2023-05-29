#include <math.h>


#include "G4SystemOfUnits.hh"
#include "G4Types.hh"

#include "G4Orb.hh"
#include "G4Sphere.hh"

#include "G4SignedDistanceField.hh"
// --------------------------------------------------------------

int main(int argc,char **argv)
{
    auto solidSdf1  = new G4SphereSDF("solidSdf1",1.0*m);
    auto solid1     = new G4Orb("solid1",1.0*m);
    auto solidSdf2  = new G4BoxSDF("solidSdf2",1*m,1*m,1*m);
    auto solidSdf3  = new G4BoxRoundSDF("solidSdf3",3*m, 3*m,3*m, 0.5*m);
    auto solidSdf4  = new G4BoxFrameSDF("solidSdf4",0.5*m, 0.3*m,0.5*m,0.1*m);
    auto solidSdf5  = new G4TorusSDF("solidSdf5",1.5*m,0.25*m);
    auto solidSdf6  = new G4TorusCappedSDF("solidSdf6",5*m,0.5*m, cos(0.5), sin(0.5));
    auto solidSdf7  = new G4LinkSDF("solidSdf7",2.5*m,0.5*m,2.5*m);

    auto solidDisplaced    = new G4DisplacedSDF("displaced",solidSdf1,new G4RotationMatrix(G4ThreeVector(0,1,0),0.0),G4ThreeVector(0*m,0*m,0.999992*m));
    auto solidScaled       = new G4ScaledSDF("scaled",solidSdf2,G4Scale3D(1,2,3));

    auto solidUnion        = new G4UnionSDF("union",solidSdf2,solidSdf1,nullptr,G4ThreeVector(1*m,0,0));
    auto solidIntersection = new G4IntersectionSDF("intersection",solidSdf2,solidSdf1,nullptr,G4ThreeVector(1*m,0,0));
    auto solidSubtraction  = new G4SubtractionSDF("subtraction",solidSdf2,solidSdf1,nullptr,G4ThreeVector(1*m,0,0));

    auto s_sdf = solidSdf1;
    auto s_g4 = solid1;

    auto iBad = 0;

    for(auto i=0;i<1000000;i++) {
        auto rx = (double) rand() / RAND_MAX;
        auto ry = (double) rand() / RAND_MAX;
        auto rz = sqrt(pow(rx, 2) + pow(ry, 2));
        auto dx = 1000 * (double) rand() / RAND_MAX;
        auto dy = 1000 * (double) rand() / RAND_MAX;
        auto dz = 1000 * (double) rand() / RAND_MAX;

        G4ThreeVector p(dx, dy, dz);
        G4ThreeVector d(rx, ry, rz);

        d = d.unit();

        auto sdf_inside = s_sdf->Inside(p);
        auto sdf_distOutDir = s_sdf->DistanceToOut(p, d);
        auto sdf_distOut = s_sdf->DistanceToOut(p);
        auto sdf_distInDir = s_sdf->DistanceToIn(p, d);
        auto sdf_distIn = s_sdf->DistanceToIn(p);

        auto inside       = s_g4->Inside(p);
        auto distOutDir = s_g4->DistanceToOut(p, d);
        auto distOut    = s_g4->DistanceToOut(p);
        auto distInDir  = s_g4->DistanceToIn(p, d);
        auto distIn     = s_g4->DistanceToIn(p);

        auto diff_inside  = inside - sdf_inside;
        auto diff_distOutDir  = distOutDir - sdf_distOutDir;
        auto diff_distOut = distOut - sdf_distOut;
        auto diff_distInDir = distInDir - sdf_distInDir;
        auto diff_distIn = distIn - sdf_distIn;

        if (diff_inside != 0 ||
            fabs(diff_distOutDir) > 5e-9 ||
            fabs(diff_distOut) > 5e-9 ||
            fabs(diff_distInDir) > 5e-9 ||
            fabs(diff_distIn) > 5e-9) {
            G4cout << "x= " << dx << " "
                   << "y= " << dy << " "
                   << "z= " << dz << " "
                   << "dx= " << rx << " "
                   << "dy= " << ry << " "
                   << "dz= " << rz << " "
                   << "inside= " << sdf_inside << " ( " << inside << " ) "
                   << "distOutDir=" << sdf_distOutDir << " ( " << distOutDir << " , " << diff_distOutDir << " ) "
                   << "distOut= " << sdf_distOut << " ( " << distOut << " , " << diff_distOut << " ) "
                   << "distInDir= " << sdf_distInDir << " ( " << distInDir << " , " << diff_distInDir << " ) "
                   << "distIn= " << sdf_distIn << " ( " << distIn << " , " << diff_distIn << " ) " << G4endl;
            iBad++;
        }
    }
    G4cout << iBad << G4endl;
}
