#include <math.h>

#include "G4Types.hh"

#include "G4Orb.hh"
#include "G4Sphere.hh"

#include "G4SignedDistanceField.hh"
// --------------------------------------------------------------

int main(int argc,char **argv)
{
    G4SphereSDF *s_sdf = new G4SphereSDF("sdf",10);
    G4Orb       *o_g4 = new G4Orb("orb",10);
    G4Sphere    *s_g4 = new G4Sphere("sphere",0,10,0,2*M_PI,0,M_PI);

    G4ThreeVector p(-10+1,0,0);
    G4ThreeVector d(1,0 ,0);

    d = d.unit();

    auto distIn = 1e99;
    auto distOut = 1e99;

    auto sdf_inside = s_sdf->Inside(p);
    auto o_inside   = o_g4->Inside(p);
    auto s_inside   = s_g4->Inside(p);

    auto sdf_distOutDir = s_sdf->DistanceToOut(p,d);
    auto sdf_distOut    = s_sdf->DistanceToOut(p);
    auto sdf_distInDir  = s_sdf->DistanceToIn(p,d);
    auto sdf_distIn     = s_sdf->DistanceToIn(p);

    auto o_distOutDir = o_g4->DistanceToOut(p,d);
    auto o_distOut    = o_g4->DistanceToOut(p);
    auto o_distInDir  = o_g4->DistanceToIn(p,d);
    auto o_distIn     = o_g4->DistanceToIn(p);

    auto s_distOutDir = s_g4->DistanceToOut(p,d);
    auto s_distOut    = s_g4->DistanceToOut(p);
    auto s_distInDir  = s_g4->DistanceToIn(p,d);
    auto s_distIn     = s_g4->DistanceToIn(p);

    G4cout << sdf_inside << " " << sdf_distOutDir << " " << sdf_distOut << " " << sdf_distInDir << " " << sdf_distIn << G4endl;
    G4cout << o_inside << " " << o_distOutDir << " " << o_distOut << " " << o_distInDir << " " << o_distIn << G4endl;
    G4cout << s_inside << " " << s_distOutDir << " " << s_distOut << " " << s_distInDir << " " << s_distIn << G4endl;

}
