#include <math.h>

#include "G4Types.hh"

#include "G4Orb.hh"
#include "G4Sphere.hh"

#include "G4SignedDistanceField.hh"
// --------------------------------------------------------------

int main(int argc,char **argv)
{
    G4SphereSDF *s_sdf = new G4SphereSDF(10);

    G4ThreeVector p(100,100,100);
    G4ThreeVector d(-1,-1 ,-1);
    d = d.unit();

    auto distIn = 1e99;
    auto distOut = 1e99;

    auto inside = s_sdf->Inside(p);
    if(inside)
        distOut = s_sdf->DistanceToOut(p,d);
    else
        distIn  = s_sdf->DistanceToIn(p,d);
    G4cout << "sphere sdf " << distIn << " " << distOut << " " << inside << G4endl;

    G4Orb *s_g4 = new G4Orb("orb",10);
    inside = s_g4->Inside(p);
    if(inside)
        distOut = s_g4->DistanceToOut(p,d);
    else
        distIn = s_g4->DistanceToIn(p,d);

    G4cout << "orb " << distIn << " " << distOut << " " << inside << G4endl;

    G4Sphere *ss = new G4Sphere("sphere",0,10,0,2*M_PI,0,M_PI);
    inside = ss->Inside(p);
    if(inside)
        distOut = ss->DistanceToOut(p,d);
    else
        distIn = ss->DistanceToIn(p,d);
    G4cout << "sphere " << distIn << " " << distOut << " " << inside << G4endl;
}
