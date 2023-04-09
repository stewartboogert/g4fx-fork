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

    auto dist = s_sdf->DistanceToIn(p,d);
    auto inside = s_sdf->Inside(p);
    G4cout << "sphere sdf " << dist << " " << inside << G4endl;

    G4Orb *s_g4 = new G4Orb("orb",10);
    dist = s_g4->DistanceToIn(p,d);
    inside = s_g4->Inside(p);
    G4cout << "orb " << dist << " " << inside << G4endl;

    G4Sphere *ss = new G4Sphere("sphere",0,10,0,2*M_PI,0,M_PI);
    dist = ss->DistanceToIn(p,d);
    inside = ss->Inside(p);
    G4cout << "sphere " << dist << " " << inside << G4endl;
}
