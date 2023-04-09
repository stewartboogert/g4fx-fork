
#include "G4SignedDistanceField.hh"
#include "G4Sphere.hh"

G4SignedDistanceField::G4SignedDistanceField() : G4VSolid("dummy") {}

G4SignedDistanceField::G4SignedDistanceField(const G4String name) : G4VSolid(name) {}

G4ThreeVector G4SignedDistanceField::RayMarch(const G4ThreeVector &p, const G4ThreeVector &d) const {

    auto p1 = p+Evaluate(p)*d.unit();

    while(Evaluate(p1) > 1e-8) {
        p1 = RayMarch(p1,d);
    }
    return p1;
}

double G4SignedDistanceField::DistanceToIn(const G4ThreeVector &p, const G4ThreeVector &d) const {
    auto p1 = RayMarch(p,d);
    return (p1-p).mag();
}

double G4SignedDistanceField::DistanceToIn(const G4ThreeVector &p) const {
    return Evaluate(p);
}

double G4SignedDistanceField::DistanceToOut(const G4ThreeVector &p, const G4ThreeVector &d,
                                            const G4bool calcNorm, G4bool* validNorm,
                                            G4ThreeVector* n) const {
    auto p1 = RayMarch(p,d);
    return (p1-p).mag();
}

double G4SignedDistanceField::DistanceToOut(const G4ThreeVector &p) const {
    return Evaluate(p);
}