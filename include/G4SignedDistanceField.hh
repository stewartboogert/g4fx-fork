//
// Created by Stewart Boogert on 09/04/2023.
//

// https://iquilezles.org/articles/distfunctions/
// Marching cubes https://en.wikipedia.org/wiki/Marching_cubes
// Marching tetraheda

#ifndef G4SIGNEDDISTANCEFIELD_HH
#define G4SIGNEDDISTANCEFIELD_HH

#include "G4ThreeVector.hh"
#include "G4VSolid.hh"


class G4SignedDistanceField : public G4VSolid {
public:
    G4SignedDistanceField();
    G4SignedDistanceField(const G4String name);
    virtual ~G4SignedDistanceField() = default;


    virtual double Evaluate(const G4ThreeVector &) const = 0;
    G4ThreeVector RayMarch(const G4ThreeVector &p, const G4ThreeVector &d) const;

    virtual G4bool CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                   G4double& pMin, G4double& pMax) const {return true;};
    virtual EInside Inside(const G4ThreeVector& p) const  {
        if(Evaluate(p) <0) {
            return EInside::kInside;
        }
        else {
            return EInside::kOutside;
        }
    }

    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const {
        return G4ThreeVector();
    }

    virtual G4double DistanceToIn(const G4ThreeVector &p, const G4ThreeVector &d) const;
    virtual G4double DistanceToIn(const G4ThreeVector &p) const;
    virtual G4double DistanceToOut(const G4ThreeVector &p, const G4ThreeVector &d,
                                   const G4bool calcNorm=false, G4bool* validNorm = nullptr,
                                   G4ThreeVector* n = nullptr) const;
    virtual G4double DistanceToOut(const G4ThreeVector &p) const;
    virtual G4GeometryType GetEntityType() const {return G4String("G4SignedDistanceField");}
    virtual std::ostream& StreamInfo(std::ostream& os) const {return os;}
    virtual void DescribeYourselfTo (G4VGraphicsScene& scene) const {}

protected:

private:
};

class G4SphereSDF : public G4SignedDistanceField {
public:
    G4SphereSDF() {fRadius = 1.0;}
    G4SphereSDF(G4double radius) {fRadius = radius;}

    void SetRadius(G4double radius) {fRadius = radius;}
    const G4double GetRadius() {return fRadius;}
    virtual double Evaluate(const G4ThreeVector &p ) const override {return p.mag() - fRadius;}
    virtual G4GeometryType GetEntityType() const override {return G4String("G4SphereSDF");}

private:
    G4double fRadius;
};

#endif //G4SIGNEDDISTANCEFIELD_HH
