//
// Created by Stewart Boogert on 09/04/2023.
//

// https://iquilezles.org/articles/distfunctions/
// Marching cubes https://en.wikipedia.org/wiki/Marching_cubes
// Marching tetraheda

#ifndef G4SIGNEDDISTANCEFIELD_HH
#define G4SIGNEDDISTANCEFIELD_HH

#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"

class G4VoxelLimits;
class G4AffineTransform;

#include "vtkImplicitFunction.h"

namespace shader {
    G4TwoVector vec2(G4double x, G4double y);
    G4ThreeVector vec3(G4double x, G4double y , G4double z);
    G4double length(G4TwoVector &v);
    G4double length(G4ThreeVector &v);
    G4TwoVector abs(G4TwoVector &v);
    G4ThreeVector abs(G4ThreeVector &v);
    G4ThreeVector max(G4TwoVector &v1, G4TwoVector &v2);
    G4ThreeVector max(G4ThreeVector &v1, G4ThreeVector &v2);
    G4ThreeVector max(G4TwoVector &v1, G4double d2);
    G4ThreeVector max(G4ThreeVector &v1, G4double &d2);
    G4ThreeVector min(G4TwoVector &v1, G4TwoVector &v2);
    G4ThreeVector min(G4ThreeVector &v1, G4ThreeVector &v2);
    G4ThreeVector min(G4TwoVector &v1, G4double d2);
    G4ThreeVector min(G4ThreeVector &v1, G4double &d2);
    G4double clamp(G4double x, G4double min, G4double max);
    G4TwoVector clamp(G4TwoVector &v, G4TwoVector min, G4TwoVector max);
    G4ThreeVector clamp(G4ThreeVector v, G4ThreeVector min, G4ThreeVector max);
    G4TwoVector clamp(G4TwoVector &v, G4double min, G4double max);
    G4ThreeVector clamp(G4ThreeVector v, G4double min, G4double max);
};


class G4SignedDistanceField : public G4VSolid {
public:
    G4SignedDistanceField();
    G4SignedDistanceField(const G4String name);
    virtual ~G4SignedDistanceField() = default;


    virtual double Evaluate(const G4ThreeVector &) const = 0;
    G4ThreeVector EvaluateGradient(const G4ThreeVector &, double epsilon=1e-5) const;
    G4ThreeVector RayMarch(const G4ThreeVector &p, const G4ThreeVector &d) const;
    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const = 0;

    virtual G4bool CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform, G4double& pMin, G4double& pMax) const;
    virtual EInside Inside(const G4ThreeVector& p) const;
    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;
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

class G4VtkSignedDistanceField : public vtkImplicitFunction {
public:
    G4VtkSignedDistanceField(G4SignedDistanceField *sdf) : fSdf(sdf) {}
    virtual double EvaluateFunction(double x[3]) {
        return fSdf->Evaluate(G4ThreeVector(x[0],x[1],x[2]));
    }
    virtual void EvaluateGradient(double x[3], double g[3]) {
        auto gdt = fSdf->EvaluateGradient(G4ThreeVector(x[0],x[1],x[2]));
        g[0] = gdt.x();
        g[1] = gdt.y();
        g[1] = gdt.z();
        return;
    }

    void CubeMarch();

protected:

private:
    G4SignedDistanceField *fSdf;
};

class G4SphereSDF : public G4SignedDistanceField {
public:
    G4SphereSDF() {fRadius = 1.0;}
    G4SphereSDF(G4double radius) {fRadius = radius;}

    void SetRadius(G4double radius) {fRadius = radius;}
    const G4double GetRadius() {return fRadius;}
    virtual double Evaluate(const G4ThreeVector &p ) const override {return p.mag() - fRadius;}
    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(-fRadius, -fRadius, -fRadius);
        bmax.set( fRadius,  fRadius,  fRadius);
    }
    virtual G4GeometryType GetEntityType() const override {return G4String("G4SphereSDF");}

private:
    G4double fRadius;
};

class G4BoxSDF : public G4SignedDistanceField {
public:
    G4BoxSDF() {
        fX = 1.0;
        fY = 1.0;
        fZ = 1.0;
        fSize = G4ThreeVector(fX,fY,fZ);
    }
    G4BoxSDF(G4double dHalfX, G4double dHalfY, G4double dHalfZ) {
        fX = 2.0*dHalfX;
        fY = 2.0*dHalfY;
        fZ = 2.0*dHalfZ;
        fSize = G4ThreeVector(fX,fY,fZ);
    }

    void SetSize(const G4ThreeVector size) {fSize = size;}
    G4ThreeVector GetSize() const {return fSize;}

    virtual double Evaluate(const G4ThreeVector &p ) const override {
        return p.mag();
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(-fX/2., -fY/2., -fZ/2.);
        bmax.set(fX/2,    fY/2,   fZ/2.);
    }

    virtual G4GeometryType GetEntityType() const override {return G4String("G4BoxSDF");}

private:
    G4double fX;
    G4double fY;
    G4double fZ;
    G4ThreeVector fSize;
};


#endif //G4SIGNEDDISTANCEFIELD_HH
