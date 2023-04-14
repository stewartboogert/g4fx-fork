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
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

class G4VoxelLimits;
class G4AffineTransform;

#include "vtkSmartPointer.h"
#include "vtkImplicitFunction.h"
class vtkPolyData;


namespace shader {
    G4TwoVector vec2(G4double x, G4double y);
    G4TwoVector vec2(G4double v);
    G4ThreeVector vec3(G4double x, G4double y , G4double z);
    G4ThreeVector vec3(G4double v);
    G4double length(const G4TwoVector &v);
    G4double length(const G4ThreeVector &v);
    G4TwoVector abs(const G4TwoVector &v);
    G4ThreeVector abs(const G4ThreeVector &v);
    G4ThreeVector max(const G4TwoVector &v1, const G4TwoVector &v2);
    G4ThreeVector max(const G4ThreeVector &v1, const G4ThreeVector &v2);
    G4ThreeVector max(const G4TwoVector &v1, G4double d2);
    G4ThreeVector max(const G4ThreeVector &v1, G4double d2);
    G4ThreeVector min(const G4TwoVector &v1, const G4TwoVector &v2);
    G4ThreeVector min(const G4ThreeVector &v1, const G4ThreeVector &v2);
    G4ThreeVector min(const G4TwoVector &v1, G4double d2);
    G4ThreeVector min(const G4ThreeVector &v1, G4double &d2);
    G4double clamp(G4double x, G4double min, G4double max);
    G4TwoVector clamp(const G4TwoVector &v, const G4TwoVector min, const G4TwoVector max);
    G4ThreeVector clamp(const G4ThreeVector v, const G4ThreeVector min, const G4ThreeVector max);
    G4TwoVector clamp(const G4TwoVector &v,  G4double min,  G4double max);
    G4ThreeVector clamp(const G4ThreeVector v, G4double min, G4double max);
};


class G4SignedDistanceField : public G4VSolid {
public:
    G4SignedDistanceField();
    G4SignedDistanceField(const G4String name);
    virtual ~G4SignedDistanceField() = default;


    virtual double Evaluate(const G4ThreeVector &) const = 0;
    G4ThreeVector EvaluateGradient(const G4ThreeVector &, double epsilon=1e-5) const;
    G4ThreeVector RayMarch(const G4ThreeVector &p, const G4ThreeVector &d) const;
    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override= 0;

    virtual G4bool CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform, G4double& pMin, G4double& pMax) const override;
    virtual EInside Inside(const G4ThreeVector& p) const override;
    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;
    virtual G4double DistanceToIn(const G4ThreeVector &p, const G4ThreeVector &d) const override ;
    virtual G4double DistanceToIn(const G4ThreeVector &p) const override;
    virtual G4double DistanceToOut(const G4ThreeVector &p, const G4ThreeVector &d,
                                   const G4bool calcNorm=false, G4bool* validNorm = nullptr,
                                   G4ThreeVector* n = nullptr) const override;
    virtual G4double DistanceToOut(const G4ThreeVector &p) const override;
    virtual G4GeometryType GetEntityType() const override {return G4String("G4SignedDistanceField");}
    virtual std::ostream& StreamInfo(std::ostream& os) const override {return os;}
    virtual void DescribeYourselfTo (G4VGraphicsScene& scene) const override;

protected:

private:
};

class G4VtkSignedDistanceField : public vtkImplicitFunction {
public:
    G4VtkSignedDistanceField(const G4SignedDistanceField *sdf) : fSdf(sdf) {}
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

    vtkSmartPointer<vtkPolyData> CubeMarch();

protected:

private:
    const G4SignedDistanceField *fSdf;
};

class G4SphereSDF : public G4SignedDistanceField {
public:
    G4SphereSDF() : G4SignedDistanceField() {fRadius = 1.0;}
    G4SphereSDF(G4String name, G4double radius) : G4SignedDistanceField(name) {fRadius = radius;}

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
    G4BoxSDF() : G4SignedDistanceField("dummy") {
        fX = 1.0;
        fY = 1.0;
        fZ = 1.0;
        fSize = G4ThreeVector(fX,fY,fZ);
    }
    G4BoxSDF(G4String name, G4double dHalfX, G4double dHalfY, G4double dHalfZ) : G4SignedDistanceField(name) {
        fX = 2.0*dHalfX;
        fY = 2.0*dHalfY;
        fZ = 2.0*dHalfZ;
        fSize = G4ThreeVector(fX,fY,fZ);
    }

    void SetSize(const G4ThreeVector size) {fSize = size;}
    G4ThreeVector GetSize() const {return fSize;}

    virtual double Evaluate(const G4ThreeVector &p ) const override {
        using namespace shader;
        using namespace std;

        auto q = abs(p) - fSize;
        auto sdf = length(max(q,0.0)) + min(max(q.x(),max(q.y(),q.z())),0.0);
        return sdf;
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(-fX, -fY, -fZ);
        bmax.set( fX,  fY,  fZ);
    }

    virtual G4GeometryType GetEntityType() const override {return G4String("G4BoxSDF");}

private:
    G4double fX;
    G4double fY;
    G4double fZ;
    G4ThreeVector fSize;
    G4double delta;
};

class G4BoxRoundSDF : public G4SignedDistanceField {
public:
    G4BoxRoundSDF() : G4SignedDistanceField("dummy") {
        fX = 1.0;
        fY = 1.0;
        fZ = 1.0;
        fSize = G4ThreeVector(fX,fY,fZ);
    }
    G4BoxRoundSDF(G4String name,
                  G4double dHalfX, G4double dHalfY, G4double dHalfZ,
                  G4double r) : G4SignedDistanceField(name) {
        fX = 2.0*dHalfX;
        fY = 2.0*dHalfY;
        fZ = 2.0*dHalfZ;
        fSize = G4ThreeVector(fX,fY,fZ);
        fR = r;
    }

    virtual double Evaluate(const G4ThreeVector &p ) const override {
        using namespace shader;
        using namespace std;

        auto q = abs(p) - fSize;
        return length(max(q,0.0)) + min(max(q.x(),max(q.y(),q.z())),0.0) - fR;
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(-fX, -fY, -fZ);
        bmax.set( fX,  fY,  fZ);
    }

    virtual G4GeometryType GetEntityType() const override {return G4String("G4BoxRoundSDF");}

private:
    G4double fX;
    G4double fY;
    G4double fZ;
    G4ThreeVector fSize;
    G4double fR;
};

class G4BoxFrameSDF : public G4SignedDistanceField {
public:
    G4BoxFrameSDF() : G4SignedDistanceField("dummy") {
        fX = 1.0;
        fY = 1.0;
        fZ = 1.0;
        fSize = G4ThreeVector(fX,fY,fZ);
    }
    G4BoxFrameSDF(G4String name,
                  G4double dHalfX, G4double dHalfY, G4double dHalfZ,
                  G4double e) : G4SignedDistanceField(name) {
        fX = 2.0*dHalfX;
        fY = 2.0*dHalfY;
        fZ = 2.0*dHalfZ;
        fSize = G4ThreeVector(fX,fY,fZ);
        fE = e;
    }

    void SetSize(const G4ThreeVector size) {fSize = size;}
    G4ThreeVector GetSize() const {return fSize;}

    virtual double Evaluate(const G4ThreeVector &p ) const override {
        using namespace shader;
        using namespace std;

        auto p1= abs(p)-fSize;
        auto q = abs(p1+vec3(fE))-vec3(fE);

        auto d = min(min(length(max(vec3(p1.x(),q.y(),q.z()),0.0))+min(max(p1.x(),max(q.y(),q.z())),0.0),
                         length(max(vec3(q.x(),p1.y(),q.z()),0.0))+min(max(q.x(),max(p1.y(),q.z())),0.0)),
                     length(max(vec3(q.x(),q.y(),p1.z()),0.0))+min(max(q.x(),max(q.y(),p1.z())),0.0));
        return d;
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(-fX, -fY, -fZ);
        bmax.set( fX,  fY,  fZ);
    }

    virtual G4GeometryType GetEntityType() const override {return G4String("G4BoxFrameSDF");}

private:
    G4double fX;
    G4double fY;
    G4double fZ;
    G4ThreeVector fSize;
    G4double fE;
};

class G4TorusSDF : public G4SignedDistanceField {
public:
    G4TorusSDF() : G4SignedDistanceField("dummy") {
    }
    G4TorusSDF(G4String name,
                  G4double dHalfX, G4double dHalfY, G4double dHalfZ,
                  G4double e) : G4SignedDistanceField(name) {
        fX = 2.0*dHalfX;
        fY = 2.0*dHalfY;
        fZ = 2.0*dHalfZ;
        fSize = G4ThreeVector(fX,fY,fZ);
        fE = e;
    }

    virtual double Evaluate(const G4ThreeVector &p ) const override {
        using namespace shader;
        using namespace std;

        return 0;
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(-fX, -fY, -fZ);
        bmax.set( fX,  fY,  fZ);
    }

    virtual G4GeometryType GetEntityType() const override {return G4String("G4TorusSDF");}

private:
    G4double fX;
    G4double fY;
    G4double fZ;
    G4ThreeVector fSize;
    G4double fE;
};



class G4DisplacedSDF : public G4SignedDistanceField {
public:
    G4DisplacedSDF(const G4String &name,
                     const G4RotationMatrix &rotation,
                     const G4ThreeVector &translation,
                     G4SignedDistanceField *sdf) :
                     G4SignedDistanceField(name),
                     fSdf(sdf) {
        fTransformation = new G4AffineTransform(rotation,translation);
    }

    G4DisplacedSDF(const G4String &name,
                   const G4Transform3D &transform,
                   G4SignedDistanceField *sdf) :
            G4SignedDistanceField(name),
            fSdf(sdf) {
        fTransformation = new G4AffineTransform(transform.getRotation().inverse(),
                                                transform.getTranslation());
    }

    G4DisplacedSDF(const G4String &name,
                   const G4AffineTransform &transform,
                   G4SignedDistanceField *sdf) :
            G4SignedDistanceField(name),
            fSdf(sdf) {
        fTransformation = new G4AffineTransform(transform);
    }

    virtual double Evaluate(const G4ThreeVector &) const override {
        return 0;
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(0, 0, 0);
        bmax.set(0, 0, 0);

    }

private:
    G4AffineTransform *fTransformation;
    G4SignedDistanceField *fSdf;
};

#endif //G4SIGNEDDISTANCEFIELD_HH
