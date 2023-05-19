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
    G4double      mix(const G4double v1, const G4double v2, const G4double a);
    G4TwoVector   mix(const G4TwoVector &v1, const G4TwoVector &v2, const G4double a);
    G4ThreeVector mix(const G4ThreeVector &v1, const G4ThreeVector &v2, const G4double a);
    G4double clamp(G4double x, G4double min, G4double max);
    G4TwoVector clamp(const G4TwoVector &v, const G4TwoVector min, const G4TwoVector max);
    G4ThreeVector clamp(const G4ThreeVector v, const G4ThreeVector min, const G4ThreeVector max);
    G4TwoVector clamp(const G4TwoVector &v,  G4double min,  G4double max);
    G4ThreeVector clamp(const G4ThreeVector v, G4double min, G4double max);
    G4TwoVector xy(const G4ThreeVector &v);
    G4TwoVector xz(const G4ThreeVector &v);
    G4TwoVector yz(const G4ThreeVector &v);
    G4double dot(const G4TwoVector &v1, const G4TwoVector &v2);
    G4double dot(const G4ThreeVector &v1, const G4ThreeVector &v2);
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
                  G4double r1,
                  G4double r2) : G4SignedDistanceField(name) {
        fR1 = r1;
        fR2 = r2;
    }

    virtual double Evaluate(const G4ThreeVector &p ) const override {
        using namespace shader;
        using namespace std;

        auto t = vec2(fR1, fR2);
        auto q = vec2(length(vec2(p.x(), p.z()))-t.x(),p.y());
        return length(q)-t.y();
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(-fR1-fR2,-fR2, -fR1-fR2);
        bmax.set( fR1+fR2, fR2,  fR1+fR2);
    }

    virtual G4GeometryType GetEntityType() const override {return G4String("G4TorusSDF");}

private:
    G4double fR1;
    G4double fR2;
};

class G4TorusCappedSDF : public G4SignedDistanceField {
public:
    G4TorusCappedSDF() : G4SignedDistanceField("dummy") {
    }

    G4TorusCappedSDF(G4String name,
               G4double r1,
               G4double r2,
               G4double cx,
               G4double cy) : G4SignedDistanceField(name) {
        fR1 = r1;
        fR2 = r2;
        fCx = cx;
        fCy = cy;
    }

    virtual double Evaluate(const G4ThreeVector &p1 ) const override {
        using namespace shader;
        using namespace std;

        auto p  = vec3(abs(p1.x()), p1.y(), p1.z());
        auto k = (fCy*p.x()>fCx*p.y()) ? dot(xy(p),vec2(fCx,fCy)) : length(xy(p));
        return sqrt( dot(p,p) + fR1*fR1 - 2.0*fR1*k ) - fR2;
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(-fR1-fR2,-fR1-fR2, -fR2);
        bmax.set( fR1+fR2, fR1+fR2, fR2 );
    }

    virtual G4GeometryType GetEntityType() const override {return G4String("G4TorusCappedSDF");}

private:
    G4double fCx;
    G4double fCy;
    G4double fR1;
    G4double fR2;

};

class G4LinkSDF : public G4SignedDistanceField {
public:
    G4LinkSDF() : G4SignedDistanceField("dummy") {
    }

    G4LinkSDF(G4String name,
                     G4double r1,
                     G4double r2,
                     G4double le) : G4SignedDistanceField(name) {
        fR1 = r1;
        fR2 = r2;
        fLe = le;
    }

    virtual double Evaluate(const G4ThreeVector &p ) const override {
        using namespace shader;
        using namespace std;

        auto q = vec3(p.x(),max(abs(p.y())-fLe,0.0), p.z());
        return length(vec2(length(xy(q))-fR1,q.z())) - fR2;
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        bmin.set(-fR1-fR2,-fR1-fR2-fLe, -fR2);
        bmax.set( fR1+fR2, fR1+fR2+fLe, fR2 );
    }

    virtual G4GeometryType GetEntityType() const override {return G4String("G4TorusCappedSDF");}

private:
    G4double fR1;
    G4double fR2;
    G4double fLe;
};

class G4DisplacedSDF : public G4SignedDistanceField {
public:
    G4DisplacedSDF(const G4String &name,
                         G4SignedDistanceField *sdf,
                         G4RotationMatrix *rotation,
                   const G4ThreeVector &translation) :
                   G4SignedDistanceField(name),
                   fSdf(sdf) {
        fTransformation = new G4AffineTransform(rotation,translation);
    }

    G4DisplacedSDF(const G4String &name,
                         G4SignedDistanceField *sdf,
                   const G4Transform3D &transform) :
            G4SignedDistanceField(name),
            fSdf(sdf) {
        fTransformation = new G4AffineTransform(transform.getRotation().inverse(),
                                                transform.getTranslation());
    }

    G4DisplacedSDF(const G4String &name,
                         G4SignedDistanceField *sdf,
                   const G4AffineTransform &transform) :
            G4SignedDistanceField(name),
            fSdf(sdf) {
        fTransformation = new G4AffineTransform(transform);
    }

    virtual double Evaluate(const G4ThreeVector &p) const override {
        return fSdf->Evaluate(fTransformation->TransformPoint(p));
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        fSdf->BoundingLimits(bmin,bmax);

        bmin = fTransformation->Inverse().TransformPoint(bmin);
        bmax = fTransformation->Inverse().TransformPoint(bmax);
    }

private:
    G4AffineTransform *fTransformation;
    G4SignedDistanceField *fSdf;
};

// G4ScaledSDF
// G4SymmetrySDF

class G4BooleanSDF : public G4SignedDistanceField {
public:
    G4BooleanSDF(const G4String& pName,
                 G4SignedDistanceField* pSolidA,
                 G4SignedDistanceField* pSolidB) : G4SignedDistanceField(pName), fPtrSolidA(pSolidA), fPtrSolidB(pSolidB) {}

    G4BooleanSDF(const G4String& pName,
                 G4SignedDistanceField* pSolidA ,
                 G4SignedDistanceField* pSolidB,
                 G4RotationMatrix* rotMatrix,
                 const G4ThreeVector& transVector) : G4SignedDistanceField(pName), fPtrSolidA(pSolidA) {
        fPtrSolidB = new G4DisplacedSDF("placedB",pSolidB,rotMatrix,transVector) ;
    };

    G4BooleanSDF(const G4String& pName,
                 G4SignedDistanceField* pSolidA ,
                 G4SignedDistanceField* pSolidB ,
                 const G4Transform3D& transform) : G4SignedDistanceField(pName), fPtrSolidA(pSolidA) {
        fPtrSolidB = new G4DisplacedSDF("placedB",pSolidB,transform) ;
    };

protected:
    G4SignedDistanceField* fPtrSolidA = nullptr;
    G4SignedDistanceField* fPtrSolidB = nullptr;
};
// G4IntersectionSDF
class G4UnionSDF : public G4BooleanSDF {
public:

    G4UnionSDF(const G4String &name,
               G4SignedDistanceField *sdf1,
               G4SignedDistanceField *sdf2) :
            G4BooleanSDF(name, sdf1, sdf2){
    }

    G4UnionSDF(const G4String &name,
               G4SignedDistanceField *sdf1,
               G4SignedDistanceField *sdf2,
               G4RotationMatrix *rotation,
               const G4ThreeVector &translation) :
            G4BooleanSDF(name, sdf1, sdf2, rotation, translation) {
    }

    G4UnionSDF(const G4String &name,
               G4SignedDistanceField *sdf1,
               G4SignedDistanceField *sdf2,
               const G4Transform3D &transform) :
            G4BooleanSDF(name, sdf1, sdf2, transform) {
    }

    G4UnionSDF(const G4String &name,
               G4SignedDistanceField *sdf1,
               G4SignedDistanceField *sdf2,
               const G4AffineTransform &transform) :
            G4BooleanSDF(name, sdf1, sdf2, transform) {
    }

    virtual double Evaluate(const G4ThreeVector &p) const override {
        using namespace std;

        auto d1 = fPtrSolidA->Evaluate(p);
        auto d2 = fPtrSolidB->Evaluate(p);
        return min(d1,d2);
    }

    virtual void BoundingLimits(G4ThreeVector &bmin, G4ThreeVector &bmax) const override {
        auto b1min = G4ThreeVector();
        auto b1max = G4ThreeVector();
        fPtrSolidA->BoundingLimits(b1min,b1max);

        auto b2min = G4ThreeVector();
        auto b2max = G4ThreeVector();
        fPtrSolidB->BoundingLimits(b2min,b2max);

        using namespace shader;

        bmin = min(b1min,b2min);
        bmax = max(b1max,b2max);
    }

private:

};
// G4SubtractionSDF
// G4MultiUnionSDF

// G4IntersectionSmoothSDF
// G4UnionSmoothSDF
// G4SubtractionSmoothSDF
// G4MultiUnionSmoothSDF

// G4ElongationSDF
// G4RoundingSDF
// G4RevolutionSDF
// G4ExtrusionSDF

// G4TwistSDF
// G4BendSDF

#endif //G4SIGNEDDISTANCEFIELD_HH
