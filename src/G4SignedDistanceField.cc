

#include "vtkImplicitFunction.h"
#include "vtkSampleFunction.h"
#include "vtkMarchingCubes.h"

#include "G4SignedDistanceField.hh"
#include "G4Sphere.hh"
#include "G4BoundingEnvelope.hh"
#include "G4AffineTransform.hh"

#include "G4VtkSurfaceMeshLoader.hh"

G4TwoVector shader::vec2(G4double x, G4double y) {return G4TwoVector(x,y);}
G4ThreeVector shader::vec3(G4double x, G4double y , G4double z) {return G4ThreeVector(x,y,z);}

G4double shader::length(G4TwoVector &v) {return v.mag();}
G4double shader::length(G4ThreeVector &v) {return v.mag();}

G4TwoVector shader::abs(G4TwoVector &v) {return G4TwoVector(fabs(v.x()), fabs(v.y()));}
G4ThreeVector shader::abs(G4ThreeVector &v) {return G4ThreeVector(fabs(v.x()), fabs(v.y()), fabs(v.z()));}

G4ThreeVector shader::max(G4TwoVector &v1, G4TwoVector &v2) {return G4TwoVector(std::max(v1.x(), v2.x()),
                                                                        std::max(v1.y(), v2.y()));}
G4ThreeVector shader::max(G4ThreeVector &v1, G4ThreeVector &v2) {return G4ThreeVector(std::max(v1.x(), v2.x()),
                                                                              std::max(v1.y(), v2.y()),
                                                                              std::max(v1.z(), v2.z()));}
G4ThreeVector shader::max(G4TwoVector &v1, G4double d2) {return G4TwoVector(std::max(v1.x(), d2),
                                                                    std::max(v1.y(), d2));}
G4ThreeVector shader::max(G4ThreeVector &v1, G4double &d2) {return G4ThreeVector(std::max(v1.x(), d2),
                                                                         std::max(v1.y(), d2),
                                                                         std::max(v1.z(), d2));}


G4ThreeVector shader::min(G4TwoVector &v1, G4TwoVector &v2) {return G4TwoVector(std::min(v1.x(), v2.x()),
                                                                        std::min(v1.y(), v2.y()));}
G4ThreeVector shader::min(G4ThreeVector &v1, G4ThreeVector &v2) {return G4ThreeVector(std::min(v1.x(), v2.x()),
                                                                              std::min(v1.y(), v2.y()),
                                                                              std::min(v1.z(), v2.z()));}
G4ThreeVector shader::min(G4TwoVector &v1, G4double d2) {return G4TwoVector(std::min(v1.x(), d2),
                                                                    std::min(v1.y(), d2));}
G4ThreeVector shader::min(G4ThreeVector &v1, G4double &d2) {return G4ThreeVector(std::min(v1.x(), d2),
                                                                         std::min(v1.y(), d2),
                                                                         std::min(v1.z(), d2));}

G4double shader::clamp(G4double x, G4double min, G4double max) {
    return x < min ? min : x > max ? max : x;
}

G4TwoVector shader::clamp(G4TwoVector &v, G4TwoVector min, G4TwoVector max) {
    return G4TwoVector(v.x() < min.x() ? min.x() : v.x() > max.x() ? max.x() : v.x(),
                       v.y() < min.y() ? min.y() : v.y() > max.y() ? max.y() : v.y());
}

G4ThreeVector shader::clamp(G4ThreeVector v, G4ThreeVector min, G4ThreeVector max) {
    return G4ThreeVector(v.x() < min.x() ? min.x() : v.x() > max.x() ? max.x() : v.x(),
                         v.y() < min.y() ? min.y() : v.y() > max.y() ? max.y() : v.y(),
                         v.z() < min.z() ? min.z() : v.z() > max.z() ? max.z() : v.z());
}

G4TwoVector shader::clamp(G4TwoVector &v, G4double min, G4double max) {
    return G4TwoVector(v.x() < min ? min : v.x() > max ? max: v.x(),
                       v.y() < min ? min : v.y() > max ? max : v.y());
}

G4ThreeVector shader::clamp(G4ThreeVector v, G4double min, G4double max) {
    return G4ThreeVector(v.x() < min ? min : v.x() > max ? max : v.x(),
                         v.y() < min ? min : v.y() > max ? max : v.y(),
                         v.z() < min ? min : v.z() > max ? max : v.z());
}


G4SignedDistanceField::G4SignedDistanceField() : G4VSolid("dummy") {}

G4SignedDistanceField::G4SignedDistanceField(const G4String name) : G4VSolid(name) {}

G4ThreeVector G4SignedDistanceField::EvaluateGradient(const G4ThreeVector &v, double epsilon) const {

    return G4ThreeVector((Evaluate(G4ThreeVector(v.x()+epsilon, v.y(), v.z())) - Evaluate(G4ThreeVector(v.x()-epsilon, v.y(), v.z())))/(2*epsilon),
                         (Evaluate(G4ThreeVector(v.x(), v.y()+epsilon, v.z())) - Evaluate(G4ThreeVector(v.x(), v.y()-epsilon, v.z())))/(2*epsilon),
                         (Evaluate(G4ThreeVector(v.x(), v.y(), v.z()+epsilon)) - Evaluate(G4ThreeVector(v.x(), v.y(), v.z()-epsilon)))/(2*epsilon));
}

G4ThreeVector G4SignedDistanceField::RayMarch(const G4ThreeVector &p, const G4ThreeVector &d) const {

    auto p1 = p+Evaluate(p)*d.unit();

    while(Evaluate(p1) > 1e-8) {
        p1 = RayMarch(p1,d);
    }
    return p1;
}

G4bool G4SignedDistanceField::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                                              const G4AffineTransform& pTransform, G4double& pMin, G4double& pMax) const
{
    G4ThreeVector bmin, bmax;

    // Get bounding box
    BoundingLimits(bmin,bmax);

    // Find extent
    G4BoundingEnvelope bbox(bmin,bmax);
    return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}

EInside G4SignedDistanceField::Inside(const G4ThreeVector& p) const  {
    if(Evaluate(p) <0) {
        return EInside::kInside;
    }
    else {
        return EInside::kOutside;
    }
}

G4ThreeVector G4SignedDistanceField::SurfaceNormal(const G4ThreeVector& p) const {
    return EvaluateGradient(p);
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

void G4VtkSignedDistanceField::CubeMarch() {
    auto resolution = 0.1;
    vtkNew<vtkSampleFunction> sampled;
    sampled->SetSampleDimensions(resolution, resolution, resolution);
    sampled->SetSampleDimensions(25,25,25);
    G4ThreeVector bmin;
    G4ThreeVector bmax;
    fSdf->BoundingLimits(bmin, bmax);
    bmin = bmin*1.1;
    bmax = bmax*1.1;
    sampled->SetModelBounds(bmin.x(), bmax.x(), bmin.y(), bmax.y(), bmin.z(), bmax.z());
    sampled->SetImplicitFunction(this);
    sampled->Update();

    vtkNew<vtkMarchingCubes> iso;
    iso->SetValue(0, 0);
    iso->SetInputConnection(sampled->GetOutputPort());
    iso->Update();

    G4VtkSurfaceMeshLoader *l = new G4VtkSurfaceMeshLoader();
    l->SetPolyData(iso->GetOutput());
    l->View();
}