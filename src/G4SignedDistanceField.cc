

#include "vtkImplicitFunction.h"
#include "vtkSampleFunction.h"
#include "vtkMarchingCubes.h"

#include "G4SignedDistanceField.hh"
#include "G4Sphere.hh"
#include "G4BoundingEnvelope.hh"
#include "G4AffineTransform.hh"

#include "G4VtkSurfaceMeshLoader.hh"

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