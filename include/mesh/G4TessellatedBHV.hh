//
// Created by Stewart Boogert on 28/03/2023.
//

#ifndef GEANT4FX_G4TESSELLATEDBHV_HH
#define GEANT4FX_G4TESSELLATEDBHV_HH

#include "G4VSolid.hh"

class G4VFacet;

class G4TessellatedBVH : public G4VSolid {
public:
    G4TessellatedBVH();
    G4TessellatedBVH(const G4String& name);
    ~G4TessellatedBVH() = default;
    void AddFacet(G4VFacet* aFacet);

    virtual G4bool CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                   G4double& pMin, G4double& pMax) const {return true;};
    virtual EInside Inside(const G4ThreeVector& p) const  {return EInside::kInside;}
    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const {return G4ThreeVector();}
    virtual G4double DistanceToIn(const G4ThreeVector& p,
                                  const G4ThreeVector& v) const {return 0;}
    virtual G4double DistanceToIn(const G4ThreeVector& p) const {return 0;}
    virtual G4double DistanceToOut(const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm=false,
                                   G4bool* validNorm = nullptr,
                                   G4ThreeVector* n = nullptr) const {return 0;}
    virtual G4double DistanceToOut(const G4ThreeVector& p) const {return 0;}
    virtual G4GeometryType  GetEntityType() const {return G4String("G4TessellatedNew");}
    virtual std::ostream& StreamInfo(std::ostream& os) const {return os;}
    virtual void DescribeYourselfTo (G4VGraphicsScene& scene) const {}

protected:

private:
    std::vector<G4VFacet*> fFacets;
};

#endif //G4TESSELLATEDBHV_HH
