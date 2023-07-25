//
// Created by Stewart Boogert on 27/03/2023.
//

#ifndef G4TESSELLATEDNEW_HH
#define G4TESSELLATEDNEW_HH

#include <vector>

#include "G4VSolid.hh"

class G4VFacet;

class G4TessellatedNew : public G4VSolid {
public:
    G4TessellatedNew();
    G4TessellatedNew(const G4String& name);
    ~G4TessellatedNew() = default;
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

#endif //G4TESSELLATEDNEW_HH
