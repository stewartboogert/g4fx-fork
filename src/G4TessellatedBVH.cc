//
// Created by Stewart Boogert on 28/03/2023.
//

#include "G4TessellatedBHV.hh"

#include "G4VFacet.hh"
#include "G4VSolid.hh"

G4TessellatedBVH::G4TessellatedBVH()  : G4VSolid("") {}


G4TessellatedBVH::G4TessellatedBVH(const G4String &name)  : G4VSolid(name) {}

void G4TessellatedBVH::AddFacet(G4VFacet* aFacet) {
    fFacets.push_back(aFacet);
}