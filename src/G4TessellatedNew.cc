//
// Created by Stewart Boogert on 27/03/2023.
//

#include "G4TessellatedNew.hh"
#include "G4VFacet.hh"

G4TessellatedNew::G4TessellatedNew()  : G4VSolid("") {}


G4TessellatedNew::G4TessellatedNew(const G4String &name)  : G4VSolid(name) {}

void G4TessellatedNew::AddFacet(G4VFacet* aFacet) {
    fFacets.push_back(aFacet);
}

