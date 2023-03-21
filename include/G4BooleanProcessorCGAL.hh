//
// Created by Stewart Boogert on 21/03/2023.
//

#ifndef G4BOOLEANPROCESSORCGAL_HH
#define G4BOOLEANPROCESSORCGAL_HH

#include "G4VBooleanProcessor.hh"
#include "G4Polyhedron.hh"
#include "G4PolyhedronArbitrary.hh"

class G4BooleanProcessorCGAL : public G4VBooleanProcessor {
public:
    G4BooleanProcessorCGAL() {}
    virtual ~G4BooleanProcessorCGAL() = default;

    virtual G4PolyhedronArbitrary* Intersection(G4Polyhedron* p1, G4Polyhedron* p2);
    virtual G4PolyhedronArbitrary* Union(G4Polyhedron* p1, G4Polyhedron* p2);
    virtual G4PolyhedronArbitrary* Subtraction(G4Polyhedron* p1, G4Polyhedron* p2);

};

#endif
