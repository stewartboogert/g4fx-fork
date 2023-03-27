//
// Created by Stewart Boogert on 24/03/2023.
//

#include "G4MeshLoader.hh"

#include "G4Types.hh"
#include "G4ios.hh"
#include "G4String.hh"
#include "G4VSurfaceMesh.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4ThreeVector.hh"

#include "vtkSTLReader.h"

G4MeshLoader::G4MeshLoader() {}

void G4MeshLoader::Load(G4String file_name) {
    vtkSTLReader *sr = vtkSTLReader::New();
    sr->SetFileName(file_name);
    sr->Update();
    pd = sr->GetOutput();

}

void G4MeshLoader::Fill(G4int meshId, G4TessellatedSolid *tess) {
    auto v = pd->GetPoints();
    auto p = pd->GetPolys();

    G4cout << v->GetNumberOfPoints() << " " << pd->GetNumberOfPolys()  << G4endl;

    for(vtkIdType i=0;i < p->GetNumberOfCells(); i++) {
        vtkNew<vtkIdList> l;
        p->GetCellAtId(i, l);
        G4TriangularFacet *tf = new G4TriangularFacet(G4ThreeVector(v->GetPoint(l->GetId(0))[0],v->GetPoint(l->GetId(0))[1],v->GetPoint(l->GetId(0))[2]),
                                                     G4ThreeVector(v->GetPoint(l->GetId(1))[0],v->GetPoint(l->GetId(1))[1],v->GetPoint(l->GetId(1))[2]),
                                                     G4ThreeVector(v->GetPoint(l->GetId(2))[0],v->GetPoint(l->GetId(2))[1],v->GetPoint(l->GetId(2))[2]),
                                                     G4FacetVertexType::ABSOLUTE);
        tess->AddFacet(tf);
    }
}