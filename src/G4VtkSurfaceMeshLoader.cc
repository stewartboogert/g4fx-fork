//
// Created by Stewart Boogert on 24/03/2023.
//

#include "G4VtkSurfaceMeshLoader.hh"

#include "G4Types.hh"
#include "G4ios.hh"
#include "G4String.hh"
#include "G4VSurfaceMesh.hh"
#include "G4TessellatedSolid.hh"
#include "G4TessellatedNew.hh"
#include "G4TriangularFacet.hh"
#include "G4ThreeVector.hh"

#include "vtkSmartPointer.h"
#include "vtkSTLReader.h"
#include "vtkOBJReader.h"
#include "vtkPLYReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkTriangle.h"
#include "vtkDelaunay3D.h"

G4VtkSurfaceMeshLoader::G4VtkSurfaceMeshLoader() {
    pd = vtkSmartPointer<vtkPolyData>::New();
}

void G4VtkSurfaceMeshLoader::Load(G4String file_name) {

    if(file_name.find("stl") !=std::string::npos) {
        vtkSmartPointer<vtkSTLReader> stlr = vtkSmartPointer<vtkSTLReader>::New();
        stlr->SetFileName(file_name);
        stlr->Update();
        pd = stlr->GetOutput();
    }
    else if(file_name.find("obj") !=std::string::npos) {
        vtkSmartPointer<vtkOBJReader> objr = vtkSmartPointer<vtkOBJReader>::New();
        objr->SetFileName(file_name);
        objr->Update();
        pd = objr->GetOutput();
    }
    else if(file_name.find("ply") !=std::string::npos) {
        vtkSmartPointer<vtkPLYReader> plyr = vtkSmartPointer<vtkPLYReader>::New();
        plyr->SetFileName(file_name);
        plyr->Update();
        pd = plyr->GetOutput();
    }
    else if(file_name.find(".vtp") != std::string::npos) {
        vtkSmartPointer<vtkXMLPolyDataReader> vtpr = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        vtpr->SetFileName(file_name);
        vtpr->Update();
        pd = vtpr->GetOutput();
    }
    else {
        G4cout << "Unknown file type" << G4endl;
    }
}

void G4VtkSurfaceMeshLoader::LoadOFF(G4String file_name) {
    std::ifstream ifstr(file_name);

    G4String dummy;
    G4int nVertex, nFace, nEdge;
    G4double x,y,z;
    G4int nFaceVert, v1, v2, v3;

    ifstr >> dummy >> nVertex >> nFace >> nEdge;

    for(G4int iVert=0; iVert < nVertex; iVert++) {
        ifstr >> x >> y >> z;
        points->InsertNextPoint(x,y,z);
    }

    for(G4int iFace=0; iFace < nFace; iFace++) {
        ifstr >> nFaceVert >> v1 >> v2 >> v3;
        vtkNew<vtkTriangle> tri;
        tri->GetPointIds()->SetId(0, v1);
        tri->GetPointIds()->SetId(1,v2);
        tri->GetPointIds()->SetId(2,v3);

        pd->InsertNextCell(tri->GetCellType(), tri->GetPointIds());
    }
}

void G4VtkSurfaceMeshLoader::Fill(G4int meshId, G4TessellatedSolid *tess) {
    auto v = pd->GetPoints();
    auto p = pd->GetPolys();

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

void G4VtkSurfaceMeshLoader::Fill(G4int meshId, G4TessellatedNew *tess) {
    auto v = pd->GetPoints();
    auto p = pd->GetPolys();

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

vtkSmartPointer<vtkUnstructuredGrid> G4VtkSurfaceMeshLoader::GetVolumeMesh() {
    vtkNew<vtkDelaunay3D> delaunay;
    delaunay->SetInputData(pd);
    delaunay->SetAlpha(0);
    delaunay->Update();
    return delaunay->GetOutput();
}