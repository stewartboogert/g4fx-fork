//
// Created by Stewart Boogert on 08/04/2023.
//

#include "vtk/G4VtkVolumeMeshLoader.hh"

#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkActor.h"
#include "vtkDataSetMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkExtractCells.h"
#include "vtkProperty.h"
#include "vtkPlane.h"
#include "vtkClipDataSet.h"
#include "vtk3DLinearGridCrinkleExtractor.h"
#include "vtkMeshQuality.h"

G4VtkVolumeMeshLoader::G4VtkVolumeMeshLoader() {
    points = vtkSmartPointer<vtkPoints>::New();
    ug     = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ug->SetPoints(points);
}

void G4VtkVolumeMeshLoader::Load(G4String file_name) {
    if(file_name.find("mesh") != std::string::npos) {

    }
    else if(file_name.find("vtk") !=std::string::npos) {
        auto ugReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        ugReader->SetFileName(file_name);
        ugReader->Update();
        ug = ugReader->GetOutput();
    }
}

void  G4VtkVolumeMeshLoader::SetUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> ugIn) {
    ug = ugIn;
}


void G4VtkVolumeMeshLoader::View() {

    // Renderer
    vtkNew<vtkRenderer> renderer;

    // Render window
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(640, 480);

    // Render window interactor
    vtkNew<vtkRenderWindowInteractor> interactor;
    interactor->SetRenderWindow(renderWindow);

    // tetrahedral cells only
    vtkNew<vtkIdList> tetIdList;
    vtkNew<vtkExtractCells> extract;

    auto ugCells = ug->GetCells();
    for(auto iCell =0; iCell < ugCells->GetNumberOfCells(); iCell++ ) {
        vtkNew<vtkIdList> ids;
        ugCells->GetCellAtId(iCell, ids);
        if(ids->GetNumberOfIds() == 4) {
            tetIdList->InsertNextId(iCell);
        }
    }
    extract->SetCellList(tetIdList);
    extract->SetInputData(ug);

    vtkNew<vtkPlane> clipPlane;
    clipPlane->SetOrigin(ug->GetCenter());
    clipPlane->SetOrigin(0,0,-1000);
    clipPlane->SetNormal(0,0,1);

    vtkNew<vtk3DLinearGridCrinkleExtractor> clipper;
    clipper->SetImplicitFunction(clipPlane);
    //clipper->SetInputConnection(extract->GetOutputPort());
    clipper->SetInputData(ug);
    clipper->Update();

    vtkNew<vtkClipDataSet> clipper2;
    clipper2->SetClipFunction(clipPlane);
    //clipper2->SetInputConnection(extract->GetOutputPort());
    clipper2->SetInputData(ug);
    clipper2->SetValue(0.0);
    clipper2->GenerateClippedOutputOn();
    clipper2->Update();

    // mapper
    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputConnection(clipper->GetOutputPort());

    auto mapper2 = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper2->SetInputConnection(clipper2->GetOutputPort());

    // create actor
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->SetVisibility(1);
    actor->GetProperty()->SetRepresentationToSurface();

    auto actor2 = vtkSmartPointer<vtkActor>::New();
    actor2->SetMapper(mapper2);
    actor2->SetVisibility(1);
    actor2->GetProperty()->SetRepresentationToSurface();

    // add to renderer
    //renderer->AddActor(actor);
    renderer->AddActor(actor2);

    // render
    renderWindow->Render();

    // start interaction
    interactor->Start();
}

void G4VtkVolumeMeshLoader::MeshQuality() {
    vtkNew<vtkMeshQuality> quality;
    quality->SetInputData(ug);
}