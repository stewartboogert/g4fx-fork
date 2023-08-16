#include <vtk/VtkPipeline.hh>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>

VtkPipeline::VtkPipeline() {}

VtkPipeline::~VtkPipeline() {}

void VtkPipeline::buildPipeline() {
    polydata->SetPoints(polydataPoints);
    polydata->SetPolys(polydataCells);

    // clean input polydata
    auto filterClean = vtkSmartPointer<vtkCleanPolyData>::New();
    filterClean->PointMergingOn();
    filterClean->AddInputData(polydata);
    polyFilters.push_back(filterClean);

    // ensure triangular mesh
    auto filterTriangle = vtkSmartPointer<vtkTriangleFilter>::New();
    filterTriangle->SetInputConnection(filterClean->GetOutputPort());
    polyFilters.push_back(filterTriangle);

    // calculate normals with a feature angle of 45 degrees
    auto filterNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
    filterNormals->SetFeatureAngle(45);
    filterNormals->SetInputConnection(filterTriangle->GetOutputPort());
    polyFilters.push_back(filterNormals);

    // mapper
    mapper->SetInputConnection(polyFilters.back()->GetOutputPort());
    mapper->SetColorModeToDirectScalars();

    // add to actor
    actor->SetMapper(mapper);
    actor->SetVisibility(1);

    // add actor to renderer

}
