#pragma once

#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"

struct VtkPipeline {
    VtkPipeline();
    ~VtkPipeline();

    void buildPipeline();

    vtkNew<vtkPoints> polydataPoints;
    vtkNew<vtkCellArray> polydataCells;
    vtkNew<vtkPolyData> polydata;
    std::vector<vtkSmartPointer<vtkPolyDataAlgorithm>> polyFilters;
    vtkNew<vtkPolyDataMapper> mapper;
    vtkNew<vtkActor> actor;
};