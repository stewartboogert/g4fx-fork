#pragma once

#include "vtkNew.h"

class vtkActor;

template<typename storageType, typename filterType, typename mapperType>
struct VtkPipeline {
    VtkPipeline() {};
    ~VtkPipeline() {};

    vtkNew<storageType> data;
    vtkNew<filterType> filter;
    vtkNew<mapperType> mapper;
    vtkNew<vtkActor> actor;
};