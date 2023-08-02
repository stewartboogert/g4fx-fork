#pragma once

#include "G4Types.hh"
#include "G4String.hh"

#include "vtkSmartPointer.h"

class vtkPoints;
class vtkUnstructuredGrid;

class G4VtkVolumeMeshLoader {
public:
    G4VtkVolumeMeshLoader();
    ~G4VtkVolumeMeshLoader() = default;

    void Load(G4String file_name);

    vtkSmartPointer<vtkUnstructuredGrid> GetUnstructuredGrid() {return ug;}
    void SetUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> ugIn);

    void View();
    void MeshQuality();

protected:

private:

    vtkSmartPointer<vtkUnstructuredGrid> ug;
    vtkSmartPointer<vtkPoints> points;

};