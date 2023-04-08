//
// Created by Stewart Boogert on 08/04/2023.
//

#ifndef G4VTKVOLUMEMESHLOADER_HH
#define G4VTKVOLUMEMESHLOADER_HH

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

    void View();

protected:

private:

    vtkSmartPointer<vtkUnstructuredGrid> ug;
    vtkSmartPointer<vtkPoints> points;

};

#endif //G4VTKVOLUMEMESHLOADER_HH
