#pragma once

#include "vtkSmartPointer.h"

class G4String;

class vtkImageData;

class G4VtkDICOMLoader {
public:
    G4VtkDICOMLoader();

    void Load(G4String dir_name);
protected:

private:
    vtkSmartPointer<vtkImageData> image;

};
