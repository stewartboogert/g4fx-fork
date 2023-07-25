//
// Created by Stewart Boogert on 10/04/2023.
//

#ifndef G4VTKDICOMLOADER_HH
#define G4VTKDICOMLOADER_HH

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

#endif //G4VTKDICOMLOADER_HH
