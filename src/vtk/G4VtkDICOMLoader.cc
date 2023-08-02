//
// Created by Stewart Boogert on 10/04/2023.
//

#include "vtk/G4VtkDICOMLoader.hh"

#include "vtkNew.h"
#include "vtkDICOMImageReader.h"
#include "vtkImageData.h"

#include "G4String.hh"

G4VtkDICOMLoader::G4VtkDICOMLoader() {}

void G4VtkDICOMLoader::Load(G4String dir_name) {
    vtkNew<vtkDICOMImageReader> dicomReader;
    dicomReader->SetDirectoryName(dir_name);
    dicomReader->Update();
    image = dicomReader->GetOutput();
}