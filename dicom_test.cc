//
// Created by Stewart Boogert on 10/04/2023.
//

#include "G4Types.hh"
#include "G4String.hh"

#include "G4VtkDICOMLoader.hh"

// --------------------------------------------------------------

int main(int argc,char **argv) {
    auto dcr = new G4VtkDICOMLoader();
    dcr->Load("./2_skull_ct/DICOMDIR");
}