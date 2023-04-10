#include "G4Types.hh"

#include "G4VtkSurfaceMeshLoader.hh"
#include "G4VtkVolumeMeshLoader.hh"
// --------------------------------------------------------------

int main(int argc,char **argv)
{
    G4VtkSurfaceMeshLoader *surface_load = new G4VtkSurfaceMeshLoader();
    surface_load->Load("T006_torus.stl");
    surface_load->View();

    G4VtkVolumeMeshLoader *volume_load = new G4VtkVolumeMeshLoader();
    volume_load->SetUnstructuredGrid(surface_load->GetVolumeMesh());
    volume_load->View();
}
