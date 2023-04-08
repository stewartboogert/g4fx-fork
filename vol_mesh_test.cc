#include "G4Types.hh"

#include "G4VtkVolumeMeshLoader.hh"
// --------------------------------------------------------------

int main(int argc,char **argv)
{
    G4VtkVolumeMeshLoader *mesh_load = new G4VtkVolumeMeshLoader();
    mesh_load->Load("T006_torus.vtk");
    mesh_load->View();
}
