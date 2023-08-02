#include "vtk/VtkViewer.hh"

int main() {

    auto vv = new VtkViewer();
    vv->setBackgroundColour(1,1,1);
    vv->start();

    return 0;
}