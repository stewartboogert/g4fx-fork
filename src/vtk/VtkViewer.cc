#include "vtk/VtkViewer.hh"

#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

VtkViewer::VtkViewer() {
    // Setup render window, renderer, and interactor
    renderWindow->SetWindowName("Debug viewer");
    renderWindow->AddRenderer(renderer);
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Setup pipelines
    // point pl

}

VtkViewer::~VtkViewer() {

}

void VtkViewer::start() {
    renderWindow->Render();
    renderWindowInteractor->Start();
}

void VtkViewer::clear() {

}

void VtkViewer::setBackgroundColour(double r, double g, double b) {
    renderer->SetBackground(r,g,b);
}

void VtkViewer::setBackgroundColour(vec3<double> colour) {
    setBackgroundColour(colour.x, colour.y, colour.z);
}

void VtkViewer::addPoint(vec3<double> p, vec4<double> colour) {
    std::cout << p.x << " " << p.y << " " << p.z << std::endl;
}