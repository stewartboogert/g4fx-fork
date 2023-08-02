#pragma once

#include <vector>

#include "vtkNew.h"
#include "common/common3d.hh"

class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;

class VtkViewer {
public:
    VtkViewer();
    ~VtkViewer();

    void start();
    void clear();

    void setBackgroundColour(double r, double g, double b);
    void setBackgroundColour(vec3<double> colour);

    // 0-d objects
    void addPoint(vec3<double> p, vec4<double> colour = vec4<double>(1,1,1,1));

    // 1-d objects
    // void addRaySegment();
    // void addCurve();

    // 2-d objects
    //void addTriangle(std::vector<vec3> v,
    //                 triangle<unsigned int> i);
    // void addQuad();

    // 2-d manifold objects
    // void addMesh(const std::string mesh_key, double **verts, double **polys);
    // void addMeshInstance(const std::string mesh_key, double **transform);

    // 3-d objects
    // void addBox();
    // void addTet();
    // void addVolume();

    // void addSDF();
    // void addBVH();


private:
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
};

