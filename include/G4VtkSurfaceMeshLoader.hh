//
// Created by Stewart Boogert on 24/03/2023.
//

#ifndef G4VTKSURFACEMESHLOADER_HH
#define G4VTKSURFACEMESHLOADER_HH

#include "G4Types.hh"
#include "G4String.hh"

class G4VSurfaceMesh;
class G4TessellatedSolid;
class G4TessellatedNew;

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"

#include <vector>
#include <map>

class G4VtkSurfaceMeshLoader {
public:
  G4VtkSurfaceMeshLoader();
  ~G4VtkSurfaceMeshLoader() = default;

  void Load(G4String file_name);
  void LoadOFF(G4String file_name);

  void Fill(G4int meshId, G4TessellatedSolid *tess);
  void Fill(G4int meshId, G4TessellatedNew *tess);

  vtkSmartPointer<vtkPolyData> GetPolyData() {return pd;}
  void SetPolyData(vtkSmartPointer<vtkPolyData> pdIn) {pd = pdIn;}

  vtkSmartPointer<vtkUnstructuredGrid> GetVolumeMesh();

  void View();

protected:

private:
  vtkSmartPointer<vtkPolyData> pd;
  vtkSmartPointer<vtkPoints> points;
};

#endif //G4VTKSURFACEMESHLOADER_HH
