//
// Created by Stewart Boogert on 24/03/2023.
//

#ifndef G4MESHLOADER_HH
#define G4MESHLOADER_HH

#include "G4Types.hh"
#include "G4String.hh"

class G4VSurfaceMesh;
class G4TessellatedSolid;

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

#include <vector>
#include <map>

class G4MeshLoader {
public:
  G4MeshLoader();
  ~G4MeshLoader() = default;

  void Load(G4String file_name);
  void Fill(G4int meshId, G4TessellatedSolid *tess);

protected:

private:
  vtkSmartPointer<vtkPolyData> pd;
};

#endif //G4MESHLOADER_HH
