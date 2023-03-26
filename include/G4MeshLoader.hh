//
// Created by Stewart Boogert on 24/03/2023.
//

#ifndef G4MESHLOADER_HH
#define G4MESHLOADER_HH

#include "G4Types.hh"
#include "G4String.hh"

class G4VSurfaceMesh;
class G4TessellatedSolid;
class aiScene;
class aiNode;
class aiMesh;

#include <vector>
#include <map>

class G4MeshLoader {
public:
  G4MeshLoader();
  ~G4MeshLoader() = default;

  void Load(G4String file_name);
  void Fill(G4String nodeName, G4int meshIndex, G4VSurfaceMesh *mesh);
  void Fill(G4String nodeName, G4int meshIndex, G4TessellatedSolid *tess);
protected:

private:
  void ProcessScene();
  void ProcessNode(aiNode* node, const aiScene* scene);

  const aiScene* scene;
  std::map<G4String, std::vector<const aiMesh*>> meshMap;

};

#endif //G4MESHLOADER_HH
