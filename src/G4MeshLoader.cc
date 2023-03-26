//
// Created by Stewart Boogert on 24/03/2023.
//

#include "G4MeshLoader.hh"

#include "G4Types.hh"
#include "G4ios.hh"
#include "G4String.hh"
#include "G4VSurfaceMesh.hh"
#include "G4TessellatedSolid.hh"

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags
#include <assimp/mesh.h>     // Post processing flags


G4MeshLoader::G4MeshLoader() {}

void G4MeshLoader::Load(G4String file_name) {
    Assimp::Importer importer;

    scene = importer.ReadFile(file_name,
                              aiProcess_CalcTangentSpace       |
                              aiProcess_Triangulate            |
                              aiProcess_JoinIdenticalVertices  |
                              aiProcess_SortByPType);

    ProcessScene();
}

void G4MeshLoader::ProcessScene() {
    ProcessNode(scene->mRootNode, scene);
}

void G4MeshLoader::ProcessNode(aiNode *node, const aiScene *scene) {

    auto node_name = G4String(node->mName.C_Str());

    if(node_name == "") {
        node_name = G4String("NoName");
    }

    G4cout << "G4MeshLoader::ProcessNode> " << node_name << G4endl;

    G4int iMesh = 0;
    // Process all the node's meshes (if any)
    for (unsigned int i = 0; i < node->mNumMeshes; i++) {
        const aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];

        if (meshMap.find(G4String(node_name)) == meshMap.end()) {
            G4cout << "G4MeshLoader::ProcessNode> add new map element " << node_name << " " << mesh << G4endl;
            meshMap.insert(std::make_pair(G4String(node_name), std::vector<const aiMesh *>()));
        }
        meshMap[G4String(node_name)].push_back(mesh);

        iMesh++;

        for(auto i : meshMap) {
            G4cout << "G4MeshLoader::ProcessNode> map " << i.first << " " << i.second[0] << " " << i.second.size() << G4endl;
        }
    }

    // Recursively process each of the node's children
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        ProcessNode(node->mChildren[i], scene);
    }
}

void G4MeshLoader::Fill(G4String nodeName, G4int meshIndex, G4VSurfaceMesh *)
{}

void G4MeshLoader::Fill(G4String nodeName, G4int meshIndex, G4TessellatedSolid *) {
    G4cout << "G4MeshLoader::Fill>" << G4endl;
    for(auto i : meshMap) {
        G4cout << "G4MeshLoader::Fill> " << i.first << " " << i.second[0] << " " << i.second.size() << G4endl;
    }

    G4cout << "G4MeshLoader::Fill> " << nodeName << " " << meshMap[nodeName].size() << G4endl;
    auto mesh = meshMap[nodeName][meshIndex];
    G4cout << "G4MeshLoader::Fill> " << mesh << " " << mesh->mNumFaces << " " << mesh->mNumVertices << G4endl;

    for(auto iVert = 0; iVert < mesh->mNumVertices; iVert++) {
        auto vert = mesh->mVertices[iVert];
        G4cout << "iVertex " << iVert << " " << vert.x << " " << vert.y << " " << vert.z << " " << G4endl;
    }

    for(auto iFace = 0; iFace < mesh->mNumFaces; iFace++) {
        auto face = mesh->mFaces[iFace];
        G4cout << iFace << " " << face.mNumIndices << " ... ";
        for(auto iIndex = 0; iIndex < face.mNumIndices; iIndex++) {
            G4cout << "iIndex " << iIndex << " " << face.mIndices[0] << " " << face.mIndices[1] << " ";
        }
        G4cout << G4endl;
    }
}
