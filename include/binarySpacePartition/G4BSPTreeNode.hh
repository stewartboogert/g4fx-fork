//
// Created by Stewart Boogert on 14/04/2023.
//


#ifndef G4BSPTREENODE_HH
#define G4BSPTREENODE_HH


#include "G4Types.hh"
#include "G4ThreeVector.hh"

#include <vector>



class G4BSPTreeNode {
public:
    G4BSPTreeNode();
    G4BSPTreeNode(G4ThreeVector normal, G4double distance);
    ~G4BSPTreeNode();
    void Set_Parent(G4BSPTreeNode *parent);
    void Add_Child(G4BSPTreeNode *child);
    void Set_Plane(G4ThreeVector normal, G4double distance);

    static G4BSPTreeNode* Cube(G4ThreeVector size);
protected:

private:
    G4BSPTreeNode *fParent;
    std::vector<G4BSPTreeNode*> fChildren;
    G4ThreeVector fNormal;
    G4double fDistance;
};


#endif //G4BSPTREENODE_HH
