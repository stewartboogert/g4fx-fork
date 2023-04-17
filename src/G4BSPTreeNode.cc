//
// Created by Stewart Boogert on 14/04/2023.
//

#include "G4BSPTreeNode.hh"

#include "G4ThreeVector.hh"

G4BSPTreeNode::G4BSPTreeNode() :
    fNormal(0,0,1),
    fDistance(0) {};

G4BSPTreeNode::G4BSPTreeNode(G4ThreeVector normal, G4double distance) :
    fNormal(normal),
    fDistance(distance)  {};

G4BSPTreeNode::~G4BSPTreeNode() {
    for(auto cn : fChildren) {
        delete cn;
    }
}

void G4BSPTreeNode::Set_Parent(G4BSPTreeNode *parent) {
    fParent = parent;
}

void G4BSPTreeNode::Add_Child(G4BSPTreeNode *child) {
    fChildren.push_back(child);
}

void G4BSPTreeNode::Set_Plane(G4ThreeVector normal, G4double distance) {
    fNormal = normal;
    fDistance = distance;
};

G4BSPTreeNode* G4BSPTreeNode::Cube(G4ThreeVector size) {
    G4BSPTreeNode *n1 = new G4BSPTreeNode(G4ThreeVector(1,0,0),size.x());
    G4BSPTreeNode *n2 = new G4BSPTreeNode(G4ThreeVector(0,1,0),size.y());
    G4BSPTreeNode *n3 = new G4BSPTreeNode(G4ThreeVector(0,0,1),size.z());
    G4BSPTreeNode *n4 = new G4BSPTreeNode(G4ThreeVector(-1,0,0),size.x());
    G4BSPTreeNode *n5 = new G4BSPTreeNode(G4ThreeVector(0,-1,0),size.y());
    G4BSPTreeNode *n6 = new G4BSPTreeNode(G4ThreeVector(0,0,-1),size.z());

    n1->Add_Child(n2);
    n2->Add_Child(n3);
    n3->Add_Child(n4);
    n4->Add_Child(n5);
    n5->Add_Child(n6);

    return n1;
}

