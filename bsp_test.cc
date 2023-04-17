#include <math.h>

#include "G4Types.hh"

#include "G4ThreeVector.hh"
#include "G4BSPTreeNode.hh"

// --------------------------------------------------------------

int main(int argc,char **argv)
{
    G4BSPTreeNode *bsp = G4BSPTreeNode::Cube(G4ThreeVector(5,5,5));
    return 0;
}

