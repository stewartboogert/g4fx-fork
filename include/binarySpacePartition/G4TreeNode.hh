//
// Created by Stewart Boogert on 27/03/2023.
//

#ifndef G4TREESTRUCTURE_HH
#define G4TREESTRUCTURE_HH

class G4TreeNode {
public:
    G4TreeNode() {};
    ~G4TreeNode() = default;
    void Set_Parent(G4TreeNode *parent) {mParent = parent;}
    void Add_Child(G4TreeNode *child) {mChildren.push_back(child);}

protected:

private:
    G4TreeNode *mParent;
    std::vector<G4TreeNode*> mChildren;
};

#endif //G4TREESTRUCTURE_HH
