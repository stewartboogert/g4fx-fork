//
// Created by Boogert, Stewart (-,DL,-) on 26/07/2023.
//

#include "occCAD/G4OCCStepLoader.hh"
#include "occCAD/G4OCCExplorer.hh"

int main() {
    auto step_reader = G4OCCStepLoader("03_Cube_with_fillet.step");

    G4OCCExplorer::Dump(step_reader.GetShapeTool());

    TDF_LabelSequence freeShapes;
    step_reader.GetShapeTool()->GetFreeShapes(freeShapes);

    std::cout << "Number of free shapes " << freeShapes.Size() << std::endl;
    auto fs1 = freeShapes.Value(1);

    std::cout << "Find by label entry" << std::endl;
    auto entry1 = G4OCCExplorer::FindLabelByEntryString(step_reader.GetMainLabel(),"0:1:1:1");
    entry1.Dump(std::cout);

    std::cout << "Find by label name" << std::endl;
    auto entry2 = G4OCCExplorer::FindLabelByNameString(step_reader.GetMainLabel(),"FusionComponent");

    for(auto l : entry2) {
        l.Dump(std::cout);
    }

    auto shape = step_reader.GetShapeTool()->GetShape(entry1);

    G4OCCExplorer::DumpShape(shape);

    return 0;
}