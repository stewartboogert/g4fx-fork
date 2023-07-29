#include <catch2/catch_all.hpp>

#include <occCAD/G4OCCStepLoader.hh>
#include <occCAD/G4OCCExplorer.hh>

TEST_CASE( "CAD file loading and exploration", "[occ]" ) {
    auto step_reader = G4OCCStepLoader("03_Cube_with_fillet.step");
    auto shape_tool = step_reader.GetShapeTool();
    auto entry1 = G4OCCExplorer::FindLabelByEntryString(step_reader.GetMainLabel(),"0:1:1:1");
    auto shape = shape_tool->GetShape(entry1);

    SECTION("Get document free shapes") {
        TDF_LabelSequence freeShapes;
        step_reader.GetShapeTool()->GetFreeShapes(freeShapes);
        REQUIRE( freeShapes.Size() == 1 );
    }

    SECTION("Find label by entry") {
        REQUIRE( G4OCCExplorer::GetNameFromLabel(entry1)->Get() == "FusionComponent" );
    }

    SECTION("Find label by name") {
        auto entry2 = G4OCCExplorer::FindLabelByNameString(step_reader.GetMainLabel(),"FusionComponent");
        auto label = entry2.First();
        REQUIRE( G4OCCExplorer::GetNameFromLabel(label)->Get() == "FusionComponent" );
    }

    SECTION("Get shape tool") {
        REQUIRE(!shape_tool.IsNull());
    }

    SECTION("Get shape") {
        REQUIRE(!shape.IsNull());
    }



}
