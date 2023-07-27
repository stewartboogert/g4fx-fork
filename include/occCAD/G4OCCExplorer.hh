//
// Created by Boogert, Stewart (-,DL,-) on 27/07/2023.
//

#ifndef GEANT4FX_G4OCCEXPLORER_HH
#define GEANT4FX_G4OCCEXPLORER_HH

#include <XCAFDoc_ShapeTool.hxx>

#include <TDataStd_Name.hxx>
#include <XCAFDoc_Location.hxx>
#include <TDataStd_TreeNode.hxx>

#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>

class G4OCCExplorer {
public:
    //G4OCCExplorer();
    //~G4OCCExplorer();

    static void Dump(Handle(XCAFDoc_ShapeTool) shapeTool);
    static opencascade::handle<TDataStd_Name> GetNameFromLabel(TDF_Label &label);
    static opencascade::handle<XCAFDoc_Location> GetLocationFromLabel(TDF_Label &label);
    static std::string GetEntryFromLabel(TDF_Label &label);

    static TDF_Label FindLabelByEntryString(TDF_Label topLabel, std::string entry);
    static TDF_LabelSequence FindLabelByNameString(TDF_Label topLabel, std::string labelName);

    static void DumpShape(TopoDS_Shape &shape);
    static void DumpFace(TopoDS_Face &face);
    static void DumpWire(TopoDS_Wire &wire);
    static void DumpEdge(TopoDS_Edge &edge);
};

#endif //GEANT4FX_G4OCCEXPLORER_HH
