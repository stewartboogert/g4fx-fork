#pragma once

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

    static void Dump(opencascade::handle<XCAFDoc_ShapeTool> shapeTool, std::ostream &ostr);
    static opencascade::handle<TDataStd_Name> GetNameFromLabel(const TDF_Label &label);
    static opencascade::handle<XCAFDoc_Location> GetLocationFromLabel(const TDF_Label &label);
    static std::string GetEntryFromLabel(const TDF_Label &label);

    static TDF_Label FindLabelByEntryString(const TDF_Label &topLabel, const std::string &entry);
    static TDF_LabelSequence FindLabelByNameString(const TDF_Label &topLabel, const std::string &labelName);

    static void DumpShape(const TopoDS_Shape &shape, std::ostream &ostr);
    static void DumpFace(const TopoDS_Face &face, std::ostream &ostr);
    static void DumpWire(const TopoDS_Wire &wire, std::ostream &ostr);
    static void DumpEdge(const TopoDS_Edge &edge, std::ostream &ostr);
};

