//
// Created by Boogert, Stewart (-,DL,-) on 27/07/2023.
//

#include "occCAD/G4OCCExplorer.hh"

#include "XCAFDoc.hxx"
#include "TDF_Tool.hxx"
#include "Standard_GUID.hxx"
#include "TDF_ChildIterator.hxx"
#include "TopExp_Explorer.hxx"
#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "Geom_Plane.hxx"
#include "Geom_CylindricalSurface.hxx"

// G4OCCExplorer::G4OCCExplorer() {}

// G4OCCExplorer::~G4OCCExplorer() {}

void G4OCCExplorer::Dump(opencascade::handle<XCAFDoc_ShapeTool> shapeTool, std::ostream &ostr) {
    /* Free shapes */

    /* Dump shape tool */
    shapeTool->Dump(ostr);
}

opencascade::handle<TDataStd_Name> G4OCCExplorer::GetNameFromLabel(const TDF_Label &label) {
    auto nameID = TDataStd_Name::GetID();
    opencascade::handle<TDataStd_Name> name;
    label.FindAttribute(nameID, name);
    return name;
}

std::string G4OCCExplorer::GetEntryFromLabel(const TDF_Label &label)
{
    TCollection_AsciiString anEntry;
    TDF_Tool::Entry(label,anEntry);
    return anEntry.ToCString();
}

opencascade::handle<XCAFDoc_Location> G4OCCExplorer::GetLocationFromLabel(const TDF_Label &label) {
    auto locationID = XCAFDoc_Location::GetID();
    opencascade::handle<XCAFDoc_Location> location;
    auto b = label.FindAttribute(locationID, location);
    if(b)
        return location;
    else
        return nullptr;
}

TDF_Label G4OCCExplorer::FindLabelByEntryString(const TDF_Label &topLabel, const std::string &entry) {
    for (TDF_ChildIterator itall (topLabel,Standard_True); itall.More(); itall.Next()) {
        auto aChild = itall.Value();
        if ( GetEntryFromLabel(aChild) == entry) {
            return aChild;
        }
    }
    return TDF_Label();
}

TDF_LabelSequence G4OCCExplorer::FindLabelByNameString(const TDF_Label &topLabel, const std::string &name) {

    TDF_LabelSequence labelSeq;
    for (TDF_ChildIterator itall (topLabel,Standard_True); itall.More(); itall.Next()) {
        auto aChild = itall.Value();
        if (GetNameFromLabel(aChild)->Get().IsEqual(name.c_str())) {
            labelSeq.Append(aChild);
        }
    }
    return labelSeq;
}

void G4OCCExplorer::DumpShape(const TopoDS_Shape &shape, std::ostream &ostr) {
    ostr << "DumpShape" << std::endl;

    auto exp = TopExp_Explorer(shape, TopAbs_FACE, TopAbs_VERTEX);

    while(exp.More()) {
        auto current = exp.Current();
        DumpFace(TopoDS::Face(current), ostr);
        exp.Next();
    }
}

void G4OCCExplorer::DumpFace(const TopoDS_Face &face, std::ostream &ostr) {
    ostr << "DumpFace" << std::endl;

    auto s = BRep_Tool::Surface(face);

    if (s->IsKind("Geom_Plane")) {
        auto ps = dynamic_cast<Geom_Plane*>(s->This());
        double a,b,c,d;
        ps->Coefficients(a,b,c,d);
        ostr << "Face: plane " << a << " " << b << " " << c << " " << d << std::endl;
    }
    else if(s->IsKind("Geom_CylindricalSurface")) {
        auto ps = dynamic_cast<Geom_CylindricalSurface*>(s->This());
        double a1, a2, a3, b1, b2, b3, c1, c2, c3 ,d;
        ps->Coefficients(a1,a2, a3, b1, b2, b3, c1, c2, c3 ,d);
        ostr << "Face: cylinder " << a1 << " " << a2 << " " << a3 << " "
                                  << b1 << " " << b2 << " " << b3 << " "
                                  << c1 << " " << c2 << " " << c3 << " " << d << std::endl;
    }
    else {
        ostr << "Face: unknown surface" << std::endl;
    }

    auto exp = TopExp_Explorer(face, TopAbs_WIRE, TopAbs_VERTEX);

    int nWire = 0;
    while(exp.More()) {
        auto current = exp.Current();
        //current.DumpJson(std::cout);
        //std::cout << std::endl;

        DumpWire(TopoDS::Wire(current), ostr);
        exp.Next();
        nWire++;
    }

    ostr << "DumpFace " << nWire << std::endl;
}

void G4OCCExplorer::DumpWire(const TopoDS_Wire &wire, std::ostream &ostr) {
    ostr << "DumpWire" << std::endl;
    auto exp = TopExp_Explorer(wire, TopAbs_EDGE, TopAbs_VERTEX);
    while(exp.More()) {
        auto current = exp.Current();

        DumpEdge(TopoDS::Edge(current), ostr);
        exp.Next();
    }
}

void G4OCCExplorer::DumpEdge(const TopoDS_Edge &edge, std::ostream &ostr) {
    ostr << "DumpEdge" << std::endl;
    auto exp = TopExp_Explorer(edge, TopAbs_VERTEX, TopAbs_SHAPE);

    double start, end;
    auto c = BRep_Tool::Curve(edge,start, end);
    ostr << "Edge: curve "  << start << " " << end << std::endl;

    if (c->IsKind("Geom_Line")) {
        ostr << "Edge: line" << std::endl;
    }
    else if (c->IsKind("Geom_Circle")) {
        ostr << "Edge: circle" << std::endl;
    }
    else if(c->IsKind("Geom_Ellipse")) {
        ostr << "Edge: ellipse" << std::endl;
    }
    else if(c->IsKind("Geom_Hyperbola")) {
        ostr << "Edge: hyperbola" << std::endl;
    }
    else if(c->IsKind("Geom_Parabola")) {
        ostr << "Edge: parabola" << std::endl;
    }
    else {
        ostr << "Edge: unknown curve" << std::endl;
    }

    while(exp.More()) {
        auto current = exp.Current();
        auto v = TopoDS::Vertex(current);

        auto pnt = BRep_Tool::Pnt(v);
        ostr << "Vertex " << pnt.X() << " " << pnt.Y() << " " << pnt.Z() << std::endl;

        exp.Next();
    }
}