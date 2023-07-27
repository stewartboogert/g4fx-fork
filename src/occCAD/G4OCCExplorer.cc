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

void G4OCCExplorer::Dump(Handle(XCAFDoc_ShapeTool) shapeTool) {
    /* Free shapes */

    /* Dump shape tool */
    shapeTool->Dump(std::cout);
}

opencascade::handle<TDataStd_Name> G4OCCExplorer::GetNameFromLabel(TDF_Label &label) {
    auto nameID = TDataStd_Name::GetID();
    Handle(TDataStd_Name) name;
    label.FindAttribute(nameID, name);
    return name;
}

std::string G4OCCExplorer::GetEntryFromLabel(TDF_Label &label)
{
    TCollection_AsciiString anEntry;
    TDF_Tool::Entry(label,anEntry);
    return anEntry.ToCString();
}

opencascade::handle<XCAFDoc_Location> G4OCCExplorer::GetLocationFromLabel(TDF_Label &label) {
    auto locationID = XCAFDoc_Location::GetID();
    Handle(XCAFDoc_Location) location;
    auto b = label.FindAttribute(locationID, location);
    if(b)
        return location;
    else
        return nullptr;
}

TDF_Label G4OCCExplorer::FindLabelByEntryString(TDF_Label topLabel, std::string entry) {
    for (TDF_ChildIterator itall (topLabel,Standard_True); itall.More(); itall.Next()) {
        auto aChild = itall.Value();
        if ( GetEntryFromLabel(aChild) == entry) {
            return aChild;
        }
    }
    return TDF_Label();
}

TDF_LabelSequence G4OCCExplorer::FindLabelByNameString(TDF_Label topLabel, std::string name) {

    TDF_LabelSequence labelSeq;
    for (TDF_ChildIterator itall (topLabel,Standard_True); itall.More(); itall.Next()) {
        auto aChild = itall.Value();
        if (GetNameFromLabel(aChild)->Get().IsEqual(name.c_str())) {
            labelSeq.Append(aChild);
        }
    }
    return labelSeq;
}

void G4OCCExplorer::DumpShape(TopoDS_Shape &shape) {
    std::cout << "DumpShape" << std::endl;

    auto exp = TopExp_Explorer(shape, TopAbs_FACE, TopAbs_VERTEX);

    while(exp.More()) {
        auto current = exp.Current();

        //current.DumpJson(std::cout);
        //std::cout << std::endl;

        DumpFace(TopoDS::Face(current));

        //std::cout << current.ShapeType() << " " << current.Closed() << " " << current.Infinite() << " " << current.Convex() << std::endl;
        exp.Next();
    }
}

void G4OCCExplorer::DumpFace(TopoDS_Face &face) {
    std::cout << "DumpFace" << std::endl;

    auto s = BRep_Tool::Surface(face);

    if (s->IsKind("Geom_Plane")) {
        auto ps = dynamic_cast<Geom_Plane*>(s->This());
        double a,b,c,d;
        ps->Coefficients(a,b,c,d);
        std::cout << "Face: plane " << a << " " << b << " " << c << " " << d << std::endl;

    }
    else if(s->IsKind("Geom_CylindricalSurface")) {
        auto ps = dynamic_cast<Geom_CylindricalSurface*>(s->This());
        double a1, a2, a3, b1, b2, b3, c1, c2, c3 ,d;
        ps->Coefficients(a1,a2, a3, b1, b2, b3, c1, c2, c3 ,d);
        std::cout << "Face: cylinder " << a1 << " " << a2 << " " << a3 << " "
                                       << b1 << " " << b2 << " " << b3 << " "
                                       << c1 << " " << c2 << " " << c3 << " " << d << std::endl;
    }
    else {
        std::cout << "Face: unknown surface" << std::endl;
    }

    auto exp = TopExp_Explorer(face, TopAbs_WIRE, TopAbs_VERTEX);

    int nWire = 0;
    while(exp.More()) {
        auto current = exp.Current();
        //current.DumpJson(std::cout);
        //std::cout << std::endl;

        DumpWire(TopoDS::Wire(current));
        exp.Next();
        nWire++;
    }

    std::cout << "DumpFace " << nWire << std::endl;


    //std::cout << "Surface " << s->get_type_name() << " " << s->get_type_descriptor() << std::endl;
    //s->DumpJson(std::cout);
    //std::cout << std::endl;

}

void G4OCCExplorer::DumpWire(TopoDS_Wire &wire) {
    std::cout << "DumpWire" << std::endl;
    auto exp = TopExp_Explorer(wire, TopAbs_EDGE, TopAbs_VERTEX);
    while(exp.More()) {
        auto current = exp.Current();
        //current.DumpJson(std::cout);
        //std::cout << std::endl;

        DumpEdge(TopoDS::Edge(current));
        exp.Next();
    }
}

void G4OCCExplorer::DumpEdge(TopoDS_Edge &edge) {
    std::cout << "DumpEdge" << std::endl;
    auto exp = TopExp_Explorer(edge, TopAbs_VERTEX, TopAbs_SHAPE);

    double start, end;
    auto c = BRep_Tool::Curve(edge,start, end);
    std::cout << "Edge: curve "  << start << " " << end << std::endl;

    if (c->IsKind("Geom_Line")) {
        std::cout << "Edge: line" << std::endl;
    }
    else if (c->IsKind("Geom_Circle")) {
        std::cout << "Edge: circle" << std::endl;
    }
    else if(c->IsKind("Geom_Ellipse")) {
        std::cout << "Edge: ellipse" << std::endl;
    }
    else if(c->IsKind("Geom_Hyperbola")) {
        std::cout << "Edge: hyperbola" << std::endl;
    }
    else if(c->IsKind("Geom_Parabola")) {
        std::cout << "Edge: parabola" << std::endl;
    }
    else {
        std::cout << "Edge: unknown curve" << std::endl;
    }

    while(exp.More()) {
        auto current = exp.Current();
        auto v = TopoDS::Vertex(current);

        auto pnt = BRep_Tool::Pnt(v);
        std::cout << "Vertex " << pnt.X() << " " << pnt.Y() << " " << pnt.Z() << std::endl;


        //current.DumpJson(std::cout);
        //std::cout << std::endl;
        exp.Next();
    }
}