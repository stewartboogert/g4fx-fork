//
// Created by Boogert, Stewart (-,DL,-) on 26/07/2023.
//

#include "occCAD/G4OCCStepLoader.hh"

#include "TCollection_ExtendedString.hxx"
#include "STEPCAFControl_Reader.hxx"
#include "Message_ProgressRange.hxx"
#include "XCAFDoc_DocumentTool.hxx"

G4OCCStepLoader::G4OCCStepLoader(const std::string &file_name) {

    /* Create application */
    app = XCAFApp_Application::GetApplication();

    /* Create new document */
    app->NewDocument(TCollection_ExtendedString("MDTV-CAF"),doc);

    /* Open and read step file */
    auto stepReader = STEPCAFControl_Reader();
    auto mpr = Message_ProgressRange();
    stepReader.ReadFile(file_name.c_str());

    /* Transfer to XCAF document */
    stepReader.Transfer(doc, mpr);

    /* Document main label */
    mainLabel = doc->Main();

    /* Create shape tool */
    shapeTool = XCAFDoc_DocumentTool::ShapeTool(mainLabel);
}

G4OCCStepLoader::~G4OCCStepLoader() {}

TDF_Label G4OCCStepLoader::GetMainLabel() {
    return mainLabel;
}

opencascade::handle<XCAFDoc_ShapeTool> G4OCCStepLoader::GetShapeTool() {
    return shapeTool;
}