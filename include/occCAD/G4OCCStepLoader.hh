#pragma once

#include <XCAFApp_Application.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFDoc_ShapeTool.hxx>

class G4OCCStepLoader {
public:
    G4OCCStepLoader(const std::string &file_name);
    ~G4OCCStepLoader();
    TDF_Label GetMainLabel();
    opencascade::handle<XCAFDoc_ShapeTool> GetShapeTool();


private:
    opencascade::handle<XCAFApp_Application> app;
    opencascade::handle<TDocStd_Document> doc;
    TDF_Label mainLabel;
    opencascade::handle<XCAFDoc_ShapeTool> shapeTool;

};