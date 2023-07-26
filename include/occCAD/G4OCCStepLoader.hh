//
// Created by Boogert, Stewart (-,DL,-) on 26/07/2023.
//

#ifndef GEANT4FX_G4OCCSTEPLOADER_HH
#define GEANT4FX_G4OCCSTEPLOADER_HH

#include <XCAFApp_Application.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFDoc_ShapeTool.hxx>

class G4OCCStepLoader {
public:
    G4OCCStepLoader(std::string file_name);
    ~G4OCCStepLoader();
    TDF_Label GetMainLabel();
    Handle(XCAFDoc_ShapeTool) GetShapeTool();


private:
    Handle(XCAFApp_Application) app;
    Handle(TDocStd_Document) doc;
    TDF_Label mainLabel;
    Handle(XCAFDoc_ShapeTool) shapeTool;

};

#endif //GEANT4FX_G4OCCSTEPLOADER_HH
