#include "GSAlignmentSelection.h"
#include "myUtility.h"
#include "DataObject/Point.h"
#include "stdlib.h"
#include "DataObject/Alignment.h"



#include <vector>
#include <string>
#include <iostream>
using namespace std;


/*
 *Input : 
   1.  fileName: The name of output file
   2.  AlVct: A vector of alignment
 *Side effects:
    store the content of alignment vector to a json file
 *output
    if succeed return true ortherwise return false
 *
 */
/*
bool WriteAligns(char* filename,std::vector<Alignment>& AlVct)
{
    Json::Value root;
    Json::StyledWriter writer;
    
    //output target
    Json::Value tgt;
    tgt["name"] = AlVct[0].getTargetName();
    tgt["FS"] = AlVct[0].getTargetSequence();
    tgt["SL"] = AlVct[0].getTargetLength();
    tgt["PSS"] = AlVct[0].getPredictedSsInfo();
    tgt["PSSC"] = AlVct[0].getPredictedSsConf();
    tgt["PSA"] = AlVct[0].getPredictedSaInfo();

    //output aligns
    for(int i = 0; i < AlVct.size(); i++)
    {
        Json::Value tpl, algn;
        //output template
        std::string fileID;

        tpl["fileID"] = AlVct[i].getTargetName()+"_"+std::to_string(i)
                        +"_"+AlVct[i].getTemplateName() ;
        tpl["name"] = AlVct[i].getTemplateName();
        tpl["RFSI"] = AlVct[i].getTemplateReferenceSequenceInfo();
        tpl["RLSI"] = AlVct[i].getTemplateRealSequenceInfo();
        tpl["TSS"] = AlVct[i].getTemplateTrueSecondaryStructure();
        tpl["SL"] = AlVct[i].getTemplateSequenceLength();
        //output alignment
        algn["TMscore"] = AlVct[i].getTmScore();
        algn["expect"] = AlVct[i].getExpectedValue();
        //algn["identities"] = AlVct[i].getIdentities();
       // algn["positives"] = AlVct[i].getPositives();
        algn["gaps"] = AlVct[i].getGaps();
        algn["trimString"] = AlVct[i].getShiftedAlign(); 
        algn["queryStart"] = AlVct[i].getQueryStart();
        algn["queryPart"] = AlVct[i].getQueryPart();
        algn["queryEnd"] = AlVct[i].getQueryEnd();
        algn["subjectStart"] = AlVct[i].getSubjectStart();
        algn["subjectPart"] = AlVct[i].getSubjectPart();
        algn["subjectEnd"] = AlVct[i].getSubjectEnd();
        root["Aligns"][i]["Template"] = tpl;
        root["Aligns"][i]["Alignment"] = algn;
    }
    
    root["Target"] = tgt;

    ofstream outfile(filename);
    //open the file
    if(!outfile)
    {   
        cerr<<"Can not open file: "<<filename<<endl;
        cerr<<strerror(errno)<<endl;
        return false;
    }
    std::string outputJson = writer.write( root );
    outfile << outputJson;
    return true;
}*/





int main()
{
    vector<Alignment> AlVct;
    string dataPath ="/Users/ming/projects/MUfold/Data/GenerateExperimentDatasetOutput/casp11_fully/"; 
    string dataName ="T0759";
    string suffix = "_blast.json";

    string n2 ="test.json";




   // dataPath ="";
   // suffix="";
   // ReadAligns((dataPath+dataName+suffix).c_str(),AlVct);
   // ReadAligns(AlVct,dataPath,n2,suffix);
   
    
    string coordName="T0759/blast_T0759_0_1LM7_A.coords";
    
    //Point * p=NULL;
    //ReadCoord((dataPath+coordName).c_str(),p);

//    for(int i =0;i<AlVct[0].getSubjectEnd();i++)
  //      cout<<i<<':'<<p[i].getX()<<','<<p[i].getY()<<','<<p[i].getZ()<<endl;

    return 1;
}
