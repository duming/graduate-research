#include "myUtility.hpp"
#include <iostream>
#include <fstream>
using namespace std;


void test()
{
    Json::Value root;  
    Json::Reader reader;
    vector<string> members;
    string filename = "test.json";
    ifstream infile(filename);
    if(!infile.is_open())
    {
        cerr<<"Can not open file: "<<filename<<endl;
        cerr<<strerror(errno)<<endl;
        return ;
    }
    
    bool parsingSuccessful = reader.parse( infile, root );

    if ( !parsingSuccessful )
    {
        std::cout << "Failed to parse configuration\n"<< reader.getFormattedErrorMessages();
        infile.close();
        return ;
                                   
    }

    root["Aligns"];


    Json::Value tmp =root[0];
    cout<<tmp.get("int",0).asString()<<endl;
    cout<<tmp.get("double",-1).asDouble()<<endl;
}

//Read json file into Alignment vector 
std::vector<std::string> ReadAligns(const char* filename, vector<Alignment>& AlVct)
{   
    Json::Value root;  
    Json::Reader reader;
    vector<string> members;
    ifstream infile(filename);
    std::vector<std::string> fileIDs;
    //open the file
    if(!infile.is_open())
    {
        cerr<<"Can not open file: "<<filename<<endl;
        cerr<<strerror(errno)<<endl;
        return fileIDs;
    }
    //parse the file
    bool parsingSuccessful = reader.parse( infile, root );

    if ( !parsingSuccessful )
    {
        std::cout << "Failed to parse configuration\n"<< reader.getFormattedErrorMessages();
        infile.close();
        return fileIDs;
                                   
    }   
    

    members = root.getMemberNames();
    Json::Value aligns =root[members[0]];
    
    AlVct.resize(aligns.size());
    for(int idx = 0; idx < aligns.size(); idx++ )
    {//asign values to Alignment instances
        Json::Value tmp =aligns[idx];

        //initialization of Target
        AlVct[idx].setTargetName(members[0]);
        AlVct[idx].setTargetSequence(tmp["targetFullSequence"].asString());
        AlVct[idx].setTargetLength(tmp["targetFullSequence"].asString().length());
        //AlVct[idx].setDifficulty(_query.getDifficulty());
        AlVct[idx].setPredictedSs(tmp["targetPredictedSecondaryStrucutre"].asString());
        //AlVct[idx].setPredictedSa(tmp["targetPredictedSecondaryStructureConfidence"].asString());
        //AlVct[idx].setProteinType(_query.getProteinType());
        //AlVct[idx].setNumAlpha(_query.getNumAlpha());
        //AlVct[idx].setNumBeta(_query.getNumBeta());
        //AlVct[idx].setNumCoil(_query.getNumCoil());


        //initialization of template

		AlVct[idx].setTemplateRealSequenceInfo(tmp.get("templateSequenceInfo", "").asString());
        AlVct[idx].setTemplateSequenceLength(atoi(tmp.get("templateSequenceLength","").asString().c_str()));

        //In the json file I got, there isnt any member named templateReferenceSequenceInfo or 
		//AlVct[idx].setTemplateReferenceSequenceInfo(tmp.get("templateReferenceSequenceInfo", "").  asString());
		//AlVct[idx].setTemplateRealSequenceInfo(tmp.get("templateRealSequenceInfo", "").asString());

        //initialization of Alignment
		AlVct[idx].setTemplateName(tmp.get("templateName", "").asString());
		if (tmp.get("isInDB", "").asString() == "1")
		{
		    AlVct[idx].isInDB = true;
		}
		else
		{
		    AlVct[idx].isInDB = false;
		}
		//AlVct[idx].setScore(atof(tmp.get("score", "").asString().c_str()));
		//AlVct[idx].setExpectedValue(atof(tmp.get("expect", "").asString().c_str()));
		//std::string identitiesString = tmp.get("identities", "").asString();
		//identitiesString = identitiesString.substr(0, identitiesString.size() - 1);
		//AlVct[idx].setIdentities(atof(identitiesString.c_str()));
		AlVct[idx].setQueryStart(atoi(tmp.get("queryStart", "").asString().c_str()));
		AlVct[idx].setQueryPart(tmp.get("queryPart", "").asString());
		AlVct[idx].setQueryEnd(atoi(tmp.get("queryEnd", "").asString().c_str()));
		AlVct[idx].setSubjectStart(atoi(tmp.get("subjectStart", "").asString().c_str()));
		AlVct[idx].setSubjectPart(tmp.get("subjectPart", "").asString());
		AlVct[idx].setSubjectEnd(atoi(tmp.get("subjectEnd", "").asString().c_str()));
		AlVct[idx].setShiftedAlign(tmp.get("shiftedAlign", "").asString());
		// AlignLength
		AlVct[idx].setAlignLength(AlVct[idx].getAlignLength());
		//Need clarification for the usage of marcro
        AlVct[idx].selectionFlag = CANDIDATE;
		//
        AlVct[idx].alignSource = tmp.get("alignSource", "").asString();
        fileIDs.push_back(tmp["fileID"].asString());
    }
    infile.close();
    return fileIDs;
}


bool ReadCoord(const char* filename, Point*& p)
{
    ifstream infile(filename);
    if(!infile)
    {
        cerr<<"can not open file:"<<filename<<endl;
        cerr<<strerror(errno)<<endl;
        return false;
    }
    
    vector<Point> pvct;
    Point tmp;
    
    double x,y,z;
    char cm1,cm2;
    int i=0;
    while(infile>>x>>cm1>>y>>cm2>>z)
    {
            tmp.setX(x);
            tmp.setY(y);
            tmp.setZ(z);
            pvct.push_back(tmp);
    }

    //allocate memory for p if p is empty
    if(p == NULL)
        p = new Point[pvct.size()];

    copy(pvct.begin(),pvct.end(),p);

    return true;
}



int countGap(std::string str)
{
    int strLen = str.length();
    int sum = 0;
    for(int i=0; i<strLen; i++)
        if( str[i] == '-')
            sum +=1;
    return sum;
}
