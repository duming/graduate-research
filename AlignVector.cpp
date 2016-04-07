#include "AlignVector.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////
//                          AlignVector part
///////////////////////////////////////////////////////////////////////////////////////////////

bool AlignVector::readAll(std::string dataPath, std::string fileName)
{
    std::vector<std::string> fileIDs;
    fileIDs = ReadAligns((dataPath+fileName).c_str(),AlVct);
    //fileIDs = readAligns((dataPath + fileName).c_str());

    vctlen = fileIDs.size();
    if(vctlen == 0)
    {
        std::cout<<"File:"<<dataPath+fileName<<" contains no aligns!"<<endl;
        return false;
    }

    info_3D.resize(vctlen);
    string secPath = AlVct[0].getTargetName()+"/AlignPDBs/";
    string suffix = ".coords";
    for(int i=0; i<vctlen; i++)
    {
        info_3D[i].readCoord((dataPath + secPath +fileIDs[i] + suffix).c_str());
        //make segement
        info_3D[i].makeSegment(AlVct[i]);
    }


    return true;
}


//read in the json file
std::vector<std::string> AlignVector::readAligns(const char* filename)
{
	Json::Value root;
	Json::Reader reader;
	std::ifstream ifs(filename);
    std::vector<std::string> coordfiles;

    //open file
    if (!ifs) {
		std::cout << "Failed to open json file" << endl;
	}
	bool parsingSuccessful = reader.parse(ifs, root);
	vector<Alignment> result;
	
    //parse the file
    if (!parsingSuccessful) {
		// report to the user the failure and their locations in the document.
		std::cout << "Failed to parse file"
				<< reader.getFormattedErrorMessages() << endl;
		return coordfiles;
	}

	std::string targetName = root["Target"]["Name"].asString();
	std::string targetFullSeq = root["Target"]["FullSeq"].asString();
	std::string predictedSA = root["Target"]["PredictedSA"].asString();
	std::string predictedSS = root["Target"]["PredictedSS"].asString();
	std::string predictedSSConf = root["Target"]["PredictedSSConf"].asString();
	int segLen = root["Target"]["SeqLen"].asInt();
    
    AlVct.resize(root["Aligns"].size());
    coordfiles.resize(root["Aligns"].size());
	for (int i = 0; i < root["Aligns"].size(); i++) 
    {
        Json::Value cur = root["Aligns"][i];
        ////////////////////////////
		//set the target information
        ////////////////////////////
		AlVct[i].setTargetName(targetName);
		AlVct[i].setTargetSequence(targetFullSeq);
		AlVct[i].setPredictedSaInfo(predictedSA);
		AlVct[i].setPredictedSsInfo(predictedSS);
		AlVct[i].setPredictedSsConf(predictedSSConf);


        //get filename for coords file
        coordfiles[i] = cur["Template"]["fileID"].asString();
        
        /////////////////////////
		// alignment information
        ////////////////////////

		//double expect = cur["Alignment"]["Expect"].asDouble();
		//AlVct[i].setExpect(expect);
        Json::Value Algnmt = cur["Alignment"];
		double gdtts = Algnmt["GDTTS"].asDouble();
		AlVct[i].setGdttsScore(gdtts);
		int gaps = Algnmt["Gaps"].asInt();
		//AlVct[i].setGaps(gaps);
		//int identities = cur["Identities"].asInt();
		//AlVct[i].setIdentities(identities);
		//int positives = cur["Positives"].asInt();
		//AlVct[i].setPositives(positives);
		int queryEnd = Algnmt["QueryEnd"].asInt();
		AlVct[i].setQueryEnd(queryEnd);
		std::string queryPart = Algnmt["QueryPart"].asString();
		AlVct[i].setQueryPart(queryPart);
		int queryStart = Algnmt["QueryStart"].asInt();
		AlVct[i].setQueryStart(queryStart);
		int subjectEnd = Algnmt["SubjectEnd"].asInt();
		AlVct[i].setSubjectEnd(subjectEnd);
		std::string subjectPart = Algnmt["SubjectPart"].asString();
		AlVct[i].setSubjectPart(subjectPart);
		int subjectStart = Algnmt["SubjectStart"].asInt();
		AlVct[i].setSubjectStart(subjectStart);
		double tmScore = Algnmt["TMscore"].asDouble();
		AlVct[i].setTmScore(tmScore);
        
        //////////////////////////////
		//set the template information
        /////////////////////////////
        Json::Value tmplt = cur["Template"];
		std::string templateName = tmplt["Name"].asString();
		AlVct[i].setTemplateName(templateName);
		std::string templateRealSeqInfo =
				tmplt["RealSeqInfo"].asString();
		AlVct[i].setTemplateRealSequenceInfo(templateRealSeqInfo);
		std::string refSeqInfo = tmplt["RefSeqInfo"].asString();
		AlVct[i].setTemplateReferenceSequenceInfo(refSeqInfo);
		int segLen = tmplt["SeqLen"].asInt();
		AlVct[i].setTemplateSequenceLength(segLen);
		std::string tempTrueSS = tmplt["TempTrueSS"].asString();
		AlVct[i].setTemplateTrueSecondaryStructure(tempTrueSS);
		//std::string fileID = tmplt["fileID"].asString();
		//AlVct[i].setFileId(fileID);
	}
	return coordfiles;
}

//////////////////////////////////////////////////////////////////////////////
//my3Dinfo part
//////////////////////////////////////////////////////////////////////////////

//read .coords  files 
bool my3Dinfo::readCoord(const char* fileName)
{
    ifstream infile(fileName);
    if(!infile)
    {
        cerr<<"can not open file:"<<fileName<<endl;
        cerr<<strerror(errno)<<endl;
        return false;
    }
    double x,y,z;
    char cm1,cm2;
    int i=0;
    std::vector<Eigen::Vector3d> tmp;
    //read in a random length file
    while(infile>>x>>cm1>>y>>cm2>>z)
    {   
        //put coordinate in matrix
        Eigen::Vector3d curr;
        curr << x , y , z;
        tmp.push_back(curr);
        i++;
    }

    //copy data to matrix
    allCrds.resize(tmp.size(),3);
    for(int i=0; i<tmp.size(); i++)
    {
        allCrds.row(i) = tmp[i];
    }
    return true;

}


/*This function have to deal with four different indexes:
 * 1. index of target : start from 1
 * 2. index of template : start from 1
 * 3. index of alignment : start from 0
 * 4. index of 3d coordinates : start from 0
 * due to the "fully-extended" operation the index in 3 and 4 are different
 * since I'm using c++ I will convert all the index to the form of starting from 0
 */
void my3Dinfo::makeSegment(Alignment& al)
{
    int qst = al.getQueryStart();
    int qed = al.getQueryEnd();
    int ql = al.getTargetLength();
    int tst = al.getSubjectStart();
    int ted = al.getSubjectEnd();
    int tl = al.getTemplateSequenceLength();
    int indx_3d=0;
    segment sg;
    int segLen;
    if(qst <tst)
    {
        segLen = qst - 1;//qst -1 +1
        // change to start from "zero": tst -segLen +1 -1
        tst = tst - segLen ;
        qst = 0;
    }
    else
    {
        segLen = tst - 1;
        qst = qst - segLen;
        tst = 0;
    }
    //get extended preffix
    prefix.crds.resize(segLen,3);
    prefix.crds = allCrds.block(indx_3d,0,segLen,3);
    prefix.qst = qst; sg.qed = qst + segLen -1;
    prefix.tst = tst; sg.ted = tst + segLen -1;

    //get actual coordinates
    qst += segLen;
    tst += segLen;
    indx_3d += segLen; 

    std::string qry_seq = al.getQueryPart();
    std::string tmplt_seq = al.getSubjectPart();
    int alignLen = qry_seq.length();
    int indx = 0;
    while(indx < alignLen-1)
    {
        //find next segment
        int len_tg = blockLength(qry_seq, indx);
        int len_tp = blockLength(tmplt_seq, indx);
        segLen = std::min(len_tg, len_tp);
        sg.crds.resize(segLen,3);
        sg.crds = allCrds.block(indx_3d, 0, segLen, 3);
        sg.qst = qst; sg.qed = qst + segLen-1;
        sg.tst = tst; sg.ted = tst + segLen-1;
        sgmts.push_back(sg);

        indx += segLen;
        indx_3d += segLen;
        qst += segLen; tst += segLen;

        //find next gap and skip
        len_tg = gapLength(qry_seq, indx);
        len_tp = gapLength(tmplt_seq, indx);
        segLen = std::max(len_tg, len_tp);
        
        //no gap at all
        if(segLen == 0)
            break;

        if(len_tg != 0)
        {//gap in target
            indx += segLen;
            tst +=segLen;
        }
        else if(len_tp !=0)
        {// gap in template
            indx +=segLen;
            indx_3d += segLen;
            qst += segLen;
        }
    }


    //get extended suffix
    segLen = std::min((tl - ted),(ql -qed));
    suffix.qst = qed;
    suffix.qed = qed + segLen -1;
    suffix.tst = tst;
    suffix.ted = tst + segLen -1;
    suffix.crds.resize(segLen,3);
    suffix.crds = allCrds.block(indx_3d,0,segLen,3);
}   


void my3Dinfo::segment:: block2array(int start, int end, double**array)
{
    int crdStart = start - qst;
    int crdEnd = end - qst;
    int colLen = crds.cols();

    int i=0;
    for(int rowIndex = crdStart; rowIndex <= crdEnd; rowIndex++,i++)
        for(int colIndex = 0; colIndex < colLen; colIndex++)
            array[i][colIndex] = crds(rowIndex,colIndex);
}


///////////////////////////////////////////////////////////////////////////////
//some utility function that should be include into utility.cpp
//For some reasons I put them here.
/////////////////////////////////////////////////////////////////////////////////

int blockLength(std::string seq, int start)
{
    int len=0;
    while(len+start < seq.length() && seq[start+len] != '-')
        len++;
    return len;
}


int gapLength(std::string seq, int start)
{
    int len=0;
    while(len+start < seq.length() && seq[start+len] =='-')
        len++;
    return len;
}

void Eigen2Array(Eigen::MatrixXd &mat, double** array)
{
    int rl = mat.rows();
    int cl = mat.cols();
    for(int i= 0; i< rl;i++)
        for(int j=0; j<cl;j++)
            array[i][j] = mat(i,j);
}

