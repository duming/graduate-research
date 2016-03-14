#include "AlignVector.hpp"

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
	int seqLen = root["Target"]["SeqLen"].asInt();
    
    AlVct.resize(root["Align"].size());
    coordfiles.resize(root["Align"].size());
	for (int i = 0; i < root["Align"].size(); i++) 
    {
		//here is to set the target information
		AlVct[i].setTargetName(targetName);
		AlVct[i].setTargetSequence(targetFullSeq);
		AlVct[i].setPredictedSaInfo(predictedSA);
		AlVct[i].setPredictedSsInfo(predictedSS);
		AlVct[i].setPredictedSsConf(predictedSSConf);


        //get filename for coords file
        coordfiles[i] = root["Align"][i]["fileID"].asString();

		//following is all the alignment information
		//double expect = root["Align"][i]["Alignment"]["Expect"].asDouble();
		//AlVct[i].setExpect(expect);
		double gdtts = root["Align"][i]["Alignment"]["GDTTS"].asDouble();
		AlVct[i].setGdttsScore(gdtts);
		int gaps = root["Align"][i]["Alignment"]["Gaps"].asInt();
		//AlVct[i].setGaps(gaps);
		//int identities = root["Align"][i]["Identities"].asInt();
		//AlVct[i].setIdentities(identities);
		//int positives = root["Align"][i]["Positives"].asInt();
		//AlVct[i].setPositives(positives);
		int queryEnd = root["Align"][i]["QueryEnd"].asInt();
		AlVct[i].setQueryEnd(queryEnd);
		std::string queryPart = root["Align"][i]["QueryPart"].asString();
		AlVct[i].setQueryPart(queryPart);
		int queryStart = root["Align"][i]["QueryStart"].asInt();
		AlVct[i].setQueryStart(queryStart);
		int subjectEnd = root["Align"][i]["SubjectEnd"].asInt();
		AlVct[i].setSubjectEnd(subjectEnd);
		std::string subjectPart = root["Align"][i]["SubjectPart"].asString();
		AlVct[i].setSubjectPart(subjectPart);
		int subjectStart = root["Align"][i]["SubjectStart"].asInt();
		AlVct[i].setSubjectStart(subjectStart);
		double tmScore = root["Align"][i]["TMscore"].asDouble();
		AlVct[i].setTmScore(tmScore);

		//following is to set the template information
		std::string templateName = root["Template"][i]["Name"].asString();
		AlVct[i].setTemplateName(templateName);
		std::string templateRealSeqInfo =
				root["Template"][i]["RealSeqInfo"].asString();
		AlVct[i].setTemplateRealSequenceInfo(templateRealSeqInfo);
		std::string refSeqInfo = root["Template"][i]["RefSeqInfo"].asString();
		AlVct[i].setTemplateReferenceSequenceInfo(refSeqInfo);
		int seqLen = root["Template"][i]["SeqLen"].asInt();
		AlVct[i].setTemplateSequenceLength(seqLen);
		std::string tempTrueSS = root["Template"][i]["TempTrueSS"].asString();
		AlVct[i].setTemplateTrueSecondaryStructure(tempTrueSS);
		std::string fileID = root["Template"][i]["fileID"].asString();
		AlVct[i].setFileId(fileID);
	}
	return coordfiles;
}
