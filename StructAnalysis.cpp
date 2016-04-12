#include "StructAnalysis.hpp"
        
#include <iostream>
using namespace std;
void StructAnalysis::test()
{
    TMscore tms;
    int alignLen = AlVct.len();
    CompareResult::Vote vote;
    CompareResult::IndexKey key;

    int i = 0;
    int j = 1;

    //for(int i=0; i < alignLen ; i++)
    {
        int sgNum = AlVct[i].info3D.segNum();
        
        std::vector<my3Dinfo::segment> &segs = AlVct[i].info3D.sgmts;
      //  for(int j=0; j < sgNum ; j++)
        {
            std::vector<Interval<CompareResult::IndexKey>> intervals;
            
            int startSource = segs[j].qst;
            int endSource = segs[j].qed;


            cout<<'{'<<startSource<<','<<endSource<<'}'<<endl;

            searchTree->findContained(
                                      startSource
                                    , endSource
                                    , intervals);


            std::cout<<"("<<i<<","<<j<<")"<<std::endl;
            //process each overlap interval 
            for(int k = 0 ; k < intervals.size() ; k++)
            {
                cout<<'['<<intervals[k].value.algnIndex
                    <<','<<intervals[k].value.sgIndex<<']'<<endl;

                if(intervals[k].value.algnIndex == 1)
                    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
            }
        }
    } 
}


ostream& operator <<(ostream& os, CompareResult::IndexKey index)
{
    os<<"["<<index.algnIndex<<","<<index.sgIndex<<"]";
    return os;
}



void StructAnalysis:: analyze(int option )
{
    TMscore tms;
    int alignLen = AlVct.len();
    CompareResult::Vote vote;
    CompareResult::IndexKey key;

    for(int i=0; i < alignLen ; i++)
    {
        int sgNum = AlVct[i].info3D.segNum();
        if(option == 1)
        {//add prefix and suffix
            
            //TODO
        }

        std::vector<my3Dinfo::segment> &segs = AlVct[i].info3D.sgmts;
        for(int j=0; j < sgNum ; j++)
        {
            std::vector<Interval<CompareResult::IndexKey>> intervals;
            
            int startSource = segs[j].qst;
            int endSource = segs[j].qed;

            searchTree->findOverlapping(
                                      startSource
                                    , endSource
                                    , intervals);


            std::cout<<"("<<i<<","<<j<<")"<<std::endl;
            //process each overlap interval 
            for(int k = 0 ; k < intervals.size() ; k++)
            {
                int algnIndex = intervals[k].value.algnIndex;
                int sgIndex = intervals[k].value.sgIndex;
                int startResult = intervals[k].start;
                int endResult = intervals[k].stop;
                //if is the same segment then skip
                if(i == algnIndex && j == sgIndex)
                    continue;

                ///////////////////////////////////////
                //calculate all infomation of the overlap
                /////////////////////////////////////////
                
                               

                int startOverlap = overlapStart(startSource, startResult);
                int endOverlap = overlapEnd(endSource, endResult);
                int overlapLen = endOverlap - startOverlap +1;

                //initialize te tmscore engine
                tms.setLength(overlapLen);
                tms.setStep(1,8);
                //input source part
                segs[j].block2array(
                                     startOverlap
                                    ,endOverlap
                                    ,tms.xa);
                
                //input result part
                AlVct[algnIndex].info3D.sgmts[sgIndex].block2array(
                                      startOverlap
                                    , endOverlap
                                    , tms.ya);




                                
                //save the result
                vote.algnIndex = algnIndex;
                vote.sgIndex = sgIndex;
                vote.segLen = AlVct[algnIndex].info3D.sgmts[sgIndex].crds.rows();
                vote.overlapLen = overlapLen;
                key.algnIndex = algnIndex;
                key.sgIndex = sgIndex;
                ///check if the distance is already calculate
                if(algnIndex < i)
                {//copy the distance from result[algnIndex][segIndex] to result[i][j] 
                    CompareResult::IndexKey searchkey;
                    CompareResult::Vote voteRet;
                    map<CompareResult::IndexKey,CompareResult:: Vote> :: iterator mit;
                    searchkey.algnIndex = i;
                    searchkey.sgIndex = j;
                    mit = algnBase[algnIndex].voteVct[sgIndex].find(searchkey);

                   
                    if(mit == algnBase[algnIndex].voteVct[sgIndex].end())
                        cout<<"not found key"<<endl;

                    vote.dist = mit->second.dist; 
                                     
                    algnBase[i].voteVct[j][key] = vote;
                    std::cout<<'['<<j<<","<<k<<":"<<vote.dist<<']';
                    continue;
                }

                //calculate the score
                double dist = tms.TMscore8_search();
                std::cout<<'['<<j<<","<<k<<":"<<dist<<']';
                vote.dist = dist;
                //insert result into map 
                algnBase[i].voteVct[j][key] = vote;

                    
                /////////////////////////////////
                //end overlap calculation
                /////////////////////////////////
                
            }
            std::cout<<std::endl<<std::endl;
        }//end seg
    }//end align
}



void StructAnalysis:: makeTree()
{
    CompareResult::IndexKey tmpKey;
    std::vector<Interval<CompareResult::IndexKey>> intervals;
    int alignLen = AlVct.len();
    
    //select all segments
    for(int i=0; i < alignLen; i++)
    {//for each alignment
        int sgNum = AlVct[i].info3D.segNum();
        std::vector<my3Dinfo::segment> &segs = AlVct[i].info3D.sgmts;

        for(int j=0 ; j < sgNum ; j++)
        {//for each segment
            //set index as data
            tmpKey.algnIndex = i;
            tmpKey.sgIndex = j;
            intervals.push_back(Interval<CompareResult::IndexKey>(
                                             segs[j].qst
                                            ,segs[j].qed
                                            ,tmpKey));
        }
    }

    //make the three
    searchTree = new IntervalTree<CompareResult::IndexKey>(intervals);
}
