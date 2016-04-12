#ifndef STRUCTANALYSIS
#define STRUCTANALYSIS

#include <map>
#include <vector>
#include "AlignVector.hpp"
#include "IntervalTree.h"
#include "TMalign.h"

//Input: starting index of two segments
//Output: starting index of the overlaping part
inline int overlapStart(int start1, int start2)
{
    return std::max(start1,start2);
}

//Input: ending index of two segments
//Output: ending index of the overlaping part
inline int overlapEnd(int end1, int end2)
{
    return std::min(end1,end2);
}



//data unit that records all overlaped segments and their similarities for each 
//segment in one alignment
class CompareResult
{
    friend class StructAnalysis;
    //keep the search result of segment that overlaps the query segment
    struct Vote
    {
        public:
        //index of the alignment in AlignmentVector
        int algnIndex;
        //index of segment in sgmts
        int sgIndex;

        int segLen;
        // length of the vote segment that overlaps the current segment
        int overlapLen;

        //the distance between the overlaping part of the two segments
        double dist;
    };
    
    public:
    struct IndexKey
    {
        int algnIndex;
        int sgIndex;
    };

    struct CompareKey
    {//use the alignment index as primary key
        bool operator () (const IndexKey &key1, const IndexKey &key2) const
        {    if(key1.algnIndex == key2.algnIndex)
                return key1.sgIndex < key2.sgIndex;
            else
                return key1.algnIndex < key2.algnIndex;
        }
    };

    static int overlapLen(int start1, int end1, int start2, int end2)
    {
        int ret = std::min(end1, end2) - std::max(start1, start2) + 1;
        return ret>0?ret:0;
    }

public:
    
public:
    std::vector<std::map<IndexKey,Vote,CompareKey>> voteVct;    
    double tureGDT;

};


class StructAnalysis
{
    public:
        StructAnalysis(AlignVector& InputAlVct): AlVct(InputAlVct)
        {
            int alignLen = InputAlVct.len();
            algnBase = NULL;
            searchTree = NULL;
            if(alignLen ==0)
            {
                cout<<"empty Alignment vector: "<<endl;
            }
            //allocate memory for array
            algnBase = new CompareResult[alignLen];
            for(int i = 0 ; i < alignLen ;i++)
                algnBase[i].voteVct.resize(AlVct[i].info3D.segNum());


            //construct search tree
            makeTree();
        }
        ~StructAnalysis()
        {
            delete [] algnBase;
        }


        //Input: 1. AlVct
        //       2. option  0: don't include extended part
        //                  !0: include extended part
        //Output: algnBase: similarity search results
        void analyze(int option = 0);

        
        void test();

    private:
        //Input Alignment vector
        //Output: searchTree an interval-tree 
        void makeTree();


    private:
        //the alignments to be analysized
        AlignVector& AlVct;
        //the result of the analysis
        CompareResult* algnBase;
        //the auxiliary search tree 
        IntervalTree<CompareResult::IndexKey>* searchTree;
};




#endif

