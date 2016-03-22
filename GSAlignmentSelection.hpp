/*
 * Author: Ming Du
 * Feb.16.2016
 *
 * The basic idea of this module comes from Dr Yi Shang's algorithm 
 * "GSAQA //Gap-filling Structure-based Alignment Quality Assessment"
 * In order to fit the naming convention of the whole MUFOLD I change it to 
 * GSAlignmentSelection here
 */


#include "AlignmentSelection/AlignmentSelectionBase.h"
#include "myUtility.hpp"
using namespace std;



typedef class GSAlignmentSelection : AlignmentSelectionBase
{
    public:
        GSAlignmentSelection();
        ~GSAlignmentSelection();
        void makeSelection(vector<Alignment>& AlignVector);
    private:
        void GSAS_all();
        void computeScores();
        double getScore();
        void createScoringDM();
        
}GSAS;
