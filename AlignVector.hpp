#ifndef ALIGNVECTOR
#define ALIGNVECTOR

#include "DataObject/DataObjects.h"

#include <stdio.h>
#include "Eigen/Dense"
#include<iostream>
#include "myUtility.hpp"

//define the 3d coordinate data structure
//each 3d structure stored in a n-by-3 matrix
//detail fo the 3D coordinate file:
//
#define DEFAULT_MATRIX_LENGTH 10

////////////////////////////////////////////
//////////////////


/*
 *Input: 1. AA sequence 
 *       2. start index
 *Output: the length from start to the end of the block
 */
int blockLength(std::string seq, int start);

/*
 *Input: 1. AA sequence
 *       2. start index
 *Output: the length from start to the end of the gap
 */
int gapLength(std::string seq, int start);



//Input: Eigen matrix
//Output: copy matrix to arry
//
//
void Eigen2Array(Eigen::MatrixXd &mat, double** arry);
//////////////////////////
////////////////////////////////////////////

class my3Dinfo
{   
    public:
        
        /*
         * input: filename of .coords file
         * output: 1. return ture if succeed false otherwise
         *         2. read file into "info_3D"
         */
        bool readCoord(const char* fileName);

       

        /*Cut 3D information in allCrds into segements corresponding to alignment info
         *Input: 1.Alignment al
         *       2.Matrix allCrds
         *Output: vector of segment sgmts
         */
        void makeSegment(Alignment& al);


        int segNum()
        {   return sgmts.size();    }



       
        Eigen::MatrixXd operator ()(int startRow, int endRow)
        {
        }


        class segment
        {
            public:
            //3d coordinates of the segment
            Eigen::MatrixXd crds;
            
            //This function takes query index start end and return the corresponds part in
            //template 3D structure
            //Input: 1.start 2.end
            //Output: copy data from allCrds to array 
            //Notice that this function do NOT allocate memory for array
            void block2array(int start, int end, double** array);


            /////following variables start from zero
            
            //start position in query
            int qst;
            //end position in query
            int qed;
            //start point in template
            int tst;
            //end point in template
            int ted;
            
            /////end start from zero

            //length of the segment
            inline int len()
            {   return crds.rows();}
        };

    public:
        //store the 3D infomations in the form of segment
        std::vector<segment> sgmts;
        
        //there may not have a prefix or suffix, and then the rows() of the segment will be zero
        //the fully extended part at the begining
        segment prefix;
        //the fully extended part at the end
        segment suffix;
        
        //store the 3D informations in the form of one n-by-3 matrix
        Eigen::MatrixXd allCrds;
};



class AlignVector
{
    public: 
        AlignVector()
        {
        }
        ~AlignVector()
        {
        }
        /*
         * input:   path that the json and .coords file are stored
         * output:  return ture if read succeed otherwise return false
         *          store result in AlVct
         */
        bool readAll(std::string dataPath,std::string fileName);


        /*
         * input: filename of Json file
         * output: 1. read all content in the Json file into AlVct
         *         2. return a list of .coords file name that correspond to
         *            each alignment
         */
        std::vector<std::string> readAligns(const char* fileName);
        
        class Align3D
        {
            public:
                Align3D(Alignment &algn, my3Dinfo &info)
                    : alignment(algn)
                    , info3D(info)
                { }
                Align3D(const Align3D& otherAlign)
                    : alignment(otherAlign.alignment)
                    , info3D(otherAlign.info3D)
                { }
            public:
                Alignment & alignment;
                my3Dinfo & info3D;
        };


        Align3D operator [](int i)
        {
            Align3D ret(AlVct[i], info_3D[i]);
            return ret;
        }

        int len()
        {   return AlVct.size();    }  

    public:
        std::vector<Alignment> AlVct;
        int vctlen;
        std::vector<my3Dinfo> info_3D;         
};

#endif 


