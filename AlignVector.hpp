#include "DataObject/DataObject.h"

#include <stdio.h>
#include "Point.h"
#include "Eigen/Dense"
#include<iostream>

//define the 3d coordinate data structure
typedef Eigen::Matrix<double,Dynamic,3> COORDS;

class AlignResult
{
    public: 
        AlignResult(std::vector<Alignment>&AlVct)
        {
            this->AlVct = AlVct;
            vctlen = AlVct.size();
        }
        ~AlignResult();
       
        /*
         * input:   path that the json and .coords file are stored
         * output:  return ture if read succeed otherwise return false
         *          store result in AlVct
         */
        bool readAll(std::string dataPath,std::string filename);


        /*
         * input: filename of Json file
         * output: 1. read all content in the Json file into AlVct
         *         2. return a list of .coords file name that correspond to
         *            each alignment
         */
        std::vector<std::string> readAligns(const char* filename);
        

        /*
         * input: filename of .coords file
         * output: 1. return ture if succeed false otherwise
         *         2. read file into "   "
         */
        bool readCoord(const char* filename);

    private:
        std::vector<Alignment> &AlVct;
        int vctlen;
        
};