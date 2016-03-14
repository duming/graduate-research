#include "GSAlignmentSelection.h"
#include "myUtility.h"
#include "DataObject/Point.h"
#include "stdlib.h"
#include "DataObject/Alignment.h"
#include "AlignVector.hpp"


#include <vector>
#include <string>
#include <iostream>
using namespace std;



int main()
{

    vector<Alignment> AlVct;
    string dataPath ="/Users/ming/projects/MUfold/Data/GenerateExperimentDatasetOutput/casp11_fully/"; 
    string dataName ="T0759";
    string suffix = "_blast.json";

    string n2 ="test.json";    
    string coordName="T0759/blast_T0759_0_1LM7_A.coords";
    
    char* filename = "";
    vector<string> fileIDs;

    AlignVector Algns;

    fileIDs = Algns.readAligns(filename);
    



    return 1;
}


/*
multi
class base
{
    public:
        base(const char* a)
        {
            setA(a);
        }
        //~base();
        void setA(const char* a)
        {
            A=a;
        }
        string getA()
        {
            return A;
        }
        virtual string getB()
        {
            string s("aaa");
            cout<<s<<endl;
            return s;
        }

    private:
        string A;
};

class derive: public base
{
    public:
        derive(const char* b):base("A form derived")
        {
            setB(b);
        }
        void setB(const char* b)
        {
            B=b;
        }

        string getB()
        {
            string s("ssss");
            cout<<s<<endl;
            return s;
        }
    private:
        string B;
};

*/

/*
void test2()
{
    base* bptr;
    base bs("base");
    derive dv("dddd");

    cout<<bs.getA()<<endl;
    dv.getB();
    
    bs = dv;

    cout<<bs.getA()<<endl;
    bs.getB();
    bptr=&dv;
    bptr->getB();
}
*/
