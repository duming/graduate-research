#include "GSAlignmentSelection.hpp"
#include "myUtility.hpp"
#include "DataObject/Point.h"
#include "stdlib.h"
#include "DataObject/Alignment.h"
#include "AlignVector.hpp"
#include "TMalign.h"
#include "basic_fun.hpp"

#include <vector>
#include <string>
#include <iostream>
using namespace std;



int main()
{

    string dataPath ="/Users/ming/projects/MUfold/Data/fakeTestData/";
    //GenerateExperimentDatasetOutput/casp11_04302014_fully/"; 
    string fileName = "T0759_blast.json";
    AlignVector alv;
    
    alv.readAll(dataPath,fileName);
    
    
    double **st1;
    double **st2;
    double t[3],u[3][3],Rcomm;
    int sglen = alv.info_3D[0].sgmts[1].len();
    NewArray(&st1,sglen,3);
    NewArray(&st2,sglen,3);
    

    double score; 
    TMscore tms(sglen,sglen);
    Eigen2Array(alv.info_3D[0].sgmts[1].crds, st1);
    Eigen2Array(alv.info_3D[0].sgmts[1].crds, st2);
    score = tms.TMscore8_search(st1,st2,sglen,t,u,1,8,&Rcomm);
    cout<<score<<endl;

    /*
    //alv.info_3D[0].makeSegment(alv.AlVct[0]);
    //test make segment
    for(int i=0;i<alv.vctlen;i++)
    {
        //cout<<alv.info_3D[i].sgmts.size()<<endl;
        //cout<<alv.AlVct[i].getQueryPart()<<endl;
        //cout<<alv.AlVct[i].getSubjectPart()<<endl;
        //cout<<alv.info_3D[i].allCrds<<endl;
        int totalLen = alv.info_3D[i].allCrds.rows();
        //cout<< alv.info_3D[i].sgmts.back().crds.rows()<<endl;
        int sum = 0;
        for(int j=0;j<alv.info_3D[i].sgmts.size();j++)
        {   
            my3Dinfo::segment& sg = alv.info_3D[i].sgmts[j];
            //check total length of all segments
            sum += sg.crds.rows();

            int sgl1, sgl2;
            sgl1 = sg.ted - sg.tst+1;
            sgl2 = sg.qed -sg.qst+1;
            if(sgl1 != sg.crds.rows() || sgl2 != sgl1)
                cout<<j<<":segment length error"<<sgl1<<"/"<<sgl2<<"/"<<sg.crds.rows()<<endl;

            //cout<<sg.qst<<"/"<<sg.qed<<"/"<<sg.tst<<"/"<<sg.ted<<endl;
            //cout<<&sg.crds<<endl;
            cout<<sg.crds.rows()<<endl;
            //cout<<sg.crds<<endl;
            //cout<<endl;
        }
        int gapNum = countGap(alv.AlVct[i].getSubjectPart());
        if(!(sum + gapNum == totalLen))
            cout<<"total length not equal"<<sum<<"/"<<totalLen<<endl;

        cout<<"++++++++++++++++"<<i<<"++++++++++++++++++"<<endl;
    }
    cout<<alv.vctlen<<endl;
    */
    /*
    //test gapLength()
    string testStr = "IVDPVSNLRLPVEEAYKRGLVGIEFKEKLLSAE------------RAVTGYNDPETGNIISLFQAMNKELIEKGH";
    for(int i=0; i<testStr.length(); i++)
    {
        cout<<gapLength(testStr,i)<<" "<<testStr[i]<<endl;
    }
    cout<<endl<<testStr<<endl;
    */

    /*
    //test file format 
    for(int i;i<alv.vctlen;i++)
    {
        string qstr = alv.AlVct[i].getQueryPart();
        string sstr = alv.AlVct[i].getSubjectPart();

       // cout<<qstr<<endl;
       // cout<<sstr<<endl;

        int tgl = alv.AlVct[i].getTargetLength();
        int tpl = alv.AlVct[i].getTemplateSequenceLength();
        
        int tgst = alv.AlVct[i].getQueryStart();
        int tpst = alv.AlVct[i].getSubjectStart();

        //int shift = alv.info_3D[i].shift;
        int mtl = alv.info_3D[i].allCrds.rows();

        int tge = alv.AlVct[i].getQueryEnd();
        int tpe = alv.AlVct[i].getSubjectEnd();


        //cout<<tgl<<"/"<<tpl<<endl;
        //cout<<tgst<<"/"<<tpst<<endl;
        //cout<<tge<<"/"<<tpe<<endl;
        //cout<<mtl<<"/"<<std::min(tgst,tpst) + (tge-tgst) +
        //cout<<(tge-tgst+1)<<"/"<<(tpe -tpst+1)<<"/"<<qstr.length()<<endl;
        int pred = (tge-tgst+1)+ std::min(tgst-1,tpst-1) + std::min((tgl -tge),(tpl-tpe));
        if(mtl!= pred)
        {
            cout<<qstr<<endl<<sstr<<endl;
            cout<<pred<<"/"<<mtl<<endl;
        }
    }*/
    

    /*
    //test info_3D operator overload
    cout<<"###################################"<<endl;
    cout<<alv.info_3D[0](0,0);
    
    cout<<"###################################"<<endl;
    int alignEnd = alv.AlVct[0].getSubjectEnd();
    cout<<alv.info_3D[0](alignEnd,alignEnd);
    */


    return 1;
}


/*
polymorphism test 
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
