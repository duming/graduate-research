/*
===============================================================================
   Implementation of TM-align in C/C++   

   This program is written by Jianyi Yang at
   Yang Zhang lab
   Center for Computational Medicine and Bioinformatics 
   University of Michigan 
   100 Washtenaw Avenue, Ann Arbor, MI 48109-2218 
                                                       
           
   Please report bugs and questions to yangji@umich.edu or zhng@umich.edu
===============================================================================
*/
#define MAXLEN 10000                        //maximum length of filenames
char version[20];                          //version 
 
 
//global variables

#include "basic_fun.hpp"
#include "Kabsch.h"
#include "TMalign.h"





void TMscore::parameter_set4search(int xlen, int ylen)
{
    this->xlen = xlen;
    this->ylen = ylen;
	//parameter initilization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
	D0_MIN=0.5; 
	dcu0=4.25;                       //update 3.85-->4.25
 
	Lnorm=getmin(xlen, ylen);        //normaliz TMscore by this in searching
    if(Lnorm<=19)                    //update 15-->19
    {
        d0=0.168;                   //update 0.5-->0.168
    }
    else
    {
        d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    }
	D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search    


	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;


    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void TMscore::parameter_set4final(double len)
{
	D0_MIN=0.5; 
 
	Lnorm=len;            //normaliz TMscore by this in searching
    if(Lnorm<=21)         
    {
        d0=0.5;          
    }
    else
    {
        d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    }
    if(d0<D0_MIN) d0=D0_MIN;   

	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  

}


void TMscore::parameter_set4scale(int len, double d_s)
{
 
	d0=d_s;          
	Lnorm=len;            //normaliz TMscore by this in searching

	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  

}




TMscore:: TMscore()
{

	//parameter_set4search(MAXVCTLEN, MAXVCTLEN);          //please set parameters in the function
    //int simplify_step     = 40;               //for similified search engine
    //int score_sum_method  = 8;                //for scoring method, whether only sum over pairs with dis<score_d8
        
	int i;
    int *invmap0          = new int[MAXVCTLEN+1]; 
    int *invmap           = new int[MAXVCTLEN+1]; 
    double TM, TMmax=-1;
	for(i=0; i<MAXVCTLEN; i++)
	{
		invmap0[i]=-1;
	}	

    //////////////////////////////////////
    //  allocate memory
    ////////////////////////////////////
    //------get length first------>
	//MAXVCTLEN=getmin(MAXVCTLEN, MAXVCTLEN);

    //------allocate memory for x and y------>
	NewArray(&xa, MAXVCTLEN, 3);
    seqx   = new char[MAXVCTLEN+1];
	secx   = new int[MAXVCTLEN];
	xresno = new int[MAXVCTLEN];

	NewArray(&ya, MAXVCTLEN, 3);
	seqy    = new char[MAXVCTLEN+1];
	yresno  = new int[MAXVCTLEN];
	secy    = new int[MAXVCTLEN];

    
    
    //------load data------>  
    
    
    //------allocate memory for other temporary varialbes------>
 	NewArray(&r1, MAXVCTLEN, 3);
	NewArray(&r2, MAXVCTLEN, 3);
	NewArray(&xtm, MAXVCTLEN, 3);
	NewArray(&ytm, MAXVCTLEN, 3);
	NewArray(&xt, MAXVCTLEN, 3);

	NewArray(&score, MAXVCTLEN+1, MAXVCTLEN+1);
	NewArray(&path, MAXVCTLEN+1, MAXVCTLEN+1);
	NewArray(&val, MAXVCTLEN+1, MAXVCTLEN+1); 
}

TMscore::~TMscore()
{
    ///////////////
    //  free memory
    ////////////////
	DeleteArray(&path, MAXVCTLEN+1);
	DeleteArray(&val, MAXVCTLEN+1);
	DeleteArray(&score, MAXVCTLEN+1);
    DeleteArray(&xa, MAXVCTLEN);
	DeleteArray(&xt, MAXVCTLEN);
	DeleteArray(&ya, MAXVCTLEN);
	DeleteArray(&r1, MAXVCTLEN);
	DeleteArray(&r2, MAXVCTLEN);
	DeleteArray(&xtm, MAXVCTLEN);
	DeleteArray(&ytm, MAXVCTLEN);
    
    
    delete [] seqx;
    delete [] seqy;
    delete [] secx;
    delete [] secy;
    delete [] xresno;
    delete [] yresno;
 
}


//     1, collect those residues with dis<d;
//     2, calculate TMscore
int TMscore::score_fun8( double **xa, 
                double **ya, 
                int n_ali,
                double d,
                int i_ali[], 
                double *score1,
                int score_sum_method
              )
{
    double score_sum=0, di;
    double d_tmp=d*d;
	double d02=d0*d0;
	double score_d8_cut = score_d8*score_d8;
    
    int i, n_cut, inc=0;

    while(1)
    {
		n_cut=0;
        score_sum=0;
        for(i=0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if(di<d_tmp)
            {
                i_ali[n_cut]=i;
                n_cut++;
            }
            if(score_sum_method==8)
            {				
                if(di<=score_d8_cut)
                {					
                    score_sum += 1/(1+di/d02);
                }                
            }
            else
            {
				score_sum += 1/(1+di/d02);
            }
        }
        //there are not enough feasible pairs, reliefe the threshold 		
        if(n_cut<3 && n_ali>3)
        {
			inc++;
			double dinc=(d+inc*0.5);
			d_tmp = dinc * dinc;
        }
        else
        {
            break;
        }

    }  

    *score1=score_sum/Lnorm;

    return n_cut;
}

// TMscore search engine
// input:   two aligned vector sets: x, y
//          scale parameter d0
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8
//                                  
//          
// output:  the best rotaion matrix t0, u0 that results in highest TMscore
//double TMscore::TMscore8_search( double **xtm, 
//                        double **ytm,
//                        double t0[3],
//                        double u0[3][3],
//                        int simplify_step,
//                        int score_sum_method,
//                        double *Rcomm
//                       )
double TMscore::TMscore8_search() 
{   
    int i, m;
    double score_max, score, rmsd;    
    const int kmax=Lali;    
    int k_ali[kmax], ka, k;
    double t[3];
    double u[3][3];
	double d;
	

	//iterative parameters
	int n_it=20;            //maximum number of iterations
    const int n_init_max=6; //maximum number of different fragment length 
    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   
    int L_ini_min=4;
    if(Lali<4) L_ini_min=Lali;   
    int n_init=0, i_init;      
    for(i=0; i<n_init_max-1; i++)
    {
        n_init++;
        L_ini[i]=(int) (Lali/pow(2.0, (double) i));
        if(L_ini[i]<=L_ini_min)
        {
            L_ini[i]=L_ini_min;
            break;
        }
    }
    if(i==n_init_max-1)
    {
        n_init++;
        L_ini[i]=L_ini_min;
    }
    
    score_max=-1;
    //find the maximum score starting from local structures superposition
    int i_ali[kmax], n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment
    for(i_init=0; i_init<n_init; i_init++)
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag;
      
        i=0;   
        while(1)
        {
            //extract the fragment starting from position i 
            ka=0;
            for(k=0; k<L_frag; k++)
            {
				int kk=k+i;
                r1[k][0]=xtm[kk][0];  
                r1[k][1]=xtm[kk][1]; 
                r1[k][2]=xtm[kk][2];   
                
                r2[k][0]=ytm[kk][0];  
                r2[k][1]=ytm[kk][1]; 
                r2[k][2]=ytm[kk][2];
                
                k_ali[ka]=kk;
                ka++;
            }
            
            //extract rotation matrix based on the fragment
            Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
            if(i_init==0)
            {
                Rcomm=sqrt(rmsd/Lali);
            }
            do_rotation(xtm, xt, Lali, t, u);
            
            //get subsegment of this fragment
            d=d0_search-1;
            n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, score_sum_method);
            if(score>score_max)
            {
                score_max=score;
                
                //save the rotation matrix
                for(k=0; k<3; k++)
                {
                    t0[k]=t[k];
                    u0[k][0]=u[k][0];
                    u0[k][1]=u[k][1];
                    u0[k][2]=u[k][2];
                }
            }
            
            //try to extend the alignment iteratively            
            d=d0_search+1;
            for(int it=0; it<n_it; it++)            
            {
                ka=0;
                for(k=0; k<n_cut; k++)
                {
                    m=i_ali[k];
                    r1[k][0]=xtm[m][0];  
                    r1[k][1]=xtm[m][1]; 
                    r1[k][2]=xtm[m][2];
                    
                    r2[k][0]=ytm[m][0];  
                    r2[k][1]=ytm[m][1]; 
                    r2[k][2]=ytm[m][2];
                    
                    k_ali[ka]=m;
                    ka++;
                } 
                //extract rotation matrix based on the fragment                
                Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
                do_rotation(xtm, xt, Lali, t, u);
                n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, score_sum_method);
                if(score>score_max)
                {
                    score_max=score;

                    //save the rotation matrix
                    for(k=0; k<3; k++)
                    {
                        t0[k]=t[k];
                        u0[k][0]=u[k][0];
                        u0[k][1]=u[k][1];
                        u0[k][2]=u[k][2];
                    }                     
                }
                
                //check if it converges                 
			
                if(n_cut==ka)
                {				
                    for(k=0; k<n_cut; k++)
                    {
                        if(i_ali[k]!=k_ali[k])
						{
							break;
						}
                    }
                    if(k==n_cut)
                    {						
                        break; //stop iteration
                    }
                }                                                               
            } //for iteration            

			if(i<iL_max)
			{
				i=i+simplify_step; //shift the fragment		
				if(i>iL_max) i=iL_max;  //do this to use the last missed fragment
			}
			else if(i>=iL_max)
			{
				break;
			}
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}


