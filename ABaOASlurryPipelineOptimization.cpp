/*
***	 BI-ATTEMPTED BASE OPTIMIZATION ALGORITHM (ABaOA)
***	 DEVELOPED BY: MEHTAP KOSE ULUKOK (https://orcid.org/0000-0003-4335-483X) AND BURHAN YILDIZ (https://orcid.org/0000-0002-0144-3562) AS A PART OF THE RESEARCH IN (doi: doi: 10.1007/s11269-023-03517-w)

THIS C PROGRAM SOLVES SLURRY PIPELINE OPTIMIZATION PROBLEM AS DEFINED IN YILDIZ ET AL (2014) (See Ref. in above article)
USING "Bi-Attempted Base Optimization Algorithm" 
The solution is saved in 4 different files.
	BestIndividuals.txt
	PenaltyFitness.txt
	Population.txt
	BestFitness.txt
*/
///////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define maxIteration 50
#define N 10000 //number of population 1000
#define m 3//problem dimension
#define n 3//problem dimension

///initialization of coefficients
#define Ss 4.74 
#define Pw 1000 
#define d 0.000045
#define CE 0.1 //electricity unit price $/Kwh
#define PI 3.1416 
#define Kur 2.00 //exchange rate
#define SigmaD 0.01
#define SigmaCV 0.01
#define SigmaD2 0.05
#define SigmaCV2 0.005
#define PenaltyCoef 1000

double population [N][2][m][n]; 
double populationAdd [2][m][n];
double populationSub [2][m][n];
double bestIndividual [2][m][n];

double W[m][n];
double c[n];
double p[m];
double pcap[m]={634.2,317.1,158.55};
double ccap[n]={317.1,317.1,317.1};
double Vl=0;

double popFitness [N];
double popFitnessAdd;
double popFitnessSub;
double bestFitness=99999999999.0;

//pipe lenght
double L [m][n]={{400000,589000,105000},{583000,901000,988000},{180000,501000,585000}}; 
double H [m][n];
double Q [m][n];
double Ro [m][n];
double CV [m][n];
double FCw [m][n];
double Vm [m][n];

//output files
FILE * fpFitness;
FILE * fpIndividual;
FILE * fpBestIndividual;
FILE * fpBest;
FILE *fp;
FILE *fpPenalty;

void initialization() 
{
       int indv, i,j;
       double y;
       srand(time(NULL));
       for (indv=0; indv<N;indv++)
       {
           for(i=0;i<m;i++)
           {
               for(j=0;j<n;j++)
               {      
                      population[indv][0][i][j] = 0+((double) (rand()%10))/10;  
                      y = 0.3+((double) (rand()%400))/1000; 
                      if (y>0.6998)
                      {                          
                          y=0.7;
                      }
                      population[indv][1][i][j] = y;
               }
           }     
       }
      
}
void initializationBase(int k)
{
  int i,j;
  //initialize base matrices           
           for(i=0;i<m;i++)
           {
               for(j=0;j<n;j++)
               {      
                     populationAdd[0][i][j]= population[k][0][i][j];  
                     populationAdd[1][i][j]= population[k][1][i][j];
                     populationSub[0][i][j]= population[k][0][i][j];  
                     populationSub[1][i][j]= population[k][1][i][j];
                    
               }
           }                 
 }

void display(int indv) 
{
       int i,j;
      fprintf(fp,"\n Individual %d D Matrix\n",indv);
       for(i=0;i<m;i++)
       {
         for(j=0;j<n;j++)
         {      
            fprintf(fp,"%f ",population[indv][0][i][j]);
          }  
          fprintf(fp,"\n");  
       } 
       fprintf(fp,"\n Individual %d CW Matrix\n",indv); 
       for(i=0;i<m;i++)
       {
         for(j=0;j<n;j++)
         {      
            fprintf(fp,"%f ",population[indv][1][i][j]); 
          } 
          fprintf(fp,"\n");   
       }  
}
void showpopulation()
{
  int mku;
  printf("\nThe Population Size is:%d\n",N);
  printf("Enter individual index to display:");
  scanf("%d",&mku);
  while((mku>=0)&(mku<N))
  {
     printf("\ndisplay is calling!!!\n");
     display(mku); 
     printf("\ndisplay is called!!!\n");
     printf("Enter individual index [Enter -1 to Stop]:");
     scanf("%d",&mku);
  }
  getchar();
  printf("Finished...");
}

void showpopulationBase()
{
  int mku=0;
  
  while(mku<N)
  {
     display(mku); 
     mku=mku+1;
  }
 
}

///////////////////////////////////////////////
void fitnessPop(int indv)
{
   int i,j;
   double C1,C2;
   double tot_pcap=0;
   double tot_ccap=0;
  
  ///Calculation of Fitness C2  
            
          C2=0;
           for(i=0;i<m;i++)
           {
               for(j=0;j<n;j++)
               {      
                   C2=C2+(210.89*pow(population[indv][0][i][j],1.3744))*L[i][j];   
                   //display each C2                  
               }
           }
           C2=C2/1000000;
            popFitness [indv]=C2;               
       
   ///Calculation of Fitness C1
      
           C1=0;
          for(i=0;i<m;i++)
           {
           for(j=0;j<n;j++)
            {      
              //Calculate Roij

              CV[i][j]=(population[indv][1][i][j])/(population[indv][1][i][j]+4.74*(1-population[indv][1][i][j]));  
              
              Ro[i][j]=Pw*(CV[i][j]*Ss+(1-CV[i][j])); 

              if(population[indv][1][i][j]<0.3)
                {
                    FCw[i][j]=1.097;         
                }
              else if(population[indv][1][i][j]<0.45)
                {
                    FCw[i][j]=0.2067*population[indv][1][i][j]+1.035;         
                }
              else if(population[indv][1][i][j]<0.55)
                {
                    FCw[i][j]=1.52*population[indv][1][i][j]+0.444;         
                }
              else if(population[indv][1][i][j]<0.7)
                {
                    FCw[i][j]=6.1*population[indv][1][i][j]-2.075;         
                }

            Vm[i][j]= 2966.45*(FCw[i][j]*pow(d,0.75)*pow(4.74,0.5)*pow(population[indv][0][i][j],0.5));
            Q[i][j]=(PI*(pow(population[indv][0][i][j],2)/4))*Vm[i][j];   
            H[i][j]=0.0039* pow(CV[i][j],0.803)*pow(population[indv][0][i][j],-1.25)*pow(Vm[i][j],1.77)*L[i][j];
            C1=(C1+(CE*7884*((Ro[i][j]*Q[i][j]*H[i][j])/101.94)));
          } 
        }         
            
           // getchar();
            C1=C1/1000000;
            popFitness [indv]=popFitness[indv]+C1; 
                    
   //// Constarints:
      
          for(i=0;i<m;i++)
           {
           p[i]=0;
           for(j=0;j<n;j++)
              {                  
                 W[i][j]=CV[i][j]*Pw*Ss*(PI/4)*pow(population[indv][0][i][j],2)*Vm[i][j];                                
                 p[i]=p[i]+W[i][j];                
              }                    
            }
                
          for(j=0;j<n;j++)
           {
           c[j]=0;
           for(i=0;i<m;i++)
              {                  
                 W[i][j]=CV[i][j]*Pw*Ss*(PI/4)*pow(population[indv][0][i][j],2)*Vm[i][j];                                
                 c[j]=c[j]+W[i][j];                
              }                    
            }
          
        
        for(i=0;i<m;i++)
        {                  
         tot_pcap=tot_pcap+pcap[i];                 
         }   
         
         for(j=0;j<n;j++)
        {                  
         tot_ccap+=ccap[j];                 
         }       
     //Constraint Vl for each individual 
  
        Vl=0;     
        if(tot_pcap>tot_ccap)
        {
          for(i=0;i<m;i++)
           {
           if( p[i]<0 )                 
             {
                Vl=Vl+PenaltyCoef*pow(p[i],2);              
             }
            if( p[i]>pcap[i])                 
             {
                Vl=Vl+PenaltyCoef*pow((p[i]-pcap[i]),2);              
             }  
           }
           
           for(j=0;j<n;j++)
           {
           if( c[j]<0.99*ccap[j] )                 
             {
                Vl=Vl+PenaltyCoef*pow((0.99*ccap[j]-c[j]),2);              
             }
            if( c[j]>ccap[j])                 
             {
                Vl=Vl+PenaltyCoef*pow((c[j]-ccap[j]),2);              
             }  
           }//end of for j
        }  //end of if          
        else
        {
          for(j=0;j<n;j++)
           {
           if( c[j]<0 )                 
             {
                Vl=Vl+PenaltyCoef*pow(c[j],2);              
             }
            if( c[j]>ccap[j])                 
             {
                Vl=Vl+PenaltyCoef*pow((c[j]-ccap[j]),2);              
             }  
           }
           
           for(i=0;i<m;i++)
           {
           if( p[i]<0.99*pcap[i] )                 
             {
                Vl=Vl+PenaltyCoef*pow((0.99*pcap[i]-p[i]),2);              
             }
            if( p[i]>pcap[i])                 
             {
                Vl=Vl+PenaltyCoef*pow((p[i]-pcap[i]),2);              
             }  
           }
           
        }//end of else

   
       popFitness[indv]=popFitness[indv]+Vl; 
      fprintf(fpPenalty,"\nSolution %d Penalties: %f, fitness %f\n", indv, Vl,popFitness[indv]);
           
  
}//end of function

void fitnessPopAdd()
{
  int i,j;
   double C1,C2;
   double tot_pcap=0;
   double tot_ccap=0;
   
  ///Calculation of Fitness C2         
          C2=0;
           for(i=0;i<m;i++)
           {
               for(j=0;j<n;j++)
               {      
                   C2=C2+(210.89*pow(populationAdd[0][i][j],1.3744))*L[i][j];   
                   //display each C2                  
               }
           }
           C2=C2/1000000;
            popFitnessAdd=C2;               
      
   ///Calculation of Fitness C1
      
           C1=0;
          for(i=0;i<m;i++)
           {
           for(j=0;j<n;j++)
            {      
              //Calculate Roij

              CV[i][j]=(populationAdd[1][i][j])/(populationAdd[1][i][j]+4.74*(1-populationAdd[1][i][j]));  
              
              Ro[i][j]=Pw*(CV[i][j]*Ss+(1-CV[i][j])); 

              if(populationAdd[1][i][j]<0.3)
                {
                    FCw[i][j]=1.097;         
                }
              else if(populationAdd[1][i][j]<0.45)
                {
                    FCw[i][j]=0.2067*populationAdd[1][i][j]+1.035;         
                }
              else if(populationAdd[1][i][j]<0.55)
                {
                    FCw[i][j]=1.52*populationAdd[1][i][j]+0.444;         
                }
              else if(populationAdd[1][i][j]<0.7)
                {
                    FCw[i][j]=6.1*populationAdd[1][i][j]+2.075;         
                }
            Vm[i][j]= 2966.45*(FCw[i][j]*pow(d,0.75)*pow(4.74,0.5)*pow(populationAdd[0][i][j],0.5));
            Q[i][j]=(PI*(pow(populationAdd[0][i][j],2)/4))*Vm[i][j];   
            H[i][j]=0.0039* pow(CV[i][j],0.803)*pow(populationAdd[0][i][j],-1.25)*pow(Vm[i][j],1.77)*L[i][j];
            C1=(C1+(CE*7884*((Ro[i][j]*Q[i][j]*H[i][j])/101.94)));
          } 

            C1=C1/1000000;
            popFitnessAdd =popFitnessAdd+C1; 
     }// end of indv
            
   //// Constarints:
     
          for(i=0;i<m;i++)
           {
           p[i]=0;
           for(j=0;j<n;j++)
              {                  
                 W[i][j]=CV[i][j]*Pw*Ss*(PI/4)*pow(populationAdd[0][i][j],2)*Vm[i][j];                                
                 p[i]=p[i]+W[i][j];                
              }                    
            }
          
          for(j=0;j<n;j++)
           {
           c[j]=0;
           for(i=0;i<m;i++)
              {                  
                 W[i][j]=CV[i][j]*Pw*Ss*(PI/4)*pow(populationAdd[0][i][j],2)*Vm[i][j];                                
                 c[j]=c[j]+W[i][j];                
              }                    
            }          
        
        for(i=0;i<m;i++)
        {                  
         tot_pcap=tot_pcap+pcap[i];                 
         }   
         
         for(j=0;j<n;j++)
        {                  
         tot_ccap+=ccap[j];                 
         }       
     //Constraint Vl for each individual 
       
        Vl=0;     
        if(tot_pcap>tot_ccap)
        {
          for(i=0;i<m;i++)
           {
           if( p[i]<0 )                 
             {
                Vl=Vl+PenaltyCoef*pow(p[i],2);              
             }
            if( p[i]>pcap[i])                 
             {
                Vl=Vl+PenaltyCoef*pow((p[i]-pcap[i]),2);              
             }  
           }
           
           for(j=0;j<n;j++)
           {
           if( c[j]<0.99*ccap[j] )                 
             {
                Vl=Vl+PenaltyCoef*pow((0.99*ccap[j]-c[j]),2);              
             }
            if( c[j]>ccap[j])                 
             {
                Vl=Vl+PenaltyCoef*pow((c[j]-ccap[j]),2);              
             }  
           }
            
        }  
        else
        {
          for(j=0;j<n;j++)
           {
           if( c[j]<0 )                 
             {
                Vl=Vl+PenaltyCoef*pow(c[j],2);              
             }
            if( c[j]>ccap[j])                 
             {
                Vl=Vl+PenaltyCoef*pow((c[j]-ccap[j]),2);              
             }  
           }
           
           for(i=0;i<m;i++)
           {
           if( p[i]<0.99*pcap[i] )                 
             {
                Vl=Vl+PenaltyCoef*pow((0.99*pcap[i]-p[i]),2);              
             }
            if( p[i]>pcap[i])                 
             {
                Vl=Vl+PenaltyCoef*pow((p[i]-pcap[i]),2);              
             }  
           }
           
        }//end of else
          
       popFitnessAdd=popFitnessAdd+Vl;    
    
}

void fitnessPopSub()
{
  int i,j;
   double C1,C2;
   double tot_pcap=0;
   double tot_ccap=0;
   
  ///Calculation of Fitness C2  
          
          C2=0;
           for(i=0;i<m;i++)
           {
               for(j=0;j<n;j++)
               {      
                   C2=C2+(210.89*pow(populationSub[0][i][j],1.3744))*L[i][j];   
                   //display each C2                  
               }
           }
           C2=C2/1000000;
            popFitnessSub=C2;               
      
   ///Calculation of Fitness C1
      
           C1=0;
          for(i=0;i<m;i++)
           {
           for(j=0;j<n;j++)
            {      
              //Calculate Roij

              CV[i][j]=(populationSub[1][i][j])/(populationSub[1][i][j]+4.74*(1-populationSub[1][i][j]));             
              Ro[i][j]=Pw*(CV[i][j]*Ss+(1-CV[i][j]));  
              if(populationSub[1][i][j]<0.3)
                {
                    FCw[i][j]=1.097;         
                }
              else if(populationSub[1][i][j]<0.45)
                {
                    FCw[i][j]=0.2067*populationSub[1][i][j]+1.035;         
                }
              else if(populationSub[1][i][j]<0.55)
                {
                    FCw[i][j]=1.52*populationSub[1][i][j]+0.444;         
                }
              else if(populationSub[1][i][j]<0.7)
                {
                    FCw[i][j]=6.1*populationSub[1][i][j]+2.075;         
                }
          
            Vm[i][j]= 2966.45*(FCw[i][j]*pow(d,0.75)*pow(4.74,0.5)*pow(populationSub[0][i][j],0.5));         
            Q[i][j]=(PI*(pow(populationSub[0][i][j],2)/4))*Vm[i][j];             
            H[i][j]=0.0039* pow(CV[i][j],0.803)*pow(populationSub[0][i][j],-1.25)*pow(Vm[i][j],1.77)*L[i][j];         
            C1=(C1+(CE*7884*((Ro[i][j]*Q[i][j]*H[i][j])/101.94)));
          } 
     
            C1=C1/1000000;

            popFitnessSub =popFitnessSub+C1; 
       
     }// end of indv
        
   //// Constarints:
     
          for(i=0;i<m;i++)
           {
           p[i]=0;
           for(j=0;j<n;j++)
              {                  
                 W[i][j]=CV[i][j]*Pw*Ss*(PI/4)*pow(populationSub[0][i][j],2)*Vm[i][j];                                
                 p[i]=p[i]+W[i][j];                
              }                    
            }
          
          for(j=0;j<n;j++)
           {
           c[j]=0;
           for(i=0;i<m;i++)
              {                  
                 W[i][j]=CV[i][j]*Pw*Ss*(PI/4)*pow(populationSub[0][i][j],2)*Vm[i][j];                                
                 c[j]=c[j]+W[i][j];                
              }                    
            }
          
        
        for(i=0;i<m;i++)
        {                  
         tot_pcap=tot_pcap+pcap[i];                 
         }   
         
         for(j=0;j<n;j++)
        {                  
         tot_ccap+=ccap[j];                 
         }       
     //Constraint Vl for each individual 
       
        Vl=0;     
        if(tot_pcap>tot_ccap)
        {
          for(i=0;i<m;i++)
           {
           if( p[i]<0 )                 
             {
                Vl=Vl+PenaltyCoef*pow(p[i],2);              
             }
            if( p[i]>pcap[i])                 
             {
                Vl=Vl+PenaltyCoef*pow((p[i]-pcap[i]),2);              
             }  
           }
           
           for(j=0;j<n;j++)
           {
           if( c[j]<0.99*ccap[j] )                 
             {
                Vl=Vl+PenaltyCoef*pow((0.99*ccap[j]-c[j]),2);              
             }
            if( c[j]>ccap[j])                 
             {
                Vl=Vl+PenaltyCoef*pow((c[j]-ccap[j]),2);              
             }  
           }
            
        }  
        else
        {
          for(j=0;j<n;j++)
           {
           if( c[j]<0 )                 
             {
                Vl=Vl+PenaltyCoef*pow(c[j],2);              
             }
            if( c[j]>ccap[j])                 
             {
                Vl=Vl+PenaltyCoef*pow((c[j]-ccap[j]),2);              
             }  
           }
           
           for(i=0;i<m;i++)
           {
           if( p[i]<0.99*pcap[i] )                 
             {
                Vl=Vl+PenaltyCoef*pow((0.99*pcap[i]-p[i]),2);              
             }
            if( p[i]>pcap[i])                 
             {
                Vl=Vl+PenaltyCoef*pow((p[i]-pcap[i]),2);              
             }  
           }
           
        }//end of else
       popFitnessSub=popFitnessSub+Vl;       
  
}
//////
void updateD2(int i, int j)
{         
   //update indv and calculate fitness   
//   SigmaD=0.01, [0.1 0.9]
//   SigmaCW=0.001 [0.3 0.7]
long double t1,t2;
t1=populationAdd [0][i][j]+SigmaD2;
t2=populationSub [0][i][j]-SigmaD2;

if (t1<0.09)
{
	t1=0;
}
else if ((t1>=0.09)&&(t1<0.11))
{
	t1=0.10;
}
else if ((t1>=0.11)&&(t1<0.135))
{
	t1=0.12;
}
else if ((t1>=0.135)&&(t1<0.175))
{
	t1=0.15;
}
else if ((t1>=0.175)&&(t1<0.225))
{
	t1=0.2;
}
else if ((t1>=0.225)&&(t1<0.275))
{
	t1=0.25;
}
else if ((t1>=0.275)&&(t1<0.325))
{
	t1=0.3;
}
else if ((t1>=0.325)&&(t1<0.375))
{
	t1=0.35;
}
else if ((t1>=0.375)&&(t1<0.425))
{
	t1=0.4;
}
else if ((t1>=0.425)&&(t1<0.475))
{
	t1=0.45;
}
else if ((t1>=0.475)&&(t1<0.525))
{
	t1=0.5;
}
else if ((t1>=0.525)&&(t1<0.575))
{
	t1=0.55;
}
else if ((t1>=0.575)&&(t1<0.625))
{
	t1=0.6;
}
else if ((t1>=0.625)&&(t1<0.675))
{
	t1=0.65;
}
else if ((t1>=0.675)&&(t1<0.725))
{
	t1=0.7;
}
else if ((t1>=0.725)&&(t1<0.775))
{
	t1=0.75;
}
else if ((t1>=0.775)&&(t1<0.825))
{
	t1=0.8;
}
else if ((t1>=0.825)&&(t1<0.875))
{
	t1=0.85;
}
else if ((t1>=0.875)&&(t1<0.925))
{
	t1=0.9;
}
else if ((t1>=0.925)&&(t1<0.975))
{
	t1=0.95;
}
else if ((t1>=0.975)&&(t1<=1)) 
{
	t1=1;
}

   populationAdd [0][i][j]=t1;

if (t2<0.09)
{
	t2=0;
}
else if ((t2>=0.09)&&(t2<0.11))
{
	t2=0.10;
}
else if ((t2>=0.11)&&(t2<0.135))
{
	t2=0.12;
}
else if ((t2>=0.135)&&(t2<0.175))
{
	t2=0.15;
}
else if ((t2>=0.175)&&(t2<0.225))
{
	t2=0.2;
}
else if ((t2>=0.225)&&(t2<0.275))
{
	t2=0.25;
}
else if ((t2>=0.275)&&(t2<0.325))
{
	t2=0.3;
}
else if ((t2>=0.325)&&(t2<0.375))
{
	t2=0.35;
}
else if ((t2>=0.375)&&(t2<0.425))
{
	t2=0.4;
}
else if ((t2>=0.425)&&(t2<0.475))
{
	t2=0.45;
}
else if ((t2>=0.475)&&(t2<0.525))
{
	t2=0.5;
}
else if ((t2>=0.525)&&(t2<0.575))
{
	t2=0.55;
}
else if ((t2>=0.575)&&(t2<0.625))
{
	t2=0.6;
}
else if ((t2>=0.625)&&(t2<0.675))
{
	t2=0.65;
}
else if ((t2>=0.675)&&(t2<0.725))
{
	t2=0.7;
}
else if ((t2>=0.725)&&(t2<0.775))
{
	t2=0.75;
}
else if ((t2>=0.775)&&(t2<0.825))
{
	t2=0.8;
}
else if ((t2>=0.825)&&(t2<0.875))
{
	t2=0.85;
}
else if ((t2>=0.875)&&(t2<0.925))
{
	t2=0.9;
}
else if ((t2>=0.925)&&(t2<0.975))
{
	t2=0.95;
}
else if ((t2>=0.975)&&(t2<=1))
{
	t2=1;
}
   
   populationSub [0][i][j]=t2;
  
   fitnessPopAdd();
   fitnessPopSub();

}
void updateCW2(int i, int j)
{         
   //update indv and calculate fitness   
//   SigmaD=0.01
//   SigmaCW=0.001  
long double t1,t2;
t1=populationAdd [1][i][j]+SigmaCV2;
t2=populationSub [1][i][j]-SigmaCV2;

if (t1<0.005)
{
	t1=0;
}
else if ((t1>=0.005)&&(t1<0.015))
{
	t1=0.01;
}
else if ((t1>=0.015)&&(t1<0.025))
{
	t1=0.02;
}
else if ((t1>=0.025)&&(t1<0.035))
{
	t1=0.03;
}
else if ((t1>=0.035)&&(t1<0.045))
{
	t1=0.04;
}
else if ((t1>=0.045)&&(t1<0.055))
{
	t1=0.05;
}
else if ((t1>=0.055)&&(t1<0.065))
{
	t1=0.06;
}
else if ((t1>=0.065)&&(t1<0.075))
{
	t1=0.07;
}
else if ((t1>=0.075)&&(t1<0.085))
{
	t1=0.08;
}
else if ((t1>=0.085)&&(t1<0.095))
{
	t1=0.09;
}
else if ((t1>=0.095)&&(t1<0.105))
{
	t1=0.1;
}
else if ((t1>=0.105)&&(t1<0.115))
{
	t1=0.11;
}
else if ((t1>=0.115)&&(t1<0.125))
{
	t1=0.12;
}
else if ((t1>=0.125)&&(t1<0.135))
{
	t1=0.13;
}
else if ((t1>=0.135)&&(t1<0.145))
{
	t1=0.14;
}
else if ((t1>=0.145)&&(t1<0.155))
{
	t1=0.15;
}
else if ((t1>=0.155)&&(t1<0.165))
{
	t1=0.16;
}
else if ((t1>=0.165)&&(t1<0.175))
{
	t1=0.17;
}
else if ((t1>=0.175)&&(t1<0.185))
{
	t1=0.18;
}
else if ((t1>=0.185)&&(t1<0.195))
{
	t1=0.19;
}
else if ((t1>=0.195)&&(t1<0.205))
{
	t1=0.2;
}
else if ((t1>=0.205)&&(t1<0.215))
{
	t1=0.21;
}
else if ((t1>=0.215)&&(t1<0.225))
{
	t1=0.22;
}
else if ((t1>=0.225)&&(t1<0.235))
{
	t1=0.23;
}
else if ((t1>=0.235)&&(t1<0.245))
{
	t1=0.24;
}
else if ((t1>=0.245)&&(t1<0.255))
{
	t1=0.25;
}
else if ((t1>=0.255)&&(t1<0.265))
{
	t1=0.26;
}
else if ((t1>=0.265)&&(t1<0.275))
{
	t1=0.27;
}
else if ((t1>=0.275)&&(t1<0.285))
{
	t1=0.28;
}
else if ((t1>=0.285)&&(t1<0.295))
{
	t1=0.29;
}
else if ((t1>=0.295)&&(t1<0.305))
{
	t1=0.30;
}
else if ((t1>=0.305)&&(t1<0.315))
{
	t1=0.31;
}
else if ((t1>=0.315)&&(t1<0.325))
{
	t1=0.32;
}
else if ((t1>=0.325)&&(t1<0.335))
{
	t1=0.33;
}
else if ((t1>=0.335)&&(t1<0.345))
{
	t1=0.34;
}
else if ((t1>=0.345)&&(t1<0.355))
{
	t1=0.35;
}
else if ((t1>=0.355)&&(t1<0.365))
{
	t1=0.36;
}
else if ((t1>=0.365)&&(t1<0.375))
{
	t1=0.37;
}
else if ((t1>=0.375)&&(t1<0.385))
{
	t1=0.38;
}
else if ((t1>=0.385)&&(t1<0.395))
{
	t1=0.39;
}
else if ((t1>=0.395)&&(t1<0.405))
{
	t1=0.40;
}
else if ((t1>=0.405)&&(t1<0.415))
{
	t1=0.41;
}
else if ((t1>=0.415)&&(t1<0.425))
{
	t1=0.42;
}
else if ((t1>=0.425)&&(t1<0.435))
{
	t1=0.43;
}
else if ((t1>=0.435)&&(t1<0.445))
{
	t1=0.44;
}
else if ((t1>=0.445)&&(t1<0.455))
{
	t1=0.45;
}
else if ((t1>=0.455)&&(t1<0.465))
{
	t1=0.46;
}
else if ((t1>=0.465)&&(t1<0.475))
{
	t1=0.47;
}
else if ((t1>=0.475)&&(t1<0.485))
{
	t1=0.48;
}
else if ((t1>=0.485)&&(t1<0.495))
{
	t1=0.49;
}
else if ((t1>=0.495)&&(t1<0.505))
{
	t1=0.5;
}
else if ((t1>=0.505)&&(t1<0.515))
{
	t1=0.51;
}
else if ((t1>=0.515)&&(t1<0.525))
{
	t1=0.52;
}
else if ((t1>=0.525)&&(t1<0.535))
{
	t1=0.53;
}
else if ((t1>=0.535)&&(t1<0.545))
{
	t1=0.54;
}
else if ((t1>=0.545)&&(t1<0.555))
{
	t1=0.55;
}
else if ((t1>=0.555)&&(t1<0.565))
{
	t1=0.56;
}
else if ((t1>=0.565)&&(t1<0.575))
{
	t1=0.57;
}
else if ((t1>=0.575)&&(t1<0.585))
{
	t1=0.58;
}
else if ((t1>=0.585)&&(t1<0.595))
{
	t1=0.59;
}
else if ((t1>=0.595)&&(t1<0.605))
{
	t1=0.6;
}
else if ((t1>=0.605)&&(t1<0.615))
{
	t1=0.61;
}
else if ((t1>=0.615)&&(t1<0.625))
{
	t1=0.62;
}
else if ((t1>=0.625)&&(t1<0.635))
{
	t1=0.63;
}
else if ((t1>=0.635)&&(t1<0.645))
{
	t1=0.64;
}
else if ((t1>=0.645)&&(t1<0.655))
{
	t1=0.65;
}
else if ((t1>=0.655)&&(t1<0.665))
{
	t1=0.66;
}
else if ((t1>=0.665)&&(t1<0.675))
{
	t1=0.67;
}
else if ((t1>=0.675)&&(t1<0.685))
{
	t1=0.68;
}
else if ((t1>=0.685)&&(t1<0.695))
{
	t1=0.69;
}
else if ((t1>=0.695)&&(t1<0.7))
{
	t1=0.6975;
}

     populationAdd [1][i][j]=t1;


if (t2<0.005)
{
	t2=0;
}
else if ((t2>=0.005)&&(t2<0.015))
{
	t2=0.01;
}
else if ((t2>=0.015)&&(t2<0.025))
{
	t2=0.02;
}
else if ((t2>=0.025)&&(t2<0.035))
{
	t2=0.03;
}
else if ((t2>=0.035)&&(t2<0.045))
{
	t2=0.04;
}
else if ((t2>=0.045)&&(t2<0.055))
{
	t2=0.05;
}
else if ((t2>=0.055)&&(t2<0.065))
{
	t2=0.06;
}
else if ((t2>=0.065)&&(t2<0.075))
{
	t2=0.07;
}
else if ((t2>=0.075)&&(t2<0.085))
{
	t2=0.08;
}
else if ((t2>=0.085)&&(t2<0.095))
{
	t2=0.09;
}
else if ((t2>=0.095)&&(t2<0.105))
{
	t2=0.1;
}
else if ((t2>=0.105)&&(t2<0.115))
{
	t2=0.11;
}
else if ((t2>=0.115)&&(t2<0.125))
{
	t2=0.12;
}
else if ((t2>=0.125)&&(t2<0.135))
{
	t2=0.13;
}
else if ((t2>=0.135)&&(t2<0.145))
{
	t2=0.14;
}
else if ((t2>=0.145)&&(t2<0.155))
{
	t2=0.15;
}
else if ((t2>=0.155)&&(t2<0.165))
{
	t2=0.16;
}
else if ((t2>=0.165)&&(t2<0.175))
{
	t2=0.17;
}
else if ((t2>=0.175)&&(t2<0.185))
{
	t2=0.18;
}
else if ((t2>=0.185)&&(t2<0.195))
{
	t2=0.19;
}
else if ((t2>=0.195)&&(t2<0.205))
{
	t2=0.2;
}
else if ((t2>=0.205)&&(t2<0.215))
{
	t2=0.21;
}
else if ((t2>=0.215)&&(t2<0.225))
{
	t2=0.22;
}
else if ((t2>=0.225)&&(t2<0.235))
{
	t2=0.23;
}
else if ((t2>=0.235)&&(t2<0.245))
{
	t2=0.24;
}
else if ((t2>=0.245)&&(t2<0.255))
{
	t2=0.25;
}
else if ((t2>=0.255)&&(t2<0.265))
{
	t2=0.26;
}
else if ((t2>=0.265)&&(t2<0.275))
{
	t2=0.27;
}
else if ((t2>=0.275)&&(t2<0.285))
{
	t2=0.28;
}
else if ((t2>=0.285)&&(t2<0.295))
{
	t2=0.29;
}
else if ((t2>=0.295)&&(t2<0.305))
{
	t2=0.30;
}
else if ((t2>=0.305)&&(t2<0.315))
{
	t2=0.31;
}
else if ((t2>=0.315)&&(t2<0.325))
{
	t2=0.32;
}
else if ((t1>=0.325)&&(t1<0.335))
{
	t1=0.33;
}
else if ((t2>=0.335)&&(t2<0.345))
{
	t2=0.34;
}
else if ((t2>=0.345)&&(t2<0.355))
{
	t2=0.35;
}
else if ((t2>=0.355)&&(t2<0.365))
{
	t2=0.36;
}
else if ((t2>=0.365)&&(t2<0.375))
{
	t2=0.37;
}
else if ((t2>=0.375)&&(t2<0.385))
{
	t2=0.38;
}
else if ((t2>=0.385)&&(t2<0.395))
{
	t2=0.39;
}
else if ((t2>=0.395)&&(t2<0.405))
{
	t2=0.40;
}
else if ((t2>=0.405)&&(t2<0.415))
{
	t2=0.41;
}
else if ((t2>=0.415)&&(t2<0.425))
{
	t2=0.42;
}
else if ((t2>=0.425)&&(t2<0.435))
{
	t2=0.43;
}
else if ((t2>=0.435)&&(t2<0.445))
{
	t2=0.44;
}
else if ((t2>=0.445)&&(t2<0.455))
{
	t2=0.45;
}
else if ((t2>=0.455)&&(t2<0.465))
{
	t2=0.46;
}
else if ((t2>=0.465)&&(t2<0.475))
{
	t2=0.47;
}
else if ((t2>=0.475)&&(t2<0.485))
{
	t2=0.48;
}
else if ((t2>=0.485)&&(t2<0.495))
{
	t2=0.49;
}
else if ((t2>=0.495)&&(t2<0.505))
{
	t2=0.5;
}
else if ((t2>=0.505)&&(t2<0.515))
{
	t2=0.51;
}
else if ((t2>=0.515)&&(t2<0.525))
{
	t2=0.52;
}
else if ((t2>=0.525)&&(t2<0.535))
{
	t2=0.53;
}
else if ((t2>=0.535)&&(t2<0.545))
{
	t2=0.54;
}
else if ((t2>=0.545)&&(t2<0.555))
{
	t2=0.55;
}
else if ((t2>=0.555)&&(t2<0.565))
{
	t2=0.56;
}
else if ((t2>=0.565)&&(t2<0.575))
{
	t2=0.57;
}
else if ((t2>=0.575)&&(t2<0.585))
{
	t2=0.58;
}
else if ((t2>=0.585)&&(t2<0.595))
{
	t2=0.59;
}
else if ((t2>=0.595)&&(t2<0.605))
{
	t2=0.6;
}
else if ((t2>=0.605)&&(t2<0.615))
{
	t2=0.61;
}
else if ((t1>=0.615)&&(t1<0.625))
{
	t2=0.62;
}
else if ((t2>=0.625)&&(t2<0.635))
{
	t2=0.63;
}
else if ((t2>=0.635)&&(t2<0.645))
{
	t2=0.64;
}
else if ((t2>=0.645)&&(t2<0.655))
{
	t2=0.65;
}
else if ((t2>=0.655)&&(t2<0.665))
{
	t2=0.66;
}
else if ((t2>=0.665)&&(t2<0.675))
{
	t2=0.67;
}
else if ((t2>=0.675)&&(t2<0.685))
{
	t2=0.68;
}
else if ((t2>=0.685)&&(t2<0.695))
{
	t2=0.69;
}
else if ((t2>=0.695)&&(t2<0.7))
{
	t2=0.6975;
}
 
     populationSub [1][i][j]=t2;

   fitnessPopAdd();
   fitnessPopSub();
   
}
//////
void updateD(int i, int j)
{   
long double t1,t2;      
   //update indv and calculate fitness   
//   SigmaD=0.01, [0.1 0.9]
//   SigmaCW=0.001 [0.3 0.7]
t1=populationAdd [0][i][j]+SigmaD;
t2=populationSub [0][i][j]-SigmaD;

if (t1<0.09)
{
	t1=0;
}
else if ((t1>=0.09)&&(t1<0.11))
{
	t1=0.10;
}
else if ((t1>=0.11)&&(t1<0.135))
{
	t1=0.12;
}
else if ((t1>=0.135)&&(t1<0.175))
{
	t1=0.15;
}
else if ((t1>=0.175)&&(t1<0.225))
{
	t1=0.2;
}
else if ((t1>=0.225)&&(t1<0.275))
{
	t1=0.25;
}
else if ((t1>=0.275)&&(t1<0.325))
{
	t1=0.3;
}
else if ((t1>=0.325)&&(t1<0.375))
{
	t1=0.35;
}
else if ((t1>=0.375)&&(t1<0.425))
{
	t1=0.4;
}
else if ((t1>=0.425)&&(t1<0.475))
{
	t1=0.45;
}
else if ((t1>=0.475)&&(t1<0.525))
{
	t1=0.5;
}
else if ((t1>=0.525)&&(t1<0.575))
{
	t1=0.55;
}
else if ((t1>=0.575)&&(t1<0.625))
{
	t1=0.6;
}
else if ((t1>=0.625)&&(t1<0.675))
{
	t1=0.65;
}
else if ((t1>=0.675)&&(t1<0.725))
{
	t1=0.7;
}
else if ((t1>=0.725)&&(t1<0.775))
{
	t1=0.75;
}
else if ((t1>=0.775)&&(t1<0.825))
{
	t1=0.8;
}
else if ((t1>=0.825)&&(t1<0.875))
{
	t1=0.85;
}
else if ((t1>=0.875)&&(t1<0.925))
{
	t1=0.9;
}
else if ((t1>=0.925)&&(t1<0.975))
{
	t1=0.95;
}
else if ((t1>=0.975)&&(t1<=1))//else if ((t1>=0.975)&&(t1<=1))
{
	t1=1;
}

   populationAdd [0][i][j]=t1;


if (t2<0.09)
{
	t2=0;
}
else if ((t2>=0.09)&&(t2<0.11))
{
	t2=0.10;
}
else if ((t2>=0.11)&&(t2<0.135))
{
	t2=0.12;
}
else if ((t2>=0.135)&&(t2<0.175))
{
	t2=0.15;
}
else if ((t2>=0.175)&&(t2<0.225))
{
	t2=0.2;
}
else if ((t2>=0.225)&&(t2<0.275))
{
	t2=0.25;
}
else if ((t2>=0.275)&&(t2<0.325))
{
	t2=0.3;
}
else if ((t2>=0.325)&&(t2<0.375))
{
	t2=0.35;
}
else if ((t2>=0.375)&&(t2<0.425))
{
	t2=0.4;
}
else if ((t2>=0.425)&&(t2<0.475))
{
	t2=0.45;
}
else if ((t2>=0.475)&&(t2<0.525))
{
	t2=0.5;
}
else if ((t2>=0.525)&&(t2<0.575))
{
	t2=0.55;
}
else if ((t2>=0.575)&&(t2<0.625))
{
	t2=0.6;
}
else if ((t2>=0.625)&&(t2<0.675))
{
	t2=0.65;
}
else if ((t2>=0.675)&&(t2<0.725))
{
	t2=0.7;
}
else if ((t2>=0.725)&&(t2<0.775))
{
	t2=0.75;
}
else if ((t2>=0.775)&&(t2<0.825))
{
	t2=0.8;
}
else if ((t2>=0.825)&&(t2<0.875))
{
	t2=0.85;
}
else if ((t2>=0.875)&&(t2<0.925))
{
	t2=0.9;
}
else if ((t2>=0.925)&&(t2<0.975))
{
	t2=0.95;
}
else if ((t2>=0.975)&&(t2<=1))//else if ((t1>=0.975)&&(t1<=1))
{
	t2=1;
}
// if (populationSub [0][i][j]-SigmaD>=0) //  if (populationSub [0][i][j]-SigmaD>=0.1)
   populationSub [0][i][j]=t2;//populationSub [0][i][j]-SigmaD;
  
   fitnessPopAdd();
   fitnessPopSub();

}
void updateCW(int i, int j)
{         
   //update indv and calculate fitness   
//   SigmaD=0.01
//   SigmaCW=0.001    
long double t1,t2;
t1=populationAdd [1][i][j]+SigmaCV;
t2=populationSub [1][i][j]-SigmaCV;
if (t1<0.005)
{
	t1=0;
}
else if ((t1>=0.005)&&(t1<0.015))
{
	t1=0.01;
}
else if ((t1>=0.015)&&(t1<0.025))
{
	t1=0.02;
}
else if ((t1>=0.025)&&(t1<0.035))
{
	t1=0.03;
}
else if ((t1>=0.035)&&(t1<0.045))
{
	t1=0.04;
}
else if ((t1>=0.045)&&(t1<0.055))
{
	t1=0.05;
}
else if ((t1>=0.055)&&(t1<0.065))
{
	t1=0.06;
}
else if ((t1>=0.065)&&(t1<0.075))
{
	t1=0.07;
}
else if ((t1>=0.075)&&(t1<0.085))
{
	t1=0.08;
}
else if ((t1>=0.085)&&(t1<0.095))
{
	t1=0.09;
}
else if ((t1>=0.095)&&(t1<0.105))
{
	t1=0.1;
}
else if ((t1>=0.105)&&(t1<0.115))
{
	t1=0.11;
}
else if ((t1>=0.115)&&(t1<0.125))
{
	t1=0.12;
}
else if ((t1>=0.125)&&(t1<0.135))
{
	t1=0.13;
}
else if ((t1>=0.135)&&(t1<0.145))
{
	t1=0.14;
}
else if ((t1>=0.145)&&(t1<0.155))
{
	t1=0.15;
}
else if ((t1>=0.155)&&(t1<0.165))
{
	t1=0.16;
}
else if ((t1>=0.165)&&(t1<0.175))
{
	t1=0.17;
}
else if ((t1>=0.175)&&(t1<0.185))
{
	t1=0.18;
}
else if ((t1>=0.185)&&(t1<0.195))
{
	t1=0.19;
}
else if ((t1>=0.195)&&(t1<0.205))
{
	t1=0.2;
}
else if ((t1>=0.205)&&(t1<0.215))
{
	t1=0.21;
}
else if ((t1>=0.215)&&(t1<0.225))
{
	t1=0.22;
}
else if ((t1>=0.225)&&(t1<0.235))
{
	t1=0.23;
}
else if ((t1>=0.235)&&(t1<0.245))
{
	t1=0.24;
}
else if ((t1>=0.245)&&(t1<0.255))
{
	t1=0.25;
}
else if ((t1>=0.255)&&(t1<0.265))
{
	t1=0.26;
}
else if ((t1>=0.265)&&(t1<0.275))
{
	t1=0.27;
}
else if ((t1>=0.275)&&(t1<0.285))
{
	t1=0.28;
}
else if ((t1>=0.285)&&(t1<0.295))
{
	t1=0.29;
}
else if ((t1>=0.295)&&(t1<0.305))
{
	t1=0.30;
}
else if ((t1>=0.305)&&(t1<0.315))
{
	t1=0.31;
}
else if ((t1>=0.315)&&(t1<0.325))
{
	t1=0.32;
}
else if ((t1>=0.325)&&(t1<0.335))
{
	t1=0.33;
}
else if ((t1>=0.335)&&(t1<0.345))
{
	t1=0.34;
}
else if ((t1>=0.345)&&(t1<0.355))
{
	t1=0.35;
}
else if ((t1>=0.355)&&(t1<0.365))
{
	t1=0.36;
}
else if ((t1>=0.365)&&(t1<0.375))
{
	t1=0.37;
}
else if ((t1>=0.375)&&(t1<0.385))
{
	t1=0.38;
}
else if ((t1>=0.385)&&(t1<0.395))
{
	t1=0.39;
}
else if ((t1>=0.395)&&(t1<0.405))
{
	t1=0.40;
}
else if ((t1>=0.405)&&(t1<0.415))
{
	t1=0.41;
}
else if ((t1>=0.415)&&(t1<0.425))
{
	t1=0.42;
}
else if ((t1>=0.425)&&(t1<0.435))
{
	t1=0.43;
}
else if ((t1>=0.435)&&(t1<0.445))
{
	t1=0.44;
}
else if ((t1>=0.445)&&(t1<0.455))
{
	t1=0.45;
}
else if ((t1>=0.455)&&(t1<0.465))
{
	t1=0.46;
}
else if ((t1>=0.465)&&(t1<0.475))
{
	t1=0.47;
}
else if ((t1>=0.475)&&(t1<0.485))
{
	t1=0.48;
}
else if ((t1>=0.485)&&(t1<0.495))
{
	t1=0.49;
}
else if ((t1>=0.495)&&(t1<0.505))
{
	t1=0.5;
}
else if ((t1>=0.505)&&(t1<0.515))
{
	t1=0.51;
}
else if ((t1>=0.515)&&(t1<0.525))
{
	t1=0.52;
}
else if ((t1>=0.525)&&(t1<0.535))
{
	t1=0.53;
}
else if ((t1>=0.535)&&(t1<0.545))
{
	t1=0.54;
}
else if ((t1>=0.545)&&(t1<0.555))
{
	t1=0.55;
}
else if ((t1>=0.555)&&(t1<0.565))
{
	t1=0.56;
}
else if ((t1>=0.565)&&(t1<0.575))
{
	t1=0.57;
}
else if ((t1>=0.575)&&(t1<0.585))
{
	t1=0.58;
}
else if ((t1>=0.585)&&(t1<0.595))
{
	t1=0.59;
}
else if ((t1>=0.595)&&(t1<0.605))
{
	t1=0.6;
}
else if ((t1>=0.605)&&(t1<0.615))
{
	t1=0.61;
}
else if ((t1>=0.615)&&(t1<0.625))
{
	t1=0.62;
}
else if ((t1>=0.625)&&(t1<0.635))
{
	t1=0.63;
}
else if ((t1>=0.635)&&(t1<0.645))
{
	t1=0.64;
}
else if ((t1>=0.645)&&(t1<0.655))
{
	t1=0.65;
}
else if ((t1>=0.655)&&(t1<0.665))
{
	t1=0.66;
}
else if ((t1>=0.665)&&(t1<0.675))
{
	t1=0.67;
}
else if ((t1>=0.675)&&(t1<0.685))
{
	t1=0.68;
}
else if ((t1>=0.685)&&(t1<0.695))
{
	t1=0.69;
}
else if ((t1>=0.695)&&(t1<0.7))
{
	t1=0.6975;
}


 //  if (populationAdd [1][i][j]+SigmaCV<=0.7)
     populationAdd [1][i][j]=t1;//populationAdd [1][i][j]+SigmaCV;

if (t2<0.005)
{
	t2=0;
}
else if ((t2>=0.005)&&(t2<0.015))
{
	t2=0.01;
}
else if ((t2>=0.015)&&(t2<0.025))
{
	t2=0.02;
}
else if ((t2>=0.025)&&(t2<0.035))
{
	t2=0.03;
}
else if ((t2>=0.035)&&(t2<0.045))
{
	t2=0.04;
}
else if ((t2>=0.045)&&(t2<0.055))
{
	t2=0.05;
}
else if ((t2>=0.055)&&(t2<0.065))
{
	t2=0.06;
}
else if ((t2>=0.065)&&(t2<0.075))
{
	t2=0.07;
}
else if ((t2>=0.075)&&(t2<0.085))
{
	t2=0.08;
}
else if ((t2>=0.085)&&(t2<0.095))
{
	t2=0.09;
}
else if ((t2>=0.095)&&(t2<0.105))
{
	t2=0.1;
}
else if ((t2>=0.105)&&(t2<0.115))
{
	t2=0.11;
}
else if ((t2>=0.115)&&(t2<0.125))
{
	t2=0.12;
}
else if ((t2>=0.125)&&(t2<0.135))
{
	t2=0.13;
}
else if ((t2>=0.135)&&(t2<0.145))
{
	t2=0.14;
}
else if ((t2>=0.145)&&(t2<0.155))
{
	t2=0.15;
}
else if ((t2>=0.155)&&(t2<0.165))
{
	t2=0.16;
}
else if ((t2>=0.165)&&(t2<0.175))
{
	t2=0.17;
}
else if ((t2>=0.175)&&(t2<0.185))
{
	t2=0.18;
}
else if ((t2>=0.185)&&(t2<0.195))
{
	t2=0.19;
}
else if ((t2>=0.195)&&(t2<0.205))
{
	t2=0.2;
}
else if ((t2>=0.205)&&(t2<0.215))
{
	t2=0.21;
}
else if ((t2>=0.215)&&(t2<0.225))
{
	t2=0.22;
}
else if ((t2>=0.225)&&(t2<0.235))
{
	t2=0.23;
}
else if ((t2>=0.235)&&(t2<0.245))
{
	t2=0.24;
}
else if ((t2>=0.245)&&(t2<0.255))
{
	t2=0.25;
}
else if ((t2>=0.255)&&(t2<0.265))
{
	t2=0.26;
}
else if ((t2>=0.265)&&(t2<0.275))
{
	t2=0.27;
}
else if ((t2>=0.275)&&(t2<0.285))
{
	t2=0.28;
}
else if ((t2>=0.285)&&(t2<0.295))
{
	t2=0.29;
}
else if ((t2>=0.295)&&(t2<0.305))
{
	t2=0.30;
}
else if ((t2>=0.305)&&(t2<0.315))
{
	t2=0.31;
}
else if ((t2>=0.315)&&(t2<0.325))
{
	t2=0.32;
}
else if ((t1>=0.325)&&(t1<0.335))
{
	t1=0.33;
}
else if ((t2>=0.335)&&(t2<0.345))
{
	t2=0.34;
}
else if ((t2>=0.345)&&(t2<0.355))
{
	t2=0.35;
}
else if ((t2>=0.355)&&(t2<0.365))
{
	t2=0.36;
}
else if ((t2>=0.365)&&(t2<0.375))
{
	t2=0.37;
}
else if ((t2>=0.375)&&(t2<0.385))
{
	t2=0.38;
}
else if ((t2>=0.385)&&(t2<0.395))
{
	t2=0.39;
}
else if ((t2>=0.395)&&(t2<0.405))
{
	t2=0.40;
}
else if ((t2>=0.405)&&(t2<0.415))
{
	t2=0.41;
}
else if ((t2>=0.415)&&(t2<0.425))
{
	t2=0.42;
}
else if ((t2>=0.425)&&(t2<0.435))
{
	t2=0.43;
}
else if ((t2>=0.435)&&(t2<0.445))
{
	t2=0.44;
}
else if ((t2>=0.445)&&(t2<0.455))
{
	t2=0.45;
}
else if ((t2>=0.455)&&(t2<0.465))
{
	t2=0.46;
}
else if ((t2>=0.465)&&(t2<0.475))
{
	t2=0.47;
}
else if ((t2>=0.475)&&(t2<0.485))
{
	t2=0.48;
}
else if ((t2>=0.485)&&(t2<0.495))
{
	t2=0.49;
}
else if ((t2>=0.495)&&(t2<0.505))
{
	t2=0.5;
}
else if ((t2>=0.505)&&(t2<0.515))
{
	t2=0.51;
}
else if ((t2>=0.515)&&(t2<0.525))
{
	t2=0.52;
}
else if ((t2>=0.525)&&(t2<0.535))
{
	t2=0.53;
}
else if ((t2>=0.535)&&(t2<0.545))
{
	t2=0.54;
}
else if ((t2>=0.545)&&(t2<0.555))
{
	t2=0.55;
}
else if ((t2>=0.555)&&(t2<0.565))
{
	t2=0.56;
}
else if ((t2>=0.565)&&(t2<0.575))
{
	t2=0.57;
}
else if ((t2>=0.575)&&(t2<0.585))
{
	t2=0.58;
}
else if ((t2>=0.585)&&(t2<0.595))
{
	t2=0.59;
}
else if ((t2>=0.595)&&(t2<0.605))
{
	t2=0.6;
}
else if ((t2>=0.605)&&(t2<0.615))
{
	t2=0.61;
}
else if ((t1>=0.615)&&(t1<0.625))
{
	t2=0.62;
}
else if ((t2>=0.625)&&(t2<0.635))
{
	t2=0.63;
}
else if ((t2>=0.635)&&(t2<0.645))
{
	t2=0.64;
}
else if ((t2>=0.645)&&(t2<0.655))
{
	t2=0.65;
}
else if ((t2>=0.655)&&(t2<0.665))
{
	t2=0.66;
}
else if ((t2>=0.665)&&(t2<0.675))
{
	t2=0.67;
}
else if ((t2>=0.675)&&(t2<0.685))
{
	t2=0.68;
}
else if ((t2>=0.685)&&(t2<0.695))
{
	t2=0.69;
}
else if ((t2>=0.695)&&(t2<0.7))
{
	t2=0.6975;
}
 
     populationSub [1][i][j]=t2;
   fitnessPopAdd();
   fitnessPopSub();
   
}

//////////////////////////////////////////////////////
void initializeBest(int z)
{
  int i,j;
  //initialize base matrices  
       for(i=0;i<m;i++)
           {
               for(j=0;j<n;j++)
               {      
                     bestIndividual[0][i][j]= population[z][0][i][j];  
                     bestIndividual[1][i][j]= population[z][1][i][j];
                     
                    
               }
           }               
          
 }
///////////////////////////////////////////////////////////////////////////////////////
void PrintBest()
{
  int i,j;
  
  fprintf(fpBestIndividual,"\nBest Individual D matrix\n");
       for(i=0;i<m;i++)
           {
               for(j=0;j<n;j++)
               {      
                    fprintf(fpBestIndividual,"%f\t",bestIndividual[0][i][j]);                                                               
               }
               fprintf(fpBestIndividual,"\n ");
           }               
  fprintf(fpBestIndividual,"\nBest Individual CW matrix\n");
           for(i=0;i<m;i++)
           {
               for(j=0;j<n;j++)
               {                                               
                     fprintf(fpBestIndividual,"%f\t",bestIndividual[1][i][j]);                   
               }
               fprintf(fpBestIndividual,"\n ");
           }  
                 
 }
/////////////////////////////////
void select3()
{ 
  float random, popFitnessAdd1,popFitnessSub1;
  int  index,individual,x,y,option;
  srand(time(NULL));
   for (index=0; index<maxIteration;index++)
   {
       
       for (individual=0;individual<N;individual++)
       {
           initializationBase(individual);
           for(x=0;x<m;x++)
           {
            for(y=0;y<n;y++)
            {              
                  updateD(x,y);
                  popFitnessAdd1=popFitnessAdd;  
                  popFitnessSub1=popFitnessSub;
                  updateD2(x,y);
                  
                  if ((popFitnessAdd1 < popFitnessSub1)&&(popFitnessAdd1 <= popFitnessSub)&&(popFitnessAdd1 <= popFitnessAdd))
                     option=1;
                  else if ((popFitnessSub1 < popFitnessAdd1)&&(popFitnessSub1 <= popFitnessSub)&&(popFitnessSub1 <= popFitnessAdd))
                     option=2;
                  else if ((popFitnessAdd < popFitnessSub)&&(popFitnessAdd <= popFitnessSub1)&&(popFitnessAdd <= popFitnessAdd1))
                     option=3;
                  else if ((popFitnessSub < popFitnessAdd)&&(popFitnessSub <= popFitnessSub)&&(popFitnessSub <= popFitnessAdd1))
                     option=4;
                  switch(option)
                  {
                     case 1:                          
                          initializationBase(individual);
                          updateD(x,y);
                          if(popFitness [individual]>popFitnessAdd)
                          {
                            population[individual][0][x][y]=populationAdd[0][x][y];
                            fitnessPop(individual);                
                          }
                          break;
                     case 2:                          
                          initializationBase(individual);
                          updateD(x,y);
                          if(popFitness [individual]>popFitnessSub)
                          {
                            population[individual][0][x][y]=populationSub[0][x][y];
                            fitnessPop(individual);                
                          }
                          break;
                     case 3:
                          if(popFitness [individual]>popFitnessAdd)
                          {
                            population[individual][0][x][y]=populationAdd[0][x][y];
                            fitnessPop(individual);                
                          }
                          break;
                     case 4:
                          if(popFitness [individual]>popFitnessSub)
                          {
                            population[individual][0][x][y]=populationSub[0][x][y];
                            fitnessPop(individual);                
                          }
                          break;                                
                  }// end of switch case                                    
             
                  updateCW(x,y);
                  popFitnessAdd1=popFitnessAdd;  
                  popFitnessSub1=popFitnessSub;
                  updateCW2(x,y);
                  
                  if ((popFitnessAdd1 < popFitnessSub1)&&(popFitnessAdd1 <= popFitnessSub)&&(popFitnessAdd1 <= popFitnessAdd))
                     option=1;
                  else if ((popFitnessSub1 < popFitnessAdd1)&&(popFitnessSub1 <= popFitnessSub)&&(popFitnessSub1 <= popFitnessAdd))
                     option=2;
                  else if ((popFitnessAdd < popFitnessSub)&&(popFitnessAdd <= popFitnessSub1)&&(popFitnessAdd <= popFitnessAdd1))
                     option=3;
                  else if ((popFitnessSub < popFitnessAdd)&&(popFitnessSub <= popFitnessSub)&&(popFitnessSub <= popFitnessAdd1))
                     option=4;
                  switch(option)
                  {
                     case 1:                          
                          initializationBase(individual);
                          updateD(x,y);
                          if(popFitness [individual]>popFitnessAdd)
                          {
                            population[individual][1][x][y]=populationAdd[1][x][y];
                            fitnessPop(individual);                
                          }
                          break;
                     case 2:                          
                          initializationBase(individual);
                          updateD(x,y);
                          if(popFitness [individual]>popFitnessSub)
                          {
                            population[individual][1][x][y]=populationSub[1][x][y];
                            fitnessPop(individual);                
                          }
                          break;
                     case 3:
                          if(popFitness [individual]>popFitnessAdd)
                          {
                            population[individual][1][x][y]=populationAdd[1][x][y];
                            fitnessPop(individual);                
                          }
                          break;
                     case 4:
                          if(popFitness [individual]>popFitnessSub)
                          {
                            population[individual][1][x][y]=populationSub[1][x][y];
                            fitnessPop(individual);                
                          }
                          break;                                
                  }// end of switch case 
                  
              
      /////////////////               
            //update best
             if (bestFitness > popFitness [individual])
             {
                 bestFitness = popFitness [individual];   
                 initializeBest(individual);                        
             }
    ///////////////////////                                 
             }//end for n
           }//end for m                            
       }//end for indv
       
     fprintf(fpFitness,"\niteration:%d\n",index);
     printf("\niteration:%d finished, Best Fitness %lf\n",index,bestFitness);
     
     fprintf(fpBest,"\niteration:%d; ",index);
     fprintf(fpBest,"Best Fitness:%f; ",bestFitness);
     fprintf(fpBestIndividual,"\niteration:%d\n ",index);
     PrintBest();
   } //end for max iteration         
    printf("Update Finished....");
    getchar();
}
////////////
int main() 
{

    int choice, indv1;
         clock_t st, end;
    double run_time;
    
    printf("\nThis program solves Slurry Pipeline  Problem using Bi-Attempted Base Optimization Algorithm\n");
    printf("\nYou may get better result when you run the program several times\n");
	printf("\nThe Population Size is:%d",N);
    printf("\nThe number of iterations is:%d\n",maxIteration);
    initialization() ;
    
    do
    {
    printf("\n1-Bi-Attempted Base Optimization:");
    printf("\n2-Exit:");
    printf("\nEnter your choice:");
    scanf("%d",&choice);
    switch(choice)
    {
     
     case 1:
  			
              fpBestIndividual=fopen("BestIndividuals.txt", "w");
              fprintf(fpBestIndividual,"\nBest Optimization Started...:");
              
              fpPenalty=fopen("PenaltyFitness.txt", "w");
              fprintf(fpPenalty,"\nPenaltyFitness Started...:");
              
              fp=fopen("Population.txt", "w");
              fprintf(fp,"\n Population Started...:");
            st=clock();  
            initialization();
            
            fpBest= fopen("BestFitness.txt", "w");
            fprintf(fpBest,"\nBase Optimization Started...:");
            fprintf(fpBest,"\nPopulation size: %d, iteration: %d, SigmaD: %f,SigmaD2: %f, SigmaCV: %f, SigmaCV2: %f\n",N,maxIteration,SigmaD,SigmaD2,SigmaCV,SigmaCV2);
           for (indv1=0;indv1<N;indv1++)
            {
                fitnessPop(indv1);                
                display(indv1);
                fprintf(fp,"\n Population %d Fitness %f:\n",indv1,popFitness[indv1]);
            }
            printf("Update Started....");
            select3();
            end=clock();
               run_time = (double)(end - st) / CLOCKS_PER_SEC;
            printf("\n%f\n",run_time);
            fprintf(fpBest,"\n%f\n second",run_time);
            printf("\nBase Optimization Finished...:");
            fclose(fpIndividual);
            fclose(fpFitness);
            fclose(fpBest);
            fclose(fp);
            fclose(fpPenalty); 
            fclose(fpBestIndividual);
            printf("File closed...:");
            break;        
     case 2:     
            printf("\nExit...:");
            getchar();
            break;
     }
   } while (choice<2);
  getchar();
  return 0;
}

