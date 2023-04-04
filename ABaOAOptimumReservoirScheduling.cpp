/*
***	 BI-ATTEMPTED BASE OPTIMIZATION ALGORITHM (ABaOA)
***	 DEVELOPED BY: MEHTAP KOSE ULUKOK (https://orcid.org/0000-0003-4335-483X) AND BURHAN YILDIZ (https://orcid.org/0000-0002-0144-3562) AS A PART OF THE RESEARCH IN (doi)

THIS C PROGRAM SOLVES 4-RESERVOIR PROBLEM AS DEFINED IN HEIDARI ET AL (1971) (See Ref. in above article)
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
#define maxIteration 5000//
#define N 5000 //number of population size 
#define n 5 //number of reservoirs+1
#define t 12 // time period
#define c1 100 //coofficients for penalty
#define c2 100
#define c3 40 //must

#define SigmaD 0.1
#define SigmaD2 0.001

//minimum and the maximum reservior and storage values
const double Rmin[n]={0.0,0.0,0.0,0.0,0.0};
const double Rmax[n]={3.0,4.0,4.0,7.0,0.0};
const double Smin[n]={0.0,0.0,0.0,0.0,0.0};
const double Smax[n]={10.0,10.0,10.0,15.0,0.0};

//bi(t) is the unit return kept as b[t][i]
double b[t][n]={{1.1,1.4,1.0,1.0,1.6},{1.0,1.1,1.0,1.2,1.7},{1.0,1.0,1.2,1.8,1.8},{1.2,1.0,1.8,2.5,1.9},{1.8,1.2,2.5,2.2,2.0},{2.5,1.8,2.2,2.0,2.0},{2.2,2.5,2.0,1.8,2.0},{2.0,2.2,1.8,2.2,1.9},{1.8,2.0,2.2,1.8,1.8},{2.2,1.8,1.8,1.4,1.7},{1.8,2.2,1.4,1.1,1.6},{1.4,1.8,1.1,1.0,1.6}};
double d[n]={5,5,5,7,0};//ending min storages
double I[N][t][n];
double S[N][t+1][n];
double SA[t+1][n];
double SS[t+1][n];
double IA[t][n];
double IS[t][n];

double population[N][t][n]; 
double populationAdd[t][n];
double populationSub[t][n];
double bestIndividual[t][n];
double bestS[t][n]={{5,5,5,5,0},{0}};
double bestI[t][n]={{2,3,0,0,0},{0}};

double popFitness [N];
double popFitnessAdd;
double popFitnessSub;
double bestFitness=0.0;

//List of output files
FILE * fpFitness;
FILE * fpIndividual;
FILE * fpBestIndividual;
FILE * fpBest;
FILE *fp;
FILE *fpPenalty;

void initialization() 
{
       int indv,i,j;
       double y;
       srand(time(NULL));
       for (indv=0; indv<N;indv++)
       {
           for(i=0;i<t;i++)
           {
           		for (j=0;j<n-1;j++)
           		{
					y = ((double) (rand()%11))/10;
                	population[indv][i][j] = (Rmax[j]-Rmin[j])*y+Rmin[j];
                	if (i==0)
                		S[indv][i][j]=5;
				}
			population[indv][i][4] = -1;	
           }     
       } 
}
///////////////
void initializationBase(int k)
{
  int i,j;
    for(i=0;i<t;i++)
    {
        for(j=0;j<n;j++)
        {      
            populationAdd[i][j]= population[k][i][j];  
            populationSub[i][j]= population[k][i][j];                             
        }
    } 
	for(i=0;i<=t;i++)
    {
        for(j=0;j<n;j++)
        {      
            SA[i][j]=S[k][i][j];
			SS[i][j]=S[k][i][j];                    
        }
    }                
 }

void display(int indv) 
{
       int i,j;
    
       fprintf(fp,"\n Individual %d R Matrix txn \n",indv);
       for(i=0;i<t;i++)
       {
         for(j=0;j<n;j++)
         {      
            fprintf(fp,"%f ",population[indv][i][j]);
          }   
          fprintf(fp,"\n");  
       }      
}
///////////////////////////////////
void showpopulationBase()
{
  int mku=0;
  
  while(mku<N)
  {
     display(mku); 
     mku=mku+1;
  }
}

void fitnessPop(int indv)
{
   int i,j;
   double penalty=0,C1=0,C2=0;

 	for(i=0;i<n-1;i++)
 	{
 		for(j=0;j<t;j++)	
 		{
 			C1=C1+population[indv][j][i]*b[j][i];	 		
		}
	}
	for(j=0;j<t;j++)	
 	{
 		C2=C2+population[indv][j][3]*b[j][4];		
	}			
    popFitness[indv]=C1+C2; 
    //continuity equation of the system
    for(j=1;j<=t;j++)	
 	{
 		S[indv][j][0]=S[indv][j-1][0]+2-population[indv][j-1][0];	 		
 		S[indv][j][1]=S[indv][j-1][1]+3-population[indv][j-1][1];	 		
 		S[indv][j][2]=S[indv][j-1][2]+population[indv][j-1][1]-population[indv][j-1][2];	 		
		S[indv][j][3]=S[indv][j-1][3]+population[indv][j-1][0]+population[indv][j-1][2]-population[indv][j-1][3];	 		
	}
   
   //penalty test
   for (i=0;i<n-1;i++) 
   {
   		if (S[indv][t][i]<=d[i])
   			penalty+=c3*pow((S[indv][t][i]-d[i]),2);		
   }
   C1=0;
   C2=0;
   for (i=0;i<n-1;i++)
   {
   		for(j=0;j<=t;j++)
   		{
   			if (S[indv][j][i]>Smax[i])
			   C1+=c1*pow((Smax[i]-S[indv][j][i]),2);
			if (S[indv][j][i]<Smin[i])
			   C2+=c2*pow((Smin[i]-S[indv][j][i]),2);	
		}
   }
   popFitness[indv]-=(C1+C2+penalty);  
     
}
////////////////////////////////////////
void fitnessPopAdd()
{
     int i,j;
   double penalty=0,C1=0,C2=0;

 	for(i=0;i<n-1;i++)
 	{
 		for(j=0;j<t;j++)	
 		{
 			C1+=populationAdd[j][i]*b[j][i];	 		
		}
	}
	for(j=0;j<t;j++)	
 	{
 		C2+=populationAdd[j][3]*b[j][4];		
	}			
    popFitnessAdd=C1+C2; 
    //continuity equation of the system
    for(j=1;j<=t;j++)	
 	{
 		SA[j][0]=SA[j-1][0]+2-populationAdd[j-1][0];	 		
		SA[j][1]=SA[j-1][1]+3-populationAdd[j-1][1];	 		
		SA[j][2]=SA[j-1][2]+populationAdd[j-1][1]-populationAdd[j-1][2];	 		
		SA[j][3]=SA[j-1][3]+populationAdd[j-1][0]+populationAdd[j-1][2]-populationAdd[j-1][3];	 		
	} 
   //penalty test
   for (i=0;i<n-1;i++) 
   {
   		if (SA[t][i]<=d[i])
   		
   			penalty+=40*pow((SA[t][i]-d[i]),2);
		
		else if (SA[t][i]>d[i])
   		
   			penalty+=0;
		
   }
   C1=0;
   C2=0;
   for (i=0;i<n-1;i++)
   {
   		for(j=0;j<=t;j++)
   		{
   			if (SA[j][i]>Smax[i])
			   C1+=c1*pow((Smax[i]-SA[j][i]),2);
			if (SA[j][i]<Smin[i])
			   C2+=c2*pow((Smin[i]-SA[j][i]),2);	
		}
   }
   popFitnessAdd-=(C1+C2+penalty); 

}

void fitnessPopSub()
{
   int i,j;
   double penalty=0,C1=0,C2=0;

 	for(i=0;i<n-1;i++)
 	{
 		for(j=0;j<t;j++)	
 		{
 			C1+=populationSub[j][i]*b[j][i];	 		
		}
	}
	for(j=0;j<t;j++)	
 	{
 		C2+=populationSub[j][3]*b[j][4];		
	}			
    popFitnessSub=C1+C2; 
    //continuity equation of the system
    for(j=1;j<=t;j++)	
 	{
 		SS[j][0]=SS[j-1][0]+2-populationSub[j-1][0];	 		
		SS[j][1]=SS[j-1][1]+3-populationSub[j-1][1];	 		
		SS[j][2]=SS[j-1][2]+populationSub[j-1][1]-populationSub[j-1][2];	 		
		SS[j][3]=SS[j-1][3]+populationSub[j-1][0]+populationSub[j-1][2]-populationSub[j-1][3];	 		
	}
   
   //penalty test
   for (i=0;i<n-1;i++) 
   {
   		if (SS[t][i]<=d[i])
   		
   			penalty+=40*pow((SS[t][i]-d[i]),2);
		
		else if (SS[t][i]>d[i])
   		
   			penalty+=0;
		
   }
   C1=0;
   C2=0;
   for (i=0;i<n-1;i++)
   {
   		for(j=0;j<=t;j++)
   		{
   			if (SS[j][i]>Smax[i])
			   C1+=c1*pow((Smax[i]-SS[j][i]),2);
			if (SS[j][i]<Smin[i])
			   C2+=c2*pow((Smin[i]-SS[j][i]),2);	
		}
   }
   popFitnessSub-=(C1+C2+penalty);        
     
}
///////////////////////////////
void updateD2(int i, int j)
{         

double t1,t2,y;
t1=populationAdd [i][j]+SigmaD2;
t2=populationSub [i][j]-SigmaD2;
if (t1>Rmax[j])
{
	y = ((double) (rand()%11))/10;
    t1= (Rmax[j]-Rmin[j])*y+Rmin[j];
}
	
if (t2<Rmin[j])
{
	y = ((double) (rand()%11))/10;
    t2= (Rmax[j]-Rmin[j])*y+Rmin[j];
}	
	
   populationAdd [i][j]=t1;//populationAdd [0][i][j]+SigmaD2;
   populationSub [i][j]=t2;//populationSub [0][i][j]-SigmaD2;
  
   fitnessPopAdd();
   fitnessPopSub();
  
}
////////////////////////////////
void updateD(int i, int j)
{   
 double t1,t2,y;      

t1=populationAdd [i][j]+SigmaD;
t2=populationSub [i][j]-SigmaD;
if (t1>Rmax[j])
{
	y = ((double) (rand()%11))/10;
    t1= (Rmax[j]-Rmin[j])*y+Rmin[j];
}
	
if (t2<Rmin[j])
{
	y = ((double) (rand()%11))/10;
    t2= (Rmax[j]-Rmin[j])*y+Rmin[j];
}	

   populationAdd [i][j]=t1;
   populationSub [i][j]=t2;//populationSub [0][i][j]-SigmaD;
  
   fitnessPopAdd();
   fitnessPopSub();

}

//////////////////////////////////////////////////////
void initializeBest(int z)
{
  int i,j;
  //initialize base matrices
      
       for(i=0;i<=t;i++)
           {
               for(j=0;j<n-1;j++)
               {      
                     bestIndividual[i][j]= population[z][i][j];  
                     bestS[i][j]=S[z][i][j];
                     bestI[i][j]=I[z][i][j];
               }
           }               
          
 }
///////////////////////////////////////////////////////////////////////////////////////
void PrintBest(int ind)
{
  int i,j;
  
  fprintf(fpBestIndividual,"\niteration:%d; ",ind);
     fprintf(fpBestIndividual,"Best Fitness:%f; ",bestFitness);
     fprintf(fpBestIndividual,"\niteration:%d\n ",ind);
  
  fprintf(fpBestIndividual,"\nBest Individual R matrix\n");
       for(i=0;i<t;i++)
           {
               for(j=0;j<n;j++)
               {      
                    fprintf(fpBestIndividual,"%f\t",bestIndividual[i][j]);                                                               
               }
               fprintf(fpBestIndividual,"\n ");
           }   
	
	fprintf(fpBestIndividual,"\nBest Individual S matrix\n");
       for(i=0;i<=t;i++)
           {
               for(j=0;j<n;j++)
               {      
                    fprintf(fpBestIndividual,"%f\t",bestS[i][j]);                                                               
               }
               fprintf(fpBestIndividual,"\n ");
           }                            
 }

/////////////////////////////////
void select3()
{ 
  double random, popFitnessAdd1,popFitnessSub1;
  int  index,individual,x,y,option;
 
  srand(time(NULL));
   for (index=0; index<maxIteration;index++)
   {
       
       for (individual=0;individual<N;individual++)
       {
           initializationBase(individual);
           for(y=0;y<n-1;y++)
           {
            for(x=0;x<t;x++)
            {
                       
                  updateD(x,y);
                  popFitnessAdd1=popFitnessAdd;  
                  popFitnessSub1=popFitnessSub;
                  //newly added 26.09.2022
                  initializationBase(individual);
                  updateD2(x,y);
                  
                  if ((popFitnessAdd1 > popFitnessSub1)&&(popFitnessAdd1 >= popFitnessSub)&&(popFitnessAdd1 >= popFitnessAdd))
                     option=1;
                  else if ((popFitnessSub1 > popFitnessAdd1)&&(popFitnessSub1 >= popFitnessSub)&&(popFitnessSub1 >= popFitnessAdd))
                     option=2;
                  else if ((popFitnessAdd > popFitnessSub)&&(popFitnessAdd >= popFitnessSub1)&&(popFitnessAdd >= popFitnessAdd1))
                     option=3;
                  else if ((popFitnessSub > popFitnessAdd)&&(popFitnessSub >= popFitnessSub)&&(popFitnessSub >= popFitnessAdd1))
                     option=4;
                  switch(option)
                  {
                     case 1:                          
                          initializationBase(individual);
                          updateD(x,y);
                          if(popFitness [individual]<popFitnessAdd)
                          {
                            population[individual][x][y]=populationAdd[x][y];
                            fitnessPop(individual);                
                          }
                          break;
                     case 2:                          
                          initializationBase(individual);
                          updateD(x,y);
                          if(popFitness [individual]<popFitnessSub)
                          {
                            population[individual][x][y]=populationSub[x][y];
                            fitnessPop(individual);                
                          }
                          break;
                     case 3:
                          if(popFitness [individual]<popFitnessAdd)
                          {
                            population[individual][x][y]=populationAdd[x][y];
                            fitnessPop(individual);                
                          }
                          break;
                     case 4:
                          if(popFitness [individual]<popFitnessSub)
                          {
                            population[individual][x][y]=populationSub[x][y];
                            fitnessPop(individual);                
                          }
                          break;                                
                  }// end of switch case                                    
                            
      /////////////////               
            //update best
             if (bestFitness < popFitness [individual])
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

     PrintBest(index);
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
    initialization() ;
    printf("\nThis program solves 4-Reservior Problem using Bi-Attempted Base Optimization Algorithm\n");
    printf("\nYou may get better result when you run the program several times\n");
	printf("\nThe Population Size is:%d",N);
    printf("\nThe number of iterations is:%d\n",maxIteration);
    
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
            fprintf(fpBest,"\nPopulation size: %d, iteration: %d, SigmaD: %f,SigmaD2: %f, c1=%d, c2=%d, c3=%d \n",N,maxIteration,SigmaD,SigmaD2,c1,c2,c3);
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
            //getchar();
            break;
     }
   } while (choice<2);
  getchar();
  return 0;
}

