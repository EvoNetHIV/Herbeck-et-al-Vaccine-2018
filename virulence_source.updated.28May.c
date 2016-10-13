/************************************************************************/
/*********************    HIV NETWORK SIMULATOR    **********************/
/************************************************************************/
// modified 9-26-14 jtm 5:57 pm
/*  Simulates changes in virulence within a population */
/*  Assume N individuals infected or not infected with viruses with different infection setpoint viral loads */
/*  Assume probality of infection is proportional to viral load */
/*  Version 1: March 26, 2010 (Just got started) */
/*  Versoin 2: May 12, 2010 -- Replac edynamic model with phenomenological model */
/*                             Try to increase population size so that susceptibles */
/*                             are roughly constant */
/*  Version 3: May 19, 2010 -- Added heritabilities expressed in terms of VL setpoint */
/*                             Create version that reads parameters from a file */
/*                             Other minor changes */
/*  Version 4: Sept 2012    -- Added simple sexual network structure into this framework */
/*  Version 5: Nov 2012     -- Added terms for treatment (percent treated & efficacy)  */
/*  Version 6: Dec. 6, 2012 -- Created modular version with distinct routines for easier reading */
/*  Version ?: May 27, 2016 -- Added a vaccination campaign (all values hard-coded into the C program!) */
/*  This version branches off of Network16May2013_Rev25July2013.c */

#include "virulence_header.updated.h"  // Includes names for variables and subroutines

long Flat_Viral_Load = 0; // Setting this to 1 forces VL to be the same value during the entire primary infection period
double PrimaryInfectionLevel = 1.0e7; // Average viral load during primary infection
double OrigAverageDegree;
double expected_time_treated;
long Progression_Model;
long VL_print_time; // Instructs program to print out VL distributions every VL_print_time days.
double VL_print_lower, VL_print_upper, VL_print_interval;  // Lower and upper end of the distribution together with bin size
long VL_print_count = 0;
void PrintVLdistribution(void);

long BreakUpList1[3000001], BreakUpList2[3000001], BreakUpList3[3000001], BreakUpList4[3000001]; // These arrays track indices of partnerships to be broken up
long RevisedNumLinks[3000001];

// New lines for vaccination
int Vaccinated[3000001], ResistantToVaccine[3000001];
long LengthVaccinated[3000001];
double PercentVaccinated, VaccineDuration, PercentResistantToVaccine, RR_Vaccinated;

double MutationVariance, H, h;
long OrigN;
long Start_Condom_Campaign;
double Percent_using_condums_after_condom_campaign;
double orig_prob_sex;

void ConstantInputNewSusceptibles(void);
void VaccinatePopulation(void);
void SimulateVaccineDecay(void);

double norm_rand(double rn_mean, double rn_variance);

int main()
{
 long jcount;
 printf("Starting!!!\n");
 printf("CD4_Exp_Flag=%d\n",CD4_Exp_Flag);
 /* Initialize the population */
 OpenFiles();
 GetParameters();
 orig_prob_sex = prob_sex;
 H = Heritability;
 h = sqrt(Heritability);
 r0 = log(V_peak/V0)/t_peak;
 Dead = 0; DiedAIDS = 0; DiedNat = 0; init_genrand(random_number_seed);
 cd4InitProbs();
 CreateCD4Table();
 InitializePopulation();
 AddInfecteds();
 DefineSexualContactNetwork();
 //PrintEdgeList();
 //PrintFirstAndLastContacts();
 time = 0;
 PrintVLdistribution();


 PrintHeaders();
 //printf("About to print stats\n");
 PrintStats();
 //printf("About to enter do-while loop\n");
 fflush(stdout);


 OrigAverageDegree = AverageDegree;

 /* Time Loop */
 do {
   if (VL_print_count >= VL_print_time) {
      PrintVLdistribution();
      VL_print_count = 0;
   }
   VL_print_count++;
   UpdateViralLoadsForAllInfectedIndividuals();
   SimulateInfectionProcessOnePersonToAnother();
   ConstantInputNewSusceptibles();
   NaturalDeath();
   DeathOfAIDSPatients();
   RemoveLinksToDeadPeople();
   SimulatePartnerShipDissolution();
   SimulateNewPartnershipFormation();
   SimulateVaccineDecay();
   CheckEdgeListForImpossibleValues();  // That is, do some elementary error checking
   time = time + 1;
   if(printOutput==0){printf("time %ld\n",time);}
   if (time == Start_SexReduction_Campaign) AverageDegree = AverageDegree/Reduction_Mean_Degree;
   if ((time == 25*365) || (time == 30*365) || (time == 35*365) )
   { /* Start of vaccination campaign */
      VaccinatePopulation();   
   }
   if (time == Start_Condom_Campaign) prob_sex = orig_prob_sex * (1.0-Percent_using_condums_after_condom_campaign);
   if (time == Start_Faithfulness_Campaign) {
     MinDuration = MinDuration * Increased_Duration;
	   MaxDuration = MaxDuration * Increased_Duration;
	   for (jcount = 1; jcount <= N; jcount++) {
	     Duration[jcount] = Duration[jcount] * Increased_Duration;
	   }
   }
//   if ((time == tfinal) || (StopEarly == 1)) printf("\ntime\tN\tLinks\tUnInf\tInf\tDead\tdAIDS\tdNat\tAveLogV\t\tVstd\t\tSet_Ave\t\tSet_std\t\td_ave\t\td_std\t\tpLinks\t\tV[1]\t\td_a[1]\t\tSts[1]\n");
   PrintStats();
   fflush(stdout);
 } while ( (time < tfinal) && (StopEarly == 0)); // End time loop

 if (StopEarly == 1) {
   PrintStats();
 }

 RecordStatusOfLivingHIVPatientsForRegistryFile();
// PrintFirstAndLastContacts();
 if (Flat_Viral_Load==1) printf("Constant VL during primary infection was %lf\n",PrimaryInfectionLevel);
 if (StopEarly == 1) {
   printf("Note: Program halted early -- most likely because of difficulties in finding partners using current search algorithm.  Consider increasing MaxLinks.\n");
 }

 printf("DONE (New) !!! %c%c",7,7); // print bell sounds (ascii character 7)
} // Main


void OpenFiles()
{
 if ( (ParFile = fopen("Parameters.txt","r")) == NULL) { printf("Cannot open file 'Parameters.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (OutputFile = fopen("ViralLoadOutput.txt","w")) == NULL) { printf("Cannot open file 'ViralLoadOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (HeritOutputFile = fopen("HeritabilityOutput.txt","w")) == NULL) { printf("Cannot open file 'HeritabilityOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (PatientRegistryOutputFile = fopen("PatientRegistryOutput.txt","w")) == NULL) { printf("Cannot open file 'PatientRegistryOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (VLOutputFile = fopen("VLHistogram.txt","w")) == NULL) { printf("Cannot open file 'VLOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (spVLOutputFile = fopen("spVLHistogram.txt","w")) == NULL) { printf("Cannot open file 'spVLOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
}


void PrintVLdistribution()
{
 double current_VL;
 long i;
 long Print_Array_VL[1000];
 long Print_Array_spVL[1000];
 long bin_count = 0;

 /* Print header giving the midpoint VL (on log scale) for each element in the histogram at time 0 */
 if (time == 0) {
   fprintf(VLOutputFile,"VL\t");
   fprintf(spVLOutputFile,"spVL\t");
   for (current_VL = VL_print_lower; current_VL <= VL_print_upper; current_VL = current_VL + VL_print_interval) {
       fprintf(VLOutputFile,"%3.2lf\t", current_VL+VL_print_interval/2.0);
       fprintf(spVLOutputFile,"%3.2lf\t", current_VL+VL_print_interval/2.0);
   }
   fprintf(VLOutputFile,"\n");
   fprintf(spVLOutputFile,"\n");
 }


 /* Zero out elements of the histogram */
 bin_count = 0;
 for (current_VL = VL_print_lower; current_VL <= VL_print_upper; current_VL = current_VL + VL_print_interval) {
   bin_count++;
   Print_Array_VL[bin_count] = 0;
   Print_Array_spVL[bin_count] = 0;
 }

 /* For each infected person, increment histogram cell spanning that person's viral load */
 for (i=1; i<=N; i++) {
   if (Status[i] == 1) {
     bin_count = 0;
     for (current_VL = VL_print_lower; current_VL <= VL_print_upper; current_VL = current_VL + VL_print_interval) {
       bin_count++;
       if ( (log10(V[i]) >= current_VL) && (log10(V[i]) < current_VL+VL_print_interval)) {
         //if (time - Time_Inf[i] > 50.0) {
           Print_Array_VL[bin_count] = Print_Array_VL[bin_count] + 1;
         //}
       }
       if ( (LogSetPoint[i] >= current_VL) && (LogSetPoint[i] < current_VL+VL_print_interval)) {
         Print_Array_spVL[bin_count] = Print_Array_spVL[bin_count] + 1;
       }
     }
   }
 }

 /* Print histogram to a file */
 bin_count = 0;
 fprintf(VLOutputFile,"%ld\t", time);
 fprintf(spVLOutputFile,"%ld\t", time);
 for (current_VL = VL_print_lower; current_VL <= VL_print_upper; current_VL = current_VL + VL_print_interval) {
   bin_count++;
   fprintf(VLOutputFile,"%ld\t", Print_Array_VL[bin_count]);
   fprintf(spVLOutputFile,"%ld\t", Print_Array_spVL[bin_count]);
 }
 fprintf(VLOutputFile,"\n");
 fprintf(spVLOutputFile,"\n");


}
//double PercentVaccinated, VaccineDuration, PercentResistantToVaccine,RR_Vaccinated;

void VaccinatePopulation() {
  long i;
  for (i = 1; i <=N; i++) {
    if ((Status[i] == 0) && (Vaccinated[i]==0)) { // Questions about second condition (not vaccinating those already vaccinated)
                                                  // Could get weird if vaccine compaigns occure before vaccine wears off.
                                                  // Current algorithm will work if vaccine campaigns are longer than vaccine duration
       if (genrand_real2() < PercentVaccinated) {
         Vaccinated[i] = 1;
         LengthVaccinated[i] = 0;
       }
    }
  }

}

void SimulateVaccineDecay() {
  long i;
  for (i = 1; i <=N; i++) {
    if (Vaccinated[i]==1) {
      LengthVaccinated[i] = LengthVaccinated[i] + 1;
      if (LengthVaccinated[i] > VaccineDuration) {
        Vaccinated[i] = 0;
      }
    }
  }

}
void InitializePopulation() {
  /* Initalizes the s[], Status[], Numrecepients[], Duration[], UnderCare[], and Treated[] arrays to their default values */
  long i;
  for (i=1;i<=N;i++) {
     s[i] = prog_rate; // Same for all patients in this version.  Could be made to vary from patient to patient in future studies
     Status[i] = 0;
     CD4[i] = 0;
     CD4_TimeToAIDS[i]= 0;
     NumRecipients[i] = 0;
     Vaccinated[i] = 0;
     ResistantToVaccine[i] = 0;
     LengthVaccinated[i] = 0;
     if (genrand_real2() < 0.15) { 
       Duration[i] = MinDuration;
     } else {
       Duration[i] = MinDuration + genrand_real2()*(MaxDuration-MinDuration);
     }
     UnderCare[i] = 0; // Set to zero to start.  This will be updated in AddInfecteds() below
     Treated[i] = 0; // Assume no treatment for initial population
  }
}

void AddInfecteds()
{
 long i;
 for (i=1;i<=Infected;i++) {

   if (genrand_real2() < PercentResistantToVaccine) { 
     ResistantToVaccine[i] = 1;
   } else {
     ResistantToVaccine[i] = 0;
   }

   ViralContribToLogSP0[i] = norm_rand(AverageLogSP0, sqrt(H*VarianceLogSP0));
   
   EnvirContribToLogSP0[i] = norm_rand(0, sqrt((1-H)*VarianceLogSP0));
   
   LogSetPoint[i] = ViralContribToLogSP0[i] + EnvirContribToLogSP0[i];
   
   SpvlCat[i]=SPVL_Category(LogSetPoint[i] );
   CD4[i]=initialCD4(SpvlCat[i]);
   CD4_initial_value[i]=CD4[i];
      
   //adding exponential distribution for time in each cd4 category
   //can tighten this up once it gets working....
   CD4_TimeToAIDS_exp[i][1]=log(genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[i]][1]) ;
   CD4_TimeToAIDS_exp[i][2]=log(genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[i]][2]) ;
   CD4_TimeToAIDS_exp[i][3]=log(genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[i]][3]) ;
   CD4_TimeToAIDS_exp[i][4]=log(genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[i]][4]) ;

 if(CD4_Exp_Flag==1)
 {
    double ll=0.0;
     int kk;
    for(kk=CD4[i];kk<=3;kk++)
    {
      ll+=CD4_TimeToAIDS_exp[i][kk];
    }
    CD4_TimeToAIDS[i]=ll;
 }

   //printf("%d----%f\n",i,CD4_TimeToAIDS_exp[i][1]);
   //printf("%d -- %d -- %d\n",i,SpvlCat[i], CD4[i]);
   //printf("-- %d --", CD4[i]);

   SetPoint[i] = pow(10.0,LogSetPoint[i]);
   d_acute[i] = log(V_peak/SetPoint[i])/(t_acute - t_peak);
   Status[i] = 1;
   Generation[i] = 1;
   if (SetPoint[i] > DefinitionHighSetPoint && CD4_Determines_Treatment==0) {
	   if (genrand_real2() <  percent_under_care + Increment_Prob_Under_Care_High_spVL) {
		   UnderCare[i] = 1;
	   }
	} else {
	   if (genrand_real2() < percent_under_care && CD4_Determines_Treatment==0) {
			UnderCare[i] = 1;
	   }
   }
   if (CD4[i] > DefinitionHighSetPoint_CD4 && CD4_Determines_Treatment==1) {
     if (genrand_real2() <  percent_under_care + Increment_Prob_Under_Care_High_CD4) {
		   UnderCare[i] = 1;
	   }
	} else {
	   if (genrand_real2() < percent_under_care && CD4_Determines_Treatment==1) {
			UnderCare[i] = 1;
	   }
   }

   if (Gamma_Death == 1) {
     ExpectedDelayTime = Dmax*pow(D50,Dk)/(pow(SetPoint[i],Dk) + pow(D50,Dk));
     theta = ExpectedDelayTime/shape_parameter;
     //printf("For patient %ld: ExpectedDelayTime = %lf, shape_parameter = %lf, theta = %lf\n",i,ExpectedDelayTime,shape_parameter,theta);
     RandomTimeToAIDS[i] =  GetGammaDelay(shape_parameter,theta);
     //RandomTimeToAIDS[i] =  ExpectedDelayTime;
   } else {
     RandomTimeToAIDS[i] = 1000000000.0;
   }

    Time_Inf[i] = -t_acute - genrand_real2()*(364-t_acute);
    CD4time[i]=-Time_Inf[i] ;
    CD4_treatment_delay_index[i]=-Time_Inf[i];

   //Time_Inf[1] = 0.0;  // This and next line are special cases for debugging purposes
   //RandomTimeToAIDS[1] = 365.0;
   if (time <= Time_Inf[i] + t_acute) {
     if (Flat_Viral_Load == 1 ) {
	    V[i] = PrimaryInfectionLevel;
	 } else {
       if (time <= Time_Inf[i] + t_peak) {
         V[i] = V0*exp(r0*(time-Time_Inf[i]));
       } else {
         V[i] = V0*exp(r0*t_peak)*exp(-d_acute[i]*(time - t_peak - Time_Inf[i]));
       }
	 }
   } else {
     Vss = V0*exp(r0*t_peak)*exp(-d_acute[i]*(t_acute-t_peak));
     V[i] = Vss*exp(s[i]*(time - t_acute - Time_Inf[i])/365.0);
     if ((Gamma_Death == 1) && (Progression_Model == 2)) {// Alternative model in which viral load jumps to half its peak level at AIDS
       if ((time > Time_Inf[i] + RandomTimeToAIDS[i]) && CD4_Death!=1) {
              V[i] = 0.5*V_peak;
            }
       if (CD4[i]==4 && CD4_Death==1) {
              V[i] = 1000000; //0.5*V_peak;
       }
     }
   }
   if (V[i] < 0.0) V[i] = 0.0;

   Donors_V[i] = -1.0; // Flag absence of information with -1's or zeros.
   Donors_ViralContribToLogSP0[i] = -1.0; // Flag absence of information with -1's or zeros
   Donors_EnvirContribToLogSP0[i] = -1.0; // Flag absence of information with -1's or zeros
   Donors_LogSetPoint[i] = -1.0; // Flag absence of information with -1's or zeros
   Donors_SetPoint[i] = 0.0; // Flag absence of information with -1's or zeros
   Donors_d_acute[i] = -1.0; // Flag absence of information with -1's or zeros
   Donors_Total_Time_Inf_At_Trans[i] = -999.0; // Set unknown donors for the initial infected population to -1 year
   Donors_Generation[i] = 0;
   Donors_Index[i] = 0;
   if (Gamma_Death == 1) {
     ExpectedDelayTime = Dmax*pow(D50,Dk)/(pow(Donors_SetPoint[i],Dk) + pow(D50,Dk));
     theta = ExpectedDelayTime/shape_parameter;
     //printf("For patient %ld: ExpectedDelayTime = %lf, shape_parameter = %lf, theta = %lf\n",i,ExpectedDelayTime,shape_parameter,theta);
     if(Donors_RandomTimeToAIDS_Flag==1)
     {
       Donors_RandomTimeToAIDS[i] =  GetGammaDelay(shape_parameter, theta);
     }else{
           Donors_RandomTimeToAIDS[i] =  ExpectedDelayTime;
          }
   } else {
     Donors_RandomTimeToAIDS[i] = 1000000000.0;
   }
   //printf("%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n",
   //        time,i,Generation[i],V[i],LogSetPoint[i],ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],time-Time_Inf[i]);
 }
}


double norm_rand(double rn_mean, double rn_variance)
{
  double rand1, rand2;
  rand1 = genrand_real2(); rand2 = genrand_real2();
  return(rn_mean + rn_variance*sqrt(-2.0*log(rand1))*cos(2.0*3.14159265358979*rand2));
}

void DefineSexualContactNetwork()
{
 long i,j;
 for (i=1;i<=N;i++) {
   NumLinks[i] = 0;
   for (j=1;j<=5;j++) {
     Links[i][j] = 0;
   }
 }
 TotalEdges = (long) ( N * AverageDegree / 2.0 + 0.5);
 //Links_Broken_Per_Day = N/(30);
 //Links_Broken_Per_Day = N/(3*365);
 for (i=1;i<=TotalEdges;i++) {
    AddNewLink();
 }
}


/*void PrintEdgeList()
{
 long i,j;
 for (i=1;i<=N;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
 }
 fflush(stdout);
//printf("Type 'c' to continue : ");
//scanf("%s",junk);
} */



void PrintFirstAndLastContacts()
{
if(printOutput==1)
{
 printf("At time %ld: the first (6) and last (6) contacts are: \n",time);
 if (N > 20) {
  for (i=1;i<=6;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   if (Status[i] <0 ) printf("Person %ld : dead ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
  }
  for (i=N/2-3;i<=N/2+3;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   if (Status[i] <0 ) printf("Person %ld : dead ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
  }
  for (i=N-6;i<=N;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   if (Status[i] <0 ) printf("Person %ld : dead ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
  }
 } else {
   for (i=1;i<=N;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   if (Status[i] <0 ) printf("Person %ld : dead ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
  }

 }
}
}

void CreateCD4Table()
{
   /* Look up table for CD4 counts from Cori,Fraser,Pickles TableS13 Supplement; years spent in cat.*/

//logSPVL <= 3
   CD4_lookup[1][1]=6.08;
   CD4_lookup[1][2]=3.83;
   CD4_lookup[1][3]=3.60;
   //CD4_lookup[1][4]=1.85;
   CD4_lookup[1][4]=1.0;	//15.36 vs 14.51
 //3 < logSPVL <= 3.5
   CD4_lookup[2][1]=4.69;
   CD4_lookup[2][2]=2.12;
   CD4_lookup[2][3]=4.18;
   //CD4_lookup[2][4]=0.72;
   CD4_lookup[1][4]=1.0;	//11.71 vs 11.99
 //3.5 < logSPVL <= 4
   CD4_lookup[3][1]=3.94;
   CD4_lookup[3][2]=3.08;
   CD4_lookup[3][3]=2.38;
   //CD4_lookup[3][4]=2.16;
   CD4_lookup[3][4]=1.0;	//11.56 vs 10.4
 //4 < logSPVL <= 4.5
   CD4_lookup[4][1]=2.96;
   CD4_lookup[4][2]=2.64;
   CD4_lookup[4][3]=3.81;
   //CD4_lookup[4][4]=0.80;
   CD4_lookup[4][4]=0.50;	//10.21 vs 9.91
 //4.5 < logSPVL <= 5
   CD4_lookup[5][1]=2.25;
   CD4_lookup[5][2]=2.13;
   CD4_lookup[5][3]=3.21;
   //CD4_lookup[5][4]=0.82;
   CD4_lookup[5][4]=0.50;	//8.41 vs 8.09
 //5 < logSPVL <= 5.5
   CD4_lookup[6][1]=1.47;
   CD4_lookup[6][2]=1.44;
   CD4_lookup[6][3]=2.27;
   //CD4_lookup[6][4]=0.56;
   CD4_lookup[6][4]=0.50;	//5.74 vs 5.68
 //5.5 < logSPVL <= 6
   CD4_lookup[7][1]=0.95;
   CD4_lookup[7][2]=1.42;
   //CD4_lookup[7][3]=3.08;
   CD4_lookup[7][3]=1.00;
   //CD4_lookup[7][4]=0.13;
   CD4_lookup[7][4]=0.25;	//5.58 vs 3.62
 //6 < logSPVL < 6.5
   CD4_lookup[8][1]=0.32;
   CD4_lookup[8][2]=0.44;
   CD4_lookup[8][3]=0.68;
   //CD4_lookup[8][4]=0.07;
   CD4_lookup[8][4]=0.25;	//1.51 vs 1.69
  //logSPVL > 6.5
   CD4_lookup[9][1]=0.30;
   CD4_lookup[9][2]=0.35;
   CD4_lookup[9][3]=0.37;
   //CD4_lookup[9][4]=0.67;
   CD4_lookup[9][4]=0.25;	//1.69 vs 1.27

 }

int SPVL_Category(double SPVL_param)
{
   int SPVL_cat;
   if(SPVL_param<=3.0){SPVL_cat = 1;}
   else{if(SPVL_param<=3.5){SPVL_cat = 2;}
   else{if(SPVL_param<=4.0){SPVL_cat = 3;}
   else{if(SPVL_param<=4.5){SPVL_cat = 4;}
   else{if(SPVL_param<=5.0){SPVL_cat = 5;}
   else{if(SPVL_param<=5.5){SPVL_cat = 6;}
   else{if(SPVL_param<=6.0){SPVL_cat = 7;}
   else{if(SPVL_param<=6.5){SPVL_cat = 8;}
   else{SPVL_cat = 9;}}}}}}}}
   return(SPVL_cat);

}


void cd4InitProbs(){

 //From Pickles CD4 report supplement, table3
   cd4_initProbs[1][1]=0.0;
   cd4_initProbs[1][2]=0.12;
   cd4_initProbs[1][3]=0.88;

   cd4_initProbs[2][1]=0.01;
   cd4_initProbs[2][2]=0.12;
   cd4_initProbs[2][3]=0.87;

   cd4_initProbs[3][1]=0.03;
   cd4_initProbs[3][2]=0.12;
   cd4_initProbs[3][3]=0.85;

   cd4_initProbs[4][1]=0.03;
   cd4_initProbs[4][2]=0.19;
   cd4_initProbs[4][3]=0.78;

   cd4_initProbs[5][1]=0.05;
   cd4_initProbs[5][2]=0.21;
   cd4_initProbs[5][3]=0.73;

   cd4_initProbs[6][1]=0.04;
   cd4_initProbs[6][2]=0.25;
   cd4_initProbs[6][3]=0.71;

   cd4_initProbs[7][1]=0.09;
   cd4_initProbs[7][2]=0.27;
   cd4_initProbs[7][3]=0.64;

   cd4_initProbs[8][1]=1.0;//having all 3 probs =1, forces into box 3
   cd4_initProbs[8][2]=1.0;
   cd4_initProbs[8][3]=1.0;

   cd4_initProbs[9][1]=1.0;
   cd4_initProbs[9][2]=1.0;
   cd4_initProbs[9][3]=1.0;

}


int initialCD4(int spvl_level){
  //can only be 1,2,3 initially
  int tempCat;
  double tempProb;
  tempProb=genrand_real2();

   if(tempProb<cd4_initProbs[spvl_level][1]){tempCat=3;}
   else{if(tempProb<cd4_initProbs[spvl_level][2]){tempCat=2;}
   else{tempCat=1;}}
  // printf("%d",tempCat);
   return(tempCat);
}




void UpdateViralLoadsForAllInfectedIndividuals()
{

   long i;
   double tempTime,TimeParam;
   for (i=1;i<=N;i++) {
     if (Status[i] == 1) {
       // CD4 update -----------------------------
       if(CD4[i]<5)
       {
         TimeParam = (CD4time[i])/365.0;
         if(CD4_Exp_Flag==0)
         {

          tempTime = CD4_lookup[SpvlCat[i]][CD4[i]];
         }else{
                tempTime = CD4_TimeToAIDS_exp[i][CD4[i]];
                //printf("tempTime---%f",tempTime);
               }

         if( TimeParam > tempTime )
          {

            CD4[i]=CD4[i]+1;
            CD4time[i]=0;

           if( CD4[i] == TreatmentThreshold_CD4){
                  CD4_treatment_delay_index[i]=time;}

            if(CD4[i]==4 && CD4time[i]==0 && CD4_Exp_Flag==0)
                {
                  CD4_TimeToAIDS[i]=(long)(time-Time_Inf[i]);
                  }
          }else{CD4time[i]=CD4time[i]+1;}



       }
     //end of CD4 update -------------------------
       if (time <= Time_Inf[i] + t_acute) {
	     if (Flat_Viral_Load == 1 ) {
	       V[i] = PrimaryInfectionLevel;
	     } else {
           if (time <= Time_Inf[i] + t_peak) {
             V[i] = V0*exp(r0*(time-Time_Inf[i]));
           } else {
             V[i] = V0*exp(r0*t_peak)*exp(-d_acute[i]*(time - t_peak - Time_Inf[i]));
           }
           //if (i==1) printf("time <= Time_Inf[1] + t_acute, V[i] = %lf (Time_Inf[1] = %lf) \n",V[i],Time_Inf[i]);
		 }
       } else {
          Vss = V0*exp(r0*t_peak)*exp(-d_acute[i]*(t_acute-t_peak));
          V[i] = Vss*exp(s[i]*(time - t_acute - Time_Inf[i])/365.0);
          if ((Gamma_Death == 1) && (Progression_Model == 2)) { // Alternative model in which viral load jumps to half its peak level at AIDS
            if ((time > Time_Inf[i] + RandomTimeToAIDS[i]) && CD4_Death!=1) {
              V[i] = 1000000; //0.1*V_peak;
            }
            if (CD4[i]==4 && CD4_Death==1) {
              V[i] = 1000000; //0.1*V_peak;
            }
          }
           //if(i==1) printf("time > Time_Inf[1] + t_acute, so V[i] = %lf (time = %ld, Time_Inf[1] = %lf, t_acute = %lf, Vss = %lf) \n",
           //V[i],time,Time_Inf[1],t_acute,Vss);
       }
       if ( (time > Start_Treatment_Campaign) && ((Treated[i] == 0) && (UnderCare[i] == 1)) ) { // Note only do this if Treated == 0, so as to not allow backsliding
          //treatment based on VL (not cd4)
          if (CD4_Determines_Treatment==0 && (time - Time_Inf[i] > TreatmentTimeInPatients) && (V[i] > TreatmentThreshold)) { /* Two ways to get under treatment: by time_treated (old method) or by VL (new method) */
                Treated[i] = 1;
          }
          //treatment based on cd4 count
          if (CD4_Determines_Treatment==1 &&
              (CD4[i] >= TreatmentThreshold_CD4 ) &&
                (time - Time_Inf[i] > TreatmentTimeInPatients) &&
                   (time-CD4_treatment_delay_index[i]) >= 365) { /* Two ways to get under treatment: by time_treated (old method) or by VL (new method) */
                Treated[i] = 1;
          }
	   }
       if (Treated[i] == 1) {
          V[i] = VL_After_Treatment;
          CD4[i]= CD4_After_Treatment;
       }
       if (V[i] < 0.0) {
         V[i] = 0.0;
     }
   }
 }
}


void SimulateInfectionProcessOnePersonToAnother()
{


  long Infector, Recipient; // Potential infector and Recipient
  long i,j;
  long random_person;
  double Specific_Prob_Sex;
  for (i=1;i<=N;i++ ) {
     HadSexAlready[i] = 0; // Refers to whether person had sex today
  }
  for (random_person=1; random_person<=N; random_person++)
  {
     i = genrand_real2()*N + 1; // Select person i at random
   	 if (NumLinks[i] >= 1)
     {
       random_partner = genrand_real2()*NumLinks[i] + 1;  // Person i selects one of his or her sex partners each day.  (Assumes coital dilution)
       j = Links[i][random_partner]; // That person is person j
       if ( ((Status[i] + Status[j] == 1) && (max(Status[i],Status[j]) == 1)) && (HadSexAlready[i] + HadSexAlready[j] == 0))
       {
          if (Status[i] == 1)
          {
             Infector = i;
             Recipient = j;
           } else {
             Infector = j;
             Recipient = i;
           }
          if (Duration[i] + Duration[j] <= 0.5*365) {
             Specific_Prob_Sex = 1*prob_sex;
          } else
          {
             if (Duration[i] + Duration[j] <= 2.5*365.0) {
                Specific_Prob_Sex = 0.05*prob_sex;
             } else {
                Specific_Prob_Sex = 0.03*prob_sex;
             }
          }
          if (genrand_real2() < Specific_Prob_Sex) {
	           HadSexAlready[Infector] = 1;
	           HadSexAlready[Recipient] = 1;
             if (Inf_Rate_Function == 1) {
                if (V[Infector] > VL_After_Treatment)
                  Inf_Prob = InfRateBaseline*pow(log10(V[Infector]),InfRateExponent);
                else
                  Inf_Prob = 0.0; // power function gives negative values if V < 1.0 (Shouldn't be much of a problem in real-life)
		         } else
             {
               if (V[Infector] > VL_After_Treatment) {
                 beta = HillCoeffInfRate;
                Inf_Prob = MaxInfRate * pow(V[Infector],beta) / (pow(V[Infector],beta) + pow(VHalfMaxInfRate,beta));
			        } else {
		  		      Inf_Prob = 0.0;
		          }
            }
            // New code to have a reduced probability of transmission in vaccinated people (assuming donor's virus is sensitive)
            if ((Vaccinated[Recipient] == 1) && (ResistantToVaccine[Infector] == 0)) {
               Inf_Prob = RR_Vaccinated * Inf_Prob;
            }
            randnum = genrand_real2();
            if (randnum < Inf_Prob) {

		          Time_Inf[Recipient] = time;
	   	        V[Recipient] = V0;

			        //ViralContribToLogSP0[Recipient] = ViralContribToLogSP0[Infector] + norm_rand(0,((time-Time_Inf[Infector])/365)*MutationVariance); 
			        // Patient inherits previous patient's virus + mutational deviation scaled by time since infection
      		    ViralContribToLogSP0[Recipient] = ViralContribToLogSP0[Infector] + norm_rand(0, MutationVariance); 
      		    // Patient inherits previous patient's virus + mutational deviation

        	    EnvirContribToLogSP0[Recipient] = norm_rand(0, sqrt((1-H)*VarianceLogSP0)); 
        	   // Environmental component is independent
        	
        	   LogSetPoint[Recipient] = ViralContribToLogSP0[Recipient] + EnvirContribToLogSP0[Recipient];
        	
        	   SetPoint[Recipient] = pow(10.0,LogSetPoint[Recipient]);

             if (SetPoint[Recipient] > DefinitionHighSetPoint && CD4_Determines_Treatment==0) {
                   if (genrand_real2() <  percent_under_care + Increment_Prob_Under_Care_High_spVL) {
              		   UnderCare[Recipient] = 1;
              	   }
              	} else {
              	   if (genrand_real2() < percent_under_care && CD4_Determines_Treatment==0) {
              			UnderCare[Recipient] = 1;
              	   }
                 }
                 if (CD4[Recipient] > DefinitionHighSetPoint_CD4 && CD4_Determines_Treatment==1) {
                   if (genrand_real2() <  percent_under_care + Increment_Prob_Under_Care_High_CD4) {
              		   UnderCare[Recipient] = 1;
              	   }
              	} else {
              	   if (genrand_real2() < percent_under_care && CD4_Determines_Treatment==1) {
              			UnderCare[Recipient] = 1;
              	   }
                 }

            d_acute[Recipient] = log(V_peak/SetPoint[Recipient])/(t_acute - t_peak);
             if (Gamma_Death == 1) {
               ExpectedDelayTime = Dmax*pow(D50,Dk)/(pow(SetPoint[Recipient],Dk) + pow(D50,Dk));
	           theta = ExpectedDelayTime/shape_parameter;
	           //printf("For newly infected patient %ld: ExpectedDelayTime = %lf, shape_parameter = %lf, theta = %lf\n",Recipient,ExpectedDelayTime,shape_parameter,theta);
               RandomTimeToAIDS[Recipient] = t_acute + GetGammaDelay(shape_parameter,theta);
               //RandomTimeToAIDS[Recipient] =  ExpectedDelayTime;
             } else {
               RandomTimeToAIDS[Recipient] = 1000000000.0;
             }
             
             ResistantToVaccine[Recipient] = ResistantToVaccine[Infector]; // Recipient inherits infectors vaccinne resistance status
             
             Generation[Recipient] = Generation[Infector] + 1;
             Status[Recipient] = 1;
             SpvlCat[Recipient]=SPVL_Category(LogSetPoint[Recipient] );
             CD4[Recipient]=initialCD4(SpvlCat[Recipient]);
             CD4_initial_value[Recipient]=CD4[Recipient];
             CD4time[Recipient]=0;
             CD4_treatment_delay_index[Recipient]=time;
             //adding exponential distribution for time in each cd4 category
             //can tighten this up once it gets working....
             CD4_TimeToAIDS_exp[Recipient][1]=log( genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[Recipient]][1]) ;
             CD4_TimeToAIDS_exp[Recipient][2]=log( genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[Recipient]][2]) ;
             CD4_TimeToAIDS_exp[Recipient][3]=log( genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[Recipient]][3]) ;
             CD4_TimeToAIDS_exp[Recipient][4]=log( genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[Recipient]][4]) ;

             if(CD4_Exp_Flag==1)
             {
              double ll=0;
              int kk;
              for(kk=CD4[Recipient];kk<=3;kk++)
              {
                ll+=CD4_TimeToAIDS_exp[Recipient][kk];
              }
              CD4_TimeToAIDS[Recipient]=ll;
             }

             Donors_ViralContribToLogSP0[Recipient] = ViralContribToLogSP0[Infector];
             Donors_EnvirContribToLogSP0[Recipient] = EnvirContribToLogSP0[Infector];
             Donors_d_acute[Recipient] = d_acute[Infector];
             Donors_RandomTimeToAIDS[Recipient] = RandomTimeToAIDS[Infector];
             Donors_Total_Time_Inf_At_Trans[Recipient] = time - Time_Inf[Infector];
             Donors_V[Recipient] = V[Infector];
             Donors_Generation[Recipient] = Generation[Infector];
             Donors_SetPoint[Recipient] = SetPoint[Infector];
             Donors_LogSetPoint[Recipient] = LogSetPoint[Infector];
             Donors_Index[Recipient] = Infector;
             NumRecipients[Infector] = NumRecipients[Infector] + 1;
          } // Infection event occurred
	    } // Couple had sex
       } // Status[i] + Status[j] == 1 (i.e., opporuntunity for new infection) & they hadn't already had sex
	 } // Person i has at least one partner
   } // for i = 1 to N


}

void ConstantInputNewSusceptibles()
{
   double NewSusceptiblesPerDayDouble;
   long NewSusceptiblesPerDay;
   long i;
   newbirths = 0;
   if (MaxNReached ==0) {
     NewSusceptiblesPerDayDouble = ((double) OrigN ) * NaturalDeathRate;
     NewSusceptiblesPerDay = (long) (NewSusceptiblesPerDayDouble + 0.5);
     if (time == 10 && printOutput==1) {
        printf("At time %ld, NewSusceptiblesPerDay = %ld\n", time, NewSusceptiblesPerDay);
     }
     for (i=1;i<=NewSusceptiblesPerDay;i++) {
       if (N < MaxN) {
          N++;
          newbirths++;
          s[N] = prog_rate; // Same for all patients in this version.  Could be made to vary from patient to patient in futures studies
          Status[N] = 0;
          NumRecipients[N] = 0;
 if (genrand_real2() < 0.15) {
Duration[N] = MinDuration;
} else {
             Duration[N] = MinDuration + genrand_real2()*(MaxDuration-MinDuration);
 }
          NumLinks[N] = 0;
          for (j=1;j<=MaxLinks;j++) Links[N][j] = 0;
		  UnderCare[N] = 0;
          Treated[N] = 0;
       } else {
         if (MaxNReached ==0) {
            printf("Warning: MaxN reached.  No more births or reincarnation events allowed\n");
            MaxNReached = 1;
            tMaxNReached = time;
         } // MaxNReached
       }  // else part of N < MaxN
     } // for i = 1 to N
   } // MaxNReached == 0
}


void BirthOfNewSusceptibles()
{
   newbirths = 0;
   if ((BirthRate > 0.0) && (MaxNReached ==0)) {
     currN = N;
     // printf("At time %ld, but before births : NumLinks[100],NumLinks[360] = %ld\n",time,NumLinks[100],NumLinks[360]);
     for (i=1;i<=currN;i++) {
       if (Status[i] == 0) {
         randnum = genrand_real2();
         //printf("Testing to see if person %ld will give birth: comparing randnum (= %lf) to BirthRate (= %lf)\n",i,randnum,BirthRate);
         if (randnum < BirthRate) {
            //printf("Potential birth for person %ld, N = %ld, MaxN = %ld\n",i,N,MaxN);
           if (N < MaxN) {
             N++;
             newbirths++;
             s[N] = prog_rate; // Same for all patients in this version.  Could be made to vary from patient to patient in futures studies
             Status[N] = 0;
             NumRecipients[N] = 0;
             if (genrand_real2() < 0.15) {
 Duration[N] = MinDuration;
 } else {
               Duration[N] = MinDuration + genrand_real2()*(MaxDuration-MinDuration); }
             NumLinks[N] = 0;
             for (j=1;j<=MaxLinks;j++) Links[N][j] = 0;
		     UnderCare[N] = 0;
             Treated[N] = 0;
		   } else {
             if (MaxNReached ==0) {
               printf("Warning: MaxN reached.  No more births or reincarnation events allowed\n");
               MaxNReached = 1;
               tMaxNReached = time;
             }
           } // MaxN Reached
         }  // randnum < birthrate
       } // Status == 1
     } // for i = 1 to N
   }


}


void NaturalDeath()
{

   long i;
   for (i=1;i<=N;i++) {
     if (Status[i] >= 0) {
	    randnum = genrand_real2();
		if (randnum < NaturalDeathRate) {
           Status[i] = -1;
      fprintf(PatientRegistryOutputFile,"%ld\t%ld\t%9.6e\t%ld\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%9.6e\t%9.6e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%d\t",
                         time,i,Duration[i],NumRecipients[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],RandomTimeToAIDS[i],Donors_RandomTimeToAIDS[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i],CD4[i],CD4_TimeToAIDS[i],
                         CD4_TimeToAIDS_exp[i][1],CD4_TimeToAIDS_exp[i][2],CD4_TimeToAIDS_exp[i][3],CD4_TimeToAIDS_exp[i][4],CD4_initial_value[i]);
   fprintf(PatientRegistryOutputFile,"DiedNatural\n");
   V[i] = 0.0;
		} // Person died of natural causes
     } // Status == 1 (i.e., person is alive)
   } // i = 1 to N
}


void DeathOfAIDSPatients()
{
   long i;
   NewAIDSDeaths = 0;
    if (CD4_Death==0) {
       if (Gamma_Death == 1) {
         for (i=1;i<=N;i++) {
           if (Status[i] == 1&& Treated[i]!=1) {
             if (time - Time_Inf[i] > RandomTimeToAIDS[i] + t_acute) {  //  "+ t_acute" for AIDS length
               if ((time-Time_Inf[i]) > t_acute) { // Only allow death in patients not acutely infected
     fprintf(PatientRegistryOutputFile,"%ld\t%ld\t%9.6e\t%ld\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%9.6e\t%9.6e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%d\t",
                         time,i,Duration[i],NumRecipients[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],RandomTimeToAIDS[i],Donors_RandomTimeToAIDS[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i],CD4[i],CD4_TimeToAIDS[i],
                         CD4_TimeToAIDS_exp[i][1],CD4_TimeToAIDS_exp[i][2],CD4_TimeToAIDS_exp[i][3],CD4_TimeToAIDS_exp[i][4],CD4_initial_value[i]);
                  V[i] = 0.0;
                 Status[i] = -2;
                 NewAIDSDeaths++;
                 fprintf(PatientRegistryOutputFile,"DiedAIDS\n");
               } // endif not acutely infected
             } // endif random AIDSTime
           } // endif HIV+ HIV+? not HIV-?
         } // endfor i = 1 to N
       } else {
       /* Gamma_Death not 1, means viral load dependent daily probability of dying */
         for (i=1;i<=N;i++) {
           if (Status[i] == 1 && Treated[i]!=1) {
             randnum = genrand_real2();
             if (V[i] > 0.0) {
               yearly_rate = death_rate_constant*jmpow(log10(V[i]),death_rate_exponent);
               // yearly_rate = 0.000003*jmpow(log10(V[i]),6.45); // CASCADE viral load and death function?
             } else {
               yearly_rate = 0.0;
             }
             daily_rate = yearly_rate/365.0;
             if (randnum < daily_rate) {
               if ((time-Time_Inf[i]) > t_acute) { // Only allow death in patients not acutely infected
     fprintf(PatientRegistryOutputFile,"%ld\t%ld\t%9.6e\t%ld\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%9.6e\t%9.6e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%d\t",
                         time,i,Duration[i],NumRecipients[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],RandomTimeToAIDS[i],Donors_RandomTimeToAIDS[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i],CD4[i],CD4_TimeToAIDS[i],
                         CD4_TimeToAIDS_exp[i][1],CD4_TimeToAIDS_exp[i][2],CD4_TimeToAIDS_exp[i][3],CD4_TimeToAIDS_exp[i][4],CD4_initial_value[i]);
                  V[i] = 0.0;
                 Status[i] = -2;
                 NewAIDSDeaths++;
                 fprintf(PatientRegistryOutputFile,"DiedAIDS\n");
               } // endif not acutely infected
             } // endif CASCADE condition
           } // status == 1
         } // endfor i = 1 to N
       } // Death functions
    } // CD4_Death == 0
    if (CD4_Death==1) {
         for (i=1;i<=N;i++) {
            if (CD4[i] == 5 && Status[i]==1 && Treated[i]!=1) {
     fprintf(PatientRegistryOutputFile,"%ld\t%ld\t%9.6e\t%ld\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%9.6e\t%9.6e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%d\t",
                         time,i,Duration[i],NumRecipients[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],RandomTimeToAIDS[i],Donors_RandomTimeToAIDS[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i],CD4[i],CD4_TimeToAIDS[i],
                         CD4_TimeToAIDS_exp[i][1],CD4_TimeToAIDS_exp[i][2],CD4_TimeToAIDS_exp[i][3],CD4_TimeToAIDS_exp[i][4],CD4_initial_value[i]);
                 V[i] = 0.0;
                  Status[i] = -2;
                  NewAIDSDeaths++;
                  fprintf(PatientRegistryOutputFile,"DiedAIDS\n");
            }
         }
    }//CD4_Death==1
}


void RemoveLinksToDeadPeople()
{
   long i, j;
   /* Remove links to dead people */
   for (i = 1; i <= N; i++) {
      if ((Status[i] < 0) && (NumLinks[i] > 0)) {
        LinksToRemove = NumLinks[i];
        for (j = 1; j <= LinksToRemove; j++) {
           Widowed_partner = Links[i][j];
            /* Need to find which of the widowed partner's links points back to i */
           Link_back = 0;
           for (k = 1; k <= NumLinks[Widowed_partner]; k++) {
               if (Links[Widowed_partner][k] == i) {
                Link_back = k;
                //printf(" --- It Does!!! (so setting Link_back = %ld)",Link_back);
              }
              //printf("\n");
           }
           if (Link_back == 0) {
             printf("Trouble with unwinding links from the widow (%ld) of dead person %ld (NumLinks[%ld] = %ld: \n",Widowed_partner,i,i,NumLinks[i]);
             printf("   Couldn't find any link backs from the Widow's list!\n"); fflush(stdout); exit(1);
           }
           Links[Widowed_partner][Link_back] = Links[Widowed_partner][NumLinks[Widowed_partner]]; // Replace link to deceased with the last person on widowed persons contact list
           Links[Widowed_partner][NumLinks[Widowed_partner]] = 0; // Remove last person (as to not to count the last person twice)

           NumLinks[Widowed_partner] = NumLinks[Widowed_partner] - 1;

        } // j = 1 to NumLinks[i]
        NumLinks[i] = 0;
      }
   }
}


void SimulatePartnerShipDissolution()
{
  //printf("Inside SimulatePartnershipDissolution:\n"); fflush(stdout);
  long person1, person2, relationship, j;
  long counter;
  breakups = 0;
  //printf("Inside SimulatePartnershipDissolution:\n"); fflush(stdout);
  for (person1 = 1; person1 <= N; person1++) {
	//printf("Testing if person %ld breaks up with anyone\n",person1); fflush(stdout);
    for (index_to_person2 = 1; index_to_person2 <= NumLinks[person1]; index_to_person2++) {
	   person2 = Links[person1][index_to_person2];
	   //printf("  Testing if person %ld breaks up with %ld\n",person1,person2); fflush(stdout);

       if (person2 >= person1) { // If person2 < person1, this relationship has already been tested.
		 // Test if persons 1 and 2 break up
         PartnershipDuration = 2.0/(1.0/Duration[person1] + 1.0/Duration[person2]); // This version assumes a harmonic mean
         ProbPartnersBreakUp = 1.0/PartnershipDuration; // Approximate method for determining probability of breakup
         if ((genrand_real2() < ProbPartnersBreakUp) || ((Status[person1] < 0) || (Status[person2] < 0)) ) { // Partnership broke up
		   //if (time <= 3) printf("Targeting the relationship between %ld and %ld ", person1,person2);
	       //if (time <= 3) printf("(index %ld from %ld's %ld sexual contacts)", index_to_person2, person1, NumLinks[person1]);

           /* Find the index in person2's list that points back to person 1 */
           found_link = 0;
           for (j=1;j<=NumLinks[person2];j++) {
             if (Links[person2][j] == person1) {
               found_link = 1;
               index_to_person1 = j;
             }
           }
	       //if (time <= 3 && printOutput==1) printf(" (index %ld from %ld's %ld sexual contacts) for breakup\n", index_to_person1, person2, NumLinks[person2]);

           if (found_link == 0 && printOutput==1) {
             printf("Big trouble in breakups : couldn't find the corresponding link for person1 (%ld) in person2's (%ld) contact list\n",person1,person2);
             printf("   NumLinks[%ld] = %ld, Links[%ld][1] = %ld, NumLinks[%ld] = %ld, Links[%ld][1] = %ld, Links[%ld][2] = %ld, Status[%ld] = %d, Status[%ld] = %d\n",
                person1,NumLinks[person1],person1,Links[person1][1],person2,NumLinks[person2],person2,Links[person2][1],person2,Links[person2][2],person1,Status[person1],person2,Status[person2]);
             fflush(stdout); exit(1);
           }
           breakups++;

		   // Add this couple to the list of who broke up with whom.  But don't remove this link now.  (Doing so would change the upper limit of the for loop above)
		   if (breakups < 3000000) {
			 BreakUpList1[breakups] = person1;
		     BreakUpList2[breakups] = index_to_person2;
			 BreakUpList3[breakups] = person2;
		     BreakUpList4[breakups] = index_to_person1;
   	        } else {
		 	  printf("Error: BreakUpList1 and BreakupList2 arrays are not large enough to track of all of the breakups.  Type ctrl-C' to exit\n");  StopEarly = 1;
	       }
		 } // Persons 1 and 2 broke up
	   } // Relationship not already tested
	 } // partners of person1
   } // person1 from 1 to N

  // Flag links for deletion
  if (time <= 3 && printOutput==1) printf("Removing %ld partnerships from the list of Linked persons\n",breakups);   fflush(stdout);

  for (relationship = 1; relationship <= breakups; relationship++) {
	 person1 = BreakUpList1[relationship];
	 index_to_person2 = BreakUpList2[relationship];
	 person2 = BreakUpList3[relationship];
	 index_to_person1 = BreakUpList4[relationship];

	 /* Remove the link from person1 */
      Links[person1][index_to_person2] = -1; // Flag this link for deletion

	  /* Remove the corresponding link index from person 2 */
      Links[person2][index_to_person1] = -1; // Flag this link for deletion
  }

  // Skip over any links flagged for deletion.
  // For example, if person1 was originally linked to 1,2,3,4, but broke up with 1 and 3, the new list would be 2, 4, 3, 4 (with last two being orphans)
  //
   for (person1 = 1; person1 <= N; person1++) {
     RevisedNumLinks[person1] = 0;
	 counter = 0;
  	 for (index_to_person2 = 1; index_to_person2 <= NumLinks[person1]; index_to_person2++) {
		if (Links[person1][index_to_person2] > 0) {
		   counter++;
		   Links[person1][counter] = Links[person1][index_to_person2];
		}
		RevisedNumLinks[person1] = counter;
	 }
	 if (Status[person1] < 0) RevisedNumLinks[person1] = 0;
  }

  // Zero out any expired links (e.g., getting rid of 3 and 4 above), then update the number of links for each person
  for (person1 = 1; person1 <= N; person1++) {
      // Zero out any inactive links
      for (index_to_person2=RevisedNumLinks[person1]+1; index_to_person2<=NumLinks[person1]; index_to_person2++) {
		 Links[person1][index_to_person2] = 0;
	  }
	  NumLinks[person1] = RevisedNumLinks[person1];
   }
}


void SimulateNewPartnershipFormation()
{
  long i;
    /* Simulate process by which singles, widows, widowers, and jilted lovers look for new contacts */
   Active = 0.0;
   NLinks = 0.0;
   for (i=1;i<=N;i++) {
      if (Status[i] >= 0) {
	    Active = Active + 1.0;
        NLinks = NLinks + (double) NumLinks[i];
	  }
   }
   ExpectedLinks = (long) ( Active* AverageDegree/2.0 + 0.5);

   dpLinksToAdd = ExpectedLinks - NLinks/2.0;
   LinksToAdd = (long) (dpLinksToAdd + 0.5);
   //printf("About to add %ld new links (NLinks = %lf, Active = %lf, ExpectedLinks = %lf, dpLinksToAdd = %lf, newbirths = %ld, NewAIDSDeaths = %ld, Breakups = %ld)\n",
   //                   LinksToAdd,       NLinks,       Active,       ExpectedLinks,       dpLinksToAdd,       newbirths,      NewAIDSDeaths,       breakups);

   for (i=1; i<=LinksToAdd; i++ ) {
	 AddNewLink();
   }
}

void CheckEdgeListForImpossibleValues()
{
   /* Do some elementary error checking before moving onto the next day*/
   for (i=1;i<=N;i++) {
     if ((NumLinks[i] < 0) || (NumLinks[i] > MaxLinks) ) {
       printf("Trouble: NumLinks[%ld] = %ld\n",i,NumLinks[i]);
       fflush(stdout); exit(1);
     }
     for (j=1; j<=NumLinks[i];j++) {
        if ((Links[i][j] <= 0) || (Links[i][j] > N)) {
          printf("Trouble at end of day %ld: Links[%ld][%ld]= %ld (NumLinks[%ld] = %ld)\n",time,i,j,Links[i][j],i,NumLinks[i]);
          fflush(stdout); exit(1);
        }
     }
   }
}

void RecordStatusOfLivingHIVPatientsForRegistryFile()
{

 for (i=1;i<=N;i++) {
    if (Status[i] == 1) {
   fprintf(PatientRegistryOutputFile,"%ld\t%ld\t%9.6e\t%ld\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%9.6e\t%9.6e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%d\tAlive(HIV+)\n",
                         time,i,Duration[i],NumRecipients[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],RandomTimeToAIDS[i],Donors_RandomTimeToAIDS[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i],CD4[i],CD4_TimeToAIDS[i],
                         CD4_TimeToAIDS_exp[i][1],CD4_TimeToAIDS_exp[i][2],CD4_TimeToAIDS_exp[i][3],CD4_TimeToAIDS_exp[i][4],CD4_initial_value[i]);

  }
 }

}


long max(long a, long b)
{
  if (a > b) return (a);
       else  return (b);
}


long GetPoisson(long num, double pois_prob)
{
   long numfound, i;
   numfound = 0;
   for (i = 1; i <= num; i++) {
     if (genrand_real2() < pois_prob) {
       numfound++;
     }
   }
   return(numfound);
}

void AddNewLink(void)
{ /* Note: While this is a clunky algorithm, it seems to work */
  long person1, person2; // Indices used to define pairings when building the sexual network
  long j, AcceptsNewPartner;
  double Prob_Same_Duration_Group;
   /* Select the first person */
  escape_counter = 0;
  do {
    person1 = genrand_real2() * N + 1;
    escape_counter++;
	AcceptsNewPartner = 0;
	if (genrand_real2() < exp(-NumLinks[person1])) {  // Probability of entering a new partnership decreases expontentially with number of existing partners; P(0) = 1, P(1) = 0.36, P(2) = 0.13, P(3) = 0.05,...
	  AcceptsNewPartner = 1;
	}

    if (escape_counter == 2897) {
      printf("Warning: Having difficulty finding links (escape_counter = %ld).  Print stats on next few attempts before giving up\n",escape_counter);
    }
    if ((escape_counter >= 2897) && (escape_counter < 2900)) {
      printf("  Testing whether person %ld (Status = %d, NumLinks = %ld, MaxLinks = %ld, Links[%ld][1] = %ld) is a suitable new link\n",person1,Status[person1],NumLinks[person1],MaxLinks,person1,Links[person1][1]);
    }
    if (escape_counter >= 2900) {
      printf("Warning: Gave up trying to find new links (escape_counter = %ld).  Type Ctrl-C to quit : \n",escape_counter);
      failed_to_find_a_new_link = 1;
    }
  } while ( ( ((NumLinks[person1] >= MaxLinks) || (AcceptsNewPartner == 0)) || (Status[person1] <0) ) && (failed_to_find_a_new_link == 0) );
  // Logic is to go back and loop through with a new person1 if the current person1 is unsuitable (e.g., person1 is dead or already has the maximum number of partners)
  //   However, we exit this while loop if there were many, many failed attempts to find a suitable person (in which case we set failed_to_find_a_new_link = 1)

  /* Select the second person */
  escape_counter = 0;
  do {
    AlreadyLinked = 0;
    person2 = genrand_real2() * N + 1;
    for (j=1; j<=NumLinks[person2]; j++) {
       if (Links[person2][j] == person1) AlreadyLinked = 1;
    }
	if (person2 == person1) AlreadyLinked = 1; // If the same person, count this as already linked for accounting purposes

	AcceptsNewPartner = 0;
	if (genrand_real2() < exp(-NumLinks[person2])) {  // Probability of entering a new partnership decreases exponentially with number of existing partners; P(0) = 1, P(1) = 0.36, P(2) = 0.13, P(3) = 0.05,...
        /*printf("Duration[%ld] = %lf\n",person1,Duration[person1]);
        printf("Duration[%ld] = %lf\n",person2,Duration[person2]);
        printf("MinDuration = %lf\n",MinDuration);
        printf("MaxDuration = %lf\n",MaxDuration);
        printf("log10(Duration[%ld]) = %lf\n",person1,log10(Duration[person1]));
        printf("log10(Duration[%ld]) = %lf\n",person2,log10(Duration[person2]));
        printf("log10(MinDuration) = %lf\n",log10(MinDuration));
        printf("log10(MaxDuration) = %lf\n",log10(MaxDuration)); */

        if (MaxDuration == MinDuration) {
            Prob_Same_Duration_Group = 1.0;
        } else {
           Prob_Same_Duration_Group = 1.0 - fabs(log10(Duration[person1]+epsilon)-log10(Duration[person2]+epsilon)) / (log10(MaxDuration+epsilon)-log10(MinDuration+epsilon));
        }
        // printf("Prob_Same_Duration_Group = %lf\n",Prob_Same_Duration_Group);

        if (genrand_real2() < Prob_Same_Duration_Group ) { // Persons less likely to form partnerships with persons with different durations
                                          // epsilon is a safey in case one of the Durations is zero. (Cannot take log10 of zero!)
 	      AcceptsNewPartner = 1;
          //printf("Partnership formed between %ld and %ld\n",person1,person2);
        }
	}

    escape_counter++;
    if (escape_counter >= 2900) {
      printf("Warning: Unable to find links (escape_counter = %ld)\nType Ctrl-C to quit : ",escape_counter);
      failed_to_find_a_new_link = 1;
     }
  } while ( (((NumLinks[person2] >= MaxLinks) || (AcceptsNewPartner == 0)) || ((Status[person2] < 0) || (AlreadyLinked == 1)) ) && (failed_to_find_a_new_link == 0) );

  if (failed_to_find_a_new_link == 0) {
    NumLinks[person1] = NumLinks[person1] + 1;
    NumLinks[person2] = NumLinks[person2] + 1;
    Links[person1][NumLinks[person1]] = person2;
    Links[person2][NumLinks[person2]] = person1;
    if (time <= -3 && printOutput==1) printf("New link established between persons %ld and %ld\n",person1,person2);
  }
  if (failed_to_find_a_new_link == 1) {
     StopEarly = 1;
  }

}

/* The following function calculates a^b (i.e., "a" raised to the power of "b"). */
double jmpow(double a, double b) {
  return(exp(b*log(a)));
}


/*************** Function GetGammaDelay ******************************/
/* This function returns a random number from the gamma distribution */
/* with a scale parameter of theta and a shape parameter of  k       */
/* Inputs:						  	     */
/*          theta -- the shape parameter			     */
/*          k     -- the scale parameter			     */
/* Ouput: A random number (delay time) 				     */
/*        with mean k*theta and variance k*theta^2   		     */
/* Dependencies: This function calls the MT random number generator  */
/*********************************************************************/

double GetGammaDelay(double k, double theta)
{
  double lnUsum, Final_Value;
  double V3m2, V3m1, V3m;
  double k_int, delta, epsilon, epsilon_m, nu_m, v0;
  long i,m,step, condition_met;
  step = 1;
  condition_met = 0;
  k_int = (int) k;
  delta = k - k_int;
  v0 = exp(1.0)/(exp(1.0) + delta);
  do {
    m = 1;
    V3m = genrand_real2();
    V3m1 = genrand_real2();
    V3m2 = genrand_real2();
    //printf("At step %ld: V3m-2 = %lf, V3m-1 = %lf, V3m = %lf\n",step,V3m2,V3m1,V3m);
    if (V3m2 <= v0) {
      epsilon_m = pow(V3m1,1/delta);
      nu_m = V3m*pow(epsilon_m,delta-1);
    } else {
      epsilon_m = 1.0 - log(V3m1);
      nu_m = V3m*exp(-epsilon_m);
    }
    if (nu_m > pow(epsilon_m,delta-1)*exp(-epsilon_m)) {
      m = m + 1;
    } else {
      condition_met = 1;
    }
    step++;
  } while (condition_met == 0);
  epsilon = epsilon_m;
  lnUsum = 0.0;
  for (i = 1; i<= k_int; i++) {
    lnUsum = lnUsum + log(genrand_real2());
  }
  Final_Value = theta*(epsilon - lnUsum);
  return(Final_Value);
}


void GetParameters(void)
{
 double Nfloat;
 fscanf(ParFile,"%lf %s", &Nfloat, descript);
   N = (long) (Nfloat + 0.5);
   fprintf(OutputFile,"# %ld %s\n",N, descript);
   fprintf(PatientRegistryOutputFile,"# %ld %s\n",N, descript);
   fprintf(HeritOutputFile,"# %ld %s\n",N, descript);
   if (N > MaxN) {
     printf("Error: N (%ld) exceeds MaxN (%ld).  Hit ctrl-C to quit!\n",N, MaxN);
     scanf("%s",junk);
   }
   OrigN = N;
 fscanf(ParFile,"%ld %s",&Infected, descript);
   fprintf(OutputFile,"# %ld %s\n",Infected, descript);
   fprintf(PatientRegistryOutputFile,"# %ld %s\n",Infected, descript);
   fprintf(HeritOutputFile,"# %ld %s\n",Infected, descript);

// New lines for vaccine module 
 fscanf(ParFile,"%lf %s",&PercentVaccinated, descript); // New code for vaccination campaign
   fprintf(OutputFile,"# %lf %s\n",PercentVaccinated, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",PercentVaccinated, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",PercentVaccinated, descript);

 fscanf(ParFile,"%lf %s",&RR_Vaccinated, descript); // New code for vaccination campaign
   fprintf(OutputFile,"# %lf %s\n",RR_Vaccinated, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",RR_Vaccinated, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",RR_Vaccinated, descript);

 fscanf(ParFile,"%lf %s",&VaccineDuration, descript); // New code for vaccination campaign
   fprintf(OutputFile,"# %lf %s\n",VaccineDuration, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",VaccineDuration, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",VaccineDuration, descript);

 fscanf(ParFile,"%lf %s",&PercentResistantToVaccine, descript); // New code for vaccination campaign
   fprintf(OutputFile,"# %lf %s\n",PercentResistantToVaccine, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",PercentResistantToVaccine, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",PercentResistantToVaccine, descript);



 fscanf(ParFile,"%lf %s",&AverageDegree, descript);
   fprintf(OutputFile,"# %lf %s\n",AverageDegree, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",AverageDegree, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",AverageDegree, descript);

// fscanf(ParFile,"%lf %s",&Percent_Links_Broken_Per_Day, descript);
//   fprintf(OutputFile,"# %lf %s\n",Percent_Links_Broken_Per_Day, descript);
//   fprintf(PatientRegistryOutputFile,"# %lf %s\n",Percent_Links_Broken_Per_Day, descript);
//   fprintf(HeritOutputFile,"# %lf %s\n",Percent_Links_Broken_Per_Day, descript);

 fscanf(ParFile,"%lf %s",&MinDuration, descript);
   fprintf(OutputFile,"# %lf %s\n",MinDuration, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",MinDuration, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",MinDuration, descript);

 fscanf(ParFile,"%lf %s",&MaxDuration, descript);
	fprintf(OutputFile,"# %lf %s\n",MaxDuration, descript);
	fprintf(PatientRegistryOutputFile,"# %lf %s\n",MaxDuration, descript);
	fprintf(HeritOutputFile,"# %lf %s\n",MaxDuration, descript);

 fscanf(ParFile,"%lf %s",&prob_sex, descript);
	fprintf(OutputFile,"# %lf %s\n",prob_sex, descript);
	fprintf(PatientRegistryOutputFile,"# %lf %s\n",prob_sex, descript);
	fprintf(HeritOutputFile,"# %lf %s\n",prob_sex, descript);

 fscanf(ParFile,"%ld %s",&MaxLinks, descript);
   fprintf(OutputFile,"# %ld %s\n",MaxLinks, descript);
   fprintf(PatientRegistryOutputFile,"# %ld %s\n",MaxLinks, descript);
   fprintf(HeritOutputFile,"# %ld %s\n",MaxLinks, descript);

 fscanf(ParFile,"%lf %s", &Start_Treatment_Campaign, descript);
   fprintf(OutputFile,"# %lf %s\n", Start_Treatment_Campaign, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n", Start_Treatment_Campaign, descript);
   fprintf(HeritOutputFile,"# %lf %s\n", Start_Treatment_Campaign, descript);

 fscanf(ParFile,"%lf %s", &percent_under_care, descript);
   fprintf(OutputFile,"# %lf %s\n", percent_under_care, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n", percent_under_care, descript);
   fprintf(HeritOutputFile,"# %lf %s\n", percent_under_care, descript);

 fscanf(ParFile,"%lf %s", &Increment_Prob_Under_Care_High_spVL, descript);
	fprintf(OutputFile,"# %lf %s\n", Increment_Prob_Under_Care_High_spVL, descript);
	fprintf(PatientRegistryOutputFile,"# %lf %s\n", Increment_Prob_Under_Care_High_spVL, descript);
	fprintf(HeritOutputFile,"# %lf %s\n", Increment_Prob_Under_Care_High_spVL, descript);

 fscanf(ParFile,"%lf %s", &DefinitionHighSetPoint, descript);
	fprintf(OutputFile,"# %lf %s\n", DefinitionHighSetPoint, descript);
	fprintf(PatientRegistryOutputFile,"# %lf %s\n", DefinitionHighSetPoint, descript);
	fprintf(HeritOutputFile,"# %lf %s\n", DefinitionHighSetPoint, descript);

 fscanf(ParFile,"%lf %s", &TreatmentTimeInPatients, descript);
   fprintf(OutputFile,"# %lf %s\n", TreatmentTimeInPatients, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n", TreatmentTimeInPatients, descript);
   fprintf(HeritOutputFile,"# %lf %s\n", TreatmentTimeInPatients, descript);

 fscanf(ParFile,"%lf %s", &TreatmentThreshold, descript);
   fprintf(OutputFile,"# %lf %s\n", TreatmentThreshold, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n", TreatmentThreshold, descript);
   fprintf(HeritOutputFile,"# %lf %s\n", TreatmentThreshold, descript);

 fscanf(ParFile,"%lf %s", &VL_After_Treatment, descript);
   fprintf(OutputFile,"# %lf %s\n", VL_After_Treatment, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n", VL_After_Treatment, descript);
   fprintf(HeritOutputFile,"# %lf %s\n", VL_After_Treatment, descript);

 fscanf(ParFile,"%d %s", &CD4_Determines_Treatment, descript);
   fprintf(OutputFile,"# %d %s\n", CD4_Determines_Treatment, descript);
   fprintf(PatientRegistryOutputFile,"# %d %s\n", CD4_Determines_Treatment, descript);
   fprintf(HeritOutputFile,"# %d %s\n", CD4_Determines_Treatment, descript);

 fscanf(ParFile,"%d %s", &TreatmentThreshold_CD4, descript);
   fprintf(OutputFile,"# %d %s\n",TreatmentThreshold_CD4, descript);
   fprintf(PatientRegistryOutputFile,"# %d %s\n", TreatmentThreshold_CD4, descript);
   fprintf(HeritOutputFile,"# %d %s\n", TreatmentThreshold_CD4, descript);

 fscanf(ParFile,"%d %s", &CD4_After_Treatment, descript);
   fprintf(OutputFile,"# %d %s\n", CD4_After_Treatment, descript);
   fprintf(PatientRegistryOutputFile,"# %d %s\n", CD4_After_Treatment, descript);
   fprintf(HeritOutputFile,"# %d %s\n", CD4_After_Treatment, descript);

fscanf(ParFile,"%lf %s", &Start_SexReduction_Campaign, descript);
   fprintf(OutputFile,"# %lf %s\n", Start_SexReduction_Campaign, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n", Start_SexReduction_Campaign, descript);
   fprintf(HeritOutputFile,"# %lf %s\n", Start_SexReduction_Campaign, descript);

 fscanf(ParFile,"%lf %s", &Reduction_Mean_Degree, descript);
   fprintf(OutputFile,"# %lf %s\n", Reduction_Mean_Degree, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n", Reduction_Mean_Degree, descript);
   fprintf(HeritOutputFile,"# %lf %s\n", Reduction_Mean_Degree, descript);

  fscanf(ParFile,"%ld %s", &Start_Condom_Campaign, descript);
    fprintf(OutputFile,"# %ld %s\n", Start_Condom_Campaign, descript);
    fprintf(PatientRegistryOutputFile,"# %ld %s\n", Start_Condom_Campaign, descript);
    fprintf(HeritOutputFile,"# %ld %s\n", Start_Condom_Campaign, descript);

  fscanf(ParFile,"%lf %s", &Percent_using_condums_after_condom_campaign, descript);
    fprintf(OutputFile,"# %lf %s\n", Percent_using_condums_after_condom_campaign, descript);
    fprintf(PatientRegistryOutputFile,"# %lf %s\n", Percent_using_condums_after_condom_campaign, descript);
    fprintf(HeritOutputFile,"# %lf %s\n", Percent_using_condums_after_condom_campaign, descript);

  fscanf(ParFile,"%lf %s", &Start_Faithfulness_Campaign, descript);
   fprintf(OutputFile,"# %lf %s\n", Start_Faithfulness_Campaign, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n", Start_Faithfulness_Campaign, descript);
   fprintf(HeritOutputFile,"# %lf %s\n", Start_Faithfulness_Campaign, descript);

 fscanf(ParFile,"%lf %s", &Increased_Duration, descript);
   fprintf(OutputFile,"# %lf %s\n", Increased_Duration, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n", Increased_Duration, descript);
   fprintf(HeritOutputFile,"# %lf %s\n", Increased_Duration, descript);

 fscanf(ParFile,"%ld %s", &Inf_Rate_Function, descript);
   fprintf(OutputFile,"# %ld %s\n", Inf_Rate_Function, descript);
   fprintf(PatientRegistryOutputFile,"# %ld %s\n", Inf_Rate_Function, descript);
   fprintf(HeritOutputFile,"# %ld %s\n", Inf_Rate_Function, descript);

 fscanf(ParFile,"%lf %s",&InfRateBaseline, descript);
   fprintf(OutputFile,"# %4.2e %s\n",InfRateBaseline, descript);
   fprintf(PatientRegistryOutputFile,"# %4.2e %s\n",InfRateBaseline, descript);
   fprintf(HeritOutputFile,"# %4.2e %s\n",InfRateBaseline, descript);

 fscanf(ParFile,"%lf %s",&InfRateExponent, descript);
   fprintf(OutputFile,"# %lf %s\n",InfRateExponent, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",InfRateExponent, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",InfRateExponent, descript);

 fscanf(ParFile,"%lf %s",&MaxInfRate, descript);
   fprintf(OutputFile,"# %lf %s\n",MaxInfRate, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",MaxInfRate, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",MaxInfRate, descript);

 fscanf(ParFile,"%lf %s",&VHalfMaxInfRate, descript);
   fprintf(OutputFile,"# %lf %s\n",VHalfMaxInfRate, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",VHalfMaxInfRate, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",VHalfMaxInfRate, descript);

 fscanf(ParFile,"%lf %s",&HillCoeffInfRate, descript);
   fprintf(OutputFile,"# %lf %s\n",HillCoeffInfRate, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",HillCoeffInfRate, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",HillCoeffInfRate, descript);

 fscanf(ParFile,"%lf %s",&BirthRate, descript);
   fprintf(OutputFile,"# %lf %s\n",BirthRate, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",BirthRate, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",BirthRate, descript);

 fscanf(ParFile,"%lf %s",&NaturalDeathRate, descript);
   fprintf(OutputFile,"# %lf %s\n",NaturalDeathRate, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",NaturalDeathRate, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",NaturalDeathRate, descript);

 fscanf(ParFile,"%ld %s",&Progression_Model, descript);
	fprintf(OutputFile,"# %ld %s\n",Progression_Model, descript);
	fprintf(PatientRegistryOutputFile,"# %ld %s\n",Progression_Model, descript);
	fprintf(HeritOutputFile,"# %ld %s\n",Progression_Model, descript);

 fscanf(ParFile,"%ld %s",&CD4_Death, descript);
  fprintf(OutputFile,"# %ld CD4_Death (%s)\n",CD4_Death, descript);
  fprintf(PatientRegistryOutputFile,"# %ld CD4_Death (%s)\n",CD4_Death, descript);
	fprintf(HeritOutputFile,"# %ld CD4_Death (%s)\n",CD4_Death, descript);

 fscanf(ParFile,"%ld %s",&Gamma_Death, descript);
  fprintf(OutputFile,"# %ld Gamma_Death (%s)\n",Gamma_Death, descript);
	fprintf(PatientRegistryOutputFile,"# %ld Gamma_Death (%s)\n",Gamma_Death, descript);
	fprintf(HeritOutputFile,"# %ld Gamma_Death (%s)\n",Gamma_Death, descript);

 fscanf(ParFile,"%lf %s",&shape_parameter, descript);
	fprintf(OutputFile,"# %lf shape_parameter (%s)\n",shape_parameter, descript);
	fprintf(PatientRegistryOutputFile,"# %lf shape_parameter(%s)\n",shape_parameter, descript);
	fprintf(HeritOutputFile,"# %lf shape_parameter (%s)\n",shape_parameter, descript);

 fscanf(ParFile,"%lf %s",&death_rate_constant, descript);
   fprintf(OutputFile,"# %lf death_rate_constant (%s)\n",death_rate_constant, descript);
   fprintf(PatientRegistryOutputFile,"# %lf death_rate_constant (%s)\n",death_rate_constant, descript);
   fprintf(HeritOutputFile,"# %lf death_rate_constant (%s)\n",death_rate_constant, descript);

 fscanf(ParFile,"%lf %s",&death_rate_exponent, descript);
   fprintf(OutputFile,"# %lf death_rate_exponent (%s)\n",D50, descript);
   fprintf(PatientRegistryOutputFile,"# %lf death_rate_exponent (%s)\n",D50, descript);
   fprintf(HeritOutputFile,"# %lf death_rate_exponent (%s)\n",D50, descript);

 fscanf(ParFile,"%lf %s",&Dmax, descript);
   fprintf(OutputFile,"# %lf Dmax (%s)\n",Dmax, descript);
   fprintf(PatientRegistryOutputFile,"# %lf Dmax (%s)\n",Dmax, descript);
   fprintf(HeritOutputFile,"# %lf Dmax (%s)\n",Dmax, descript);

 fscanf(ParFile,"%lf %s",&D50, descript);
   fprintf(OutputFile,"# %lf D50 (%s)\n",D50, descript);
   fprintf(PatientRegistryOutputFile,"# %lf D50 (%s)\n",D50, descript);
   fprintf(HeritOutputFile,"# %lf D50 (%s)\n",D50, descript);

 fscanf(ParFile,"%lf %s",&Dk, descript);
   fprintf(OutputFile,"# %lf Dk (%s)\n",Dk, descript);
   fprintf(PatientRegistryOutputFile,"# %lf Dk (%s)\n",Dk, descript);
   fprintf(HeritOutputFile,"# %lf Dk (%s)\n",Dk, descript);

 fscanf(ParFile,"%lf %s",&V0, descript);
   fprintf(OutputFile,"# %lf %s\n",V0, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",V0, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",V0, descript);

 fscanf(ParFile,"%lf %s",&V_peak, descript);
   fprintf(OutputFile,"# %lf %s\n",V_peak, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",V_peak, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",V_peak, descript);

 fscanf(ParFile,"%lf %s",&t_peak, descript);
   fprintf(OutputFile,"# %lf %s\n",t_peak, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",t_peak, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",t_peak, descript);

 fscanf(ParFile,"%lf %s",&t_acute, descript);
   fprintf(OutputFile,"# %lf %s\n",t_acute, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",t_acute, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",t_acute, descript);

 fscanf(ParFile,"%lf %s",&AverageLogSP0, descript);
   fprintf(OutputFile,"# %lf %s\n",AverageLogSP0, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",AverageLogSP0, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",AverageLogSP0, descript);

 fscanf(ParFile,"%lf %s",&VarianceLogSP0, descript);
//   printf("# %lf %s\n",VarianceLogSP0, descript);
   fprintf(OutputFile,"# %lf %s\n",VarianceLogSP0, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",VarianceLogSP0, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",VarianceLogSP0, descript);

 fscanf(ParFile,"%lf %s",&MutationVariance, descript);
//   printf("# %lf %s\n",MutationVariance, descript);
   fprintf(OutputFile,"# %lf %s\n",MutationVariance, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",MutationVariance, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",MutationVariance, descript);

 fscanf(ParFile,"%ld %s",&tfinal, descript);
//   printf("# %ld %s\n",tfinal, descript);
   fprintf(OutputFile,"# %ld %s\n",tfinal, descript);
   fprintf(PatientRegistryOutputFile,"# %ld %s\n",tfinal, descript);
   fprintf(HeritOutputFile,"# %ld %s\n",tfinal, descript);

 fscanf(ParFile,"%lf %s",&prog_rate, descript);
//   printf("# %lf %s\n",prog_rate, descript);
   fprintf(OutputFile,"# %lf %s\n",prog_rate, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",prog_rate, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",prog_rate, descript);

 fscanf(ParFile,"%lf %s",&Heritability, descript);
//   printf("# %lf %s\n",Heritability, descript);
   fprintf(OutputFile,"# %lf %s\n",Heritability, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",Heritability, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",Heritability, descript);

 fscanf(ParFile,"%ld %s",&tprintHerit, descript);
//   printf("# %ld %s\n",tprintHerit, descript);
   fprintf(OutputFile,"# %ld %s\n",tprintHerit, descript);
   fprintf(PatientRegistryOutputFile,"# %ld %s\n",tprintHerit, descript);
   fprintf(HeritOutputFile,"# %ld %s\n",tprintHerit, descript);

 fscanf(ParFile,"%ld %s",&VL_print_time, descript);
//   printf("# %ld %s\n",VL_print_time, descript);
   fprintf(OutputFile,"# %ld %s\n",VL_print_time, descript);
   fprintf(PatientRegistryOutputFile,"# %ld %s\n",VL_print_time, descript);
   fprintf(HeritOutputFile,"# %ld %s\n",VL_print_time, descript);

 fscanf(ParFile,"%lf %s",&VL_print_lower, descript);
//   printf("# %lf %s\n",VL_print_lower, descript);
   fprintf(OutputFile,"# %lf %s\n",VL_print_lower, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",VL_print_lower, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",VL_print_lower, descript);

 fscanf(ParFile,"%lf %s",&VL_print_upper, descript);
//   printf("# %lf %s\n",VL_print_upper, descript);
   fprintf(OutputFile,"# %lf %s\n",VL_print_upper, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",VL_print_upper, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",VL_print_upper, descript);

 fscanf(ParFile,"%lf %s",&VL_print_interval, descript);
//   printf("# %lf %s\n",VL_print_interval, descript);
   fprintf(OutputFile,"# %lf %s\n",VL_print_interval, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",VL_print_interval, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",VL_print_interval, descript);

 fscanf(ParFile,"%ld %s",&random_number_seed, descript);
//   printf("# %ld %s\n",random_number_seed, descript);
   fprintf(OutputFile,"# %ld %s\n",random_number_seed, descript);
   fprintf(PatientRegistryOutputFile,"# %ld %s\n",random_number_seed, descript);
   fprintf(HeritOutputFile,"# %ld %s\n",random_number_seed, descript);

 fscanf(ParFile,"%d %s",&printOutput, descript);

 fscanf(ParFile,"%d %s",&CD4_Exp_Flag, descript);
   fprintf(OutputFile,"# %d %s\n",CD4_Exp_Flag, descript);
   fprintf(PatientRegistryOutputFile,"# %d %s\n",CD4_Exp_Flag, descript);
   fprintf(HeritOutputFile,"# %d %s\n",CD4_Exp_Flag, descript);

  fscanf(ParFile,"%d %s",&Donors_RandomTimeToAIDS_Flag, descript);
   fprintf(OutputFile,"# %d %s\n",Donors_RandomTimeToAIDS_Flag, descript);
   fprintf(PatientRegistryOutputFile,"# %d %s\n",Donors_RandomTimeToAIDS_Flag, descript);
   fprintf(HeritOutputFile,"# %d %s\n",Donors_RandomTimeToAIDS_Flag, descript);

  fscanf(ParFile,"%lf %s",&Increment_Prob_Under_Care_High_CD4, descript);
   fprintf(OutputFile,"# %lf %s\n",Increment_Prob_Under_Care_High_CD4, descript);
   fprintf(PatientRegistryOutputFile,"# %lf %s\n",Increment_Prob_Under_Care_High_CD4, descript);
   fprintf(HeritOutputFile,"# %lf %s\n",Increment_Prob_Under_Care_High_CD4, descript);

  fscanf(ParFile,"%d %s",&DefinitionHighSetPoint_CD4, descript);
   fprintf(OutputFile,"# %d %s\n",DefinitionHighSetPoint_CD4, descript);
   fprintf(PatientRegistryOutputFile,"# %d %s\n",DefinitionHighSetPoint_CD4, descript);
   fprintf(HeritOutputFile,"# %d %s\n",DefinitionHighSetPoint_CD4, descript);


}


void PrintHeaders()
{

 fprintf(OutputFile,"\n\ntime\tN\tSusc\tInfected\tVaccinated\tVaccinatedInfected\tResistInfected\tInfectedTreated\tDead\tdAIDS\tdNat\tAveLogV\t\tVstd\t\tSet_Ave\t\tSet_std\t\td_ave\t\td_std\t\tLinks\t\tV1\t\tda1\t\tSts1\n");
 fprintf(HeritOutputFile,"time\tIndex\tDuration\tDonorsIndex\tGen\tDonorsGen\tTimeInf\tVL\tTimeToAIDS\tDonorsTimeToAIDS\tDonorsVLatTrans\tLogSP0\tDonors_LogSP0\tVirContLogSP0\tDonors_VirContLogSP0\tEnvContLogSP0\tDonors_EnvContLogSP0\td_acute\tDonors_d_acute\tTotal_Time_Inf\tDonors_Total_Time_Inf_At_Trans\n");
 //printf("Initial Infected Population\nTime\tIndex\tGen\tVL\t\tLogSetPoint\tVContrSP0\tEnvContrSP0\td_acute\tTotal_Time_Inf\n");
 fprintf(PatientRegistryOutputFile,"time\tIndex\tDuration\tNumRecipients\tDonorsIndex\tGen\tDonorsGen\tTimeInf\tVL\tTimeToAIDS\tDonors_TimeToAIDS\tDonorsVLatTrans\tLogSP0\tDonors_LogSP0\tVirContrLogSP0\tDonors_VirContLogSP0\tEnvContLogSP0\td_acute\tDonors_d_acute\tTotal_Time_Inf\tDonors_Total_Time_Inf_At_Trans\tCD4\tCD4_TimeToAIDS\tCD4_TimeToAIDS_exp_1\tCD4_TimeToAIDS_exp_2\tCD4_TimeToAIDS_exp_3\tCD4_TimeToAIDS_exp_4\tCD4_initial_value\tStatus\n");
 if(printOutput==1)
 {
  printf("time\tN\tLinks\tUnInf\tInf\tDead\tdAIDS\tdNat\tAveLogV\t\tVstd\t\tSet_Ave\t\tSet_std\t\td_ave\t\td_std\t\tpLinks\t\tV[1]\t\td_a[1]\t\tSts[1]\n");
 }
}


void PrintStats(void)
{
  long num_notprimary;
  long vaccinated, vaccinated_infected, resistant_to_vaccine;
  long InfectedTreated = 0; // Number of HIV+ persons receiving treatment
   //printf("Inside PrintStats: N = %ld\n",N);
   if ((labs(time - tprintHerit) < 0.01) || (StopEarly == 1) ) {  //  print data for heritabilities graphs
     for (i=1; i<= N; i++) {
       if ((Status[i] == 1) && (Generation[i] > 1)) {
         fprintf(HeritOutputFile,"%ld\t%ld\t%6.5e\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n",
           time,i,Duration[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],RandomTimeToAIDS[i],Donors_RandomTimeToAIDS[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],Donors_EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i]);
       }
     }
   }

   Vsum = 0.0; d_sum = 0.0; G_sum = 0.0; set_sum = 0.0; Susceptible = 0; Dead = 0; DiedNat = 0; DiedAIDS = 0;Infected = 0; SumLinks = 0.0; AveLinks = 0.0; num_notprimary = 0;
   vaccinated = 0; vaccinated_infected = 0;resistant_to_vaccine = 0;
   for (i =1 ; i<= N;i++) {
     if (Status[i] == -2) DiedAIDS++;
     if (Status[i] == -1) DiedNat++;
     if ((Status[i] == -1)|| (Status[i] == -2)) Dead++;
     if (Status[i] == 0) {
        Susceptible++;
        if (Vaccinated[i] == 1) vaccinated++;
     }
     if (Status[i] == 1) {
       Infected++;
       G_sum = G_sum + 1.0*Generation[i];
       d_sum = d_sum + d_acute[i];
       set_sum = set_sum + LogSetPoint[i];
	     if (Treated[i] == 1) {
		     InfectedTreated++;
	     }
       if (Vaccinated[i] == 1) vaccinated_infected++; 
       if (ResistantToVaccine[i] == 1) resistant_to_vaccine++; 
     }
     if ( (Status[i] == 1) && (time - Time_Inf[i] > 50.0) )
     {
       Vsum = Vsum + log10(V[i]);
       num_notprimary++;
     }

     SumLinks = SumLinks + (double) NumLinks[i];
   }
   AveLinks = SumLinks / ((double) (Susceptible+Infected));
   if (Infected >= 1) {
     Vave = Vsum/(1.0*num_notprimary);
     d_ave = d_sum/(1.0*Infected);
     G_ave = G_sum/(1.0*Infected);
     set_ave = set_sum/(1.0*Infected);
   } else {
     Vave = 0.0;
     d_ave = 0.0;
     G_ave = 0.0;
     set_ave = 0.0;
   }
   Vstdsum = 0.0;  d_stdsum = 0.0; G_stdsum = 0.0; set_stdsum = 0.0; vcount = 0.0;
   for (i=1;i<=N;i++) {
     if (V[i] > epsilon) {
       vcount = vcount + 1.0;
        d_stdsum = d_stdsum + (d_acute[i] - d_ave)*(d_acute[i]-d_ave);
       set_stdsum = set_stdsum + (LogSetPoint[i] - set_ave)*(LogSetPoint[i]-set_ave);
       G_stdsum = G_stdsum + ((1.0*Generation[i]) - G_ave)* ((1.0*Generation[i]) - G_ave);
     }
     if ((Status[i] == 1) && ( time - Time_Inf[i] > 50.0)) {
       Vstdsum = Vstdsum + (log10(V[i]) - Vave)*(log10(V[i])-Vave);
     }

   }
   if (num_notprimary >= 1) Vstd = sqrt(Vstdsum/num_notprimary);
   else Vstd = 0.0;
   d_std = sqrt(d_stdsum/vcount); G_std = sqrt(G_stdsum/vcount); set_std = sqrt(set_stdsum/vcount);
   printcount++;
   if (printcount >= 20 && printOutput==1) {
    printf("\n\ntime\tN\tLinks\tUnInf\tInf\tInfTrt\tDead\tdAIDS\tdNat\tAveLogV\t\tVstd\t\tSet_Ave\t\tSet_std\t\td_ave\t\td_std\t\tpLinks\t\tV[1]\t\td_a[1]\t\tSts[1]\n");
     printcount = 0;
   }
   if(printOutput==1)
   {
     printf("%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t",time,N,(long)(SumLinks+0.5),Susceptible,Infected,InfectedTreated,Dead,DiedAIDS,DiedNat);
     printf("%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t",Vave,Vstd,set_ave,set_std,d_ave,d_std);
     printf("%3.2e\t%3.2e\t%3.2e\t%d\n",AveLinks,V[1],d_acute[1],Status[1]);
   }
   // fprintf(OutputFile,"\n\ntime\tN\tSusc\tInfected\tVaccinated\tVaccinatedInfected\tResistVaccine\tInfectedTreated\tDead\tdAIDS\tdNat\tAveLogV\t\tVstd\t\tSet_Ave\t\tSet_std\t\td_ave\t\td_std\t\tLinks\t\tV1\t\tda1\t\tSts1\n");
   //   long vaccinated, vaccinated_infected, resistant_to_vaccine;

   fprintf(OutputFile,"%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t",time,N,Susceptible,Infected,vaccinated, vaccinated_infected, resistant_to_vaccine,
                                                                 InfectedTreated,Dead,DiedAIDS,DiedNat);
   fprintf(OutputFile,"%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t",Vave,Vstd,set_ave,set_std,d_ave,d_std);
   fprintf(OutputFile,"%3.2e\t%3.2e\t%3.2e\t%d\n",AveLinks,V[1],d_acute[1],Status[1]);
   fflush(stdout);
}

double nrand(void)
{
  return(sqrt(-2.0*log(genrand_real2()))*cos(2.0*3.14159265358979*genrand_real2()));
}

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/


/* Period parameters */
#define N 624
#define M397 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long SSS)
{
    mt[0]= SSS & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M397;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M397] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M397-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M397-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */
