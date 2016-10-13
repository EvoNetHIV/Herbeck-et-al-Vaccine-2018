/*
 *  NetworkHeader.h
 *  
 *
 *  Created by John Mittler on 12/6/12.
 *  Copyright 2012 University of Washington. All rights reserved.
 *
 */


#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*Subroutines */
void GetParameters(void), PrintStats(void);
long printcount;
int printOutput;
double jmpow(double a, double b);
double GetGammaDelay(double k, double theta);
void AddNewLink(void);
void AddLinkForSinglePerson(long person1);
long GetPoisson(long num, double pois_prob);
long max(long a, long b);

void OpenFiles(void);
void GetParameters(void);
void InitializePopulation(void);
void AddInfecteds(void);
void DefineSexualContactNetwork(void);
void PrintEdgeList(void);
void PrintHeaders(void);
void PrintStats(void);
void UpdateViralLoadsForAllInfectedIndividuals(void);
void SimulateInfectionProcessOnePersonToAnother(void);
void BirthOfNewSusceptibles(void);
void NaturalDeath(void);
void DeathOfAIDSPatients(void);
void RemoveLinksToDeadPeople(void);
void SimulatePartnerShipDissolution(void);
void SimulateNewPartnershipFormation(void);
void CheckEdgeListForImpossibleValues(void);  // Does some elementary error checking
void RecordStatusOfLivingHIVPatientsForRegistryFile(void); 
void PrintFirstAndLastContacts(void);

/* Files used by this program */
FILE *ParFile; // Contains parameters
FILE *OutputFile; // Records average and standard deviation of viral loads (and other parameters) at each time step
FILE *HeritOutputFile; // Records relationship between donors and recipients at a single timepoint (tprintHerit)
FILE *PatientRegistryOutputFile; // Records characteristics of each patient who died and state of living patients at tfinal
FILE *VLOutputFile; // Gives a histogram of log viral loads at different times
FILE *spVLOutputFile; // Gives a histogram of (log-transformed) set point viral loads at different times

/*Epidemiological parameters */
double Inf_Prob; // Probability of infection (calculated from quantities below)
double MaxInfRate, VHalfMaxInfRate, HillCoeffInfRate, beta; // Maximum probability of infection,  viral load at which infection prob is 0.5 of the max
double BaselineDeathRate, MaxAIDSDeathRate, VHalfMaxAIDSDeathRate, HillCoeffAIDSDeathRate, alpha; // Rate of inputs and deaths of people
double yearly_rate, daily_rate; // rates of progression to AIDS using John's new formula
double prog_rate; // prog_rate is rate of increase of viral load per year
long Inf_Rate_Function; // 1 = power function from Lingappa et al, 0 = hill function  
double InfRateBaseline, InfRateExponent; // Prob(Infmission) = InfRateBaseline * [log(VL)] ^ InfRateExponent
double BirthRate; // Rate of entry of new susceptibles into the population.
double NaturalDeathRate; // Rate at which people normally die (of causes other than AIDS)
long MaxN = 2000000; // NOTE: This is the extent of the remaining arrays minus 1
long MaxNReached = 0; // Flag to indicate if population was about to exceed the maximum array size
double tMaxNReached = 0.0; // Time that births or reincarnations were halted
long NewAIDSDeaths = 0; // Number of deaths in each time period (only used when reincarnation is allowed)
long Widowed_partner;
long newbirths; // Internal counter to count number of births each day
double Active; // Internal -- number of sexually active people
 
/* Virological parameters */
double V[2000001], Donors_V[2000001], s[2000001];
double V0, r0, V_peak, t_peak, t_acute, Vss;
double d_acute[2000001], Donors_d_acute[2000001];
double death_rate_constant, death_rate_exponent; // Death rate parameters for daily probability death function (used only when Gamma_Death !=1)
long CD4_Death; // Flag indicates whether to use Fraser et al.s unpublished CD4 progression rate function (1=use CD4 death function)
long Gamma_Death; // Flag indicates whether to use Frasier's gamma delay term (1=use gamma delay death function)
double RandomTimeToAIDS[2000001]; 
double Donors_RandomTimeToAIDS[2000001];
int Donors_RandomTimeToAIDS_Flag;

double ExpectedDelayTime; // Assume gamma distributed delays
double Dmax, Dk, D50; // Parameters governing Fraser's gamma distributed times to AIDS
double shape_parameter,theta; // These parameters control the mean and variance of the gamma delay

/* Infection status parameters */
long N, Infected, Susceptible; // Initial number of people, Initial number infected, Number who could be infected
long currN;  // Value of N used to calculate the number of births that will occur
long Dead, DiedAIDS, DiedNat; // Cumulative Number Died of AIDS or natural causes
double prob_sex; // Probability that partners will have sex each day
double Time_Inf[2000001], Donors_Total_Time_Inf_At_Trans[2000001]; // Time infected
int Status[2000001]; // Status of each person (Infected=1, Suscp = 0, Died = -1, DiedAIDS = -2)
long Generation[2000001], Donors_Generation[2000001], Donors_Index[2000001], NumRecipients[2000001];
int HadSexAlready[2000001];

/* Parameters governing heritability of VL setpoint */
double AverageLogSP0, VarianceLogSP0, Heritability; // Primary inputs
double SetPoint[20000010], Donors_SetPoint[20000010]; // Patient specific values (drawn from random number generator)
double LogSetPoint[20000010], Donors_LogSetPoint[20000010]; // Patient specific values (drawn from random number generator)
double ViralContribToLogSP0[2000001],EnvirContribToLogSP0[2000001]; // Viral and environmental deviations that contribute to final setpoint
double Donors_ViralContribToLogSP0[2000001], Donors_EnvirContribToLogSP0[2000001]; // Viral and environmental deviations that contribute to final setpoint
int CD4[2000001];
int CD4_initial_value[2000001];
double CD4time[2000001];
int SpvlCat[2000001];
double CD4_TimeToAIDS[2000001];
double CD4_TimeToAIDS_exp[2000001][5];
int CD4_Exp_Flag;
double Increment_Prob_Under_Care_High_CD4;
int DefinitionHighSetPoint_CD4;
long CD4_treatment_delay_index[2000001];


/* Look up table for CD4 counts */
double CD4_lookup[10][5];//cd4 as function of SPVL and time since infection
void CD4_Category();
int SPVL_Category(double SPVL_param);
double cd4_initProbs[10][4];
void cd4InitProbs();
void CreateCD4Table();
int initialCD4(int spvl_level);
int TreatmentThreshold_CD4; // VL threshold that triggers initiation of therapy in patients under care.
int CD4_After_Treatment; // Viral load after treatment starts
int CD4_Determines_Treatment; //if 0, vl determines treament, if 1, then cd4

/* Statistical results (and internal params needed to get those statistics) */
double Vsum, Vave, Vstd, Vstdsum; // Statistics of viral load
double set_ave, set_sum, set_std, set_stdsum; // Statistics on logSetPoint
double d_ave, d_sum, d_std, d_stdsum; // Statistics on decay rate of infected cells after peak (our way of setting setpoints)
double G_ave, G_sum, G_std, G_stdsum; // Statistics on viral generation
double vcount; // Needed for averaging above quantities
double AveLinks,SumLinks;  // Average number of links per person

/* Treatment Parameters */
long UnderCare[2000001]; // Pre-defined list of patients who would get treated after infection (e.g., they live near a clinic)
long Treated[2000001]; // Patients currently being treated
double Start_Treatment_Campaign; // Time (day of epidemic) when public health authorities start new treatment campaign
double percent_under_care; // Percent of people who get treated after being infected
double Increment_Prob_Under_Care_High_spVL;  // Increase in probabity of going under care for folks with high setpoint VLs
double DefinitionHighSetPoint;  // VL Threshold above which someone is more likely to come under care
double TreatmentTimeInPatients; // Time (within person) after infection that people treated if under care.  Assumed to be the same for all patients
double TreatmentThreshold; // VL threshold that triggers initiation of therapy in patients under care.
double VL_After_Treatment; // Viral load after treatment starts

double Start_SexReduction_Campaign; // # Day of epidemic when public health authorities start treatment campaign
double Reduction_Mean_Degree; //  Reduction in mean degree resulting from this campaign

double Start_Faithfulness_Campaign; // Day of epidemic when public health authorities start campaign to reduce partnership turnover
double Increased_Duration; // Increase in partnership durations resulting from this campaign

/* Time-keeping variables */
long time, tfinal, tprintHerit;

/* Stopping paramters */
int StopEarly = 0;   // This variable is set to 1 when the program is unable to find a new link

/* Internal variables */
long i, j, k, count, daycount, deadcount;
double epsilon = 1e-10;
char junk[12],descript[30];

/* Variables required for random number generation*/
double nrand(void); // Normal random numbers with mean 0 and stdev 1
double randnum, normalrand; // internal variables
long random_number_seed; // This is the randum number seed
void init_genrand(unsigned long SSS); // All of the rest are existing MT functions
long genrand_int31(void);
unsigned long genrand_int32(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

/* Sexual network parameters */
long NumLinks[2000001], Links[2000001][10]; // Number of sexual links for each person, Identities of those links (up to 5)
long TotalEdges; // Total number of sexual links in the entire population (note that one link links two individuals
double AverageDegree; // Average number of sexual partners person
long MaxLinks; // Maximum number of links per person (MaxLinks = 1 implies monogamy)
long NewLinks; // Temporary variable giving the number of links for each new person entering the population
long Link_back; // Internal variable to find links back from widow's of the dead person back to the dead person
double dp_TotalEdges; // double precision version of TotalEdges (for reading in only, immediately converted to long)
double Duration[2000001]; // Length that person i tends to stay in a relationship
double MinDuration, MaxDuration; // Mins and Maxs used to draw a uniform random number for each person
double PartnershipDuration;  // Duration of each partnership (mean of the underlygin partnership duration of each partner)
double ProbPartnersBreakUp;  // Prob of partnerships breaking up 
long PartnersOfPerson1; // Internal parameter: index of each person's partners 
long LinksToRemove; // Temporary variable indicating number of links per dead person that should be purged
long AlreadyLinked; // Bookkeeping variable to prevent same persons from linking up more than once (but still counting this as two links)
long failed_to_find_a_new_link = 0; // Set to one if algorithm is unable to find a new pair (indicates some kind of programming problem)
double tr; // Turnover rate = Links_Broken_Per_Day/TotalEdges;
long escape_counter; // Used to prevent endless loops (in case of programming error)
long random_partner; // Which of person i's partners (if more than one) does person i have sex with today
long index_to_person2, found_link, index_to_person1; // Used internally to add and break links
long link_break; // Used internally to count up links that need to broken
long breakups; // Internal parameter: how many people broke up each day
double NLinks; // Number of people linked after deaths, births, and breakups
double dpLinksToAdd; // Number of links to add after accounting for deaths, births, and breakups
long LinksToAdd; // Number of links to add after accounting for deaths, births, and breakups
double ExpectedLinks; // Expected number of links that you would have had had their been now births, deaths, or breakups 
long HadDifficutlyFindingLinks = 0; // Flag set when the program had difficulty finding new links

