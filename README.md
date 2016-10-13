## Overview

We use two agent-based, stochastic models of HIV epidemiology and evolution to study the potential impact of HIV evolution in response to vaccination.

The models used in this study are based on models that have been described previously (Herbeck, 2014, 2016; Handcock, 2008; Krivitsky, 2012), with modifications to enable simulation of vaccination programs and viral resistance and sensitivity to the vaccine response. These are individual-based stochastic, dynamic models that track inter-host and intra-host dynamics. Agents are endowed with a large number of individual attributes, including demographics, clinical features and behavioral characteristics including sexual role preference (for MSM). HIV-infected individuals possess numerous viral characteristics such as a “set point viral load” (SPVL) and a viral load at each time step over the course of infection.

### Description of the models

#### Heterosexual epidemic model (ViralDynamics)

We simulated a heterosexual HIV epidemic that was calibrated to reproduce incidence and prevalence trajectories based on data from South Africa. We chose to calibrate the model in this way for several reasons: first, because South African has one of the largest HIV epidemics; second, because the large-scale HIV vaccine trial HVTN 702 is located in South Africa, and is based upon the RV144 trial design and results; and third, so that the epidemic output was directly comparable to the 12 HIV epidemic models described in Eaton et al. (29). We first set empirical parameter values (i.e., viral load trajectories, CD4 progression rates, set point viral load distribution), then subsequently calibrated the assumption-based parameter values to fit expected epidemic dynamics. Parameters and initial values are listed in the Supplementary Table 1.

For each individual, upon HIV-1 infection with an initial viral population size, viral load increased exponentially until a peak viremia at the midpoint of the length of acute infection. Viral load then declined exponentially at an individual-specific decay rate until it reached set point viral load (SPVL). After reaching set point, the viral load increased in a log linear manner until the onset of AIDS. Viral load upon onset of AIDS is defined as the same for all individuals, and is independent of SPVL. To model the relationship between SPVL and disease progression, the model related individual SPVL to the starting CD4 count category (individual CD4 count immediately after infection) and to disease progression based on waiting times in four CD4 count categories: CD4≥500; 500>CD4≥350; 350>CD4≥200; CD4<200 (AIDS).

Transmission rates were assumed to follow available data from serodiscordant heterosexual partners; the probability that an HIV-infected person will transmit to an HIV-negative person is determined by a Hill function that follows from Fraser (54). The individual SPVL is determined by both viral and environmental factors, with the viral component being heritable across transmissions. For each individual in a simulation, the model maintains a list of sexual partnerships and viral transmission pairs.

For any one epidemic simulation, an initial set of contacts was formed by randomly choosing pairs from the population, prior to the simulation of viral transmission. The probability of each person entering into this link was set to e-Li, where Li is the number contacts for person i. If this probability was not met (or if Li > MaxLinks), another partner was selected at random from the population. This process was repeated until the total number of links equals N*M/2, where N is the population and M is the mean degree. The probability of a connection between individuals i and j dissolving is set to 2/(Duration[i] * Duration[j]), where Duration[x] is the expected time that person x stays in a relationship. 

The viral and CD4-based progression functions are embedded in a host population without demographic (sex, age) heterogeneity, but with behavioral heterogeneity. Each individual is assigned a relational duration propensity; when two individuals partner, their relationship lasts for the mean of these individual effects. The daily probability of sex varied with relationship duration. We did not include a change in behavior patterns over calendar time; e.g., we did not include a reduction in the average sexual contact rate, across individuals, as the epidemic progresses. Supplementary Table 1 lists parameters and parameter settings.

We ran the heterosexual epidemic simulation for 40 years, in order to mimic HIV epidemics from 1990 to 2030 in South Africa. We calibrated epidemic simulations to begin with 2% HIV incidence (per 100 person years), which rose quickly to ~3%, followed by a decline to ~1.4% and stabilization by 2010. Prevalence rose to ~12%, followed by a decline to ~9% with ART coverage at approximately 30%. 

#### Men-who-have-sex-with-men (MSM) epidemic model

In the MSM HIV epidemic model, the dynamics of sexual relationships were governed by network models based in the separable-temporal exponential-family random graph model (STERGM) framework.

Each simulation started with a population of HIV-uninfected and infected individuals at time zero. Dyads formed and dissolved partnerships to match target network parameters. Within ongoing partnerships, we modeled individual coital acts, as well as decisions about condom use per act. For each serodiscordant coital act that occurred, the probability of transmission from an HIV-infected person to an HIV-negative person was based primarily on viral load but also on disease stage (above and beyond viral load), condom use, and antiretroviral use; this is described by a complementary-log log link function that follows from Hughes (2012).

Individual viral load trajectories were modeled with five phases in the absence of treatment:  an initial increase from infection to acute peak; a two-phase decay to set point viral load; a long, chronic phase with slow viral load increase, and a final AIDS stage. Individuals varied in their set point viral load, with individual values set to have both a heritable component from their infector, and a random (environmental) component.

The MSM model varied in several key ways from the source model, as befits the aims of this work. We simplified from three separate network models (for main, casual, and one-time partners) in those models to just one overall network per model, since we were not focused on understanding the evolutionary response to vaccination stratified by partner and exposure types. Four features already mentioned—variation in SPVL, post-acute biphasic decay in viral load, administration of vaccine, and the existence and heritability of vaccine resistance—were new to this model. Finally, as a result of these changes the model contains numerous novel parameters, which were derived from the data described in Supplementary Table 2.

#### Vaccine programs

The vaccine parameters related to duration and the effect on transmission/infection are identical in both of our models, and are as described in the main text. 20 replicate simulations were performed for each combination of parameters for both models. For the MSM simulations, the vaccination program was initiated 10 years after equilibrium prevalence of ~25% was reached, and repeated additionally at five year intervals (i.e., at years 5 and 10 after initial vaccine rollout). Across viral transmissions there was 100% heritability of the vaccine-related genotype, and no within-host viral evolution took place that affected the genotype. As described in the main text, the duration of vaccine efficacy (VE) was three years, with VE reduced immediately to 0% at the end of this period; while VE in RV144 declined from ~60% to ~30% over 3.5 years, we chose to employ a simplified flat VE trajectory, as recent RV144-related modeling work by Gilbert et al. has shown that the shape of the temporal VE profile has a negligible, and statistically insignificant, impact on model results (58).

#### Antiretroviral therapy

Antiretroviral therapy (ART) became available in 2012 in the heterosexual South African model, and was available at all time points in the MSM model. In both models, individuals became eligible for treatment at CD4<350 and population coverage was set at approximately 30% and 40%, for the HET and MSM models, respectively. Identical ART scenarios were included in the counterfactual epidemic scenarios (no vaccine and vaccine with no viral variation) and the epidemic scenarios that included viral resistance to the vaccine effect.

### Authors of EvoNetHIV and ViralDynamics

* Steve Goodreau <goodreau@uw.edu>
* Josh Herbeck <herbeck@uw.edu>
* John Mittler <jmittler@uw.edu>
* James Murphy <jtm6@uw.edu>

