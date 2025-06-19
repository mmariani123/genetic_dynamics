#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(){

void transpose_1d_array(int arr[], int size) {
    int temp;
    for (int i = 0; i < size / 2; i++) {
        temp = arr[i];
        arr[i] = arr[size - 1 - i];
        arr[size - 1 - i] = temp;
    }
}

int main() {
    int arr[] = {1, 2, 3, 4, 5};
    int size = sizeof(arr) / sizeof(arr[0]);

    printf("Original array: ");
    for (int i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    transpose_1d_array(arr, size);

    printf("Transposed array: ");
    for (int i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    return 0;
}

// matrix dimensions so that we dont have to pass them as
// parametersmat1[R1][C1] and mat2[R2][C2]
#define R1 2 // number of rows in Matrix-1
#define C1 2 // number of columns in Matrix-1
#define R2 2 // number of rows in Matrix-2
#define C2 3 // number of columns in Matrix-2

void multiplyMatrix(int m1[][C1], int m2[][C2])
{
    int result[R1][C2];

    printf("Resultant Matrix is:\n");

    for (int i = 0; i < R1; i++) {
        for (int j = 0; j < C2; j++) {
            result[i][j] = 0;

            for (int k = 0; k < R2; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }

            printf("%d\t", result[i][j]);
        }

        printf("\n");
    }
}

//# -*- coding: utf-8 -*-
//#! /usr/bin/env python3


//"""
//Created on Sat Jan 11 12:24:02 2025

//@author: michael p. mariani PhD, 2025

//"""
//import math
//import numpy as np
//from pytictoc import TicToc

//# -- Simulate CTMC model of within-vector parasite dynamics using modified
//# Gillespie Algorithm (outlined in Supplementary Materials), with
//# continuous rupture function.

//#ported over from the
//#nParasiteGroups_ContRuptFunc.m
//#MATLAB / Octave file

//#close all
//#clear

//SEED = rng('shuffle'); # set seed

srand(time(NULL));   // Initialization, should only be called once.
int seed = rand();      // Returns a pseudo-random integer between 0 and RAND_MAX.

//# Set parameters equal to values defined in Create_Parameter_Set.m:
//#par = Create_Parameter_Set;
//student = {
//    "name": "Alice",
//    "age": 20,
//    "major": "Computer Science"
//}

/*
#The numpy equivalent of repmat(a, m, n) is tile(a, (m, n)).
#From the "Create_Parameter_Set.m" file:
par.N        = 2 #% Number of different starting genotypes
par.a        = 24/(20/60) #% failure rate of male gametocytes (per day)
par.b        = 24/(25/60) #% failure rate of female gametocytes (per day)
par.r        = 0.08 #% fertilization of male and female gametes (per day)
par.mu_z     = 1 #% death rate of zygotes (per day)
par.sigma_z  = 24/19 #% transformation rate of zygotes (per day)
par.sigma_e  = 0.6 #% transformation rate of ookinetes (per day)
par.mu_e     = 1.4 #% death rate of ookinetes (per day)
par.mu_o     = 0 #% death rate of oocysts (per day) **changed from paper**
par.n        = 3e3 #% number of sporozoites per oocyst
par.p        = 0.2 #% proportion of sporozoites that make it to salivary gland
par.alpha    = 0.39 #% fraction of male gametes that are viable
par.beta     = 1 #%0.96 % fraction of female gametes that are viable
par.eta      = 1 #% number of female gametes per female gametocyte
par.rho      = 8 #% number of male gametes per male gametocytes
par.k        = 1/7
par.t0       = 10
par.max_bias = 0.1 #%0.5 %0.1
par.NumSim   = 1 #%10000

params = {
    "N" = par.N;
    "a" = par.a;
    "b" = par.b;
    "r" = repmat(par.r,N,N);
    "mu_z" = par.mu_z;
    "sigma_z" = par.sigma_z;
    "sigma_e" = par.sigma_e;
    "mu_e" = par.mu_e;
    "mu_o" = par.mu_o;
    "n" = par.n;
    "p" = par.p;
    "alpha" = par.alpha;
    "beta" = par.beta;
    "eta" = par.eta;
    "rho" = par.rho;
    "k" = par.k;
    "max_bias" = par.max_bias;
    "t0" = par.t0;
}
*/

unsigned int sN        = 2          //#% Number of different starting genotypes
float a                = 24/(20/60) //#% failure rate of male gametocytes (per day)
float b                = 24/(25/60) //#% failure rate of female gametocytes (per day)
float r                = 0.08       //#% fertilization of male and female gametes (per day)
unsigned short int muZ = 1          //#% death rate of zygotes (per day)
float sigmaZ           = 24/19      //#% transformation rate of zygotes (per day)
float sigmaE           = 0.6        //#% transformation rate of ookinetes (per day)
float muE              = 1.4        //#% death rate of ookinetes (per day)
unsigned short int mu0 = 0          //#% //death rate of oocysts (per day) **changed from paper**
unsigned short n       = 3e3        //#% number of sporozoites per oocyst
float p                = 0.2        //#% proportion of sporozoites that make it to salivary gland
float alpha            = 0.39       //#% fraction of male gametes that are viable
float beta             = 1          //#%0.96 % fraction of female gametes that are viable
unsigned short eta     = 1          //#% number of female gametes per female gametocyte
unsigned short rho     = 8          //#% number of male gametes per male gametocytes
float k                = 1/7
unsigned short t0      = 10
float maxBias          = 0.1        //#%0.5 %0.1
unsigned short NumSim  = 1          //#%10000

//#example:
//#arr2 = np.array([[1, 2], [3, 4]])
//#result2 = np.tile(arr2, (2, 3))
//#print(result2)

//#percent_male = repmat(0.25,N,1);
//#% proportion of gametocytes that are male of each genotype
//percent_male = tile(0.25,N,1)
float percentMale[sN] = {0.25, 0.25};

//#% Define vector of fitness bias for each genotype:
//#    (this can be modified to consider >3 genotypes)
//#% For N>3, must define vector of fitness biases.

//#if N == 1
//#    bias = 0;
//#elseif N == 2
//#    bias = [0,max_bias];
//#elseif N == 3
//#    bias = [0, 0.1, 0.5];
//#end

//if N==1
//    bias = 0
//elif
//    N==2
//    bias = [0,max_bias]
//elif N==3
//    bias = [0,0.1,0.5]

if (sN==1) {
    float bias[1] = {0.0};
} else if (sN==2) {
    float bias[2] = {0.0, maxBias};
} else {
    float bias[3] = {0.0, 0.1, 0.5};
}

//#MATLAB uses the apostrophe operator ( ' )
//#to perform a complex conjugate transpose,
//#and the dot-apostrophe operator ( . ' )
//#to transpose without conjugation.
//#For matrices containing all real elements,
//#the two operators return the same result. produce the same scalar result.

//#bias = bias';
//bias = np.transpose(bias)
bias = transpose_1d_array(bias);

//#.* means matrix product, if you don't write .
//#will Matlab product the numbers on the same position.

//#a_vec = repmat(a,N,1).*(1 + bias);
//a_vec = np.matmul(tile(a,(N,1)), (1 + bias))
float aMat1[sN] = {a,a};
float aMat2[1] = {1+bias};
float aVec = aMat1*aMat2

//#b_vec = repmat(b,N,1).*(1 + bias);
//b_vec = np.matmul(tile(b,(N,1)), (1 + bias))
float bMat_1[sN] = {b,b};
float bMat_2[1] = {1+bias};
float bVec = bMat1*bMat2

//#alpha_vec = repmat(alpha,N,1).*(1 - bias);
//alpha_vec = np.matmul(tile(alpha,(N,1),(1-bias))
float alpha_mat_1[sN] = {alpha,alpha};
float alpha_mat_2[1] = {1+bias};
float alpha_vec = alphaMat1*alphaMat_2

//#beta_vec = repmat(beta,N,1).*(1 - bias);
//beta_vec = np.matmul(tile(beta,(N,1),[1 - bias])
float betaMat_1[sN] = {beta,beta};
float betaMat_2[1] = {1+bias};
float betaVec = betaMat1*betaMat2

//#mu_z_vec = repmat(mu_z,N,1).*(1 + bias);
//mu_z_vec = np.matmul(tile(z,(N,1)),(1 + bias))
float muz_mat_1[sN] = {muZ,muZ};
float muz_mat_2[1] = {1+bias};
float muz_vec = muz_mat_1*muz_mat_2

//#% zygote mortality rate for zygotes whose parents are of the same genotype

//#mu_e_vec = repmat(mu_e,N,1).*(1 + bias);
//mu_e_vec = np.matmul(tile(mu_e,(N,1)),(1+bias)
float mueMat1[sN] = {muE,muE};
float mueMat2[1] = {1+bias};
float mueVec = mueMat1*mueMat2

//#% ookinete mortality rate for ookinetes " ...

//#mu_o_vec = repmat(mu_o,N,1).*(1 + bias);
//mu_o_vec = np.matmul(mu_o,(N,1),(1+bias))
float muoMat_1[sN] = {muO,muO};
float muoMat_2[1] = {1+bias};
float muoVec = muoMat_1*muoMat_2

//#% oocyst mortality rate oocysts " ...

//#% Create symmetric matrices of biased parameters.  The (i,j) entry is the
//#% parameter value for a parasite that has a genotype (i) female, and
//#% genotype (j) male parent.

//#np.full(shape=5, fill_value=None)
//#np.full(shape=5, fill_value=np.nan)

//#np.full(shape=(N,N), fill_value=np.nan)

//#mu_z_mat = NaN(N,N);
//#mu_e_mat = NaN(N,N);
//#mu_o_mat = NaN(N,N);

//#check = np.full(shape=(2,2), fill_value=np.nan)
//#check.shape

//mu_z_mat = np.full(shape=(2,2), fill_value=np.nan)
//mu_e_mat = np.full(shape=(2,2), fill_value=np.nan)
//mu_o_mat = np.full(shape=(2,2), fill_value=np.nan)

int muzMat[sN][sN] = {NAN};
int mueMat[sN][sN] = {NAN};
int muoMat[sN][sN] = {NAN};

bool checkNan = is.nan(muzMat[0][0]);
printf(checkNan)

//#for index, value in enumerate(my_list):
//#    print(index, value)

//for i in range(N):

//   for j in range(N):

//        mu_z_mat[i,j] = np.mean([mu_z_vec(i),mu_z_vec(j)])
//        mu_e_mat[i,j] = np.mean([mu_e_vec(i),mu_e_vec(j)])
//        mu_o_mat[i,j] = np.mean([mu_o_vec(i),mu_o_vec(j)])

for (int i = 0; i < sN; i++) {

    for(int j = 0; j < sN; j++){

        muzMat[i][j] = (muzvec[i]+muzVec[j])/2
        mueMat[i][j] = (mueVec[i]+mueVec[j])/2
        muoMat[i][j] = (muoVec[i]+muoVec[j])/2

    }

}

//#% --------------------------------------------------------------------%
//#% ---- Run NumSim number of simulations for each value of G0 in G0_vec as
//#% defined in CreateParameterSet.m:

//#% Number of simulations %
//NumSim = par.NumSim;

//#% Max number of iterations for each simulation %
unsigned short maxIter = 1e3;

#% Vector of Initial gametocyte densities %
G0_vec = list(range(150,450+50,50))

unsigned short lowDensity = 150;
unsigned short highDensity = 450;
unsigned short densityStep = 50;
unsigned short densityArraySize = (highDensity-lowDensity)/densityStep;
unsigned short densityArray[densityArraySize] = {};
unsigned short densityCounter = 0;
for(int i = lowDensity; i <= highDensity; i += densityStep) {

    densityArray[counter] = i
    densityCounter+=1;

}

//#% Max time (in days)
//Tfinal = 21
unsigned short tFinal = 21;

//#% Time-step to increment by if the time to next time-step is greater than
//#% dt.
//dt = 0.1
float dt = 0.1;

#% Fraction of G0 of each subtype (Nx1 vector):
#strainprop = repmat(1/N,N,1);
//strainprop = tile(1/N,(N,1))
float strainProp[sN] = {1/sN, 1/sN};

//#TransitionIDs = NaN(NumSim,maxiter);
//TransiitonIDs = np.full(shape=(NumSim,maxiter), fill_value=np.nan)
float transitionsIDs[NumSim][maxIter] = {NaN}

lengthG0Vec = sizeof(g0Vec)
for ii = 0; ii < lengthG0Vec, ii++){

    //g0 = round(strainProp*g0Vec[ii])
    //#% Vector of initial gametocyte densities.
    //#G0(i) is the # of genotype i gametocytes.
    g0 = round(strainProp*g0Vec[ii])

    //#Gf0 = ceil((1 - percent_male).*G0.*eta.*beta)
    //gF0 = math.ceil((1 - percentMale)*G0*eta*beta)
    //#% Female gamete densities of each genotype
    gF0 = math.ceil((1 - percentMale)*G0*eta*beta)

    //#Gm0 = ceil(percent_male.*G0.*rho.*alpha)
    //gM0 = math.ceil(percent_male*G0*rho*alpha)
    //#% Male gamete densities of each genotype
    gM0 = math.ceil(percent_male*G0*rho*alpha)

    //#Time_data = NaN(maxiter,NumSim); % Initialize Time_data vector.
    //#Time_data = NaN(maxiter,NumSim); % Initialize Time_data vector.
    //timeData = np.full(shape=(maxiter,numSim), fill_value=np.nan)
    timeData = np.full(shape=(maxiter,numSim), fill_value=np.nan)

    //#Male_data = zeros(N,maxiter,NumSim);
    //maleData = np.zeros((N, 2, 3))
    //#% Each row is a different genotype;
    //#each column is a different time step;
    //#each layer is a different simulation
    maleData = np.zeros((N, 2, 3))

    //#Female_data = zeros(N,maxiter,NumSim);
    //femaleData = np.zeros(N,maxiter,NumSim)
    //#% Same as for Male_data
    femaleData = np.zeros(N,maxiter,NumSim)

    //#Zygote_data = zeros(N,N,maxiter,NumSim);
    //zygoteData = zeros((N,N,maxIter,numSim)
    //#% 4D array: Zygote(i,j,k,l) is the # of Zygotes
    //#with female parent i, male parent j, at timestep k,
    //#in simulation l
    zygoteData = zeros((N,N,maxIter,numSim)

    //#Ookinete_data = zeros(N,N,maxiter,NumSim)
    //#% "
    //ookineteData = np.zeros(N,N,maxIter,numSim)
    //#Oocyst_data = zeros(N,N,maxiter,NumSim)
    //#% "
    ookineteData = np.zeros(N,N,maxIter,numSim)

    //#Sporozoite_data = zeros(N,N,maxiter,NumSim)
    //#% "
    sporozoiteData = np.zeros(N,N,maxIter,numSim)

    //#Burst_time = zeros(N,N,maxiter,NumSim)
    //#% "
    //burstTime = np.zeros(N,N,maxIter,numSim)
    burstTime = np.zeros(N,N,maxIter,numSim)

    //#from pytictoc import TicToc
    //#t = TicToc()
    //#t.tic()  # Start timer
    //# Code to time
    //#t.toc()  # Print elapsed time

    //#tic

    //t = TicToc()
    //t.tic()

    for j in range(NumSim):

        Male_data[:,1,j] = Gm0
        #% Set initial conditions for male and female gamete data.

        Female_data[:,1,j] = Gf0
        #% "

        #burst_count[:,:,j] = zeros(N,N); % Initialize burst_count 3-D array.
        burst_count[:,:,j] = np.zeros((N,N))

        t[1] = 0
        time = t[1]

        for i in range(2:maxiter):

            #y1 = rand;
            #y2 = rand;

            y1=numpy.random.rand
            y2=numpy.random.rand

            #% State variables
            m = Male_data[:,i-1,j]
            #% Column vector where each row is a different genotype,
            #and the elements are the # of Male gametes at the previous
            #time-step (t_{i-1}) in simulation j.  So, m(k) is the number
            #of genotype-k males at time t_{i-1} in sim j.
            f = Female_data[:,i-1,j]
            z = Zygote_data[:,:,i-1,j]
            #% z(k,l) = # zygotes with female-k parent and male-l parent
            #at time t_{i-1} in sim j.
            e = Ookinete_data[:,:,i-1,j]
            o = Oocyst_data[:,:,i-1,j]
            s = Sporozoite_data[:,:,i-1,j]

            #% Break if states m through o are extinct:
            state_variables = [m,\
                               f,\
                               np.reshape(z,(N*N,1)),\
                               np.reshape(e,(N*N,1)),\
                               np.reshape(o,(N*N,1))]

            if np.sum(state_variables(:)) == 0
                LastIter[j] = i-1
                break

            #% Assign data values at current time-step, i,
            #to values at previous time-step, i-1.
            Male_data[:,i,j] = m
            Female_data[:,i,j] = f
            Zygote_data[:,:,i,j] = z
            Ookinete_data[:,:,i,j] = e
            Oocyst_data[:,:,i,j] = o
            Sporozoite_data[:,:,i,j] = s

            #% Define bursting function.
            bursting = 1/(1 + (math.exp(t0 - t(i-1))))
            #% probability of bursting at time 'time'
            #% Create vectors/matrices of all possible transitions,
            #arranged by type: death, mating, maturation, bursting.
            Male_death = a_vec*m;
            Female_death = b_vec*f;
            Zygote_death = mu_z_mat*z;
            Ookinete_death = mu_e_mat*e;
            Oocyst_death = mu_o_mat*o;
            #Mating = r.*(f*m');
            #* is fine for numpy arrays in python for
            #element-wise multiplication
            Mating= r*(f*np.transpose(m))
            Zyg_maturation = sigma_z*z;
            Ook_maturation = sigma_e*e;
            Ooc_Bursting = k*bursting*o;

            #% Calculate all possible transitions:
            transitions = \
                [Male_death,\
                Female_death,\
                np.reshape(Zygote_death,(N*N,1)),\
                np.reshape(Ookinete_death,(N*N,1)),\
                np.reshape(Oocyst_death,(N*N,1)),\
                np.reshape(Mating,(N*N,1)),\
                np.reshape(Zyg_maturation,(N*N,1)),\
                np.reshape(Ook_maturation,(N*N,1)),\
                np.reshape(Ooc_Bursting,(N*N,1))]

            #Create vector of indices marking different transition
            #types in the transition vector:

            transition_type_indices = \
                [0,\
                cumsum([length(Male_death),\
                length(Female_death),\
                length(np.reshape(Zygote_death,(N*N,1))),\
                length(np.reshape(Ookinete_death,(N*N,1))),\
                length(np.reshape(Oocyst_death,(N*N,1))),\
                length(np.reshape(Mating,(N*N,1)))\
                length(np.reshape(Zyg_maturation,(N*N,1))),\
                length(np.reshape(Ook_maturation,(N*N,1))),\
                length(np.reshape(Ooc_Bursting,(N*N,1))])]

            #% Calculate vector of transition probabilities %
            #total_transitions = sum(transitions);
            total_transitions = np.sum(transitions)
            trans_probs = transitions/total_transitions

            #% Calculate time of next event.  If time > Tfinal, break!
            stoch_time_step = -np.log(y1)/(total_transitions)
            time1 = stoch_time_step + t(i-1) #% time until next event;
            time2 = t[i-1] + dt

            if min(time1,time2) > Tfinal
                LastIter[j] = i-1
                break

            elif stoch_time_step > dt
                t(i) = time2

            elif stoch_time_step <= dt
                t(i) = time1

                #% ---------------------- %
                #% Determine next event and update population matrices %

                #k = find( X ) returns a vector containing
                #the linear indices of each nonzero element in array X
                #he default for direction is 'first', which finds the first n
                #indices corresponding to nonzero elements.

                #choose_transition = find(np.cumsum(trans_probs)>y2,1,'first');
                choose_transition = np.where(np.cumsum(trans_probs) != 0)[0]

                #TransitionIDs(j,i) = choose_transition;
                TransitionIDs[j,i] = choose_transition

                if choose_transition <= transition_type_indices[2]
                    #% Male gamete death 1
                    CT = choose_transition - transition_type_indices[1]
                    Male_data[CT,i,j] = m[CT] - 1

                elif choose_transition > transition_type_indices[2]
                    && choose_transition <= transition_type_indices[3]
                    #% Female gamete death 2
                    CT = choose_transition - transition_type_indices[2]
                    Female_data[CT,i,j] = f[CT] - 1

                elif choose_transition > transition_type_indices[3]
                    && choose_transition <= transition_type_indices[4]
                    #% Zygote Death 3
                    CT_mat = np.reshape(transition_type_indices(3)+1:\
                                        transition_type_indices(4),N,N)
                    #[ix1,ix2] = find(CT_mat == choose_transition);
                    [ix1,ix2] = where(CT_mat == choose_transition)

                    Zygote_data[ix1,ix2,i,j] = z[ix1,ix2] - 1

                elif choose_transition > transition_type_indices[4] \
                    && choose_transition <= transition_type_indices[5]
                    #% Ookinete Death 4
                    CT_mat = np.reshape(transition_type_indices[4]+1:\
                                        transition_type_indices[5],N,N)
                    #[ix1,ix2] = find(CT_mat == choose_transition);
                    [ix1,ix2] = where(CT_mat == choose_transition)

                    Ookinete_data[ix1,ix2,i,j] = e[ix1,ix2] - 1

                elif choose_transition > transition_type_indices[5] \
                    && choose_transition <= transition_type_indices[6]
                    #% Oocyst Death 5
                    CT_mat = np.reshape(transition_type_indices[5]+1:\
                                        transition_type_indices[6],N,N)
                    #[ix1,ix2] = find(CT_mat == choose_transition);
                    [ix1,ix2] = where(CT_mat == choose_transition)

                    Oocyst_data[ix1,ix2,i,j] = o[ix1,ix2] - 1

                elif choose_transition > transition_type_indices[6] \
                    && choose_transition <= transition_type_indices[7]
                    #% Mating 6
                    CT_mat = np.reshape(transition_type_indices[6]+1:\
                                        transition_type_indices[7],N,N)
                    #[ix1,ix2] = find(CT_mat == choose_transition);
                    [ix1,ix2] = where(CT_mat == choose_transition)

                    Male_data[ix2,i,j] = m[ix2] - 1
                    Female_data[ix1,i,j] = f[ix1] - 1;
                    Zygote_data[ix1,ix2,i,j] = z[ix1,ix2] + 1;

                elif choose_transition > transition_type_indices[7] \
                    && choose_transition <= transition_type_indices[8]
                    #% Zygote maturation 7
                    CT_mat = np.reshape(transition_type_indices[7]+1:\
                                        transition_type_indices(8),N,N);
                    #[ix1,ix2] = find(CT_mat == choose_transition);
                    [ix1,ix2] = where(CT_mat == choose_transition)

                    Zygote_data[ix1,ix2,i,j] = z[ix1,ix2] - 1;
                    Ookinete_data[ix1,ix2,i,j] = e[ix1,ix2] + 1;

                elif choose_transition > transition_type_indices[8] \
                    && choose_transition <= transition_type_indices[9]
                    #% Ookinete maturation 8
                    CT_mat = np.reshape(transition_type_indices[8]+1:\
                                        transition_type_indices[9],N,N)
                    #[ix1,ix2] = find(CT_mat == choose_transition);
                    [ix1,ix2] = where(CT_mat == choose_transition)

                    Ookinete_data[ix1,ix2,i,j] = e[ix1,ix2] - 1
                    Oocyst_data[ix1,ix2,i,j] = o[ix1,ix2] + 1

                elif choose_transition > transition_type_indices[9] &&\
                    choose_transition <= transition_type_indices[10]
                    #% Oocyst bursting 9
                    CT_mat = np.reshape(transition_type_indices[9]+1:\
                                        transition_type_indices[10],N,N);
                    #[ix1,ix2] = find(CT_mat == choose_transition);
                    [ix1,ix2] = find(CT_mat == choose_transition)
                    Oocyst_data[ix1,ix2,i,j] = o[ix1,ix2] - 1
                    Sporozoite_data[ix1,ix2,i,j] = s[ix1,ix2] + my_binornd(poissrnd(n),p)
                    Burst_time[ix1,ix2,i,j] = time1
                    burst_count[ix1,ix2,j] = burst_count[ix1,ix2,j] + 1

    }

    //#Time_data(1:length(t),j) = t;
    //#clear t
    //t.toc()

}

//#%% Uncomment to Overwrite Data in Directory RawData/ContRuptFunc ------------ %
//#filename = 'RawData/ContRuptFunc/nParasiteGroupsData_G0' +
//#    num2str(G0_vec(ii)) + '_max_bias' + num2str(max_bias) +
//#    '_NumStrains' + num2str(N) + '.mat'

//filename1 = 'G://My Drive/mariani_systems/malaria_ctmc_project/' +\
//'MS1ReviewCode\MS1ReviewCode/MM_RawData/'

//#filename2 = 'nParasiteGroupsData_G0' +\
//#num2str(G0_vec(ii)) + '_max_bias' + num2str(max_bias) +\
//#'_NumStrains' + num2str(N) + '.pkl'

//#save(filename,
//#'Time_data',
//#'Male_data',
//#'Female_data',
//#'Zygote_data',
//#'Ookinete_data',
//#'Oocyst_data',
//#'Sporozoite_data',
//#'Burst_time',
//#'burst_count',
//#'TransitionIDs',
//#'SEED')

//import pickle

//data = {"name": "John", "age": 30, "city": "New York"}

//# Save the data to a file
//with open(filename1+"test_output.pkl", "wb") as f:
//    pickle.dump(data, f)

//## Load the data from the file
//#with open("data.pkl", "rb") as f:
//#    loaded_data = pickle.load(f)

return 0;

}

