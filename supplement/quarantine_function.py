'''
=================================================================================================================================
Quarantine
12/28/2020

01/03/2021
- modified escape rate to accept max parameter, scale based on
  max
- added a bit more plotting code
- test ran quarantine scenarios

01/10/2021
- added number of quarantine compartments as a parameter
- added "close contact" compartment
- separated "escaped people" and "people who completed quarantine"
- same variance across different scenarios with different length of quarantine
- collinear escape rate function
- more scenarios

01/23/2021
- make escape rate plateau after 14 days
- only focus on 7, 10, and 14 days for quarantine
- leave an option to replace QuarantineEndRate with QuarantineRate
- make max excape rate easier to understand by indicating number of people left in quarantine after certain days (escape.py)
- explore risk ratios between 7-day, 10-day, and 14-day quarantine with compliance and in-compliance assumptions (riskratio.py)

02/02/2021
- allow escaped people to isolate after the onset of symptoms
- allow the linear scalar to change across different quarantine lengths

02/21/2021
- change non-compliance function (two)

06/26/2021
- a thorough check
- probability of infection per contact 0.20 -> 0.10
- update escape function
- EscapeIsolationRate     = 1/10 -> 1/2
  ComplyIsolationRate     = 1/10 -> 1/2

06/27/2021
- EscapeIsolationRate     = 1/10 -> 1/5
  ComplyIsolationRate     = 1/10 -> 1/5
  ref: https://academic.oup.com/ofid/article/8/2/ofab023/6104793

7/7/2021
- the EscapeFunctionSlope issue has been solved

7/30/2021
- updated the 102 compliance model function

8/5/2021
- updated the 103 compliance model function

10/18/2021
- updated the parameter values for Delta variant 
=================================================================================================================================
'''

### Packages
import math
import copy
import numpy as np 
import matplotlib.pyplot as pl 
from scipy.integrate import solve_ivp

### OED Model
def quarantine_model(t,x,param):

    ##Parameters
    NumQuarantineCmp        = param[0]
    NumRow                  = param[1] 
    ProbInfection           = param[2]  #probability of infection per contact
    BetaDefault             = param[3]
    BetaOneEsc              = param[4]
    BetaTwoEsc              = param[5]
    LatentRate              = param[6]
    ProAsym                 = param[7]  #proportion of people become asymptomatic after latency
    AsymRecoverRate         = param[8]
    OnsetRate               = param[9]
    SymptomaticRate         = param[10]
    TraceRate               = param[11]
    ClsCntBackRate          = param[12]
    QuarantineRate          = param[13]
    MaxEscapeRate           = param[14] #escape rate = LinearScalar*EscapeFunction*MaxEscapeRate*(j/NumQuarantineCmp)
    LinearScalar            = param[15] 
    QuarantineEndRate       = param[16]
    IsolationRate           = param[17]
    EscapeIsolationRate     = param[18]
    ComplyIsolationRate     = param[19]
    EscapeFunctionNum       = param[20] #0 means default, other numbers mean other functions specified outside of the ODE model
    QuarantineLength        = param[21] #this parameter is global, but just in case

    ##Make state variables into a matrix
    numrow = NumRow
    numcol = 1 + NumQuarantineCmp + 2
    x = x.reshape(numrow, numcol)

    ##Infectious & susceptible population
    InfectiousDefault   = x[6, 0] + x[7, 0] + x[8, 0] + x[9, 0] + x[10, 0] \
                        + x[6, -2] + x[7, -2] + x[8, -2] + x[9, -2] + x[10, -2]
    InfectiousEsc       = x[6, -1] + x[7, -1] + x[8, -1] + x[9, -1] + x[10, -1]
    SusceptibleDefault  = x[0, 0] + x[1, 0] + x[1, -2]
    SusceptibleEsc      = x[1, -1] 

    ##Equations
    dxdt = np.zeros((numrow, numcol))
    for i in range(numrow):
        for j in range(numcol):

            #Create EscapeFunctions
            #EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
            if EscapeFunctionNum == 101:
                EscapeFunction = 1
            elif EscapeFunctionNum == 102 and j != 0: #specify that j can't be 0
                Slope1 = (    0                                                             )
                Slope2 = (    1                                                             ) ### modified on 7/30/2021 -- (QuarantineLength) / (QuarantineLength - 3)
                ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(3)

                EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
                                 +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )
            elif EscapeFunctionNum == 103 and j != 0: #specify that j can't be 0
                Slope1 = (    0                                                             ) ### modified on 8/5/2021
                Slope2 = (    1                                                             ) ### modified on 8/5/2021
                ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(QuarantineLength-2) 

                EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
                                +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )
            else: 
                EscapeFunction = 0 #added on Feb 21th 


            #First row -- susceptible people without recent close contacts 
            if i == 0 and j == 0: 
                dxdt[i, j] = - (BetaDefault/ProbInfection)*x[i, j]*InfectiousDefault \
                             - (BetaOneEsc /ProbInfection)*x[i, j]*InfectiousEsc \
                             + ClsCntBackRate * dxdt[1, 0]
            elif i == 0 and 1 <= j <= NumQuarantineCmp + 2:
                dxdt[i, j] = 0

            #Second row -- susceptible people with recent close contacts
            elif i == 1 and j == 0: 
                dxdt[i, j] = + (BetaDefault/ProbInfection)*(1 - ProbInfection)*x[0,  0]*InfectiousDefault \
                             + (BetaOneEsc /ProbInfection)*(1 - ProbInfection)*x[0,  0]*InfectiousEsc \
                             + (BetaDefault/ProbInfection)*(1 - ProbInfection)*x[1, -2]*InfectiousDefault \
                             + (BetaOneEsc /ProbInfection)*(1 - ProbInfection)*x[1, -2]*InfectiousEsc \
                             + (BetaOneEsc /ProbInfection)*(1 - ProbInfection)*x[1, -1]*InfectiousDefault \
                             + (BetaTwoEsc /ProbInfection)*(1 - ProbInfection)*x[1, -1]*InfectiousEsc \
                             - BetaDefault*x[i, j]*InfectiousDefault \
                             - BetaOneEsc *x[i, j]*InfectiousEsc \
                             - ClsCntBackRate * dxdt[i, j] \
                             - TraceRate * x[i, j]
            elif i == 1 and j == 1: 
                dxdt[i, j] = + TraceRate * x[i, j-1] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j] #escape rate function 
                EscapedS = 0
                EscapedS += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
            elif i == 1 and 2 <= j <= NumQuarantineCmp - 1: 
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                EscapedS += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
            elif i == 1 and j == NumQuarantineCmp: 
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             - QuarantineEndRate * x[i, j]
            elif i == 1 and j == NumQuarantineCmp + 1: #those who waited until the end of the Q
                dxdt[i, j] = + QuarantineEndRate * x[i, j-1] \
                             - (BetaDefault/ProbInfection)*x[i, j]*InfectiousDefault \
                             - (BetaOneEsc/ProbInfection)*x[i, j]*InfectiousEsc
            elif i == 1 and j == NumQuarantineCmp + 2: #those who escaped
                dxdt[i, j] = + EscapedS \
                             - (BetaOneEsc/ProbInfection)*x[i, j]*InfectiousDefault \
                             - (BetaTwoEsc/ProbInfection)*x[i, j]*InfectiousEsc
            
            #Latent rows -- L1
            elif i == 2 and j == 0:
                dxdt[i, j] = + (BetaDefault)*x[0, 0]*InfectiousDefault \
                             + (BetaOneEsc) *x[0, 0]*InfectiousEsc \
                             + (BetaDefault)*x[1, -2]*InfectiousDefault \
                             + (BetaOneEsc) *x[1, -2]*InfectiousEsc \
                             + (BetaOneEsc) *x[1, -1]*InfectiousDefault \
                             + (BetaTwoEsc) *x[1, -1]*InfectiousEsc \
                             + BetaDefault  *x[i-1, j]*InfectiousDefault \
                             + BetaOneEsc   *x[i-1, j]*InfectiousEsc \
                             - LatentRate * x[i, j] \
                             - TraceRate * x[i, j]
            elif i == 2 and j == 1:
                dxdt[i, j] = + TraceRate * x[i, j-1] \
                             - LatentRate * x[i, j] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                EscapedL1 = 0
                EscapedL1 += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
            elif i == 2 and 2 <= j <= NumQuarantineCmp - 1: 
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             - LatentRate * x[i, j] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                EscapedL1 += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
            elif i == 2 and j == NumQuarantineCmp:
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             - LatentRate * x[i, j] \
                             - QuarantineEndRate * x[i, j]
            elif i == 2 and j == NumQuarantineCmp + 1: #those who waited until the end of the Q
                dxdt[i, j] = + QuarantineEndRate * x[i, j-1] \
                             - LatentRate * x[i, j] 
            elif i == 2 and j == NumQuarantineCmp + 2: #those who escaped
                dxdt[i, j] = + EscapedL1 \
                             - LatentRate * x[i, j]

            #Latent rows -- L2 L3 L4
            elif 3 <= i <= 5 and j == 0:
                dxdt[i, j] = + LatentRate * x[i-1, j] \
                             - TraceRate * x[i, j] \
                             - LatentRate * x[i, j]
            elif 3 <= i <= 5 and j == 1:
                dxdt[i, j] = + LatentRate * x[i-1, j] \
                             + TraceRate * x[i, 0] \
                             - LatentRate * x[i, j] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                EscapedL2L3L4 = 0
                EscapedL2L3L4 += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
            elif (3 <= i <= 5 and 2 <= j <= NumQuarantineCmp - 1): 
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             + LatentRate * x[i-1, j] \
                             - LatentRate * x[i, j] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                EscapedL2L3L4 += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]            
            elif 3 <= i <= 5 and j == NumQuarantineCmp:
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             + LatentRate * x[i-1, j] \
                             - LatentRate * x[i, j] \
                             - QuarantineEndRate * x[i, j]
            elif 3 <= i <= 5 and j == NumQuarantineCmp + 1:
                dxdt[i, j] = + QuarantineEndRate * x[i, j-1] \
                             + LatentRate * x[i-1, j] \
                             - LatentRate * x[i, j]
            elif 3 <= i <= 5 and j == NumQuarantineCmp + 2:
                dxdt[i, j] = + EscapedL2L3L4 \
                             + LatentRate * x[i-1, j] \
                             - LatentRate * x[i, j]

            #Asymptomatic
            elif i == 6 and j == 0:
                dxdt[i, j] = + LatentRate * ProAsym * x[i-1, j] \
                             - TraceRate * x[i, j] \
                             - AsymRecoverRate * x[i, j]
                AsymRecovered = 0
                AsymRecovered += AsymRecoverRate * x[i, j]
            elif i == 6 and j == 1:
                dxdt[i, j] = + TraceRate * x[i, 0] \
                             + LatentRate * ProAsym * x[i-1, j] \
                             - AsymRecoverRate * x[i, j] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                EscapedAsym = 0
                EscapedAsym += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                AsymRecovered += AsymRecoverRate * x[i, j]
            elif (i == 6 and 2 <= j <= NumQuarantineCmp - 1): 
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             + LatentRate * ProAsym * x[i-1, j] \
                             - AsymRecoverRate * x[i, j] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                EscapedAsym += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]  
                AsymRecovered += AsymRecoverRate * x[i, j]
            elif i == 6 and j == NumQuarantineCmp:
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             + LatentRate * ProAsym * x[i-1, j] \
                             - AsymRecoverRate * x[i, j] \
                             - QuarantineEndRate * x[i, j]
                AsymRecovered += AsymRecoverRate * x[i, j]
            elif i == 6 and j == NumQuarantineCmp + 1:
                dxdt[i, j] = + QuarantineEndRate * x[i, j-1] \
                             + LatentRate * ProAsym * x[i-1, j] \
                             - AsymRecoverRate * x[i, j]
                AsymRecovered += AsymRecoverRate * x[i, j]
            elif i == 6 and j == NumQuarantineCmp + 2:
                dxdt[i, j] = + EscapedAsym \
                             + LatentRate * ProAsym * x[i-1, j] \
                             - AsymRecoverRate * x[i, j]
                AsymRecovered += AsymRecoverRate * x[i, j]

            #Pre-symptomatic
            elif i == 7 and j == 0:
                dxdt[i, j] = + LatentRate * (1 - ProAsym) * x[i-2, j] \
                             - TraceRate * x[i, j] \
                             - OnsetRate * x[i, j]
            elif i == 7 and j == 1:
                dxdt[i, j] = + TraceRate * x[i, 0] \
                             + LatentRate * (1 - ProAsym) * x[i-2, j] \
                             - OnsetRate * x[i, j] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                EscapedPresym = 0
                EscapedPresym += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                IsolationPeople = 0
                IsolationPeople += OnsetRate * x[i, j]
            elif (i == 7 and 2 <= j <= NumQuarantineCmp - 1): 
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             + LatentRate * (1 - ProAsym) * x[i-2, j] \
                             - OnsetRate * x[i, j] \
                             - QuarantineRate * x[i, j] \
                             - EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j]
                EscapedPresym += EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp) * x[i, j] 
                IsolationPeople += OnsetRate * x[i, j]
            elif i == 7 and j == NumQuarantineCmp:
                dxdt[i, j] = + QuarantineRate * x[i, j-1] \
                             + LatentRate * (1 - ProAsym) * x[i-2, j] \
                             - OnsetRate * x[i, j] \
                             - QuarantineEndRate * x[i, j]
                IsolationPeople += OnsetRate * x[i, j]
            elif i == 7 and j == NumQuarantineCmp + 1:
                dxdt[i, j] = + QuarantineEndRate * x[i, j-1] \
                             + LatentRate * (1 - ProAsym) * x[i-2, j] \
                             - OnsetRate * x[i, j]
            elif i == 7 and j == NumQuarantineCmp + 2:
                dxdt[i, j] = + EscapedPresym \
                             + LatentRate * (1 - ProAsym) * x[i-2, j] \
                             - OnsetRate * x[i, j]

            #Symptomatic 1 2 3 -- only 3 isolation compartments
            elif i == 8 and j == 0:
                dxdt[i, j] = + OnsetRate * x[i-1, 0] \
                             - SymptomaticRate * x[i, j] \
                             - TraceRate * x[i, j]
            elif i == 8 and j == 1:
                dxdt[i, j] = + IsolationPeople \
                             + TraceRate * x[i, j-1] \
                             + ComplyIsolationRate * x[i, (NumQuarantineCmp + 1)] \
                             + EscapeIsolationRate * x[i, (NumQuarantineCmp + 2)] \
                             - SymptomaticRate * x[i, j] \
                             - IsolationRate * x[i, j]
            elif i == 8 and 2 <= j <= NumQuarantineCmp:
                dxdt[i, j] = 0
            elif i == 8 and j == NumQuarantineCmp + 1:
                dxdt[i, j] = + OnsetRate*x[i-1, j] \
                             + IsolationRate * x[i, 1] \
                             - SymptomaticRate * x[i, j] \
                             - ComplyIsolationRate * x[i, j]
            elif i == 8 and j == NumQuarantineCmp + 2: #didn't consider escaped people during isolation 
                dxdt[i, j] = + OnsetRate*x[i-1, j] \
                             - SymptomaticRate * x[i, j] \
                             - EscapeIsolationRate * x[i, j]

            elif 9 <= i <= 10 and j == 0:
                dxdt[i, j] = + SymptomaticRate * x[i-1, j] \
                             - SymptomaticRate * x[i, j] \
                             - TraceRate * x[i, j]
                SymRecovered = 0
                SymRecovered += SymptomaticRate * x[i, j] * int(i == 10)
            elif 9 <= i <= 10 and j == 1:
                dxdt[i, j] = + SymptomaticRate * x[i-1, j] \
                             + TraceRate * x[i, j-1] \
                             + ComplyIsolationRate * x[i, (NumQuarantineCmp + 1)] \
                             + EscapeIsolationRate * x[i, (NumQuarantineCmp + 2)] \
                             - SymptomaticRate * x[i, j] \
                             - IsolationRate * x[i, j]
                SymRecovered += SymptomaticRate * x[i, j] * int(i == 10)
            elif 9 <= i <= 10 and 2 <= j <= NumQuarantineCmp:
                dxdt[i, j] = 0
            elif 9 <= i <= 10 and j == NumQuarantineCmp + 1:
                dxdt[i, j] = + SymptomaticRate*x[i-1, j] \
                             + IsolationRate * x[i, 1] \
                             - SymptomaticRate * x[i, j] \
                             - ComplyIsolationRate * x[i, j]
                SymRecovered += SymptomaticRate * x[i, j] * int(i == 10)
            elif 9 <= i <= 10 and j == NumQuarantineCmp + 2:
                dxdt[i, j] = + SymptomaticRate*x[i-1, j] \
                             - SymptomaticRate * x[i, j] \
                             - EscapeIsolationRate * x[i, j]
                SymRecovered += SymptomaticRate * x[i, j] * int(i == 10)

            #Recovered -- one compartment
            if i == 11 and 0 <= j <= NumQuarantineCmp + 2:
                if j == 0: 
                    dxdt[i, j] = + AsymRecovered + SymRecovered
                elif j != 0:
                    dxdt[i, j] = 0 

    ##Return the equations in one-dimension 
    return dxdt.flatten()



### Running the Model

##Default Parameters
#Parameters to vary 
SimulationDays          = 300
L1Proportion            = 0.001
RProportion             = 0.57*0.796 + 0.09*0.307   #due to vaccination campaign 
MaxEscapeRateDay14      = 1         #The max escape rate on day 14
DelayVariance           = 112/64    #The variance is fixed
QuarantineLength        = 14        #Length of quarantine
PerfectCompliance       = 0         
NOQuarantine            = 0
NewQuarantineEndRate    = 0
EscapeFunctionSlope     = 1         #1 means the same slope set by MaxEscapeRateDay14
EscapeFunctionNum       = 0

#Fixed parameters
NumQuarantineCmp        = int(round(QuarantineLength**2/DelayVariance))
NumRow                  = 12 
ProbInfection           = 0.1693      #probability of infection per contact - modified from 0.2 to 0.1
BetaDefault             = (5.08/12)
BetaOneEsc              = (5.08/12)
BetaTwoEsc              = (5.08/12)
LatentRate              = 4/3.71
ProAsym                 = 0.30      #proportion of people become asymptomatic after latency
AsymRecoverRate         = 1/8
OnsetRate               = 1/2.09
SymptomaticRate         = 3/9.91
TraceRate               = int(NOQuarantine == 0)*(1/1)                      #we assume almost perfect tracing 
ClsCntBackRate          = 1/7           #assumed
QuarantineRate          = (QuarantineLength/DelayVariance)
MaxEscapeRate           = (MaxEscapeRateDay14*(QuarantineLength/14))*int(QuarantineLength < 14) + MaxEscapeRateDay14*int(QuarantineLength >= 14)  #escape rate = EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp)
LinearScalar            = 1 #changed on 7/7/2021, the previous was "int(PerfectCompliance == 0)*1*EscapeFunctionSlope"
QuarantineEndRate       = QuarantineRate*int(NewQuarantineEndRate == 0) + NewQuarantineEndRate              #QuarantineEndRate equals to QuarantineRate
IsolationRate           = 1/10
EscapeIsolationRate     = 1/5              #assume it to be 1/5, the rate at which symptomatic people who escaped from quarantine go tested and begin isolation
ComplyIsolationRate     = 1/5              #also 1/5, the rate at which symptomatic people who perfectly finished quarantine go tested and begin isolation
EscapeFunctionNum       = EscapeFunctionNum #0 means default, other numbers mean other functions specified outside of the ODE model
QuarantineLength        = QuarantineLength  #this parameter is global, but just in case

#Default parameter vector
param_default=  [
                NumQuarantineCmp,
                NumRow,
                ProbInfection,
                BetaDefault,
                BetaOneEsc,
                BetaTwoEsc,
                LatentRate,
                ProAsym,
                AsymRecoverRate,
                OnsetRate,
                SymptomaticRate,
                TraceRate,
                ClsCntBackRate,
                QuarantineRate,
                MaxEscapeRate,
                LinearScalar,
                QuarantineEndRate,
                IsolationRate,
                EscapeIsolationRate,
                ComplyIsolationRate,
                EscapeFunctionNum,
                QuarantineLength
                ]

#Parameter Vector Generator
def parameter_vector(QuarantineLength, MaxEscapeRateDay14, EscapeFunctionSlope, PerfectCompliance = 0, NOQuarantine = 0, NewQuarantineEndRate = 0, EscapeFunctionNum = 0):
    param_define    =   copy.deepcopy(param_default)
    param_define[0] =   int(round(QuarantineLength**2/DelayVariance))
    param_define[-11]=  int(NOQuarantine == 0)*param_define[-11]
    param_define[-9]=       (QuarantineLength/DelayVariance)
    param_define[-8]=       ((MaxEscapeRateDay14*(QuarantineLength/14))*int(QuarantineLength < 14) + MaxEscapeRateDay14*int(QuarantineLength >= 14))
    param_define[-7]=   int(PerfectCompliance == 0)*(   1*int(QuarantineLength == 14)                   + \
                                                        EscapeFunctionSlope*int(QuarantineLength == 10) + \
                                                        EscapeFunctionSlope**2*int(QuarantineLength == 7)  ) #modified on 7/7/2021, deleted unneccessary parts on 8/5/2021
    param_define[-6]=       (param_define[-6]*int(NewQuarantineEndRate == 0) + NewQuarantineEndRate)
    param_define[-2]=       (EscapeFunctionNum)
    param_define[-1]=       (QuarantineLength)
    
    return param_define

# ##Time range
# tSpan   = (0, SimulationDays)
# outT    = np.linspace(0, SimulationDays, (SimulationDays)*10 + 1) #time range

##Scenario Machine Function
def scenario_machine(
    SimulationDays          = 300, 
    L1Proportion            = 0.001, 
    # RProportion             = 0.57*0.796 + 0.09*0.307, #commented out on Feb 21th 
    MaxEscapeRateDay14      = 1,            #The max escape rate on day 14
    DelayVariance           = 112/64,       #The variance is fixed
    QuarantineLength        = 14,           #Length of quarantine
    PerfectCompliance       = 0,         
    NOQuarantine            = 0,
    NewQuarantineEndRate    = 0,
    EscapeFunctionSlope     = 1,
    EscapeFunctionNum       = 0
    ):

    #Generate parameter vector
    NewParam = parameter_vector(QuarantineLength, MaxEscapeRateDay14, EscapeFunctionSlope, 
                                PerfectCompliance, NOQuarantine, NewQuarantineEndRate, EscapeFunctionNum)

    #Generate time range
    tSpan   = (0, SimulationDays)
    outT    = np.linspace(0, SimulationDays, (SimulationDays)*10 + 1)

    #Generate initial states (optional)
    NumQuarantineCmp = NewParam[0]
    x0 = np.zeros(NumRow * (NumQuarantineCmp + 3))
    x0[(NumQuarantineCmp + 3)*2] = L1Proportion        #modified on Feb 21th 2023
    x0[0] = 1 - L1Proportion #modified on Feb 21th 2023

    #Output
    return NewParam, tSpan, outT, x0  



# ##Scenarios to simulate (initial states, parameters)

# # Scenario Default
# NewParam, tSpan, outT, x0 = scenario_machine(
#                                             SimulationDays, 
#                                             L1Proportion,
#                                             MaxEscapeRateDay14,
#                                             DelayVariance, 
#                                             QuarantineLength,
#                                             PerfectCompliance,
#                                             NOQuarantine,
#                                             NewQuarantineEndRate,
#                                             EscapeFunctionSlope
#                                             )

# NumQuarantineCmp = NewParam[0]
# NumRow  = NewParam[1]
# out     = solve_ivp(quarantine_model, tSpan, x0, args=(NewParam,), dense_output=True)
# x       = out.sol(outT).T
# S       = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L       = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym    = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre     = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym     = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R       = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)

# pl.figure(1)
# pl.plot(outT, S,    color='deepskyblue', linestyle='-', label='Susceptible')
# pl.plot(outT, L,    color='gold', linestyle='--', label='Latent')
# pl.plot(outT, Asym, color='orange', linestyle='-.', label='Asymptomatic')
# pl.plot(outT, Pre,  color='tomato', linestyle=':', label='Presymptomatic')
# pl.plot(outT, Sym,  color='maroon', linestyle='-', label='Symptomatic')
# pl.plot(outT, R,    color='olive', linestyle='-', label='Recovered')
# pl.xlabel('Days')
# pl.ylabel('Proportion')
# pl.title('14-Day Quarantine Overall')
# pl.legend(loc='upper right')
# pl.savefig('14_day_quarantine_overall.png')



# #Scenario 1 -- 14 days quarantine
# QuarantineLength = 14
# param1= parameter_vector(param_default, QuarantineLength, MaxEscapeRateDay14)
# NumQuarantineCmp = param1[0]
# x01 = np.zeros(NumRow * (NumQuarantineCmp + 3))
# x01[(NumQuarantineCmp + 3)*2] = L1Proportion            
# x01[0] = 1 - L1Proportion              

# #Scenario 2 -- 14 days quarantine with perfect compliance
# QuarantineLength = 14
# PerfectCompliance= 1
# param2= parameter_vector(param_default, QuarantineLength, MaxEscapeRateDay14, PerfectCompliance)
# NumQuarantineCmp = param2[0]
# x02 = np.zeros(NumRow * (NumQuarantineCmp + 3))
# x02[(NumQuarantineCmp + 3)*2] = L1Proportion    
# x02[0] = 1 - x02[(NumQuarantineCmp + 3)*2]               

# #Scenario 3 -- 0 day quarantine
# QuarantineLength = 14
# PerfectCompliance= 0
# NOQuarantine = 1
# param3= parameter_vector(param_default, QuarantineLength, MaxEscapeRateDay14, PerfectCompliance, NOQuarantine)
# NumQuarantineCmp = param3[0]
# x03 = np.zeros(NumRow * (NumQuarantineCmp + 3))
# x03[(NumQuarantineCmp + 3)*2] = L1Proportion             
# x03[0] = 1 - x03[(NumQuarantineCmp + 3)*2]                

# #Scenario 4 -- 7 days quarantine
# QuarantineLength = 7
# param4= parameter_vector(param_default, QuarantineLength, MaxEscapeRateDay14)
# NumQuarantineCmp = param4[0]
# x04 = np.zeros(NumRow * (NumQuarantineCmp + 3))
# x04[(NumQuarantineCmp + 3)*2] = L1Proportion            
# x04[0] = 1 - x04[(NumQuarantineCmp + 3)*2]                

# #Scenario 5 -- 21 days quarantine
# QuarantineLength = 21
# param5= parameter_vector(param_default, QuarantineLength, MaxEscapeRateDay14)
# NumQuarantineCmp = param5[0]
# x05 = np.zeros(NumRow * (NumQuarantineCmp + 3))
# x05[(NumQuarantineCmp + 3)*2] = L1Proportion                
# x05[0] = 1 - x05[(NumQuarantineCmp + 3)*2]  

# #Scenario 6 -- 3 day quarantine
# QuarantineLength = 3
# param6= parameter_vector(param_default, QuarantineLength, MaxEscapeRateDay14)
# NumQuarantineCmp = param6[0]
# x06 = np.zeros(NumRow * (NumQuarantineCmp + 3))
# x06[(NumQuarantineCmp + 3)*2] = L1Proportion             
# x06[0] = 1 - x06[(NumQuarantineCmp + 3)*2]  

# #Scenario 7 -- 5 day quarantine
# QuarantineLength = 5
# param7= parameter_vector(param_default, QuarantineLength, MaxEscapeRateDay14)
# NumQuarantineCmp = param7[0]
# x07 = np.zeros(NumRow * (NumQuarantineCmp + 3))
# x07[(NumQuarantineCmp + 3)*2] = L1Proportion             
# x07[0] = 1 - x07[(NumQuarantineCmp + 3)*2]  

# #Scenario 8 -- MaxEscapeRateDay14 = 0.5
# QuarantineLength = 14
# MaxEscapeRateDay14 = 0.5
# param8= parameter_vector(param_default, QuarantineLength, MaxEscapeRateDay14)
# NumQuarantineCmp = param8[0]
# x08 = np.zeros(NumRow * (NumQuarantineCmp + 3))
# x08[(NumQuarantineCmp + 3)*2] = L1Proportion             
# x08[0] = 1 - x08[(NumQuarantineCmp + 3)*2]    

# #Scenario 9 -- MaxEscapeRateDay14 = 0.1
# QuarantineLength = 14
# MaxEscapeRateDay14 = 0.1
# param9= parameter_vector(param_default, QuarantineLength, MaxEscapeRateDay14)
# NumQuarantineCmp = param9[0]
# x09 = np.zeros(NumRow * (NumQuarantineCmp + 3))
# x09[(NumQuarantineCmp + 3)*2] = L1Proportion             
# x09[0] = 1 - x09[(NumQuarantineCmp + 3)*2]           



##Simulations

# #Scenario 1 -- 14 days quarantine
# out1    = solve_ivp(quarantine_model, tSpan, x01, args=(param1,), dense_output=True)
# x       = out1.sol(outT).T
# NumQuarantineCmp = param1[0]
# S1      = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L1      = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym1   = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre1    = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym1    = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R1      = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)

# #Scenario 2 -- 14 days quarantine with perfect compliance
# out2    = solve_ivp(quarantine_model, tSpan, x02, args=(param2,), dense_output=True)
# x       = out2.sol(outT).T
# NumQuarantineCmp = param2[0]
# S2      = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L2      = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym2   = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre2    = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym2    = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R2      = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)

# #Scenario 3 -- 0 day quarantine
# out3    = solve_ivp(quarantine_model, tSpan, x03, args=(param3,), dense_output=True)
# x       = out3.sol(outT).T
# NumQuarantineCmp = param3[0]
# S3      = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L3      = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym3   = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre3    = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym3    = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R3      = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)

# #Scenario 4 -- 7 days quarantine
# out4    = solve_ivp(quarantine_model, tSpan, x04, args=(param4,), dense_output=True)
# x       = out4.sol(outT).T
# NumQuarantineCmp = param4[0]
# S4      = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L4      = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym4   = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre4    = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym4    = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R4      = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)

# #Scenario 5 -- 21 days quarantine
# out5    = solve_ivp(quarantine_model, tSpan, x05, args=(param5,), dense_output=True)
# x       = out5.sol(outT).T
# NumQuarantineCmp = param5[0]
# S5      = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L5      = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym5   = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre5    = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym5    = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R5      = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)

# #Scenario 6 -- 3 days quarantine
# out6    = solve_ivp(quarantine_model, tSpan, x06, args=(param6,), dense_output=True)
# x       = out6.sol(outT).T
# NumQuarantineCmp = param6[0]
# S6      = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L6      = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym6   = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre6    = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym6    = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R6      = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)

# #Scenario 7 -- 5 days quarantine
# out7    = solve_ivp(quarantine_model, tSpan, x07, args=(param7,), dense_output=True)
# x       = out7.sol(outT).T
# NumQuarantineCmp = param7[0]
# S7      = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L7      = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym7   = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre7    = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym7    = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R7      = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)

# #Scenario 8 -- MaxEscapeRateDay14 = 0.5
# out8    = solve_ivp(quarantine_model, tSpan, x08, args=(param8,), dense_output=True)
# x       = out8.sol(outT).T
# NumQuarantineCmp = param8[0]
# S8      = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L8      = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym8   = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre8    = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym8    = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R8      = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)

# #Scenario 9 -- MaxEscapeRateDay14 = 0.5
# out9    = solve_ivp(quarantine_model, tSpan, x09, args=(param9,), dense_output=True)
# x       = out9.sol(outT).T
# NumQuarantineCmp = param9[0]
# S9      = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
# L9      = np.sum(x[:, (2*NumQuarantineCmp+6):(6*NumQuarantineCmp+18)], axis=1)
# Asym9   = np.sum(x[:, (6*NumQuarantineCmp+18):(7*NumQuarantineCmp+21)], axis=1)
# Pre9    = np.sum(x[:, (7*NumQuarantineCmp+21):(8*NumQuarantineCmp+24)], axis=1)
# Sym9    = np.sum(x[:, (8*NumQuarantineCmp+24):(11*NumQuarantineCmp+33)], axis=1)
# R9      = np.sum(x[:, (11*NumQuarantineCmp+33):(-1)], axis=1)



##Figures

# #Overall with 14 days quarantine
# pl.figure(1)
# pl.plot(outT, S1,    color='deepskyblue', linestyle='-', label='Susceptible')
# pl.plot(outT, L1,    color='gold', linestyle='--', label='Latent')
# pl.plot(outT, Asym1, color='orange', linestyle='-.', label='Asymptomatic')
# pl.plot(outT, Pre1,  color='tomato', linestyle=':', label='Presymptomatic')
# pl.plot(outT, Sym1,  color='maroon', linestyle='-', label='Symptomatic')
# pl.plot(outT, R1,    color='olive', linestyle='-', label='Recovered')
# pl.xlabel('Days')
# pl.ylabel('Proportion')
# pl.title('14-Day Quarantine Overall')
# pl.legend(loc='upper right')
# pl.savefig('14_day_quarantine.png')

# #Comparisons across different length of quarantine period
# pl.figure(2)
# pl.plot(outT, S3,   color='tomato', linestyle='-', label='0 Day Quarantine')
# pl.plot(outT, S4,   color='gold', linestyle='--', label='7 Days Quarantine')
# pl.plot(outT, S1,   color='darkolivegreen', linestyle='-.', label='14 Days Quarantine')
# pl.plot(outT, S5,   color='deepskyblue', linestyle=':', label='21 Days Quarantine')
# pl.xlabel('Days')
# pl.ylabel('Susceptible People')
# pl.title('Different Length of Quarantine Comparison')
# pl.legend(loc='upper right')
# pl.savefig('0_7_14_21_quarantine.png')

# #Perfect compliance versus incompliance
# pl.figure(3)
# pl.plot(outT, S1,    color='darkred', linestyle='-', label='Incompliance')
# pl.plot(outT, S2,    color='deepskyblue', linestyle=':', label='Perfect Compliance')
# pl.xlabel('Days')
# pl.ylabel('Susceptible People')
# pl.title('Compliance versus Incompliance')
# pl.legend(loc='lower left')
# pl.savefig('compliance_comparison.png')

# #Comparisons across different length of quarantine period
# pl.figure(4)
# pl.plot(outT, S3,   color='tomato', linestyle='-', label='0 Day Quarantine')
# pl.plot(outT, S6,   color='steelblue', linestyle='-', label='3 Day Quarantine')
# pl.plot(outT, S7,   color='pink', linestyle='-', label='5 Day Quarantine')
# pl.plot(outT, S4,   color='gold', linestyle='--', label='7 Days Quarantine')
# pl.plot(outT, S1,   color='darkolivegreen', linestyle='-.', label='14 Days Quarantine')
# pl.plot(outT, S5,   color='deepskyblue', linestyle=':', label='21 Days Quarantine')
# pl.xlabel('Days')
# pl.ylabel('Susceptible People')
# pl.title('Different Length of Quarantine Comparison')
# pl.legend(loc='upper right')
# pl.savefig('0_3_5_7_14_21_quarantine.png')

# #Perfect compliance versus different kinds of incompliance
# pl.figure(5)
# pl.plot(outT, S1,    color='darkred', linestyle='-', label='Incompliance--MaxEscapeRateDay14 = 1')
# pl.plot(outT, S8,    color='crimson', linestyle='--', label='MaxEscapeRateDay14 = 0.5')
# pl.plot(outT, S9,    color='pink', linestyle='-.', label='MaxEscapeRateDay14 = 0.1')
# pl.plot(outT, S2,    color='deepskyblue', linestyle=':', label='Perfect Compliance')
# pl.xlabel('Days')
# pl.ylabel('Susceptible People')
# pl.title('Compliance versus Different Max Escape Rate')
# pl.legend(loc='lower left')
# pl.savefig('maxescaperate_comparison.png')



# ##Diagnosis
# CompartmentTotal = np.sum(x, axis=1)
# # print(a)
# # print(a.shape)
# # print(a == 1.)
# non_ones_list = []
# for i in range(len(CompartmentTotal)):
#     if CompartmentTotal[i] != 1.0:
#         non_ones_list.append([i])
#         non_ones_list.append([CompartmentTotal[i]])
# non_ones_list = np.array(non_ones_list)
# non_ones_list = non_ones_list.reshape(int(len(non_ones_list)/2), 2)
# np.savetxt('non_ones_list.txt', non_ones_list, delimiter=',') 