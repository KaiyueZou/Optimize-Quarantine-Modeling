'''
Escape Functions Test
- 02/23/2021

- 06/26/2021
Modified the function 102

- 7/30/2021
Added more inputs into the function

'''

### Setup
import math
import copy
import numpy as np 
import matplotlib.pyplot as pl 
from scipy.integrate import solve_ivp



### Make it a function
def Escape_Plot(QuarantineLength, EscapeFunctionNum, MaxEscapeRateDay14, EscapeFunctionSlope):

    QuarantineLength        = QuarantineLength
    EscapeFunctionNum       = EscapeFunctionNum
    MaxEscapeRateDay14      = MaxEscapeRateDay14                #The max escape rate on day 14
    EscapeFunctionSlope     = EscapeFunctionSlope               #1 means the same slope set by MaxEscapeRateDay14

    DelayVariance           = 112/64    #The variance is fixed
    PerfectCompliance       = 0         
    NOQuarantine            = 0
    NewQuarantineEndRate    = 0
    NumQuarantineCmp        = int(round(QuarantineLength**2/DelayVariance)) * 5 ###modified on 8/5/2021, " * 5" is totally unnessacary for the model itself, but important for 10-day quarantine plotting
    MaxEscapeRate           = (MaxEscapeRateDay14*(QuarantineLength/14))*int(QuarantineLength < 14) + MaxEscapeRateDay14*int(QuarantineLength >= 14)  #escape rate = LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp)
    if QuarantineLength == 14:
        LinearScalar    = int(PerfectCompliance == 0)*1*1
    elif QuarantineLength == 10:
        LinearScalar    = int(PerfectCompliance == 0)*1*EscapeFunctionSlope
    elif QuarantineLength == 5:
        LinearScalar    = int(PerfectCompliance == 0)*1*(EscapeFunctionSlope**2)
    
    j_list = np.arange(1, NumQuarantineCmp + 1)
    Y = [0]*(NumQuarantineCmp + 1)

    ### Function
    for j in j_list: 
        
        if EscapeFunctionNum == 101:
            EscapeFunction = 1
        elif EscapeFunctionNum == 102: #
            Slope1 = (    0                                                             )
            Slope2 = (    1                                                             ) ### modified on 7/30/2021 -- (QuarantineLength) / (QuarantineLength - 3)
            ElbowPoint =  (NumQuarantineCmp/QuarantineLength)*(3)

            EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
                             +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )
        elif EscapeFunctionNum == 103: #
            Slope1 = (    0                                                             ) ### modified on 8/5/2021
            Slope2 = (    1                                                             ) ### modified on 8/5/2021
            ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(QuarantineLength-2) 

            EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
                             +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )
        
        Y[j] = EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp)

    ### Return
    X = (np.array(j_list,dtype=float)) / (NumQuarantineCmp/QuarantineLength)
    X = np.insert(X, 0, float(0))
    return X, Y
