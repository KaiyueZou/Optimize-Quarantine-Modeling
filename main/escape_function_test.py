'''
========================================================
Escape Functions Test
- 02/23/2021
Created

- 05/04/2021
Should I keep change???

- 6/26/2021
Modify the mistake with function 102

- 8/5/2021
Modified scenario 2 and 3
=========================================================
'''

### Setup
import math
import copy
import numpy as np 
import matplotlib.pyplot as pl 
from scipy.integrate import solve_ivp

QuarantineLength_list = [14, 10, 7]
print(QuarantineLength_list)

for QuarantineLength in QuarantineLength_list: 

    EscapeFunctionNum       = 103
    MaxEscapeRateDay14      = 1        #The max escape rate on day 14
    EscapeFunctionSlope     = 0.6      #1 means the same slope set by MaxEscapeRateDay14

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
    elif QuarantineLength == 7:
        LinearScalar    = int(PerfectCompliance == 0)*1*(EscapeFunctionSlope**2)
    
    ##what's happening here???
    # j_list = range(NumQuarantineCmp)
    j_list = np.arange(1, NumQuarantineCmp + 1)
    Y = [0]*(NumQuarantineCmp + 1)

    ### Function
    for j in j_list: 
        # j += 1
        if EscapeFunctionNum == 101:
            EscapeFunction = 1
        elif EscapeFunctionNum == 102: #
            Slope1 = (    0                                                             )
            Slope2 = (    1                                                             ) ### modified on 8/5/2021 -- (QuarantineLength) / (QuarantineLength - 3)
            ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(3)

            EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                             ) \
                             +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + Slope2 - (1/j) * Slope2 * ElbowPoint )
        elif EscapeFunctionNum == 103: #
            # XTimes = 0.3
            # Slope1 = (    XTimes/LinearScalar                                           )
            # Slope2 = (   (14)/((2)*LinearScalar) \
            #            - (   (QuarantineLength-2)*XTimes     )/(2*LinearScalar)         )
            Slope1 = (    0                                                             ) ### modified on 8/5/2021
            Slope2 = (    1                                                             ) ### modified on 8/5/2021
            ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(QuarantineLength-2) 

            EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
                             +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )

        # ### Not Included
        # elif EscapeFunctionNum == 1: #day 3 elbow, k1 = 0.2 k, k2 = (14*k-3*k1)/(14-3) = (13.4/11)*k
        #     EscapeFunction = (int( j <= (NumQuarantineCmp/QuarantineLength)*3 ) * 0.2) \
        #                     +(int( j >  (NumQuarantineCmp/QuarantineLength)*3 ) * (NumQuarantineCmp/(j)) * 0.2 * (3/QuarantineLength) \
        #                     +(int( j >  (NumQuarantineCmp/QuarantineLength)*3 ) / j) * (j - (NumQuarantineCmp/QuarantineLength)*3) * (13.4/11))
        # elif EscapeFunctionNum == 2: #
        #     XTimes = 0.5
        #     Slope1 = (    XTimes/LinearScalar                                           )
        #     Slope2 = (   (14)/((QuarantineLength-3)*LinearScalar) \
        #                - (3*XTimes        )/((QuarantineLength-3)*LinearScalar)         )
        #     ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(3)

        #     EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
        #                      +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )
        # elif EscapeFunctionNum == 3: #LinearScalar accounted
        #     XTimes = 0.5
        #     Slope1 = (    XTimes/LinearScalar                                           )
        #     Slope2 = (   (14)/((QuarantineLength-3)*1) \
        #                - (3*XTimes        )/((QuarantineLength-3)*LinearScalar)         )
        #     ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(3)

        #     EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
        #                      +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )
        # elif EscapeFunctionNum == 5: #LinearScalar accounted
        #     XTimes = 0.3
        #     Slope1 = (    XTimes/LinearScalar                                           )
        #     Slope2 = (   (14)/((2)*1) \
        #                - (   (QuarantineLength-2)*XTimes     )/(2*LinearScalar)         )
        #     ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(QuarantineLength-2)

        #     EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
        #                      +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )



        # #The rest of this part is just a copy of the other forms of the function, not being used anymore
        # elif EscapeFunctionNum == 6: #day 7 elbow, k1 = 0.2 k, k2 = (14*k-7*k1)/(14-7) = 1.8*k
        #     EscapeFunction = (int( j <= (NumQuarantineCmp/QuarantineLength)*7 ) * 0.2) \
        #                     +(int( j >  (NumQuarantineCmp/QuarantineLength)*7 ) * (NumQuarantineCmp/j) * 0.2 * (7/QuarantineLength) \
        #                     +(int( j >  (NumQuarantineCmp/QuarantineLength)*7 ) / j) * (j - (NumQuarantineCmp/QuarantineLength)*7) * 1.8)
        #             #k1 for first stage, constant for second stage, k2 for second stage 
        # elif EscapeFunctionNum == 7: #2 days before quarantine ends, the same slope (1.5 times the Day14 Slope)
        #     XTimes = 1.5
        #     Slope1 = (   (QuarantineLength)/(QuarantineLength-2) \
        #                - (2*XTimes        )/((QuarantineLength-2)*LinearScalar)         )
        #     Slope2 = (    XTimes/LinearScalar                                           )
        #     ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(QuarantineLength-2)

        #     EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
        #                      +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )
        # elif EscapeFunctionNum == 8: 
        #     XTimes = 1.5
        #     Slope1 = (   (14)/(QuarantineLength-2)                                      \
        #                - (2*XTimes        )/((QuarantineLength-2)*LinearScalar)         )
        #     Slope2 = (    XTimes/LinearScalar                                           )
        #     ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(QuarantineLength-2)

        #     EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
        #                      +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )
        # elif EscapeFunctionNum == 9: 
        #     XTimes = 2.5
        #     Slope1 = (   (14)/(QuarantineLength-2)/LinearScalar                         \
        #                - (2*XTimes        )/((QuarantineLength-2)*LinearScalar)         )
        #     Slope2 = (    XTimes/LinearScalar                                           )
        #     ElbowPoint = (NumQuarantineCmp/QuarantineLength)*(QuarantineLength-2)

        #     EscapeFunction = +  int( j <=  ElbowPoint) * (   Slope1                                                         ) \
        #                      +  int( j > ElbowPoint  ) * (   (1/j) * Slope1 * ElbowPoint + (Slope2 * (j-ElbowPoint) / j)    )
        
         
        Y[j] = EscapeFunction*LinearScalar*MaxEscapeRate*(j/NumQuarantineCmp)

    ### Plot
    X = (np.array(j_list,dtype=float)) / (NumQuarantineCmp/QuarantineLength)
    # print(type(X))
    # print('This is the length of X before insertation: ')
    # print(len(X))
    X = np.insert(X, 0, float(0))
    # print('This is the length of X after insertation: ')
    # print(len(X))
    pl.figure(1)
    pl.plot(X, Y)

    print('The Max Escape Rate in the Plot:')
    print(Y[-1])
    print('The Assigned MaxEscapeRate:')
    print(MaxEscapeRate)
    print(Y)

### Save
pl.savefig('escape_function_curve.png')
