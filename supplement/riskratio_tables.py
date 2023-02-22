'''
===============================================================================================================
Risk Ratio (attack rate ratio) Matrix
1/25/2021

02/02/2021 - 02/07/2021
- make simulation days to 2000 days to get the equilibrium state
- add perfect compliance scenario
- add risk difference 
- consider varying escape rate slope
- change matrix structures
- add sensitivity analysis: varying tracing rate

02/17/2021
- use available data from last time to create curves in Temp.py
    -- with attack rate and risk difference at the same time
- --(not done)-- in-compliance pattern from literature or not 
- actual scalar value for 7-day quarantine
- --(not done)-- indicate how many people are traced and quarantined
- --(not done)-- wider range of trace rate (1/10 - 1/30)

7/7/2021
- the EscapeFunctionSlope issue is not resolved, still mismatches with the quarantine_function.py

8/5/2021
- updated how EscapeFunctionSlope is created

10/18/2021
- updated RProportion for the Delta variant
===============================================================================================================
'''



### Packages and Quarantine Model Function
import math
import copy
import numpy as np 
import matplotlib
import matplotlib.pyplot as pl 
from scipy.integrate import solve_ivp
from quarantine_function import quarantine_model, scenario_machine
# import seaborn as sns

# ### Matrix
EscapeFunctionNum_list = [101, 102, 103]
for EscapeFunctionNum in EscapeFunctionNum_list:

    Timer                   = 0
    QuarantineLengthMatrix  = [7, 10, 14]
    SimulationDays          = 2000
    MaxEscapeRateDay14Matrix= np.arange(0.0, 1.1,  0.1)
    EscapeFunctionSlopeMatrix=np.arange(0.6, 1.01, 0.1)
    # MaxEscapeRateDay14Matrix= np.arange(0.0, 1.1, 0.1) #default
    # EscapeFunctionSlopeMatrix=np.arange(0.7, 1.0, 0.02)#default
    # # SimulationDays          = 500
    # # MaxEscapeRateDay14Matrix= np.arange(0.0, 0.4, 0.1)
    # # EscapeFunctionSlopeMatrix=np.arange(0.7, 1.0, 0.1)
    # # SevenDayQuaARMatrix     = np.random.random((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    # # TenDayQuaARMatrix       = np.random.random((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    # # FourteenDayQuaARMatrix  = np.random.random((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    # # AttackRate7_14Matrix    = np.random.random((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    # # AttackRate10_14Matrix   = np.random.random((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    # # ARDiff7_14Matrix        = np.random.random((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    # # ARDiff10_14Matrix       = np.random.random((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    SevenDayQuaARMatrix     = np.zeros((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    TenDayQuaARMatrix       = np.zeros((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    FourteenDayQuaARMatrix  = np.zeros((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    AttackRate7_14Matrix    = np.zeros((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    AttackRate10_14Matrix   = np.zeros((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    ARDiff7_14Matrix        = np.zeros((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))
    ARDiff10_14Matrix       = np.zeros((len(MaxEscapeRateDay14Matrix), len(EscapeFunctionSlopeMatrix)))

    ### For loop
    for h in range(len(QuarantineLengthMatrix)):
        for i in range(len(MaxEscapeRateDay14Matrix)):
            NewSusceptibleList = []

            for j in range(len(EscapeFunctionSlopeMatrix)):
                SimulationDays          = SimulationDays
                L1Proportion            = 0.001
                RProportion             = 0.57*0.796 + 0.09*0.307
                MaxEscapeRateDay14      = MaxEscapeRateDay14Matrix[i]         
                DelayVariance           = 112/64    
                QuarantineLength        = QuarantineLengthMatrix[h]       
                PerfectCompliance       = 0         
                NOQuarantine            = 0
                NewQuarantineEndRate    = 0
                EscapeFunctionSlope     = EscapeFunctionSlopeMatrix[j] #updated on 8/5/2021 to work with modified quarantine_function
                EscapeFunctionNum       = EscapeFunctionNum #newly added, vary this variable directly

                ##Print the process and parameters
                Timer += 1
                print('Simulating', str(Timer), '/', 
                    str(len(QuarantineLengthMatrix) * len(MaxEscapeRateDay14Matrix) * len(EscapeFunctionSlopeMatrix)), '\n')
                print("SimulationDays", SimulationDays, '\n'
                    "L1Proportion", L1Proportion, '\n'
                    "RProportion", RProportion, '\n'
                    "MaxEscapeRateDay14", MaxEscapeRateDay14, '\n'
                    "DelayVariance", DelayVariance, '\n'
                    "QuarantineLength", QuarantineLength, '\n'
                    "PerfectCompliance", PerfectCompliance, '\n'
                    "NOQuarantine", NOQuarantine, '\n'
                    "NewQuarantineEndRate", NewQuarantineEndRate, '\n'
                    "EscapeFunctionSlope", EscapeFunctionSlope, '\n'
                    "EscapeFunctionNum", EscapeFunctionNum, '\n', '\n')

                ##If 14-day quarantine, skip unnecessary steps
                if QuarantineLength == 14 and j >= 1:
                    pass 
                else: #else, simulate
                    ##Get the Scenario
                    NewParam, tSpan, outT, x0 = scenario_machine(
                                                SimulationDays, 
                                                L1Proportion,
                                                RProportion, 
                                                MaxEscapeRateDay14,
                                                DelayVariance, 
                                                QuarantineLength,
                                                PerfectCompliance,
                                                NOQuarantine,
                                                NewQuarantineEndRate,
                                                EscapeFunctionSlope,
                                                EscapeFunctionNum
                                                )

                    ##Simulate and Save S
                    NumQuarantineCmp = NewParam[0]
                    NumRow  = NewParam[1]
                    out     = solve_ivp(quarantine_model, tSpan, x0, args=(NewParam,), dense_output=True)
                    x       = out.sol(outT).T
                    S       = np.sum(x[:, (0):(2*NumQuarantineCmp+6)], axis=1)
                    NewSusceptibleList.append(S)

                    ##Diagnosis for Equalibrium -- Susceptible Population
                    if (j == len(EscapeFunctionSlopeMatrix) - 1) or (QuarantineLength == 14 and j == 0):
                        pl.figure(1000 + (h + 1)*100 + (i + 1))
                        for S_object in NewSusceptibleList:
                            pl.plot(outT, S_object, color='deepskyblue', linestyle='-')
                        pl.xlabel('Days')
                        pl.ylabel('Susceptible')
                        TitleForFigure = 'Equalibrium Check ' + str(QuarantineLength) + '-Day Q ' + str(np.around(MaxEscapeRateDay14, decimals = 1)) + ' Rate.png'
                        pl.title(TitleForFigure)
                        pl.savefig(TitleForFigure, format = 'png')

                    ##Save Attack Rate
                    AttackRate = 1 - S[-1] 
                    if QuarantineLength == 7: 
                        SevenDayQuaARMatrix[-i-1, j] = AttackRate 
                    elif QuarantineLength == 10: 
                        TenDayQuaARMatrix[-i-1, j] = AttackRate 
                    elif QuarantineLength == 14: 
                        FourteenDayQuaARMatrix[-i-1, :] = AttackRate 

    ### Save the Matrix
    AttackRate7_14Matrix    = SevenDayQuaARMatrix/FourteenDayQuaARMatrix
    AttackRate10_14Matrix   = TenDayQuaARMatrix/FourteenDayQuaARMatrix
    ARDiff7_14Matrix        = SevenDayQuaARMatrix - FourteenDayQuaARMatrix
    ARDiff10_14Matrix       = TenDayQuaARMatrix - FourteenDayQuaARMatrix
    FileName1 = str(EscapeFunctionNum) + "SevenDayQuaARMatrix.csv"
    FileName2 = str(EscapeFunctionNum) + "TenDayQuaARMatrix.csv"
    FileName3 = str(EscapeFunctionNum) + "FourteenDayQuaARMatrix.csv"
    FileName4 = str(EscapeFunctionNum) + "AttackRate7_14Matrix.csv"
    FileName5 = str(EscapeFunctionNum) + "AttackRate10_14Matrix.csv"
    FileName6 = str(EscapeFunctionNum) + "ARDiff7_14Matrix.csv"
    FileName7 = str(EscapeFunctionNum) + "ARDiff10_14Matrix.csv"
    np.savetxt(FileName1, SevenDayQuaARMatrix, delimiter=",")
    np.savetxt(FileName2, TenDayQuaARMatrix, delimiter=",")
    np.savetxt(FileName3, FourteenDayQuaARMatrix, delimiter=",")
    np.savetxt(FileName4, AttackRate7_14Matrix, delimiter=",")
    np.savetxt(FileName5, AttackRate10_14Matrix, delimiter=",")
    np.savetxt(FileName6, ARDiff7_14Matrix, delimiter=",")
    np.savetxt(FileName7, ARDiff10_14Matrix, delimiter=",")





# ### Load in the matrix from outside
# SevenDayQuaARMatrix     = np.genfromtxt('SevenDayQuaARMatrix.csv', delimiter = ',')
# TenDayQuaARMatrix       = np.genfromtxt('TenDayQuaARMatrix.csv', delimiter = ',')
# FourteenDayQuaARMatrix  = np.genfromtxt('FourteenDayQuaARMatrix.csv', delimiter = ',')
# AttackRate7_14Matrix    = np.genfromtxt('AttackRate7_14Matrix.csv', delimiter = ',')
# AttackRate10_14Matrix   = np.genfromtxt('AttackRate10_14Matrix.csv', delimiter = ',')
# ARDiff7_14Matrix        = np.genfromtxt('ARDiff7_14Matrix.csv', delimiter = ',')
# ARDiff10_14Matrix       = np.genfromtxt('ARDiff10_14Matrix.csv', delimiter = ',')



# # ### Area figures
# # pl.figure(1)
# # pl.figure(num=None, figsize=(8, 4), dpi = 200, facecolor='w')


# # ##
# # pl.subplot(2, 3, 1)
# # x   = MaxEscapeRateDay14Matrix
# # y7  = SevenDayQuaARMatrix[:, -1][::-1]
# # y10 = TenDayQuaARMatrix[:, -1][::-1]
# # y14 = FourteenDayQuaARMatrix[:, -1][::-1]
# # color_for_14 = 'lightskyblue'
# # color_for_10 = 'yellowgreen'
# # color_for_7  = 'coral'

# # pl.plot(x, y7,   color=color_for_7,  lw=1.5, linestyle = ':',    label = '7-Day')
# # pl.plot(x, y10,  color=color_for_10, lw=1.5, linestyle = '--',   label = '10-Day')
# # pl.plot(x, y14,  color=color_for_14, lw=1.5, linestyle = '-',    label = '14-Day')
# # pl.fill_between(x, y7, 0,
# #                 facecolor =color_for_7, alpha = 0.2) 
# # pl.fill_between(x, y10, 0, 
# #                 facecolor =color_for_10, alpha = 0.2)
# # pl.fill_between(x, y14, 0, 
# #                 facecolor =color_for_14, alpha = 0.3)
# # pl.axis([0, 1, 0, 1])
# # pl.legend(loc='upper left')

# # ##
# # pl.subplot(1, 2, 2)
# # x   = MaxEscapeRateDay14Matrix
# # y7  = SevenDayQuaARMatrix[:, -11][::-1]
# # y10 = TenDayQuaARMatrix[:, -11][::-1]
# # y14 = FourteenDayQuaARMatrix[:, -11][::-1]
# # color_for_14 = 'lightskyblue'
# # color_for_10 = 'yellowgreen'
# # color_for_7  = 'coral'

# # pl.plot(x, y7,   color=color_for_7,  lw=1.5, linestyle = ':',    label = '7-Day')
# # pl.plot(x, y10,  color=color_for_10, lw=1.5, linestyle = '--',   label = '10-Day')
# # pl.plot(x, y14,  color=color_for_14, lw=1.5, linestyle = '-',    label = '14-Day')
# # pl.fill_between(x, y7, 0,
# #                 facecolor =color_for_7, alpha = 0.2) 
# # pl.fill_between(x, y10, 0, 
# #                 facecolor =color_for_10, alpha = 0.2)
# # pl.fill_between(x, y14, 0, 
# #                 facecolor =color_for_14, alpha = 0.3)
# # pl.axis([0, 1, 0, 1])
# # pl.legend(loc='upper left')
# # pl.savefig('areaplot.png')





# ### Compact
# ### Area figures Try 2
# pl.figure(2)
# pl.figure(num=None, dpi = 300, facecolor='w')
# pl.rcParams.update({'font.size': 25})
# color_for_14 = 'lightskyblue'
# color_for_10 = 'yellowgreen'
# color_for_7  = 'coral'
# x   = MaxEscapeRateDay14Matrix
# scalarDict =    {
#                 -1: '10-Day scalar = 1\n7-Day scalar = 1^2', 
#                 -3: '10-Day scalar = 0.96\n7-Day scalar = 0.96^2', 
#                 -5: '10-Day scalar = 0.92\n7-Day scalar = 0.92^2', 
#                 -7: '10-Day scalar = 0.88\n7-Day scalar = 0.88^2', 
#                 -9: '10-Day scalar = 0.84\n7-Day scalar = 0.84^2', 
#                 -11: '10-Day scalar = 0.80\n7-Day scalar = 0.80^2', 
#                 -13: '10-Day scalar = 0.76\n7-Day scalar = 0.76^2', 
#                 -15: '10-Day scalar = 0.72\n7-Day scalar = 0.72^2', 
#                 }

# widths  = [2, 2, 2, 2]
# heights = [2, 1.2, 0.8, 2, 1.2, 0.8]

# gs_kw = dict(width_ratios=widths, height_ratios=heights)
# fig6, f6_axes = pl.subplots(ncols=4, nrows=6, 
#                             sharex=False, sharey=False,
#                             figsize=(50, 50), 
#                             constrained_layout=True,
#                             gridspec_kw=gs_kw)

# for i, row in enumerate(f6_axes):
#     for j, cell in enumerate(row):
#         ## Area Surface Figures
#         if i in [0, 3]: 
#             # Assign y7 y10 y14
#             if i == 0:
#                 if   j == 0:
#                     ColumnIndex = -1
#                 elif j == 1:
#                     ColumnIndex = -3
#                 elif j == 2:
#                     ColumnIndex = -5
#                 elif j == 3:
#                     ColumnIndex = -7
#             elif i == 3:
#                 if j == 0:
#                     ColumnIndex = -9
#                 elif j == 1:
#                     ColumnIndex = -11
#                 elif j == 2:
#                     ColumnIndex = -13
#                 elif j == 3:
#                     ColumnIndex = -15
#             y7  = SevenDayQuaARMatrix[:, ColumnIndex][::-1]
#             y10 = TenDayQuaARMatrix[:, ColumnIndex][::-1]
#             y14 = FourteenDayQuaARMatrix[:, ColumnIndex][::-1]
#             # Draw the figure
#             cell.axis([-0.02, 1.02, 0, 0.6])
#             cell.set_aspect('auto')
#             label = scalarDict[ColumnIndex]
#             cell.plot([], [], ' ', label = label)
#             cell.plot(x, y7,   color=color_for_7,  lw=5, linestyle = ':',    label = '7-Day')
#             cell.plot(x, y10,  color=color_for_10, lw=5, linestyle = '--',   label = '10-Day')
#             cell.plot(x, y14,  color=color_for_14, lw=5, linestyle = '-',    label = '14-Day')
#             cell.fill_between(x, y7, 0,
#                             facecolor =color_for_7, alpha = 0.2) 
#             cell.fill_between(x, y10, 0, 
#                             facecolor =color_for_10, alpha = 0.2)
#             cell.fill_between(x, y14, 0, 
#                             facecolor =color_for_14, alpha = 0.3)
#             cell.legend(loc='upper left')

#         ## Stem Plots
#         elif i in [1, 4]:
#             # Assign y7 y10 y14
#             if i == 1:
#                 if   j == 0:
#                     ColumnIndex = -1
#                 elif j == 1:
#                     ColumnIndex = -3
#                 elif j == 2:
#                     ColumnIndex = -5
#                 elif j == 3:
#                     ColumnIndex = -7
#             elif i == 4:
#                 if j == 0:
#                     ColumnIndex = -9
#                 elif j == 1:
#                     ColumnIndex = -11
#                 elif j == 2:
#                     ColumnIndex = -13
#                 elif j == 3:
#                     ColumnIndex = -15
#             y7  = SevenDayQuaARMatrix[:, ColumnIndex][::-1]
#             y10 = TenDayQuaARMatrix[:, ColumnIndex][::-1]
#             y14 = FourteenDayQuaARMatrix[:, ColumnIndex][::-1]
#             # Draw the Figures
#             cell.set_aspect('auto')
#             cell.axis([-0.02, 1.02, -0.3, 0.3])
#             (markers, stemlines, baseline) =  cell.stem(x, y7 - y14)
#             pl.setp(stemlines, linestyle="-", color=color_for_7, linewidth=5 )
#             pl.setp(markers, color=color_for_7, markersize=20)
#             pl.setp(baseline, linestyle="--", color="grey", linewidth=2)

#         elif i in [2, 5]:
#             # Assign y7 y10 y14
#             if i == 2:
#                 if   j == 0:
#                     ColumnIndex = -1
#                 elif j == 1:
#                     ColumnIndex = -3
#                 elif j == 2:
#                     ColumnIndex = -5
#                 elif j == 3:
#                     ColumnIndex = -7
#             elif i == 5:
#                 if j == 0:
#                     ColumnIndex = -9
#                 elif j == 1:
#                     ColumnIndex = -11
#                 elif j == 2:
#                     ColumnIndex = -13
#                 elif j == 3:
#                     ColumnIndex = -15
#             y7  = SevenDayQuaARMatrix[:, ColumnIndex][::-1]
#             y10 = TenDayQuaARMatrix[:, ColumnIndex][::-1]
#             y14 = FourteenDayQuaARMatrix[:, ColumnIndex][::-1]
#             # Draw the Figures
#             cell.set_aspect('auto')
#             cell.axis([-0.02, 1.02, -0.2, 0.2])
#             (markers, stemlines, baseline) = cell.stem(x, y10 - y14)
#             pl.setp(stemlines, linestyle="-", color=color_for_10, linewidth=5 )
#             pl.setp(markers, color=color_for_10, markersize=20)
#             pl.setp(baseline, linestyle="--", color="grey", linewidth=2)

#         ## Assign Labels
#         if i == len(f6_axes) - 1:
#             cell.set_xlabel("Escape Rate")
#         if j == 0 and i in [0, 3]:
#             cell.set_ylabel("Attack Rate")
#         elif j == 0 and i in [1, 4]:
#             cell.set_ylabel("7-Day 14-Day Risk Difference")
#         elif j == 0 and i in [2, 5]:
#             cell.set_ylabel("10-Day 14-Day Risk Difference")

# # pl.tight_layout()
# pl.savefig('areaplot_TraceRate_025.png')





# # ### Area figures Try 1

# # ## Setup
# # pl.figure(2)
# # pl.figure(num=None, dpi = 300, facecolor='w')
# # pl.rcParams.update({'font.size': 25})
# # fig, axes2d = pl.subplots(  nrows=6, ncols=4,
# #                             sharex=False, sharey=False,
# #                             figsize=(60, 60), 
# #                             constrained_layout=True)
# # color_for_14 = 'lightskyblue'
# # color_for_10 = 'yellowgreen'
# # color_for_7  = 'coral'
# # x   = MaxEscapeRateDay14Matrix

# # ## Plotting
# # for i, row in enumerate(axes2d):
# #     for j, cell in enumerate(row):

# #         ## Area Surface Figures
# #         if i in [0, 3]: 
# #             # Assign y7 y10 y14
# #             if i == 0:
# #                 if   j == 0:
# #                     ColumnIndex = -1
# #                 elif j == 1:
# #                     ColumnIndex = -3
# #                 elif j == 2:
# #                     ColumnIndex = -5
# #                 elif j == 3:
# #                     ColumnIndex = -7
# #             elif i == 3:
# #                 if j == 0:
# #                     ColumnIndex = -9
# #                 elif j == 1:
# #                     ColumnIndex = -11
# #                 elif j == 2:
# #                     ColumnIndex = -13
# #                 elif j == 3:
# #                     ColumnIndex = -15
# #             y7  = SevenDayQuaARMatrix[:, ColumnIndex][::-1]
# #             y10 = TenDayQuaARMatrix[:, ColumnIndex][::-1]
# #             y14 = FourteenDayQuaARMatrix[:, ColumnIndex][::-1]
# #             # Draw the figure
# #             cell.axis([0, 1, 0, 1])
# #             cell.set_aspect('auto')
# #             cell.plot(x, y7,   color=color_for_7,  lw=3, linestyle = ':',    label = '7-Day')
# #             cell.plot(x, y10,  color=color_for_10, lw=3, linestyle = '--',   label = '10-Day')
# #             cell.plot(x, y14,  color=color_for_14, lw=3, linestyle = '-',    label = '14-Day')
# #             cell.fill_between(x, y7, 0,
# #                             facecolor =color_for_7, alpha = 0.2) 
# #             cell.fill_between(x, y10, 0, 
# #                             facecolor =color_for_10, alpha = 0.2)
# #             cell.fill_between(x, y14, 0, 
# #                             facecolor =color_for_14, alpha = 0.3)
# #             cell.legend(loc='upper left')

# #         ## Stem Plots
# #         elif i in [1, 4]:
# #             # Assign y7 y10 y14
# #             if i == 1:
# #                 if   j == 0:
# #                     ColumnIndex = -1
# #                 elif j == 1:
# #                     ColumnIndex = -3
# #                 elif j == 2:
# #                     ColumnIndex = -5
# #                 elif j == 3:
# #                     ColumnIndex = -7
# #             elif i == 4:
# #                 if j == 0:
# #                     ColumnIndex = -9
# #                 elif j == 1:
# #                     ColumnIndex = -11
# #                 elif j == 2:
# #                     ColumnIndex = -13
# #                 elif j == 3:
# #                     ColumnIndex = -15
# #             y7  = SevenDayQuaARMatrix[:, ColumnIndex][::-1]
# #             y10 = TenDayQuaARMatrix[:, ColumnIndex][::-1]
# #             y14 = FourteenDayQuaARMatrix[:, ColumnIndex][::-1]
# #             # Draw the Figures
# #             cell.set_aspect('auto')
# #             cell.axis([0, 1, -0.3, 0.3])
# #             (markers, stemlines, baseline) =  cell.stem(x, y7 - y14)
# #             pl.setp(stemlines, linestyle="-", color=color_for_7, linewidth=5 )
# #             pl.setp(markers, color=color_for_7, markersize=20)
# #             pl.setp(baseline, linestyle="--", color="grey", linewidth=2)

# #         elif i in [2, 5]:
# #             # Assign y7 y10 y14
# #             if i == 2:
# #                 if   j == 0:
# #                     ColumnIndex = -1
# #                 elif j == 1:
# #                     ColumnIndex = -3
# #                 elif j == 2:
# #                     ColumnIndex = -5
# #                 elif j == 3:
# #                     ColumnIndex = -7
# #             elif i == 5:
# #                 if j == 0:
# #                     ColumnIndex = -9
# #                 elif j == 1:
# #                     ColumnIndex = -11
# #                 elif j == 2:
# #                     ColumnIndex = -13
# #                 elif j == 3:
# #                     ColumnIndex = -15
# #             y7  = SevenDayQuaARMatrix[:, ColumnIndex][::-1]
# #             y10 = TenDayQuaARMatrix[:, ColumnIndex][::-1]
# #             y14 = FourteenDayQuaARMatrix[:, ColumnIndex][::-1]
# #             # Draw the Figures
# #             cell.set_aspect('auto')
# #             cell.axis([0, 1, -0.2, 0.2])
# #             (markers, stemlines, baseline) = cell.stem(x, y10 - y14)
# #             pl.setp(stemlines, linestyle="-", color=color_for_10, linewidth=5 )
# #             pl.setp(markers, color=color_for_10, markersize=20)
# #             pl.setp(baseline, linestyle="--", color="grey", linewidth=2)

# #         ## Assign Labels
# #         if i == len(axes2d) - 1:
# #             cell.set_xlabel("Escape Rate")
# #         if j == 0 and i in [0, 3]:
# #             cell.set_ylabel("Attack Rate")
# #         elif j == 0 and i in [1, 4]:
# #             cell.set_ylabel("7-Day 14-Day Risk Difference")
# #         elif j == 0 and i in [2, 5]:
# #             cell.set_ylabel("10-Day 14-Day Risk Difference")

# # # pl.tight_layout()
# # pl.savefig('areaplot3.png')





# ### Heatmaps
# AllMatrix = [
#             SevenDayQuaARMatrix,
#             TenDayQuaARMatrix,
#             FourteenDayQuaARMatrix,
#             AttackRate7_14Matrix,
#             AttackRate10_14Matrix,
#             ARDiff7_14Matrix,
#             ARDiff10_14Matrix
#             ]
# TitleMatrix=[
#             'Attack Rate with 7-Day Quarantine',
#             'Attack Rate with 10-Day Quarantine',
#             'Attack Rate with 14-Day Quarantine',
#             'Risk Ratio between 7-Day Quarantine and 14-Day Quarantine', 
#             'Risk Ratio between 10-Day Quarantine and 14-Day Quarantine',
#             'Risk Difference between 7-Day Quarantine and 14-Day Quarantine', 
#             'Risk Difference between 10-Day Quarantine and 14-Day Quarantine'
#             ]

# for k in range(len(AllMatrix)):
#     MatrixForHeatmap = AllMatrix[k]

#     fig, ax = pl.subplots()
#     if k <= 2:
#         im = ax.imshow(MatrixForHeatmap, cmap = 'autumn_r')
#     elif k <= 4:
#         im = ax.imshow(MatrixForHeatmap, cmap = 'Blues')
#     elif k <= 6:
#         im = ax.imshow(MatrixForHeatmap, cmap = 'summer')
#     # We want to show all ticks...
#     ax.set_xticks(np.arange(MatrixForHeatmap.shape[1]))
#     ax.set_yticks(np.arange(MatrixForHeatmap.shape[0]))
#     # ... and label them with the respective list entries
#     ax.set_xticklabels(np.around(EscapeFunctionSlopeMatrix, decimals = 2))
#     ax.set_yticklabels(np.flipud(np.around(MaxEscapeRateDay14Matrix, decimals = 2)))

#     # Rotate the tick labels and set their alignment.
#     pl.setp(ax.get_xticklabels(), rotation=45, ha="right",
#             rotation_mode="anchor")

#     # Loop over data dimensions and create text annotations.
#     for i in range(MatrixForHeatmap.shape[0]):
#         for j in range(MatrixForHeatmap.shape[1]):
#             if k <= 2:
#                 text = ax.text(j, i, round(MatrixForHeatmap[i, j], 3), 
#                         ha="center", va="center", color="k", fontsize = 6.9)
#             elif k <= 4:
#                 text = ax.text(j, i, round(MatrixForHeatmap[i, j], 3), 
#                         ha="center", va="center", color="darkgreen", fontsize = 6.9)
#             elif k <= 6:
#                 text = ax.text(j, i, round(MatrixForHeatmap[i, j], 3), 
#                         ha="center", va="center", color="k", fontsize = 6.9)
#     TitleForFigure = TitleMatrix[k]
#     ax.set_title(TitleForFigure)
#     pl.xlabel('Escape Rate scalar')
#     pl.ylabel('Forced Original Max Escape Rate on Day 14th')
#     TitleForFigure += '.png'
#     pl.savefig(TitleForFigure, format = 'png')

# ### Another Set of Heatmaps
# for m in range(len(EscapeFunctionSlopeMatrix)):
#     AttackRateMatrix71014   = np.zeros((len(MaxEscapeRateDay14Matrix), 3))
#     RiskRatio71014          = np.zeros((len(MaxEscapeRateDay14Matrix), 2))
#     AttackRateDiff71014     = np.zeros((len(MaxEscapeRateDay14Matrix), 2))

#     AttackRateMatrix71014[:, 0] = SevenDayQuaARMatrix[:, m]  
#     AttackRateMatrix71014[:, 1] = TenDayQuaARMatrix[:, m]    
#     AttackRateMatrix71014[:, 2] = FourteenDayQuaARMatrix[:, m]  
#     RiskRatio71014[:, 0]        = AttackRate7_14Matrix[:, m] 
#     RiskRatio71014[:, 1]        = AttackRate10_14Matrix[:, m] 
#     AttackRateDiff71014[:, 0]   = ARDiff7_14Matrix[:, m]
#     AttackRateDiff71014[:, 1]   = ARDiff10_14Matrix[:, m]

#     AllMatrix2  =   [
#                     AttackRateMatrix71014, 
#                     RiskRatio71014, 
#                     AttackRateDiff71014
#                     ]
#     for n in range(len(AllMatrix2)):
#         MatrixForHeatmap = AllMatrix2[n]

#         fig, ax = pl.subplots()
#         if n == 0:
#             im = ax.imshow(MatrixForHeatmap, cmap = 'autumn_r')
#         elif n == 1:
#             im = ax.imshow(MatrixForHeatmap, cmap = 'Blues')
#         elif n == 2:
#             im = ax.imshow(MatrixForHeatmap, cmap = 'summer')
#         # We want to show all ticks...
#         ax.set_xticks(np.arange(MatrixForHeatmap.shape[1]))
#         ax.set_yticks(np.arange(MatrixForHeatmap.shape[0]))
#         # ... and label them with the respective list entries
#         if n == 0:
#             ax.set_xticklabels(['7-Day', '10-Day', '14-Day'])
#         elif n == 1:
#             ax.set_xticklabels(['7:10', '10:14'])
#         elif n == 2:
#             ax.set_xticklabels(['7-10', '10-14'])
#         ax.set_yticklabels(np.flipud(np.around(MaxEscapeRateDay14Matrix, decimals = 2)))

#         # Rotate the tick labels and set their alignment.
#         pl.setp(ax.get_xticklabels(), rotation=45, ha="right",
#                 rotation_mode="anchor")

#         # Loop over data dimensions and create text annotations.
#         for i in range(MatrixForHeatmap.shape[0]):
#             for j in range(MatrixForHeatmap.shape[1]):
#                 if n == 0:
#                     text = ax.text(j, i, round(MatrixForHeatmap[i, j], 3), 
#                             ha="center", va="center", color="k", fontsize = 6.9)
#                 elif n == 1:
#                     text = ax.text(j, i, round(MatrixForHeatmap[i, j], 3), 
#                             ha="center", va="center", color="darkgreen", fontsize = 6.9)
#                 elif n == 2:
#                     text = ax.text(j, i, round(MatrixForHeatmap[i, j], 3), 
#                             ha="center", va="center", color="k", fontsize = 6.9)
#         TitleForFigure = 'EscSlop=' + str(np.around(EscapeFunctionSlopeMatrix[m], decimals = 3)) + '-' + str(n + 1)
#         pl.ylabel('Forced Original Max Escape Rate on Day 14th')
#         TitleForFigure = '/mnt/raid1/kybryant/1.10 Model Code/More Figures/' + TitleForFigure + '.png'
#         pl.savefig(TitleForFigure, format = 'png')



### Diagnosis
# print(SevenDayQuaARMatrix)
# print(TenDayQuaARMatrix)
# print(FourteenDayQuaARMatrix)
# print(AttackRate7_14Matrix)
# print(AttackRate10_14Matrix)