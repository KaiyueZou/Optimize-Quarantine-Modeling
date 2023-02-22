'''
======================================================================
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

02/27/2021
- simple version of figures

07/19/2021
- formatting the plots

09/25/2021
- format change

09/26/2021
- panel labels

10/19/2021
- adjusted how attack rate is calculated (minus the protected proportion)
======================================================================
'''



### Packages and Quarantine Model Function
import numpy as np 
import matplotlib.pyplot as pl 
from matplotlib.lines import Line2D
from escape_function_function import Escape_Plot



### Matrix
QuarantineLengthMatrix  = [7, 10, 14]
SimulationDays          = 2000
MaxEscapeRateDay14Matrix= np.arange(0.0, 1.1,  0.1)
EscapeFunctionSlopeMatrix=np.arange(0.6, 1.01, 0.1)
penal_labels =  np.array(   [
                            ["(A)", "(B)", "(C)"], 
                            ["(D)", "(E)", "(F)"], 
                            ["(G)", "(H)", "(I)"]
                            ])



### Load in the matrix from outside
EscapeFunctionNumVector = [101, 102, 103]
RProportion             = 0.57*0.796 + 0.09*0.307
for EscapeFunctionNum in EscapeFunctionNumVector:
    FileName1 = str(EscapeFunctionNum) + "SevenDayQuaARMatrix.csv"
    FileName2 = str(EscapeFunctionNum) + "TenDayQuaARMatrix.csv"
    FileName3 = str(EscapeFunctionNum) + "FourteenDayQuaARMatrix.csv"
    SevenDayQuaARMatrix     = np.genfromtxt(FileName1, delimiter = ',') - RProportion
    TenDayQuaARMatrix       = np.genfromtxt(FileName2, delimiter = ',') - RProportion
    FourteenDayQuaARMatrix  = np.genfromtxt(FileName3, delimiter = ',') - RProportion



    ### Compact Figures
    ##Settings
    color_for_7  = 'orange'
    color_for_10 = 'mediumseagreen'
    color_for_14 = 'dodgerblue'

    ##Figure for escape function 101 - 102 - 103
    pl.figure(1)
    pl.rcParams.update({'font.size': 7})
    pl.rcParams["font.family"] = "Times New Roman"

    X_Axis_real  = MaxEscapeRateDay14Matrix
    X_Axis  =  np.around(np.arange(0, 1.1, 0.1)/14, 3)
    ColumnIndex_Dict = { 0: -1, 1: -3, 2: -5 }
    widths  = [2, 2, 2]
    heights = [2, 1, 1]
    gs_kw = dict(width_ratios=widths, height_ratios=heights)
    Fig2, F2_axes = pl.subplots(ncols=3, nrows=3, 
                                sharex=False, sharey=False,
                                # figsize=(7.5, 5), 
                                constrained_layout=True,
                                gridspec_kw=gs_kw)
    # pl.subplots_adjust(wspace=0.20, hspace=0.20)

    for i, row in enumerate(F2_axes):
        for j, cell in enumerate(row):
            # Add panel labels
            if i != 2: 
                cell.set_title( penal_labels[i, j],
                        fontsize='large', loc='left', fontfamily="Times New Roman", fontweight='bold')


            ## row 1 - 2 - 3
            if i == 0:
                #get the values from the spreadsheet
                ColumnIndex = ColumnIndex_Dict[j]
                y7  = SevenDayQuaARMatrix[:, ColumnIndex][::-1]
                y10 = TenDayQuaARMatrix[:, ColumnIndex][::-1]
                y14 = FourteenDayQuaARMatrix[:, ColumnIndex][::-1]
                #get vertical lines
                x14_7  = -1
                x14_10 = -1
                for h in range(len(X_Axis_real)):
                    if y14[h] > y7[h] and y14[h-1] < y7[h-1]:
                        x14_7 = X_Axis_real[h]
                    if y14[h] > y10[h] and y14[h-1] < y10[h-1]:
                        x14_10= X_Axis_real[h]
                #draw the figure and adjust the axis
                if EscapeFunctionNum == 101: 
                    cell.axis([-0.03, 1.03, 0, 25])
                elif EscapeFunctionNum == 102:
                    cell.axis([-0.03, 1.03, 0, 2])
                elif EscapeFunctionNum == 103:
                    cell.axis([-0.03, 1.03, 0, 1.2])
                cell.set_aspect('auto')
                cell.plot(X_Axis_real, y7*100,   color=color_for_7,  lw=2)
                cell.plot(X_Axis_real, y10*100,  color=color_for_10, lw=2)
                cell.plot(X_Axis_real, y14*100,  color=color_for_14, lw=2)
                if x14_10 in X_Axis_real:
                    cell.axvline(x = x14_10, color = 'grey', lw=1.0, linestyle=(0, (3, 5, 1, 5)))
                if x14_7  in X_Axis_real:
                    cell.axvline(x = x14_7,  color = 'grey', lw=1.0, linestyle=(0, (1, 5)))
                LegendLines =   [
                                Line2D([0], [0], color=color_for_14, linewidth=2),
                                Line2D([0], [0], color=color_for_10, linewidth=2),
                                Line2D([0], [0], color=color_for_7, linewidth=2),
                                # Line2D([0], [0], color='grey', linewidth=3.5, linestyle=(0, (3, 5, 1, 5))),
                                # Line2D([0], [0], color='grey', linewidth=3.5, linestyle=(0, (1, 5)))
                                ]
                LegendLabels = [
                                '14-Day Quarantine', '10-Day Quarantine', '7-Day Quarantine',
                                # '14-Day Quarantine\nCross 10-Day Quarantine', '14-Day Quarantine\nCross 7-Day Quarantine'
                                ]
                cell.legend(LegendLines, LegendLabels, prop = {'size': 7}, handlelength = 2, loc='upper left')
                cell.set_xlabel('Compliance function slope\n(rate of increase in the proportion\nleaving quarantine per unit time)')
                cell.set_xticks(X_Axis_real)
                cell.set_xticklabels(X_Axis, rotation = 45, ha='right')
                cell.set_ylabel('Attack Rate(%)')
                # Figure2Title = str( - 10*(ColumnIndex + 1) ) + "% Improvement in Compliance\nDue to Shorter Quarantine"
                # cell.set_title(Figure2Title)

            elif i == 1:
                ColumnIndex = ColumnIndex_Dict[j]
                y7  = SevenDayQuaARMatrix[:, ColumnIndex][::-1]
                y10 = TenDayQuaARMatrix[:, ColumnIndex][::-1]
                y14 = FourteenDayQuaARMatrix[:, ColumnIndex][::-1]
                #get vertical lines
                x14_7  = -1
                x14_10 = -1
                for h in range(len(X_Axis_real)):
                    if y14[h] > y7[h] and y14[h-1] < y7[h-1]:
                        x14_7 = X_Axis_real[h]
                    if y14[h] > y10[h] and y14[h-1] < y10[h-1]:
                        x14_10= X_Axis_real[h]
                # Draw the Figures
                cell.set_aspect('auto')
                if EscapeFunctionNum == 101: 
                    cell.axis([-0.03, 1.03, -0.15, 0.15])
                else:
                    cell.axis([-0.03, 1.03, -0.0025, 0.0025])
                (markers, stemlines, baseline) =  cell.stem(X_Axis_real, y10 - y14)
                if x14_10 in X_Axis_real:
                    cell.axvline(x = x14_10, color = 'grey', lw=1.0, linestyle=(0, (3, 5, 1, 5)))
                pl.setp(stemlines, linestyle="-", color=color_for_10, linewidth=2 )
                pl.setp(markers, color=color_for_10, markersize=4)
                pl.setp(baseline, linestyle="--", color="red", linewidth=0.5)
                # cell.set_xlabel('Maximum -Rate of Leaving Quarantine\n(Per Day)')
                cell.set_xticks(X_Axis_real)
                cell.set_xticklabels(X_Axis, rotation = 45, ha='right')
                if j == 0:
                    cell.set_ylabel('Risk Difference\n(10 vs 14 Day Quarantine)')

            elif i == 2:
                ColumnIndex = ColumnIndex_Dict[j]
                y7  = SevenDayQuaARMatrix[:, ColumnIndex][::-1]
                y10 = TenDayQuaARMatrix[:, ColumnIndex][::-1]
                y14 = FourteenDayQuaARMatrix[:, ColumnIndex][::-1]
                #get vertical lines
                x14_7  = -1
                x14_10 = -1
                for h in range(len(X_Axis_real)):
                    if y14[h] > y7[h] and y14[h-1] < y7[h-1]:
                        x14_7 = X_Axis_real[h]
                    if y14[h] > y10[h] and y14[h-1] < y10[h-1]:
                        x14_10= X_Axis_real[h]
                # Draw the Figures
                cell.set_aspect('auto')
                if EscapeFunctionNum == 101: 
                    cell.axis([-0.03, 1.03, -0.15, 0.15])
                else:
                    cell.axis([-0.03, 1.03, -0.025, 0.025])
                (markers, stemlines, baseline) =  cell.stem(X_Axis_real, y7 - y14)
                if x14_7  in X_Axis_real:
                    cell.axvline(x = x14_7,  color = 'grey', lw=1.0, linestyle=(0, (1, 5)))
                pl.setp(stemlines, linestyle="-", color=color_for_7, linewidth=2 )
                pl.setp(markers, color=color_for_7, markersize=4)
                pl.setp(baseline, linestyle="--", color="red", linewidth=0.5)
                cell.set_xlabel('Compliance function slope\n(rate of increase in the proportion\nleaving quarantine per unit time)')
                cell.set_xticks(X_Axis_real)
                cell.set_xticklabels(X_Axis, rotation = 45, ha='right')
                if j == 0:
                    cell.set_ylabel('Risk Difference\n(7 vs 14 Day Quarantine)')
    
    Figure2Name = str(EscapeFunctionNum) + '-3X3-AreaPlot.png'
    Fig2.set_size_inches(7.5, 6)
    Fig2.savefig(Figure2Name, dpi=300)