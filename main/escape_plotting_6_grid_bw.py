'''
==============================================
New Formatted Escape Figures
7/15/2021

- 9/25/2021
  Changed the labels and captions

- 9/26/2021
  Added panel labels 

==============================================
'''

### Packages and Quarantine Model Function
import math
import copy
import numpy as np 
import matplotlib
import matplotlib.pyplot as pl 
from matplotlib.lines import Line2D
from numpy.core.fromnumeric import size
from scipy.integrate import solve_ivp
from quarantine_function import quarantine_model, scenario_machine
from escape_function_function import Escape_Plot
import matplotlib.patches as patches

### Pre-defined Vectors
EscapeFunctionMatrix    = [101, 102, 103]
MaxEscapeRateDay14Matrix= [0.11, 0.5] ### modified on 7/30/2021, only the first element matters
QuarantineLengthMatrix  = [5, 10, 14]
EscapeFunctionSlope     =  0.6

### Column 1 Data
X7list  = []
X10list = []
X14list = []
Y7list  = []
Y10list = []
Y14list = []
for EscapeFunction in EscapeFunctionMatrix:
    # if EscapeFunction == 102:
    #     MaxEscapeRateDay14 = MaxEscapeRateDay14Matrix[0]*11/14 ### modified on 7/30/2021
    MaxEscapeRateDay14 = MaxEscapeRateDay14Matrix[0]
    Escape_X7, Escape_Y7    = Escape_Plot(5,  EscapeFunction, MaxEscapeRateDay14, EscapeFunctionSlope) ### modified on 7/30/2021
    Escape_X10, Escape_Y10  = Escape_Plot(10, EscapeFunction, MaxEscapeRateDay14, EscapeFunctionSlope)
    Escape_X14, Escape_Y14  = Escape_Plot(14, EscapeFunction, MaxEscapeRateDay14, EscapeFunctionSlope)
    X7list.append(Escape_X7)
    X10list.append(Escape_X10)
    X14list.append(Escape_X14)
    Y7list.append(Escape_Y7)
    Y10list.append(Escape_Y10)
    Y14list.append(Escape_Y14)

### Column 2 Data
ResultArray = np.zeros((3, 2, 3, 141))
for h in range(len(EscapeFunctionMatrix)):
    for i in range(len(MaxEscapeRateDay14Matrix)): 
        for j in range(len(QuarantineLengthMatrix)):
            SimulationDays          = QuarantineLengthMatrix[j]
            L1Proportion            = 0
            MaxEscapeRateDay14      = MaxEscapeRateDay14Matrix[i]         
            DelayVariance           = 112/64    
            QuarantineLength        = QuarantineLengthMatrix[j]       
            PerfectCompliance       = 0         
            NOQuarantine            = 0
            NewQuarantineEndRate    = 0
            EscapeFunctionSlope     = EscapeFunctionSlope ###6/30/2021 update -- so that it corresponds to the slope in compliance model
            EscapeFunctionNum       = EscapeFunctionMatrix[h] 
            # if EscapeFunction == 102:
            #     MaxEscapeRateDay14 = MaxEscapeRateDay14Matrix[0]*(SimulationDays-3)/SimulationDays  ### modified on 7/30/2021

            NewParam, tSpan, outT, x0 = scenario_machine(
                                        SimulationDays, 
                                        L1Proportion,
                                        MaxEscapeRateDay14,
                                        DelayVariance, 
                                        QuarantineLength,
                                        PerfectCompliance,
                                        NOQuarantine,
                                        NewQuarantineEndRate,
                                        EscapeFunctionSlope,
                                        EscapeFunctionNum
                                                        )

            ## we need to reset the parameter values and initial states
            NumQuarantineCmp= NewParam[0]
            x0[0] = 0
            x0[2*NumQuarantineCmp + 7] = 1
            NewParam[6] = 0 #set the latent rate to be 0

            ## simulate
            out             = solve_ivp(quarantine_model, tSpan, x0, args=(NewParam,), dense_output=True)
            x               = out.sol(outT).T
            LeftPro         = np.sum(x[:, (2*NumQuarantineCmp + 7):(3*NumQuarantineCmp+8)], axis=1)
            ResultArray[h, i, j, 0:len(LeftPro)] = LeftPro 



### 9-Grid Figure

## color configuration
# color_for_7  = 'orange'
# color_for_10 = 'mediumseagreen'
# color_for_14 = 'dodgerblue'

## line type configuration 
ls_for_7 = 'solid'
ls_for_10 = 'dashed'
ls_for_14 = (0, (1, 1))

## dimension and spacing configuration
pl.figure(1)
pl.rcParams.update({'font.size': 9})
pl.rcParams["font.family"] = "Times New Roman"

horizontal_size = 7.5 #inch
vertical_size   = 8.5 #inch
widths      = [1, 1]
heights     = [1, 1, 1]
gs_kw       = dict(width_ratios=widths, height_ratios=heights)
fig, axs    = pl.subplots(  ncols=2, nrows=3, 
                            sharex=False, sharey=False, 
                            constrained_layout=False,
                            gridspec_kw=gs_kw   )
pl.subplots_adjust(left=0.10, bottom=0.055, right=0.98, top=0.97, wspace=0.25, hspace=0.25)
pl.xticks(np.arange(0, 15, 1))

## plotting
penal_labels =  np.array([
                ["(A)", "(B)"], 
                ["(C)", "(D)"], 
                ["(E)", "(F)"]
                        ])

for i, row in enumerate(axs):
    for j, cell in enumerate(row):
        # Add panel labels
        cell.set_title( penal_labels[i, j],
                        fontsize='large', loc='left', fontfamily="Times New Roman", fontweight='bold')
        

        # the first column
        if j == 0:
            if i != 2: 
                cell.axis([0, 14.1, -0.005, 0.115]) ### modified on 7/30/2021 -- cell.axis([0, 14.1, -0.01, 1.01])
                cell.set_xticks(np.arange(0, 15, 1))
                cell.set_yticks(np.arange(0, 0.15, 0.05))
                cell.set_aspect('auto')
            else: 
                cell.axis([0, 14.1, -0.005, 0.115]) ### modified on 7/30/2021 -- cell.axis([0, 14.1, -0.01, 1.01])
                cell.set_xticks(np.arange(0, 15, 1))
                cell.set_yticks(np.arange(0, 0.15, 0.05))
                cell.set_aspect('auto')

            cell.plot(X7list[i],  Y7list[i],   linestyle=ls_for_7,  color="black", alpha = 1.0, lw=2.0)
            cell.plot(X10list[i], Y10list[i],  linestyle=ls_for_10, color="black", alpha = 0.7, lw=2.0)
            cell.plot(X14list[i], Y14list[i],  linestyle=ls_for_14, color="black", alpha = 0.5, lw=2.0)

            LegendLines =   [
                            Line2D([0], [0], linestyle=ls_for_14, color="black", alpha = 0.5, linewidth=2),
                            Line2D([0], [0], linestyle=ls_for_10, color="black", alpha = 0.7, linewidth=2),
                            Line2D([0], [0], linestyle=ls_for_7, color="black", alpha = 1.0, linewidth=2),
                            ]
            LegendLabels = [
                            '14-Day Quarantine', '10-Day Quarantine', '5-Day Quarantine'
                            ]
            cell.legend(LegendLines, LegendLabels, prop = {'size': 9}, handlelength = 2, loc='upper left')
            
            SharedYLabel =  r'$ (\frac{persons}{persons\times day}) $'
            pl.figtext(0.035, 0.5, 'Rate of Leaving Quarantine ' + SharedYLabel, ha='center', va='center', rotation='vertical', fontfamily = "Times New Roman", fontsize = 12)
            if i == 2:
                cell.set_xlabel('Quarantine Length (Days)')
        
                
        # the other two columns
        elif j == 1:
            if i != 2:
                cell.axis([0, 14.1, -0.01*100, 1.01*100])
                cell.set_xticks(np.arange(0, 15, 1))
                cell.set_aspect('auto')
            else:
                cell.axis([0, 14.1, 0.9495*100, 1.0005*100])
                cell.set_xticks(np.arange(0, 15, 1))
                cell.set_aspect('auto')

            QuarantineLength= QuarantineLengthMatrix[0]
            cell.plot(outT[0:(QuarantineLength*10 + 1)], ResultArray[i, j-1, 0, 0:(QuarantineLength*10 + 1)]*100, 
            linestyle=ls_for_7,  color="black", alpha = 1.0, lw=2.0)
            QuarantineLength= QuarantineLengthMatrix[1]
            cell.plot(outT[0:(QuarantineLength*10 + 1)], ResultArray[i, j-1, 1, 0:(QuarantineLength*10 + 1)]*100, 
            linestyle=ls_for_10, color="black", alpha = 0.7, lw=2.0)
            QuarantineLength= QuarantineLengthMatrix[2]
            cell.plot(outT[0:(QuarantineLength*10 + 1)], ResultArray[i, j-1, 2, 0:(QuarantineLength*10 + 1)]*100, 
            linestyle=ls_for_14, color="black", alpha = 0.5, lw=2.0)

            LegendLines =   [
                            Line2D([0], [0], linestyle=ls_for_14, color="black", alpha = 0.5, linewidth=2),
                            Line2D([0], [0], linestyle=ls_for_10, color="black", alpha = 0.7, linewidth=2),
                            Line2D([0], [0], linestyle=ls_for_7,  color="black", alpha = 1.0, linewidth=2),
                            ]
            LegendLabels = [
                            '14-Day Quarantine', '10-Day Quarantine', '5-Day Quarantine'
                            ]
            if i == 0 and j == 1:
                cell.legend(LegendLines, LegendLabels, prop = {'size': 9}, handlelength = 2, loc='lower left')
            else:
                cell.legend(LegendLines, LegendLabels, prop = {'size': 9}, handlelength = 2, loc='lower left')
            
            pl.figtext(0.535, 0.5, 'Proportion of Initial Population still Remaining in Quarantine (%)', ha='center', va='center', rotation='vertical', fontfamily = "Times New Roman", fontsize = 12)
            if i == 2:
                cell.set_xlabel('Quarantine Length (Days)')



## save
fig.set_size_inches(horizontal_size, vertical_size) 
fig.savefig('6-grid-figure.png', dpi=300, )