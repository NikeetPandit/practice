
import DataFormatCallFunctions as DFF
import DataFormat2 as DF2
import OptimalEstimatorProject1_Routines as OE

DataProd = [DFF.GNI1B_posvel()] #Data Products to Fit
ID = ['C']                      #GRACE-FO identifier
import numpy as np


Hours = 5                                       #Hours of data to fit
Hours_Extrap = 28.5                               #Hours of data to extrapolate
[Ind,Exp] = OE.SelectDataProd(DataProd,ID,Hours)    #Formating data based on selection 
OE.PlotData(Ind,Exp)                                #Plotting so user can select geometry to fit
OptionToSet = 'Ellipse'                             #User selection of geometry (Line / Ellipse/ Polynomial)
P = [0,1] #The first entry is to use weight and the second entry is the weight type (0 for no and 1 for yes)
[Fit,Type,Omega,Equations,Residuals,Variances,ChiSqr,Fit_Uncertainty] = OE.LSFittoData(Ind,Exp,P,ID,OptionToSet,Hours) #Fits data, Plots, Error Analysis
[Residuals1,Variances,ChiSqr,Fit_Uncertainty] = OE.ExtrapolateAndPlot(DataProd,Ind,ID,Fit,Type,Omega,Hours,Hours_Extrap)#Extrapolates data, Plots, Error Analysis
Variances = np.trunc(Variances)
ChiSqr = np.trunc(ChiSqr)
Fit_Uncertainty = np.trunc(Fit_Uncertainty)