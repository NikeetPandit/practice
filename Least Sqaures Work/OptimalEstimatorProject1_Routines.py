''''
Optimal Estimator Project 1 - Fitting line/curve/surface by using Least Squares
by Nikeet Pandit
******************** ALL UTILITIES ARE GIVEN SUMMARIES DESCRIBING FUNCTIONALITY ******************** 

This module comprises of a series of executables that allows you to fit a geometry to a dataset using LS and WLS
(1) line, curve (polynomial-quadractic), surface (plane), and sinusoid (or ellipse in 3D)
(2) user must select this my setting OptionToSet variable which is -> Linear, Polynomial, Ellipse
(3) all descriptions of routines are provided so that can easily be manipulated 
(4) functions will work (may require some manipulation) for any dataset, but it is optimized to work with GRACE-FO Ephemeris files
(5) module calls LSSA to estimate fundamental frequency of sine waves written by Professor Ebrahim Ghaderpour
(6) module calls DataFormCallFunctions and DataFormat2 which read in the GRACE-FO data files provided by JPL which I have previously written
(7) the plane fitting option is not included as an OptionToSet because it is not necessary for this project, but can easily be if needed
(7) for a further application. It was added to show that it could easily be done for my own learning. 
(8) user must specifiy the DataProd using the DataFormatCallFunctions, the ID (GRACE-FO identifier), the Hours of the dataset to fit
(8) and the weighted scheme which is a [1,2] list. The first element denotes if weighted is used (1) or not (0). The second element 
(9) denotes the type of weighting scheme to be used. For example to not weight you would specify P as P = [0,0]. To weight
(9) and specify weighting scheme 3... you would specify P as P = [1,3]. More info is given in SUMMARY OF WeightedObsTransfer

'''

# ------- IMPORTING LIBRARIES ------- #
import numpy as np
import DataFormatCallFunctions as DFF
import DataFormat2 as DF2
import matplotlib.pyplot as plt
import LSSA as LS
# ------- IMPORTING LIBRARIES ------- #


#----SUMMARY OF SelectDataProd(DataProd,ID,Hours)----
#SELECTS DATA PRODUCTS TO FIT
#SELECTS ID OF GRACE SATELLITE (IDENTIFIER) TO CONSTRUCT DATAFRAME PROVIDED A GIVEN DATE-RANGE IN INPUTS.PY SUBMODULE
#SELECTS HOURS FOR FURTHER MANIPULATION TO GIVE USER CAPABABILITY TO DICTATE DURATION TO WHICH DATA IS TO BE FITTED TO [IN HOURS]... CAN BE MANIPULATED EASILY
#RETURNS IND/EXP VARIABLES WHICH ARE INPUTS TO MANY OF THESE FUNCTIONS 
def SelectDataProd(DataProd,ID,Hours):
    Index = int(Hours*3600-1) 
    DataBase = list(map(DF2.CreateDataFrame,DataProd,ID))
    Pos = DataBase[0].iloc[0:Index,1:4].to_numpy()
    Time = DataBase[0].iloc[:Index,0].to_numpy()
    Time = Time[:]-Time[0]
    return Time,Pos


#----SUMMARY of PlotData(DataProd,ID,Hours)----
#THIS FUNCTION SERVES AS A PLOTTING ROUTINE TO PLOT X,Y,Z POSITION DATA IN A SUBPLOT IN RESPECTIVE DIMENSIONS
#ALSO IT WILL PROVIDE ANOTHER PLOT IN 3-DIMENSIONS
#SINCE LS CURVE FITTING REQUIRES KNOWLEDGE OF GEOMETRY (OR GENERAL SURFACE) WHICH IS TO FIT TO DATA...
#THIS IS AN INCREDIBLY HELPFUL UTILITY TO DECIDE WHICH OF THE LINE/PLANE/ELLIPSE/POLYNOMIAL IS SUFFICIENT TO FIT DATA
#OR IF MORE GEOMETRY IS NEEDED TO FIT
#THIS GIVES USER ABILITY TO EASILY MANIPULATE FUNCTIONS IN THIS CODE TO FIT CURVE/LINE/SURFACE TO ANY DATASET BY MANIPULATING THESE BLOCKS OF CODE
def PlotData(Ind,Exp):
    fig, axs = plt.subplots(3)
    axs[0].plot(Ind,Exp[:,0])
    axs[0].set_title('GNSS Data in X')
    axs[1].plot(Ind,Exp[:,1],'--')
    axs[1].set_title('GNSS POD Data in Y')
    axs[2].plot(Ind,Exp[:,2],'--')
    axs[2].set_title('GNSS POD Data in Z')
    for ax in axs.flat:
          ax.set(xlabel='Time [s]', ylabel = 'Position in ECI [m]')
    for ax in axs.flat:
        ax.label_outer()
    fig.show()
    fig1 = plt.figure()
    ax1 = plt.axes(projection='3d')
    ax1.plot(Exp[:,0],Exp[:,1],Exp[:,2],label = 'LS FIT')
    ax1.set_title('GNSS POD Data in 3D [ECI Coord]')
    ax1.legend()
    ax1.set(xlabel='Data in X [m]',ylabel='Data in Y [m]',zlabel='Data in Z [m]')
    fig1.show()
    return print("Examine Geometry and Select Curve to Fit...\n\n")


#----SUMMARY OF LSFitToData(Ind,Exp,P,ID,OptiontoSet)----
#Ind -> Independent Variable, Exp->Explanatory Variable, P-> Weighted LS Method Identifier, ID->Exp. Above, OptiontoSet->Geometry Selection for LS Fit
#IND/EXP ARE OUTPUTTED OUT OF SelectDataProd(DataProd,ID,Hours) FUNCTION
#P VAR (IDENTIFIER FOR METHOD OF WLS REGRESSION) ARE DISCUSSED IN SelectCurveTypeAndFit(Ind,Exp,P,ID,OptiontoSet) 
#OptiontoSet ALLOWS FOR: (1) ELLIPSE [USES SINUSOIDAL BASE FUNCTIONS [COS AND SINE] + CONSTANT]... FUNDAMENTAL FREQUENCY PARAMETER IS EST. USING LSSA IN NESTED FUNC
#(2) LINEAR (3) POLYNOMIAL - USES QUADRACTIC BUT CAN EXTEND EASILY TO GENERAL CASE (4) PLANE FUNCTION IS FitPlaneType(Ind,Exp) AND WOULD NEED TO BE ADDED TO BE USED AS AN OPTION
#IT IS NOT USED AN OPTION CURRENTLY BECAUSE IT IS NOT NEEDED FOR THIS APPLICATION OF FITTING DATA TO POD DATA FOR GRACE-FO
#ANY OTHER FUNCTION (LINEAR) CAN EASILY BE ADDED FOLLOWING THE METHOD OF USING THE NORMAL EQUATIONS AS NEEDED
#-----
#IN FUNCTION -> CALLS NESTED FUNCTION SelectCurveTypeAndFit(Ind,Exp,P,ID,OptiontoSet) WHICH DETERMINES COEFFICIENTS ESTIMATED FROM LS PROCESS
#IT ALSO WILL RETURN OMEGA IF APPLICABLE (FOR ELLIPSE) AND WILL RETURN 999 OTHERWISE
#IT ALSO RETURNS THE TYPE (ELLIPSE/LINE/POLYNOMIAL) TO CONDITIONALLY CONSTRUCT EQUATIONS REQURIED FOR PLOTTING
#ALSO RESIDUALS ARE CALCULATED AND PLOTTED. DATA SEPERATED INTO X,Y,Z COMPONENTS AND PLOTTED IN 3D ARE GIVEN WITH RESIDUALS AND LS BEST FITS
#ALL PLOTS ARE LABELLED AS NECESSARY FOR THE APPLICATION FOR FITTING LINE/CURVE/SURFACE TO GRACE-FO DATA IN ECI COORD.
def LSFittoData(Ind,Exp,P,ID,OptiontoSet,Hours):
    #[Ind,Exp] = SelectDataProd(DataProd,ID,Hours)
    [Fit,Type,Omega] = SelectCurveTypeAndFit(Ind,Exp,P,ID,Hours,OptiontoSet)
    if Type == 'Ellipse':
        Eq0 = Fit[0,0] + Fit[1,0]*np.cos(2*np.pi*Omega*Ind) + Fit[2,0]*np.sin(2*np.pi*Omega*Ind)
        Eq1 = Fit[0,1] + Fit[1,1]*np.cos(2*np.pi*Omega*Ind) + Fit[2,1]*np.sin(2*np.pi*Omega*Ind)
        Eq2 = Fit[0,2] + Fit[1,2]*np.cos(2*np.pi*Omega*Ind) + Fit[2,2]*np.sin(2*np.pi*Omega*Ind)

        Eq0_Residual = Exp[:,0]-Eq0
        Eq1_Residaul = Exp[:,1]-Eq1
        Eq2_Residual = Exp[:,2]-Eq2

    elif Type == 'Linear':
        Eq0 = Fit[0,0]+Fit[1,0]*Ind
        Eq1 = Fit[0,1]+Fit[1,1]*Ind
        Eq2 = Fit[0,2]+Fit[1,2]*Ind

        Eq0_Residual = Exp[:,0]-Eq0
        Eq1_Residaul = Exp[:,1]-Eq1
        Eq2_Residual = Exp[:,2]-Eq2

    elif Type == 'Polynomial':
        Eq0 = Fit[0,0]+Fit[1,0]*Ind+Fit[2,0]*Ind**2
        Eq1 = Fit[0,1]+Fit[1,1]*Ind+Fit[2,1]*Ind**2
        Eq2 = Fit[0,2]+Fit[1,2]*Ind+Fit[2,2]*Ind**2

        Eq0_Residual = Exp[:,0]-Eq0
        Eq1_Residaul = Exp[:,1]-Eq1
        Eq2_Residual = Exp[:,2]-Eq2

    fig, axs = plt.subplots(3)
    axs[0].plot(Ind,Exp[:,0],'r--',linewidth = 3,label='POD DATA')
    axs[0].plot(Ind,Eq0,'b',label='LS FIT')
    axs[0].set_title('GNSS Data in X with LS Fit')
    axs[1].plot(Ind,Exp[:,1],'r--',linewidth = 3,label='POD DATA')
    axs[1].plot(Ind,Eq1,'b',label='LS FIT')
    axs[1].set_title('GNSS POD Data in Y with LS Fit')
    axs[2].plot(Ind,Exp[:,2],'r--',linewidth = 3,label='POD DATA')
    axs[2].plot(Ind,Eq2,'b',label='LS FIT')
    axs[2].set_title('GNSS POD Data in Z with LS Fit')
    axs[1].set(ylabel = 'Position in ECI [m]')

    axs[2].set(xlabel = 'Time [s]')
    for ax in axs.flat:
        ax.label_outer()
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    fig.show()
    fig1 = plt.figure()
    ax1 = plt.axes(projection='3d')
    ax1.plot(Exp[:,0],Exp[:,1],Exp[:,2],'r--',linewidth = 3,label='POD DATA')
    ax1.plot(Eq0,Eq1,Eq2,'b',label='LS FIT')
    ax1.set_title('GNSS POD Data in 3D')
    ax1.set(xlabel='Data in X [m]',ylabel='Data in Y [m]',zlabel='Data in Z [m]')
    ax1.legend()
    fig1.show()

    fig2, axs2 = plt.subplots(3)
    axs2[0].plot(Ind,Eq0_Residual)
    axs2[0].set_title('GNSS POD in X LS Fit Residual')
    axs2[1].plot(Ind,Eq1_Residaul)
    axs2[1].set_title('GNSS POD Data in Y LS Fit Residual')
    axs2[2].plot(Ind,Eq2_Residual)
    axs2[2].set_title('GNSS POD Data in Z LS Fit Residual')
    for ax in axs2.flat:
          ax.set(xlabel='Time [s]', ylabel = 'Position in ECI [m]')
    for ax in axs2.flat:
        ax.label_outer()
    fig2.show()

    Equations = [Eq0,Eq1,Eq2]
    Residuals = [Eq0_Residual,Eq1_Residaul,Eq2_Residual]

    Variances = [SigmaExp(Eq0_Residual,Type),SigmaExp(Eq1_Residaul,Type),SigmaExp(Eq2_Residual,Type)] #calculting variance 
    Index,DataBase = int(Hours*3600-1), list(map(DF2.CreateDataFrame,[DFF.GNV1B_all()],ID))
    Pos_ERR = ((DataBase[0].iloc[:Index,4:7].to_numpy()))
    ChiSqr = ChiSqrFun(Residuals,Pos_ERR,Type) #calculating Chisqr
    Fit_Uncertainty =  PropErrorToFit(Fit,Pos_ERR,Residuals) #calculating uncertainty in coefficients of fit
    
    return Fit,Type,Omega,Equations,Residuals,Variances,ChiSqr,Fit_Uncertainty


#----SUMMMARY OF SelectCurveTypeAndFit(Ind,Exp,P,ID,OptiontoSet)----
#P VAR (IDENTIFIER FOR METHOD OF WLS REGRESSION) REQUIRES TO ELEMENTS IN A LIST
#THE FIRST ELEMENT IN THE LIST DENOTES IF WLS IS USED AT ALL 
#THE SECOND ELEMENT DENOTES THE TYPE OF WLS TO BE USED (EXPERIMENTAL). 
#**THE TYPES AND METHODS/REASONINGS BEHIND THE UTILITIES ARE GIVEN IN: WeightedObsTransfer(ID,Hours,Exp,P,OptionToSet)**
#THIS UTILITY CALLS THE WeightedObsTransfer FUNCTION TO UPDATE (EQUALIZE) THE SERIES TO BE AN EQUIVALENT SERIES OF UNIT WEIGHT
#BASED ON THE OptiontoSet TYPE IT CONSTRUCTS THE NORMAL EQUATIONS AND SOLVES FOR THE COEFFICIENTS TO BE USED TO CONSTRUCT THE EQUATIONS FOR PLOTTING/ RESIDUALS ETC.
#ALSO IN THIS FUNCTION... THE LSSA IS CALLED TO ESTIMATE THE FUNDAMENTAL FREQUENCY WHICH IS NOT A PARAMETER THAT IS TO BE ESTIMATED FROM LS ALONE
#BANDWIDTH FOR LSSA IS CHOSEN BASED ON VISUALLY INSPECTING PLOTS IN X,Y,Z AND WOULD NEED TO BE UPDATED IF A DIFFERENT DATASET IS USED
#LSSA IS NOT ESTIMATED WITH WEIGHTS ALTHOUGH IT WOULD BE INTERESTING TO EXPERIMENT WITH

def SelectCurveTypeAndFit(Ind,Exp,P,ID,Hours,OptionToSet):
    OptiontoSet = OptionToSet
    if P[0] == 1: #If want weight LS... equalize observation series to be equivalent obs series of unit weight
        Exp1 = WeightedObsTransfer(ID,Hours,Exp,P,OptionToSet,Ind)
        Exp = Exp1
    if OptiontoSet == 'Polynomial':
        Col1,Col2,Col3 = np.repeat(1,len(Ind)),Ind,Ind*Ind
        A = np.stack([Col1,Col2,Col3],axis=1)
        c = np.linalg.inv(A.T @ A) @ A.T @ Exp #normal equation
        return c,OptiontoSet,999
    elif OptiontoSet == 'Linear':
        Col1,Col2 = np.repeat(1,len(Ind)),Ind
        A = np.stack([Col1,Col2],axis=1)
        c = np.linalg.inv(A.T @ A) @ A.T @ Exp
        return c,OptiontoSet,999
    elif OptiontoSet == 'Ellipse':
        Omega = 1/(np.linspace(10000,500,500)) #Setting BW for LSSA 
        maxSpectrum = []
        for x in range(0,3):
            maxSpectrum.append(np.argmax(LS.LSSA(Ind,Exp[:,x],P=1,Omega=Omega,ind=[],level=0.01,trend='linear')[0]))
        w = np.mean([x for x in Omega[maxSpectrum]])
        Col1,Col2,Col3 = np.repeat(1,len(Ind)),np.cos(Ind*2*np.pi*w),np.sin(Ind*2*np.pi*w)
        A = np.stack([Col1,Col2,Col3],axis=1)
        c = np.linalg.inv(A.T @ A) @ A.T @ Exp
        return c,OptiontoSet,w
    else:
        print("Valid Option Not Selected. Try again.")
        SelectCurveTypeAndFit()


#----SUMMMARY OF FitPlaneType(Ind,Exp)----
#FITS PLANE TO DATASET USING LS. IS NOT USED FOR THIS PROJECT, HOWEVER, WAS USING IT FOR TESTING
#SPECIFICALLY METHODS ONLINE FOR FITTING A GEOMETRY IN 3D THAT IS IMPLICITLY DEFINED IN 2D (EXAMPLE A LINE)
#WOULD USE METHOD OF SVD TO FIT A PLANE TO THE DATASET AND THEN PROJECT ALL POINTS ONTO THIS PLANE
#THEN FIT THE GEOMETRY USING LS THEN PROJECT BACK AFTER THE PARAMETERS HAVE BEEN ESTIMATED
#WAS INTERESTED TO SEE ABOUT USING LS TO FIT THE PLANE TO THE DATASET AND ALTHOUGH THE FUNCTION DOES WORK (TESTED)
#AND COMPARE IT TO THE METHOD OF MULTIPLE REGRESSION USED IN THIS CODE
#FIRSTLY FOUND IT DIFFICULT TO DETERMINE HOW TO PROJECT ALL THE POINTS TO THE PLANE BUT ALSO THIS METHOD SEEMS THAT IT WOULD INTRODUCE ADDITIONAL ERRORS
#THEN THE METHOD PRESENTED HERE - SPECIFICALLY THERE WILL BE ERROR IN FIT OF PLANE AND THEN ERROR IN THE SUBSEQUENT ESTIMATION OF THE PARAMETERS ON THIS PLANE
#WHICH SEEMS HEURISTICALLY TO BE WORSE OFF THEN FITTING AT ONE TIME USING MULTIPLE REGRESSION
def FitPlaneType(Ind,Exp): ## LS Solution to Fitting a Plane to Data
    Col1,Col2,Col3 = Exp[:,0],Exp[:,1],np.repeat(1,len(Ind))
    A = np.stack([Col1,Col2,Col3],axis=1)
    B = Exp[:,2].T
    c = np.linalg.inv(A.T @ A) @ A.T @ B
    return c


#----SUMMARY OF WeightedObsTransfer(ID,Hours,Exp,P,OptionToSet)----
#FOR THE SELECTED TIMEFRAME OF POSITION DATA OF GRACE-FO THIS FUNCTION FIRST CREATES A MATRIX OF THE ASSOCIATED FORMAL ERRORS
#THE FORMAL ERRORS ARE A THEORETICAL ESTIMATE BASED ON MODELLING WHICH ESTIMATES THE STANDARD DEVIATION OF THE MEASUREMENT
#IN LINEAR REGRESSION WITH ONE EXPLANATORY AND ONE INDEPENDANT VARIABLE... THE WEIGHT MATRIX SIMPLY DENOTES THE WEIGHTING OF EACH INDIVIDUAL MATRIX
#IN MULTIPLE LINEAR REGRESSION... USING A WEIGHTED MATRIX WHICH ALSO INCLUDES CROSS CORRELATION BETWEEN THE MEASUREMENT ERRORS FOR EACH OF THE DEPENDENT VARIABLES
#BECOMES A COMPLEX TASK. A PAPER ON THE SUBJECT IS WRITTEN BY B. LI WHERE THEY DISCUSS AND COME UP WITH AN ALGORITHM TO: PERFORM MULTIPLE LINEAR REGRESSION
#WITH CORRELATED EXPLANATORY VARIABLES AND RESPONSES. THIS METHOD SEEMED BEYOND THE SCOPE OF THIS PROJECT AND 6 DIFFERENT METHODS WERE TESTED EXPERIMENTALLY 
#WITH DECISIONS DRIVEN BASED ON LEARNED COURSE KNOWLEDGE. https://doi.org/10.1179/1752270615Y.0000000006 [PAPER LINK]

#OPTION 0 -> INVERSE OF FORMAL ERROR IS TAKEN AS WEIGHT MATRIX AND THEN MULTIPLIED BY ASSOCIATED EXPLANATORY VARIABLES TO EQUALIZE OBS TO WEIGHT 1...
    #VARIANCE OF UNIT FACTOR IS GIVEN VALUE OF 1
#OPTION 1 -> SAME AS OPTION 0, EXCEPT FOR Z WHOSE MEASUREMENTS ARE GIVEN CONSISTENT WEIGHT OF 1 (NOT ADJUSTED). THIS WAS EXPERIMENTED WITH BECAUSE THE WLS
    #WAS NOT WORKING PROPERLY FOR OPTION 0. ALSO, WHEN PLOTTING A HISTOGRAM OF THE MEASUREMENT ERRORS FOR X,Y,Z... X AND Y FOLLOWED CLOSELY A NORMAL DIST. 
    #BUT ERROR Z DID NOT. AS LS ASSUMES MEASUREMENT ERRORS ARE NORMALLY DIST. I IMAGINED PERHAPS THIS WAS CAUSING PROBLEMS IN THE WLS... SO I SET IT TO 1 FOR Z.
    #THIS WAS ALSO NOT YIELDING A SENSICAL SOLUTION OR APPROXIMATING ANY GEOMETRY BY FITTING THE GEOMETRY AS A BEST FIT TO THE DATASET
#OPTION 2 -> I ESTIMATED THE UNIT VARIANCE FACTOR K. USING THE FORMAL IN THE SLIDES. IN THE FORMULA, I ASSUME THE WEIGHTS ARE ALL GIVEN A VALUE OF 1. 
    #THIS WAS DONE BECAUSE THIS METHOD TRANSFERS THE WEIGHT TO THE OBSERVATIONS TO MAKE AN ASSOCIATED WEIGHT MATRIX OF 1
#OPTION 3 -> REPEATED OPTION 2 BUT WITHOUT INCORPERATING MEASUREMENT ERRORS FOR Z AND USING MEASUREMENTS AS IS WITHOUT ANY WEIGHTING TRANSFER
#OPTION 4 -> INSTEAD OF ASSUMING ERRORS ARE FORMAL ERRORS GIVEN IN THE DATA-FILES... ERRORS ARE ESTIMATED ARE SAID TO BE EQUAL TO THE RESIDUALS
    #FROM THE LS FIT. WITH THESE ERRORS... THE VARIANCE FACTORS ARE THEN ESTIMATED FOR THE EXPLANATORY VARIABLES. THEN THE TRANSFERED WEIGHTED
    #LS IS PERFORMED
#OPTION 5 -> I TRY TO INDUCE ARTIFICAL CORRELATION TO THE MEASUREMENT ERRORS TO SEE THE EFFECT IN THE WLS. AS I MENTIONED ABOVE, MULTIPLE REGRESSION 
    #WEIGHTED LS WITH CORRELATED EXPLANATORY VARIABLES, TO MY KNOWLEDGE, SEEMS A BIT ADVANCED FOR NOW. I REALLY DO NOT BELIEVE THIS METHOD 
    #HAS ANY MERIT, BUT IT WAS WORTH THE TRY. ESSENTIALLY I TRIED INTRODUCING THE CORRELATION OF MEASUREMENT ERRORS BY ESTIMATING THE VARIANCE FACTOR
    #FOR EACH POSITION X, Y, Z AS A WEIGHTED SUM OF EACH OTHER. WHEN LOOKING AT THE GRAPHS FOR ERRORS, I NOTICED THAT AS Y ERRORS GO UP Z ERRORS GO DOWN.
    #AS X GOES UP, Y GOES DOWN AND X GOES UP Z GOES UP. SO IN THE WEIGHTED AVERAGE OF THE THREE ERRORS I MULT EITHER 0.5 FACTOR OR -0.5 FACTOR DEP
    #IF IT THE ERROR WOULD GO UP OR DOWN. IN PRACTICE THE WEIGHTED LEAST SQUARES DID NOT WORK WELL IN THE X OR Y, BUT SOMEHOW THE WEIGHTED LEAST
    #SQUARES IN THE Z WORKED QUITE WELL FOLLOWING THIS METHOD. PERHAPS, AN ITTERATIVE SOLUTION COULD HAVE BEEN DONE IN THIS MANNER. PERHAPS
    #IT WORKED WELL IN THE Z BECAUSE THE ERRORS IN THE Z DID NOT FOLLOW A NORMAL DISTRIBUTION VERY WELL AND WEIGHTING THEM WITH X AND Y ERRORS
    #SKEWED THE ERROR TO FOLLOW A MORE NORMAL DISTRIBUTION. CHANGING THE VARIANCE FACTOR FOR THE X POSITION MEASUREMENTS FOR EXAMPLE WOULD 
    #AFFECT THE FIT IN THE Y AND THE Z, WHICH WAS UNEXPECTED BUT MAKES INTUITIVE SENSE.

def WeightedObsTransfer(ID,Hours,Exp,P,OptionToSet,Ind): #Determining Weight Matrix 
    Index,DataBase = int(Hours*3600-1), list(map(DF2.CreateDataFrame,[DFF.GNV1B_all()],ID))
    Pos_ERR = ((DataBase[0].iloc[:Index,4:7].to_numpy()))
    if P[1] == 0: #Est of Var Unit Weight is Taken as 1 and divided by SD to give Weight Mat
            P_Posx, P_posy, P_posz = 1/Pos_ERR[:,0],1/Pos_ERR[:,1],1/Pos_ERR[:,2] #Var Unit Weight == 1
    elif P[1] == 1: #Est of Var Unit Weight is Taken as 1 and divided by SD to give Weight Mat Except for Z
            P_Posx, P_posy, P_posz = 1/Pos_ERR[:,0],1/Pos_ERR[:,1],np.repeat(1,len(Pos_ERR[:,2])) #Var Unit Weight == 1 and Non-weighting to Z
    elif P[1] == 2: #Est of Var Unit Weight is estimated by sigma_hat = sqrt(p*delta**2)/n where the delta (measurement error) is assumed to equal the SD
            EstVarOfUnitWeight = EstVarUnitWeight(ID,Ind,Exp,P,OptionToSet,Pos_ERR,0,Hours)
            P_Posx, P_posy, P_posz = EstVarOfUnitWeight[0]/Pos_ERR[:,0],EstVarOfUnitWeight[1]/Pos_ERR[:,1],EstVarOfUnitWeight[2]/Pos_ERR[:,2] 
    elif P[1] == 3: #Only weight x, y (follow normal dist with ~0 mean)
            EstVarOfUnitWeight = EstVarUnitWeight(ID,Ind,Exp,P,OptionToSet,Pos_ERR,0,Hours)
            P_Posx, P_posy, P_posz = EstVarOfUnitWeight[0]/Pos_ERR[:,0],EstVarOfUnitWeight[1]/Pos_ERR[:,1],np.repeat(1,len(Pos_ERR[:,2])) 
    elif P[1] == 4: #Errors are residuals from LS best fit without weights. Then Estimate Variance Factor
            EstVarOfUnitWeight = EstVarUnitWeight(ID,Ind,Exp,P,OptionToSet,Pos_ERR,1,Hours)
            P_Posx, P_posy, P_posz = EstVarOfUnitWeight[0]/Pos_ERR[:,0],EstVarOfUnitWeight[1]/Pos_ERR[:,1],EstVarOfUnitWeight[2]/Pos_ERR[:,2] 
    elif P[1] == 5: #Do Weighted LS With POSx full error then POSy half error then POSz half error ->VarOfUnitWeight for X... then for Y... then for Z... then do weighted LS
                    #With these variance factors
            #Fixing X error and introducing artificial correlation
            P_Posx = (Pos_ERR[:,0]+-0.5*Pos_ERR[:,1]+Pos_ERR[:,2]*0.5)/3
            Pos_ERR = np.stack([P_Posx,P_Posx,P_Posx],axis=1)
            VarEst1 = EstVarUnitWeight(ID,Ind,Exp,P,OptionToSet,Pos_ERR,0,Hours)[0]

            #Fixing Y error and introducing artifical correlation 
            P_Posx = (Pos_ERR[:,1]+-0.5*Pos_ERR[:,0]+Pos_ERR[:,2]*0.5)/3
            Pos_ERR = np.stack([P_Posx,P_Posx,P_Posx],axis=1)           
            VarEst2 = EstVarUnitWeight(ID,Ind,Pos_ERR,P,OptionToSet,Pos_ERR,0,Hours)[1]

            P_Posx = (Pos_ERR[:,2]+0.5*Pos_ERR[:,1]+Pos_ERR[:,0]*0.5)/3
            Pos_ERR = np.stack([P_Posx,P_Posx,P_Posx],axis=1)          
            VarEst3 = EstVarUnitWeight(ID,Ind,Pos_ERR,P,OptionToSet,Pos_ERR,0,Hours)[2]         

            #Using Individual Variance of Unit Estimate Weight
            P_Posx, P_posy, P_posz = VarEst1/Pos_ERR[:,0],VarEst2/Pos_ERR[:,1],VarEst3/Pos_ERR[:,2]

    else:
        P_Posx,P_posy,P_posz = np.repmat(1,len(Ind)),np.repmat(1,len(Ind)),np.repmat(1,len(Ind))

    P_Pos = np.stack([P_Posx, P_posy, P_posz],axis=1)
    Xweighted,Yweighted,Zweighted = Exp[:,0]*P_Pos[:,0], Exp[:,1]*P_Pos[:,1], Exp[:,2]*P_Pos[:,2] #turns weight matrix to 1
    PosWeighted = np.stack([Xweighted,Yweighted,Zweighted],axis=1)
    return PosWeighted

#----SUMMARY OF EstVarUnitWeight(ID,Ind,Exp,P,OptionToSet,Pos_ERR,Weight_Itt)----
#THIS FUNCTION ESTIMATES THE ESTIMATED VARIANCE OF UNIT WEIGHT
#IF THE WEIGHT_ITT IS EQUAL TO 1... IT ESTIMATES THE VARIANCE FACTOR BY USING ERRORS TAKEN FROM THE LS RESIDUALS AND EST. THE FACTOR FROM THIS
#IF THE WEIGHT_ITT IS SET TO 0... IT ESTIMATES THE VARIANCE FACTOR BY USING ERRORS TAKEN FROM THE FORMAL ERRORS PROVIDED IN THE EPHEMERIS DATA FILES
#IN BOTH CASES THE WEIGHT IS SET TO 1 FOR THE MEASUREMENTS BECAUSE IN THIS LS PROCESS THE WEIGHT IS BEING TRANSFERRED TO THE OBSERVATIONS 
#TO MAKE A WEIGHT MATRIX OF 1. AS A RESULT - I BELIEVED IT TO HOLD TO KEEP THE WEIGHT MATRIX SET TO 1 IN THE VARIANCE FACTOR ESTIMATION
def EstVarUnitWeight(ID,Ind,Exp,P,OptionToSet,Pos_ERR,Weight_Itt,Hours): #Estimating Variance of Unit Weight
    if Weight_Itt == 1:
        P = [0,0]
        [Fit,Type,Omega,Equations,Residuals,Variances,ChiSqr,Fit_Uncertainty] = LSFittoData(Ind,Exp,P,ID,OptionToSet,Hours)
        Eqs = Equations
        VarUnitW_hat = []   
        Weighted1 = np.repeat(1,len(Ind))
        for x in range(0,3):
            VarUnitW_hat.append(np.sqrt(np.sum(Eqs[x]**2*Weighted1)/len(Ind)))
        return VarUnitW_hat
    else:
        Weighted1 = np.repeat(1,len(Ind))
        VarUnitW_hat = []
        for x in range(0,3):
            VarUnitW_hat.append(np.sqrt(np.sum(Weighted1*Pos_ERR[:,x]**2)/len(Ind)))
        return VarUnitW_hat


#----SUMMARY OF PlotError(Hours,DataBase,ID)----
#THIS FUNCTION PLOTS HISTOGRAM OF THE ERRORS GIVEN IN THE EPHEMERIS FILES
#IN A SUBPLOT
def PlotError(Hours,DataBase,ID):
    Index,DataBase = int(Hours*3600-1), list(map(DF2.CreateDataFrame,[DFF.GNV1B_all()],ID))
    Pos_ERR = ((DataBase[0].iloc[:Index,4:7].to_numpy()))
    fig, axs = plt.subplots(3)
    axs[0].hist(Pos_ERR[:,0],label='POS X Sigma')
    axs[1].hist(Pos_ERR[:,1],label='POS Y Sigma')
    axs[2].hist(Pos_ERR[:,2],label='POS Z Sigma')
    axs[2].set(xlabel='Model Estimate of Sigma')
    axs[1].set(ylabel = 'Occurence of Error')
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    fig.show()   

#----SUMMARY OF SigmaExp(Residual,Type)----
#THIS FUNCTION RETURNS THE VARIANCE OF THE LEAST SQUARES REGRESSION
def SigmaExp(Residual,Type):
    if Type == 'Linear': 
        Sigma = np.sqrt((1/(len(Residual)-2))*(np.sum(Residual**2)))
    else:
        Sigma = np.sqrt((1/(len(Residual)-3))*(np.sum(Residual**2)))
    return Sigma

#----SUMMARY OF ChiSqrFun(Residual_List,ID,Hours,Type)
#THIS FUNCTION CALCULATES THE CHISQR TO SEE IF IT IS ~1 AND INDEED A GOOD FIT
def ChiSqrFun(Residual_List,Pos_ERR,Type):
    ChiSqr = []
    for x in range(0,len(Residual_List)):
        if Type == 'Linear':
            ChiSqr.append(np.sum((1/(Pos_ERR[:,x]**2))*Residual_List[x]**2)/(len(Residual_List[x]-2))) # Should ~ N-n
        else:
            ChiSqr.append(np.sum((1/(Pos_ERR[:,x]**2))*Residual_List[x]**2)/(len(Residual_List[x]-3))) # Should ~ N-n
    return ChiSqr


#----SUMMARY OF PropErrorToFit(Fit,Pos_ERR,Residual_List)
#PROPAGATES ERROR FROM OBSERVABLES INTO COEFFICIENT ESTIMATES USING GENERAL LAW OF ERROR PROPAGATION
def PropErrorToFit(Fit,Pos_ERR,Residual_List):
    Pos_Err_Var = []
    for x in range(0,3):
        Pos_Err_Var.append(np.var(Pos_ERR[:,x]))
    Pos_Err_Var = np.array(Pos_Err_Var)
    Error_WithError = Fit**2 @ Pos_Err_Var**2

    Residual_array = np.array(Residual_List)
    Residual_Var = []
    for x in range(0,3):
        Residual_Var.append(np.var(Residual_array[:,x]))
    Residual_Var = np.array(Residual_Var)
    Error_WithRes = Fit**2 @ Residual_Var**2
    Fit_Uncertainty = [np.sqrt(Error_WithError), np.sqrt(Error_WithRes)]
    return Fit_Uncertainty

#----SUMMARY OF ExtrapolateAndPlot(DataProd,ID,Fit,Type,Omega,Hours_Fit,Hours_Extrap)
#EXTRAPOLATES POSITION BASED ON LS FIT
def ExtrapolateAndPlot(DataProd,Ind,ID,Fit,Type,Omega,Hours_Fit,Hours_Extrap):
    IndEnd = len(Ind)
    [Ind,Exp] = SelectDataProd(DataProd,ID,Hours_Extrap)
    Ind = Ind[IndEnd:] #Extracting Data Only After Extrapolation
    Exp = Exp[IndEnd:,:]
    if Hours_Fit >= Hours_Extrap:
        Hours_Extrap = input("Hour to extrapolate to is <= hours of data-set. Please enter new hour to extrapolate to ")
    else:
        if Type == 'Ellipse':

            Eq0 = Fit[0,0] + Fit[1,0]*np.cos(2*np.pi*Omega*Ind) + Fit[2,0]*np.sin(2*np.pi*Omega*Ind)
            Eq1 = Fit[0,1] + Fit[1,1]*np.cos(2*np.pi*Omega*Ind) + Fit[2,1]*np.sin(2*np.pi*Omega*Ind)
            Eq2 = Fit[0,2] + Fit[1,2]*np.cos(2*np.pi*Omega*Ind) + Fit[2,2]*np.sin(2*np.pi*Omega*Ind)
            print(Exp[:,0])

            Eq0_Residual = Exp[:,0]-Eq0
            Eq1_Residaul = Exp[:,1]-Eq1
            Eq2_Residual = Exp[:,2]-Eq2

        elif Type == 'Linear':
            Eq0 = Fit[0,0]+Fit[1,0]*Ind
            Eq1 = Fit[0,1]+Fit[1,1]*Ind
            Eq2 = Fit[0,2]+Fit[1,2]*Ind

            Eq0_Residual = Exp[:,0]-Eq0
            Eq1_Residaul = Exp[:,1]-Eq1
            Eq2_Residual = Exp[:,2]-Eq2

        elif Type == 'Polynomial':
            Eq0 = Fit[0,0]+Fit[1,0]*Ind+Fit[2,0]*Ind**2
            Eq1 = Fit[0,1]+Fit[1,1]*Ind+Fit[2,1]*Ind**2
            Eq2 = Fit[0,2]+Fit[1,2]*Ind+Fit[2,2]*Ind**2

            Eq0_Residual = Exp[:,0]-Eq0
            Eq1_Residaul = Exp[:,1]-Eq1
            Eq2_Residual = Exp[:,2]-Eq2

        
    fig, axs = plt.subplots(3)
    axs[0].plot(Ind,Exp[:,0],'r--',linewidth = 3,label='POD DATA - TRUE')
    axs[0].plot(Ind,Eq0,'b',label='LS FIT - EXTRAP.')
    axs[0].set_title('GNSS Data in X with LS Fit EXTRAPOLATED')
    axs[1].plot(Ind,Exp[:,1],'r--',linewidth = 3,label='POD DATA - TRUE')
    axs[1].plot(Ind,Eq1,'b',label='LS FIT - EXTRAP.')
    axs[1].set_title('GNSS POD Data in Y with LS Fit EXTRAPOLATED')
    axs[2].plot(Ind,Exp[:,2],'r--',linewidth = 3,label='POD DATA - TRUE')
    axs[2].plot(Ind,Eq2,'b',label='LS FIT - EXTRAP.')
    axs[2].set_title('GNSS POD Data in Z with LS Fit EXTRAPOLATED')
    axs[1].set(ylabel = 'Position in ECI [m]')

    axs[2].set(xlabel = 'Time [s]')
    for ax in axs.flat:
        ax.label_outer()
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    fig.show()
    fig1 = plt.figure()
    ax1 = plt.axes(projection='3d')
    ax1.plot(Exp[:,0],Exp[:,1],Exp[:,2],'r--',linewidth = 3,label='POD DATA')
    ax1.plot(Eq0,Eq1,Eq2,'b',label='LS FIT')
    ax1.set_title('GNSS POD Data in 3D EXTRAPOLATED')
    ax1.set(xlabel='Data in X [m]',ylabel='Data in Y [m]',zlabel='Data in Z [m]')
    ax1.legend()
    fig1.show()

    fig2, axs2 = plt.subplots(3)
    axs2[0].plot(Ind,Eq0_Residual)
    axs2[0].set_title('GNSS POD in X LS Fit Residual EXTRAPOLATED')
    axs2[1].plot(Ind,Eq1_Residaul)
    axs2[1].set_title('GNSS POD Data in Y LS Fit Residual EXTRAPOLATED')
    axs2[2].plot(Ind,Eq2_Residual)
    axs2[2].set_title('GNSS POD Data in Z LS Fit Residual EXTRAPOLATED')
    for ax in axs2.flat:
          ax.set(xlabel='Time [s]', ylabel = 'Position in ECI [m]')
    for ax in axs2.flat:
        ax.label_outer()
    fig2.show()

    Residuals = [Eq0_Residual,Eq1_Residaul,Eq2_Residual]

    Variances = [SigmaExp(Eq0_Residual,Type),SigmaExp(Eq1_Residaul,Type),SigmaExp(Eq2_Residual,Type)] #calculting variance 
    Index,DataBase = int(Hours_Extrap*3600-1), list(map(DF2.CreateDataFrame,[DFF.GNV1B_all()],ID))
    Pos_ERR = ((DataBase[0].iloc[:Index,4:7].to_numpy()))
    Pos_ERR = Pos_ERR[IndEnd:,:]
    ChiSqr = ChiSqrFun(Residuals,Pos_ERR,Type) #calculating Chisqr
    Fit_Uncertainty =  PropErrorToFit(Fit,Pos_ERR,Residuals) #calculating uncertainty in coefficients of fit
    
    return Residuals,Variances,ChiSqr,Fit_Uncertainty


#------------------------------ LEGACY CODE WHEN USING CALCULUS BEFORE USING NORMAL EQUATION ------------------------------#
    """"
        def SystemOfEquationDet(Ind,Exp):
        #SUMMARY: FUNCTION CALLS SELECT CURVE TYPE
        #PROVIDES LS ESTIMATE NORMAL EQUATIONS FOR MODEL CONSTANTS 
                Vars = SelectCurveType()
                [sum_y, sum_x, sum_x_sqr, sum_xy, N] = sums(Ind,Exp)
                Eq1 = sum_y-Vars[0]*N-Vars[1]*sum_x
                Eq2 = sum_xy-sum_x*Vars[0]-Vars[1]*sum_x_sqr
                Solved = sp.linsolve([Eq1,Eq2],(Vars[0],Vars[1]))
                return Solved
        def SystemOfEquationDet_Symbolic(Ind,Exp):
        #SUMMARY: FUNCTION CALLS SELECT CURVE TYPE
        #PROVIDES LS ESTIMATE NORMAL EQUATIONS FOR MODEL CONSTANTS 
                Vars = SelectCurveType()
                [sum_y, sum_x, sum_x_sqr, sum_xy, N] = sums(Ind,Exp)
                xi,yi,i,n = sp.symbols('xi yi i n')
                Exp = sp.expand(sp.summation(sp.Pow(yi-Vars[0]-Vars[1]*(xi),2),(i,1,n)))

                Eq1 = sp.diff(Exp,Vars[0])
                
                Eq1 = Eq1.subs(n,N)
                Eq1 = Eq1.subs(xi,sum_x/N)
                Eq1 = Eq1.subs(yi,sum_y/N)
                Eq2 = sp.diff(Exp,Vars[1])
            
                Eq2 = Eq2.subs(n,N)
                Eq2 = Eq2.subs(xi**2,sum_x_sqr)
                Eq2 = Eq2.subs(xi*yi,sum_xy)
                Eq2 = Eq2.subs(xi,sum_x)

                Solved = sp.linsolve([Eq1,Eq2],(Vars[0],Vars[1]))
                return Solved

        def SystemOfEquationDet_Symbolic(Ind,Exp):
        #SUMMARY: FUNCTION CALLS SELECT CURVE TYPE
        #PROVIDES LS ESTIMATE NORMAL EQUATIONS FOR MODEL CONSTANTS 
                SelectCurveTypeReturn  = SelectCurveType()   
                Order = SelectCurveTypeReturn[-1]
                Vars = list(SelectCurveTypeReturn[:-1])
                xi,yi,i,n = sp.symbols('xi yi i n')

                if Vars[0] == 'Wave':
                    Wave = 'Wave'
                    Vars = Vars[1:]
                    NewExp = list([yi,-Vars[0],-Vars[1]*sp.cos(Vars[-1]*xi),-Vars[2]*sp.sin(Vars[-1]*xi)])
                    Vars = Vars[:-1]
                else:
                    Wave = []
                    MinVarEq = []
                    for x in range(0,len(Vars)): 
                        MinVarEq.append(Vars[x]*xi**x)
                    MinVarEq.insert(0,yi)
                    NewExp = []
                    for x in range (0,len(MinVarEq)):
                        OldExp = MinVarEq[x]
                        if x == 0:
                            NewExp.append(OldExp)
                        else:
                            NewExp.append(-OldExp)

                ExpSym = sp.expand(sp.summation(sp.Pow(sum(NewExp),2),(i,1,n))) #Expression to be solved
                if Wave == 'Wave':
                    Omega = np.linspace(10000,500,500) #Setting BW for LSSA over entire interval
                    Omega = 1/Omega
                    spectrum = LS.LSSA(Ind,Exp,P=1,Omega=Omega,ind=[],level=0.01,trend='linear')[0]
                    maxSpectrum = np.argmax(spectrum)
                    maxOmegaIndex = Omega[maxSpectrum]
                    w = sp.symbols('w')
                    ExpSym = ExpSym.subs(w,maxOmegaIndex)
                Eq = []
                EqSub = []
                for x in range(0,len(Vars)):
                    Eq.append(sp.diff(ExpSym,Vars[x]))
                if Wave == 'Wave':
                    Eq = SymPySubWrapperWave(Eq,xi,yi,n,Ind,Exp)
                    Solved = sp.linsolve(list(Eq),(Vars))
                    return Solved, maxOmegaIndex

                else: 
                    for x in range(0,len(Vars)):
                        EqSub.append(SymPySubWrapper(Eq[x],xi,yi,n,Order,Ind,Exp))
                    Solved = sp.linsolve(list(EqSub),(Vars))               
                    return Solved

        def SymPySubWrapper(Eq,xi,yi,n,Order,Ind,Exp):
            [sum_y, sum_x, sum_x_sqr, sum_xy, No] = sums(Ind,Exp)
            Eq = Eq.subs(n,No)

            for x in range(Order,0,-1):
                sum_xy_N = sums_xy(Ind,Exp,x,1)
                Eq = Eq.subs(xi**Order*yi,sum_xy_N/No)

            for x in range(Order,-1,-1):
                sum_x_sqr_N = SumIndOrderN(Ind,x+2)       #IF ORDER IS LINEAR (ORDER 1), THEN EQUATION WILL REQUIRE SUM (Xi)^2 ONLY
                Eq = Eq.subs(xi**(x+2),sum_x_sqr_N/No)          #IF ORDER IS POLYNOMIAL (ORDER 2), THEN EQUATION WILL REQUIRE SUM (Xi)^2 SUB THEN (Xi)^3

            Eq = Eq.subs(xi*yi,sum_xy/No)
            Eq = Eq.subs(xi,sum_x/No)
            Eq = Eq.subs(yi,sum_y/No)  
            return Eq

        def SymPySubWrapperWave(Eq,xi,yi,n,Ind,Exp):
            [sum_y, sum_x, sum_x_sqr, sum_xy, No] = sums(Ind,Exp)
            A,B,C,OmegaSym = sp.symbols('A B C w')
            Eq1, Eq2, Eq3 = Eq[0], Eq[1], Eq[2]

            Eq1 = Eq1.subs(A*n,A*No)
            Eq1 = Eq1.subs(B*n,B)
            Eq1 = Eq1.subs(C*n,C)
            Eq1 = Eq1.subs(n*yi,sum_y)
            Eq1 = Eq1.subs(xi,sum_x)
            Eq1 = Eq1.subs(n*yi,sum_y)

            Eq2 = Eq2.subs(A*n,A)
            Eq2 = Eq2.subs(B*n,B)
            Eq2 = Eq2.subs(C*n,C)
            Eq2 = Eq2.subs(n*yi,sum_y)
            Eq2 = Eq2.subs(xi,sum_x)

            Eq3 = Eq3.subs(A*n,A)
            Eq3 = Eq3.subs(B*n,B)
            Eq3 = Eq3.subs(C*n,C)
            Eq3 = Eq3.subs(n*yi,sum_y)
            Eq3 = Eq3.subs(xi,sum_x)

            Eq = []
            Eq = [Eq1,Eq2,Eq3]
            return Eq

        def SumIndOrderN(Ind,No): 
            sum_x_sqr_N = np.sum(Ind**No)
            return sum_x_sqr_N

    """