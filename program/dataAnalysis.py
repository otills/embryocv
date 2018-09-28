import cv2
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import pyqtgraph as pg
import pandas as pd
import glob
import os
import re
from skimage.segmentation import clear_border
from skimage.morphology import disk
from PyQt5.Qt import *
import sys
#from imageAnalysis import imageAnalysis
import eggUI
#import viewOutput
import time
import pathos
import json
#import tables
import scipy
import scipy.fftpack
import pylab
from scipy import pi
import scipy.signal as signal
import peakutils
import xarray as xr
import math
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from matplotlib.backends.backend_pdf import PdfPages
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

class dataAnalysis(object):    

#==============================================================================
#       Produce summary reports for each embryo and save to savePath. 
#       Note: savePath must exist i.e. create the folder beforehand..
#==============================================================================
    def generateSummaryReports(self,savePath):
        for e in range(len(self.embryoLabels)):
            self.embryo = self.embryoLabels[e]
            self.loadXRResults()
            if self.results is not 'NoData':
                data = self.results['TimeSpecificSummaryData'].to_pandas()
                f, ax = plt.subplots(4, sharex=True)
                ax[0].set_yscale('log')
                # Size
                for i in range(13,24):
                    ax[0].fill_between(np.arange(0,len(data.ix[:,i].values)),data.ix[:,i].values, alpha = 0.3)
                ax[0].set_title('Enery at different frequencies')
                ax[1].plot(data.ix[:,0].values)
                ax[1].fill_between(np.arange(0,len(data.ix[:,1].values)),data.ix[:,1].values,data.ix[:,2].values, alpha = 0.3)
                ax[1].set_title('Min:Max and mean area')     
                ax[2].plot(data.ix[:,3].values)
                ax[2].fill_between(np.arange(0,len(data.ix[:,1].values)),data.ix[:,4].values,data.ix[:,5].values, alpha = 0.3)
                ax[2].set_title('Min:Max and mean BB')
                ax[3].plot(data.ix[:,7].values)
                ax[3].set_yscale('log')
                ax[3].fill_between(np.arange(0,len(data.ix[:,1].values)),data.ix[:,8].values,data.ix[:,9].values, alpha = 0.3)
                ax[3].set_title('Min:Max and mean distance')
                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
                pdf_pages = PdfPages(savePath + self.embryoLabels[e] + '_phenomeSummary.pdf')
                pdf_pages.savefig(f)
                pdf_pages.close()
                plt.plot()
                plt.close()
                
#==============================================================================
#    Measure heart rate for all embryos 
#==============================================================================
    def measureHeartRateForAllEmbryos(self,savePath,filtMin=False,filtMax=False,histCount=False,segRegfiltXVal=False,segRegfiltYVal=False,minXVal=False,peakThresh=False):
        if self.species == 'rbalthica':
            # If user has not specified these 'tweakable paramters' in the HR ID
            # and model fitting use these values
            if filtMin == False:
                filtMin = 1
            if filtMax == False:
                filtMax = 2
            if histCount == False:
                histCount = 6
            if segRegfiltXVal == False:
                segRegfiltXVal = 5
            if segRegfiltYVal == False:
                segRegfiltYVal = 6
            if minXVal == False:
                minXVal = 60
            if peakThresh == False:
                peakThresh = 0.4   
            for e in range(len(self.embryoLabels)):
                self.embryo = self.embryoLabels[e]
                self.loadXRResults()
                self.measureHeartRate_radix(savePath,filtMin,filtMax,histCount,segRegfiltXVal,segRegfiltYVal,minXVal,peakThresh)
            np.save(savePath + '/hrdata.npy', self.heartRate_data)
        if self.species == 'ogammarellus':
            # If user has not specified these 'tweakable paramters' in the HR ID
            # and model fitting use these values
            if filtMin == False:
                filtMin = 0.5
            if filtMax == False:
                filtMax = 4
            if histCount == False:
                histCount = 6
            if minXVal == False:
                minXVal = False
            if peakThresh == False:
                peakThresh = 0.2   
            print str(savePath)  + '_' + str(filtMin) + '_' +str(filtMax) + '_' + str(histCount) + '_' + str(peakThresh)
            for e in range(len(self.embryoLabels)):
                self.embryo = self.embryoLabels[e]
                self.loadXRResults()
                self.measureHeartRate_orchestia(savePath,filtMin,filtMax,histCount,peakThresh)
            np.save(savePath + '/hrdata.npy', self.heartRate_data)
#==============================================================================
#     Modify heart rate modelling for specific embryo
#==============================================================================
    def measureHeartRateForSpecificEmbryos(self,embryo,savePath,filtMin=False,filtMax=False,histCount=False,segRegfiltXVal=False,segRegfiltYVal=False,minXVal=False,peakThresh=False):
        if self.species == 'rbalthica':
            # If user has not specified these 'tweakable paramters' in the HR ID
            # and model fitting use these values
            if filtMin == False:
                filtMin = 1
            if filtMax == False:
                filtMax = 2
            if histCount == False:
                histCount = 6
            if segRegfiltXVal == False:
                segRegfiltXVal = 5
            if segRegfiltYVal == False:
                segRegfiltYVal = 6
            if minXVal == False:
                minXVal = 60
            if peakThresh == False:
                peakThresh = 0.4                 
            if isinstance(embryo, str):
                self.heartRate_data = np.load(savePath + '/hrdata.npy')
                self.heartRate_data = self.heartRate_data[()]
                self.embryo = embryo
                self.loadXRResults()
                self.measureHeartRate_radix(savePath,filtMin,filtMax,histCount,segRegfiltXVal,segRegfiltYVal,minXVal,peakThresh)
                np.save(savePath + '/hrdata.npy', self.heartRate_data)
            if isinstance(embryo, list):   
                self.heartRate_data = np.load(savePath + '/hrdata.npy')
                self.heartRate_data = self.heartRate_data[()]
                for e in range(len(embryo)):            
                    self.embryo = embryo[e]
                    self.loadXRResults()
                    self.measureHeartRate_radix(savePath,filtMin,filtMax,histCount,segRegfiltXVal,segRegfiltYVal,minXVal,peakThresh)
                np.save(savePath + '/hrdata.npy', self.heartRate_data)
        elif self.species == 'ogammarellus':
            if filtMin == None:
                filtMin = 0.5
            if filtMax == None:
                filtMax = 4
            if histCount == None:
                histCount = 6
            if minXVal == None:
                minXVal = False
            if peakThresh == None:
                peakThresh = 0.2  
            if isinstance(embryo, str):
                self.heartRate_data = np.load(savePath + '/hrdata.npy')
                self.heartRate_data = self.heartRate_data[()]
                self.embryo = embryo
                self.loadXRResults()
                self.measureHearRate_orchestia(self,savePath,filtMin,filtMax,histCount,peakThresh)
                np.save(savePath + '/hrdata.npy', self.heartRate_data)
            if isinstance(embryo, list):   
                self.heartRate_data = np.load(savePath + '/hrdata.npy')
                self.heartRate_data = self.heartRate_data[()]
                for e in range(len(embryo)):            
                    self.embryo = embryo[e]
                    self.loadXRResults()
                    self.measureHearRate_orchestia(self,savePath,filtMin,filtMax,histCount,peakThresh)
                np.save(savePath + '/hrdata.npy', self.heartRate_data)
                
#==============================================================================
#   Load heart rate data
#==============================================================================
    def loadHeartRateData(self,savePath):
        self.heartRate_data = np.load(savePath + '/hrdata.npy')
        self.heartRate_data = self.heartRate_data[()]    
#==============================================================================
#   Umbrella function for measuring heart rate - Radix balthica
#==============================================================================
    def measureHeartRate_orchestia(self,savePath,filtMin,filtMax,histCount,peakThresh):        
        if not hasattr(self,'heartRate_data'):
            self.heartRate_data = dict()
        if self.results is 'NoData':
            print 'No data for ' + str(self.embryo)
            
        freqres = self.results['FreqOutput_8x8'].values   
        gs, fs = np.meshgrid(range(8),range(8))
        gs = gs.flatten()
        fs = fs.flatten()   
        HRs=[]
        # Loop over time points
        for t in range(self.results['FreqOutput_8x8'].values.shape[0]):
            HRFreqs=[]  
            # Loop over blockwise signals
            for b in range(len(gs)):
                powerSpect = freqres[t,fs[b],gs[b],1,:]
                if (powerSpect.max() == 0.0) or np.isnan(powerSpect[0]):
                    HRFreqs.append(np.NaN)
                else:
                    sig = pd.rolling_mean(self.results['BlockWise_8x8'].loc[t,:,fs[b],gs[b]].values, window=3)
                    sig = pd.DataFrame(sig).interpolate(limit_direction = 'both').values.ravel()
                    # Interpolate to fill missing data
                    frameRate = sig.shape[0]/((self.results['SizePos'].to_pandas().ix[t,sig.shape[0]-1,'elapsedTime'] - self.results['SizePos'].to_pandas().ix[t,0,'elapsedTime'])/1000)
                    sampFreqs, powerSpect = signal.welch(sig,frameRate, scaling='spectrum')
                    baselineRemoved = np.log(powerSpect) - peakutils.baseline(np.log(powerSpect)- min(np.log(powerSpect)))
                    # Now ID peaks.anan
                    indexes = peakutils.indexes(baselineRemoved[0:200], peakThresh, min_dist=0)
                    outPeakFreqs = sampFreqs[indexes]
                    inds = np.where(((outPeakFreqs < filtMax) & (outPeakFreqs > filtMin)))
                    filtFreqs = outPeakFreqs[inds]
                    filtPower = baselineRemoved[inds]
                    # If more than one freq peak identified within the range 0.5-3 Hz
                    if inds[0].shape[0] > 1:
                        #print str(b)
                        # Identify freq with most power
                        maxInd = np.argmax(filtPower)
                        HRFreq = filtFreqs[maxInd]
                        HRFreqs.append(HRFreq)
                    elif inds[0].shape[0] == 1:
                        HRFreq = outPeakFreqs[inds[0][0]]
                        HRFreqs.append(HRFreq)
                    else:
                        HRFreqs.append(np.NaN)
            else:
                HRFreqs.append(np.NaN)
            HRFreqs = np.array(HRFreqs)
            # Remove NaNs
            count, freqs = np.histogram(HRFreqs[np.invert(np.isnan(HRFreqs))],bins=50)
            # Identify whether a sufficiently dominant frequency to reliably ID HR.
            inds = np.where(count > histCount)            
            #plt.hist(HRFreqs[~np.isnan(HRFreqs)],20)
            if len(inds[0]) is not 0:
            # Get max ind
                maxInds = inds[-1]
                HR = np.mean(freqs[maxInds])    
                print HR
                HRs.append(HR)
            else:
                HRs.append(np.NaN)
        # Attempt to fit a linear model to data from 24 h post HR detection..
        try:
            # Use for fitting model
            reducedHRs = HRs[23:]
            # Fit a linear model to Orchestia heart function ontogeny
            # Format data
            pdData = np.array([np.arange(len(reducedHRs))[~np.isnan(reducedHRs)], np.array(reducedHRs)[~np.isnan(reducedHRs)]]).T
            df = pd.DataFrame(pdData,columns=['Ind','HR'])
            mod = smf.ols(formula='HR ~ Ind', data=df)
            res = mod.fit()
            print(res.summary())
            
            # If any data points have residuals > 0.5, remove these.
            if np.any(np.sqrt(res.resid*res.resid) > 0.5):
                df = df.drop(df.ix[np.sqrt(res.resid*res.resid) > 0.5,:].index)
                mod = smf.ols(formula='HR ~ Ind', data=df)
                res = mod.fit()
            # Get predicted HR values
            # predictedHRs = res.predict(np.arange(24))
            predictedHRs = res.predict(pd.DataFrame(data = np.arange(24), columns = ['Ind']))
            # Get intercept and gradient
            p = mod.fit().params
            print 'Intercept: ' + str(p[0]) + '; Slope: ' + str(p[1])
                        # Show fit
            fig, ax = plt.subplots(1,1)
            plt.title(str(self.embryo))
            plt.xlabel('Time point')
            plt.ylabel('Frequency (Hz)')
            ax.plot(df.ix[:,'Ind'],df.ix[:,'HR'],'o')
            ax.plot(df.ix[:,'Ind'], p[0] + p[1] * df.ix[:,'Ind'],color='black')
            plt.ylim(0,df.ix[:,'HR'].max()*1.2)
            fig.subplots_adjust(bottom=0.3)
            lab = str('Intercept: ') + str(p[0]) + '; Slope: ' + str(p[1])
            fig.text(.1,.1,lab)
            #fig.canvas.manager.window.raise_()
            pdf_pages = PdfPages(savePath + self.embryo + '_heartModelSummary.pdf')
            pdf_pages.savefig(fig)
            pdf_pages.close()
            plt.plot()

            
            zip(np.arange(0, len(self.results['dateTime'].values[23:])),self.results['dateTime'].values,self.results['Metadata'].to_pandas().ix[23:,0,'currentFolder'],HRs[23:],predictedHRs)
            
            self.heartRate_data[str(self.embryo)] = {'Embryo': self.embryo,'PredHR':predictedHRs, 'Intercept': str(p[0]), 'Slope':str(p[1]),
                'modPars_filtMin':filtMin,'modPars_filtMax':filtMax,'modPars_histCount':histCount,'modPars_peakThresh':peakThresh}
        except:
            print '*** Heart rate not succesfully modelled for ' + str(self.embryo) + ' ***'


#==============================================================================
#   Umbrella function for measuring heart rate - Radix balthica
#==============================================================================
    def measureHeartRate_radix(self,savePath,filtMin=1,filtMax=2,histCount=6,segRegfiltXVal=5,segRegfiltYVal=6,minXVal=False,peakThresh=0.4):
        """ Attempt to extract and model heart rate from 8x8 blockwise frequency data.
        Note that some user guidance is recommended and if the heart rate cannot be succesfully
        modelled try changing some of these optional arguments:
            # Provide a folder for saving output plots and statistics.
            # The minimum acceptable HR (Hz)
            filtMin=1
            # The maximum acceptable HR (Hz)
            filtMax=2
            # The number of signals in which a frequency must occur as being dominant to be considered.
            histCount=6
            # The extent to which outlying values are filtered in the X dimension prior to fitting a segmented regression.
            segRegfiltXVal=4
            # The extent to which outlying values are filtered in the Y dimension prior to fitting a segmented regression.
            segRegfiltYVal=6
            # The threshold for identifying peaks in the frequency output
            peakThresh=0.4
            """
        #if 'heartRate_data' not in locals():
        if not hasattr(self,'heartRate_data'):
            self.heartRate_data = dict()
        if self.results is 'NoData':
            print 'No data for ' + str(self.embryo)
        else:
            self.freqRes = self.results['FreqOutput_8x8'].values
            hrvals = []
            gs, fs = np.meshgrid(range(8),range(8))
            gs = gs.flatten()
            fs = fs.flatten()
            for i in range(self.freqRes.shape[0]):
                hrvals.append(self.modelHR_radix(gs,fs,i,filtMin,filtMax,histCount,minXVal,peakThresh))
                self.test = hrvals
            #try:
            self.predHR, self.predBP, self.BP, self.S1,self.S2,self.startTime,self.endTime,self.firstHR,self.lastHR,self.HRData,self.outputSummary = self.segReg_radix(hrvals,segRegfiltXVal,segRegfiltYVal)
            # Show fit
            fig, ax = plt.subplots(1,1)
            ax.plot(self.predHR)
            ax.plot(self.HRData,'o')
            #ax.plot(sortedTimes, sortedPredicted)
            #ax.plot(times, np.array(HRData)[~np.isnan(HRData)], 'o')
            plt.title(str(self.embryo))
            plt.show()
            plt.xlabel('Time point')
            plt.ylabel('Frequency (Hz)')
            plt.ylim(0,(np.nanmax(np.array((self.predHR,self.HRData)))*1.2))
            fig.subplots_adjust(bottom=0.3)
            fig.text(.1,.1,str(self.outputSummary[0] + '\n' + self.outputSummary[1] + '\n' + self.outputSummary[2]))
            fig.canvas.manager.window.raise_()
            pdf_pages = PdfPages(savePath + self.embryo + '_heartModelSummary.pdf')
            pdf_pages.savefig(fig)
            pdf_pages.close()
            plt.plot()
            #time.sleep(0.2)
            # Add data to a dict - new key for each embryo
            self.heartRate_data[str(self.embryo)] = {'Embryo': self.embryo,'PredHR':self.predHR, 
            'PredBP':self.predBP, 'BP': self.BP, 'S1': self.S1,'S2':self.S2,'startTime':self.startTime,
            'endTime':self.endTime,'firstHR':self.firstHR,'lastHR':self.lastHR,'HRData':self.HRData, 
            'modPars_filtMin':filtMin,'modPars_filtMax':filtMax,'modPars_histCount':histCount,'modPars_segRegfiltXVal':segRegfiltXVal,
            'modPars_segRegfiltYVal':segRegfiltYVal,'modPars_minXVal':minXVal,'modPars_peakThresh':peakThresh}
            #except:
              #  print '*** Heart rate not succesfully modelled for ' + str(self.embryo) + ' ***'

#==============================================================================
#   Model HR for Radix balthica
#==============================================================================
    # Extract peaks from the freq output that may be HR.
    def modelHR_radix(self,gs,fs, t,filtMin, filtMax, histCount, minXVal,peakThresh):
        """ Called from measureHeartRate_radix(). 
        """
        HRFreqs=[]
        # If no minVal (time point) provided, or the t value is above the minimum value
        # attempt to ID the heart rate.
        if type(minXVal) is not bool and minXVal<= t or not minXVal:
            if not math.isnan(self.freqRes[t,0,0,1,0]):
                for b in range(len(gs)):
                    powerSpect = self.freqRes[t,fs[b],gs[b],1,:]
                    if (powerSpect.max() == 0.0):
                        HRFreqs.append(np.NaN)
                    else:
                        sampFreqs = self.freqRes[t,fs[b],gs[b],0,:]
                        try:
                            baselineRemoved = np.log(powerSpect) - peakutils.baseline(np.log(powerSpect)- min(np.log(powerSpect)))
                            # Now ID peaks.anan
                            indexes = peakutils.indexes(baselineRemoved, peakThresh, min_dist=1)
                            outPeakFreqs = sampFreqs[indexes]
                            # And get magnitudes (may be useful for downstream filtering).
                            outMagnitudes = baselineRemoved[indexes]
                            # Identify index(or indices)
                            # Gam
                            #inds = np.where(((outPeakFreqs < 5) & (outPeakFreqs > 0.5)))
                            # Rad
                            #inds = np.where(((outPeakFreqs < 2) & (outPeakFreqs > 1)))
                            inds = np.where(((outPeakFreqs < filtMax) & (outPeakFreqs > filtMin)))
                            filtFreqs = outPeakFreqs[inds]
                            filtPower = baselineRemoved[inds]
                            # If more than one freq peak identified within the range 0.5-3 Hz
                            if inds[0].shape[0] > 1:
                                #print str(b)
                                # Identify freq with most power
                                maxInd = np.argmax(filtPower)
                                HRFreq = filtFreqs[maxInd]
                                HRFreqs.append(HRFreq)
                            elif inds[0].shape[0] == 1:
                                HRFreq = outPeakFreqs[inds[0][0]]
                                HRFreqs.append(HRFreq)
                            else:
                                HRFreqs.append(np.NaN)
                        except:
                            HRFreqs.append(np.NaN)
                else:
                    HRFreqs.append(np.NaN)
                HRFreqs = np.array(HRFreqs)
                # Remove NaNs
                count, freqs = np.histogram(HRFreqs[np.invert(np.isnan(HRFreqs))],bins=50)
                # Identify whether a sufficiently dominant frequency to reliably ID HR.
                inds = np.where(count > histCount)
                if len(inds[0]) is not 0:
                    # Get max ind
                    maxInds = inds[-1]
                    HR = np.mean(freqs[maxInds])
                else:
                    HR = np.NaN
            else:
                HR = np.NaN
        else:
            # If current time iteration is below the min value (time) assign NaN.
            HR = np.NaN
        return HR
    
#==============================================================================
#   HR - fits a segmented regression using R.
#   Used to model HR for R. balthica
#==============================================================================
    def segReg_radix(self,HRData,filtXVal,filtYVal):
        #  Fit segmented regressionpN in R to the HR data. Note various filtering done here..
        base = importr('base')
        r = robjects.r
        seg = importr('segmented')
        # Filter data to remove lone x axis
        # ID where values exist
        filt = ~np.isnan(HRData)
        # Use histogram to ID outlier
        count, bins = np.histogram(np.array(range(0,len(HRData)))[filt])
        # Filter to remove outlying x values
        filtInd = bins[:-1][(count > 0) & (count < filtXVal)]
        filtInds = np.nonzero([(count > 0) & (count < filtXVal)])[1]
        # If any time points to filter
        HRData = np.array(HRData)
        if filtInd.shape[0] is not 0:
            #HRData[filtInd.astype(np.int)] = np.NaN
            # Add filter to remove any values within particular bin requiring filtering..
            for i in range(len(filtInds)):
                for t in range(len(HRData)):
                    if (t >= bins[filtInds][i]) & (t <= bins[filtInds[i]+1]):
                        # Debug
                        # print 'X filtered'
                        HRData[t] = np.NaN
        # Filter to remove outlying y values
        Yax, bins = np.histogram(HRData[~np.isnan(HRData)])
        filtInd = Yax < filtYVal
        filtInds = np.nonzero([Yax < filtYVal])[1]
        if filtInds.shape[0] is not 0:
            #Yax[1][1:][Yax[0] < 2]
            for i in range(len(filtInds)):
                for t in range(len(HRData)):
                    if (HRData[t] >= bins[filtInds][i]) & (HRData[t] <= bins[filtInds[i]+1]):
                        #print 'success'
                        HRData[t] = np.NaN
        # Format/convert data
        HR = robjects.FloatVector(HRData)
        Time = robjects.FloatVector(range(0,len(HRData)))
        robjects.globalenv["HR"] = HR
        robjects.globalenv["Time"] = Time
        # Run linear model
        lmHR = r.lm("HR ~ Time")
        #print(base.summary(lmHR))
        # Model for seg.z
        formula = robjects.Formula("~Time")
        # Also use mean of data containing values to get starting psi
        startingPSI = np.int(np.median(np.arange(0,len(HRData))[~np.isnan(HRData)]))
        # Run segmented model and print results
        segModHR = seg.segmented(lmHR, seg_Z=formula, psi=startingPSI, model=True)
        #print(seg.summary_segmented(segModHR, short = True))
        # print(seg.print.segmented(segModHR))
        # Plotting segmented in Python
        resultsDict = dict(zip(segModHR.names, list(segModHR)))
        breakPoint = resultsDict['psi'][1]
        
        # Get fitted values 
        predicted = np.array(seg.predict_segmented(segModHR))
        times = np.arange(0,len(HRData))[~np.isnan(HRData)]
        #plt.plot(times, np.array(HRData)[~np.isnan(HRData)], 'o')
    
        # Get slopes
        slopes = seg.slope(segModHR)
        slopes = np.array(slopes[0])
        slopeOne = slopes[0,0]
        slopeTwo = slopes[1,0]
        #print(seg.summary_segmented(segModHR))
        self.out = seg.summary_segmented(segModHR)
        self.tmp = segModHR.items
        
        # Get fitted values for all times (predicted values above skip missing data)
        # Create an R dataframe with x-values to predict from
        d = {'Time': robjects.FloatVector(np.arange(times.min(),times.max()))}
        dataf = robjects.DataFrame(d)
        # Predict values - all times (although when plotting in Matplotlib the breakpoint will appear slightly flattened)
        predictedVals = r.predict(segModHR, newdata = dataf)
        #plt.plot(np.arange(times.min(),times.max()),np.array(predictedVals))
        
        # Get fitted value for breakpoint - this must be used in matplotlib plotting for break to appear appropriately
        d = {'Time': robjects.FloatVector([breakPoint])}
        breakpointTime = robjects.DataFrame(d)
        pyBreakpointTime = np.array(breakpointTime)[0][0]
        # Predict values
        predictedBreakPointVal = np.array(r.predict(segModHR, newdata = breakpointTime))
        #plt.plot(np.append(np.arange(times.min(),times.max()),breakPoint),np.append(predictedVals,predictedBreakPointVal), 'o')
       
        # Sort to show segmented regression including appropriate breakpoint.
        timeSort = np.append(np.arange(times.min(),times.max()),breakPoint)
        sortedTimes = timeSort[np.argsort(timeSort)]
        predicted = np.append(predictedVals,predictedBreakPointVal)
        sortedPredicted = predicted[np.argsort(timeSort)]
        
        # Pad predicted values with NaNs to shift appropriately.
        startTime = np.min(np.where((~np.isnan(HRData))))
        endTime = np.max(np.where((~np.isnan(HRData))))
        #endHR = np.argmax(~np.isnan(HRData[::-1]))
        predictedHRArray = np.zeros(len(HRData))
        predictedHRArray[:] = np.NaN
        predictedHRArray[times.min():times.max()] = predictedVals

        # Get first predicted rate
        firstHR = predictedHRArray[np.where(~np.isnan(predictedHRArray))[0][0]]
        # Get last predicted rate
        lastHR = predictedHRArray[np.where(~np.isnan(predictedHRArray))[0][-1]]
        # Get rate at breakpoint
        #breakpointHR = predictedHRArray[np.int(breakPoint)]
        
        print 'Heart rate modelling for: ', str(self.embryo)
        print 'Breakpoint HR: ', str(predictedBreakPointVal[0]), ', Slope 1: ', str(slopeOne), ', Slope 2: ', str(slopeTwo)
        print 'Heart rate detected: ', str(startTime), 'Breakpoint occured: ', str(pyBreakpointTime), 'Heart rate no longer detected: ', str(endTime)
        print 'Start heart rate: ', str(firstHR), ', End heart rate: ', str(lastHR)
        
        outputSummary = [str('Heart rate modelling for: '+ str(self.embryo)),
        str('Breakpoint HR: '+ str(predictedBreakPointVal[0])+ ' Slope 1: '+ str(slopeOne)+ ' Slope 2: '+ str(slopeTwo)),
        str('Start heart rate: '+ str(firstHR)+ ' End heart rate: '+ str(lastHR))]
         
        return list(predictedHRArray), list(predictedBreakPointVal), pyBreakpointTime, slopeOne, slopeTwo, startTime, endTime, firstHR, lastHR, list(HRData),outputSummary
 



    def identifyLethalEndPoints(self, savePath, developmentalStage):
        """ 
        Identify lethal end points. Currently optimised for 24 h exposure (3 recordings/hour) and tested for three
        developmental stages of Radix balthica (E3 = trochophore, E7 = mid-hippo and E11 = late hippo), employing different strategies
        for each (drops in frequency energy, peaks in size indicative of a 
        failure in osmotic control and a combination of both).
        """
        # If developmentalStage is late hippo identify drops in frequency energy
        if developmentalStage is 'E11':
            self.lethalEndPoint_data = dict()
            for e in range(len(self.embryoLabels)):
                self.embryo = self.embryoLabels[e]
                self.loadXRResults()
                # Extract required data from TimeSpecificSummaryData component of XArray dataset
                totalLower = np.nansum(self.results['TimeSpecificSummaryData'].to_pandas().ix[:72,13:17],axis=1)
                avHigher = np.nanmean(self.results['TimeSpecificSummaryData'].to_pandas().ix[:72,18:],axis=1)
                # Normalise - this seems to work well.
                data = pd.DataFrame(totalLower/avHigher).interpolate(limit_direction = 'both').values.ravel()
                data = pd.DataFrame(data).interpolate(limit_direction = 'both').values.ravel()
                # Identify where energy levels fall below threshold - likely to require diff
                # thresholds for diff stages etc - first filter.
                inds = np.where(np.log(data) < 3)
                # If more than one drop identified - take just the first (falling below)
                if len(inds[0]) > 1:
                    lethalIndex = np.min(inds[0])
                    lethalTime = self.results['TimeSpecificSummaryData'].to_pandas().index[lethalIndex] 
                    self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': lethalTime,'LethalIndex': lethalIndex}
                else:
                    # Identify where energy levels fall below threshold - likely to require diff
                    # thresholds for diff stages etc - second filter. 
                    inds = np.where(np.log(data) < 3.3)
                    if len(inds[0]) > 1:
                        lethalIndex = np.min(inds[0])
                        lethalTime = self.results['TimeSpecificSummaryData'].to_pandas().index[lethalIndex] 
                        self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': lethalTime,'LethalIndex': lethalIndex}
                                                 
                    else:
                        lethalIndex = np.NaN
                        self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': lethalIndex}
            np.save(savePath + '/lethalEndPoints.npy', self.lethalEndPoint_data)
        
        # If developmentalStage is mid hippo first attempt to identify drops in frequency energy
        # And failing this try peaks in size indicative of a failure in osmotic control.
        if developmentalStage is 'E7':
            self.lethalEndPoint_data = dict()
            for e in range(len(self.embryoLabels)):
                self.embryo = self.embryoLabels[e]
                self.loadXRResults()
                # Extract required data from TimeSpecificSummaryData component of XArray dataset
                totalLower = np.nansum(self.results['TimeSpecificSummaryData'].to_pandas().ix[:,13:17],axis=1)
                avHigher = np.nanmean(self.results['TimeSpecificSummaryData'].to_pandas().ix[:,18:],axis=1)
                # Normalise - this seems to work well.
                #data = pd.DataFrame(totalLower/avHigher).interpolate(limit_direction = 'both').values.ravel()
                #data = pd.DataFrame(data).interpolate(limit_direction = 'both').values.ravel()
                # Identify where energy levels fall below threshold - likely to require diff
                # thresholds for diff stages etc - first filter.
                inds = np.where(np.log(totalLower/avHigher) < 3)
                # If more than one drop identified - take just the first (falling below)
                if len(inds[0]) > 1:
                    lethalIndex = np.min(inds[0])
                    lethalTime = self.results['TimeSpecificSummaryData'].to_pandas().index[lethalIndex] 
                    self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': lethalTime,'LethalIndex': lethalIndex}
                    #print str(lethalIndex)
                else:
                    self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': np.NaN}
            # If no lethal end point has been found also check for peaks in embryo size  
            # indicative of a loss of osmotic control
            for e in range(len(self.embryoLabels)):
                self.embryo = self.embryoLabels[e]
                self.loadXRResults()
                if self.lethalEndPoint_data[str(self.embryo)]['LethalTime'] is np.NaN:
                    data = pd.rolling_mean(self.results['TimeSpecificSummaryData'].to_pandas().ix[:,0],window=6)
                    if np.count_nonzero(np.isnan(data))/np.float(len(data)) > 0.7:
                        self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': np.NaN}
                    else:
                        # Interpolate to fill missing data
                        data = pd.DataFrame(data).interpolate(limit_direction = 'both').values.ravel()
                        #print np.count_nonzero(np.isnan(data))
                        # Identify peaks/loss of embryo osmotic control
                        indexes = peakutils.indexes(data-peakutils.baseline(data), thres=0.99, min_dist=1)
                        if len(indexes) == 1:
                            # Check that the peak is sufficiently distant from the baseline..
                            # Assess proportionate distance of peak from baseline.  
                            baseline = peakutils.baseline(data)
                            propDist = np.abs(baseline-data[indexes]).min()/np.mean(baseline)     
                            #print propDist
                            if propDist < 0.2:
                                self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': np.NaN}
                            else:
                                lethalTime = self.results['TimeSpecificSummaryData'].to_pandas().index[indexes[0]] 
                                lethalIndex = indexes[0]
                                self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': lethalTime,'LethalIndex': lethalIndex}
                                #print str(lethalIndex)        
                                #print 'Lethal end point identified via peak in area'
                        if len(indexes) > 1:
                            baseline = peakutils.baseline(data)
                            index = np.argmax(data[indexes])
                            # Check that the peak is sufficiently distant from the baseline..
                            # Assess proportionate distance of peak from baseline.   
                            propDist = np.abs(baseline-data[indexes][index]).min()/np.mean(baseline)    
                            print propDist
                            if propDist < 0.2:
                                #print 'No peak found'
                                 self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': np.NaN}
                            else:
                                #mode[e] = 'ar'
                                lethalTime = self.results['TimeSpecificSummaryData'].to_pandas().index[index] 
                                lethalIndex =index
                                self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': lethalTime,'LethalIndex': lethalIndex}
                                #print str(lethalIndex)        
                                #print 'Lethal end point identified via peak in area'
                        if len(indexes) < 1:
                             self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': np.NaN}
                             #print 'Cant be sure'
                            
        # If developmentalStage is trochophore identify drops in frequency energy
        if developmentalStage is 'E3':                
            self.lethalEndPoint_data = dict()
            for e in range(len(self.embryoLabels)):
                self.embryo = self.embryoLabels[e]
                self.loadXRResults()
                data = pd.rolling_mean(self.results['TimeSpecificSummaryData'].to_pandas().ix[:,0],window=2)    
                if np.count_nonzero(np.isnan(data))/np.float(len(data)) > 0.7:
                    print 'Not enough data for ' + str(self.embryo)
                    self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': np.NaN}
                else:
                    # Interpolate to fill missing data
                    data = pd.DataFrame(data).interpolate(limit_direction = 'both').values.ravel()
                    # Identify peaks/loss of embryo osmotic control
                    indexes = peakutils.indexes(data-peakutils.baseline(data), thres=0.9, min_dist=1)
                    if len(indexes) == 1:
                        # Check that the peak is sufficiently distant from the baseline..
                        # Assess proportionate distance of peak from baseline.  
                        baseline = peakutils.baseline(data)
                        propDist = np.abs(baseline-data[indexes]).min()/np.mean(baseline)     
                        if propDist < 0.1:
                            self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': np.NaN}
                        else:
                            #lethalIndex = indexes[0]-1
                            lethalIndex = indexes[0]-1
                            lethalTime = self.results['TimeSpecificSummaryData'].to_pandas().index[lethalIndex] 
                            self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': lethalTime,'LethalIndex': lethalIndex}
                    if len(indexes) > 1:
                        baseline = peakutils.baseline(data)
                        index = np.argmax(data[indexes])
                        # Check that the peak is sufficiently distant from the baseline..
                        # Assess proportionate distance of peak from baseline.   
                        propDist = np.abs(baseline-data[indexes][index]).min()/np.mean(baseline)    
                        if propDist < 0.1:
                            self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': np.NaN}
                        else:
                            lethalIndex = indexes[index]-1
                            #lethalIndex = indexes[index]-1
                            lethalTime = self.results['TimeSpecificSummaryData'].to_pandas().index[lethalIndex] 
                            self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': lethalTime,'LethalIndex': lethalIndex}
                    if len(indexes) < 1:
                        self.lethalEndPoint_data[str(self.embryo)] = {'Embryo': self.embryo,'Dev stage': str(developmentalStage),'LethalTime': np.NaN,'LethalIndex': np.NaN}
        # Visualise and save output
        for keys,values in sorted(self.lethalEndPoint_data.items()):
            print(keys)
            print(values)
        np.save(savePath + '/lethalEndPoints.npy', self.lethalEndPoint_data)
            
        
    
