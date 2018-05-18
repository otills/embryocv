import cv2
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import pathos
import scipy
import scipy.fftpack
import xarray as xr
from scipy import stats
import scipy.signal as signal
import math
import warnings
import shutil
import glob
import numpy as np

class dataIntegration(object):    
    """
    This class generates summary stats from XArray Datasets and adds these to 
    the DataSet file. Two types of stats are generated:
        
        - global (statistics relevant to the entire period of the experiment
        e.g. growth rates)
        - timeSpecific (e.g. specific to a particular time point). The length of 
        timeSpecifc stats is equal to the duration of the experiment, with missing 
        values where either no data was present or stats cannot be calculated for
        some reason.
    """
#==============================================================================

#%%
#==============================================================================
#   Apply calculatePhenomeMeasures to all embryos
#==============================================================================
    def savePhenomeMeasuresForAllEmbryos(self, savePath,ignoreMeta = False):
        # Note can add ignoreMeta = True if metadata does not need to be copied to new dataset folder.
        # This can be useful if the results folder and file structure is abnormal.
        if not ignoreMeta:
            shutil.copy2(str(glob.glob(self.parentPath + "*.npy")[0]), savePath + "phenomeMetadata.npy")
            print self.embryoLabels
        for e in range(len(self.embryoLabels)):
            self.calculatePhenomeMeasures(e, savePath)
            
#==============================================================================
#   Integrate data to i)Produce phenome measures, ii)Save as XArray dataset and iii)Produce summary reports
#==============================================================================
    def calculatePhenomeMeasures(self, embryo, savePath):
        warnings.filterwarnings('ignore','All-NaN slice encountered')
        warnings.filterwarnings('ignore','Mean of empty slice')
        #   Load results
        ts = time.time()
        self.embryo = self.embryoLabels[embryo]
        self.loadResults()
        print 'Loading ' + self.embryoLabels[embryo]
        # If no/little data/tracking...
        if np.isnan(np.nanmean(self.results.ix[:,:,'area'], axis=0).astype(np.float64)).all():
            print self.embryoLabels[embryo] + ' contains no embryo data.'
            print 'This is most likely to have occured due to it not being located succesfully.'
            print 'Consequently no results file will be produced, but a .log file will be to serve as a record.' 
            with open(savePath + self.embryoLabels[embryo] + 'dataset.HDF5_log', 'w') as f:
                f.write('No data for ' + self.embryoLabels[embryo] + ".")
            savePath + self.embryoLabels[embryo] + 'dataset.HDF5_log'
            print self.embryoLabels[embryo] + "dataset.HDF5_log saved"
        elif np.sum(~np.isnan(np.nanmean(self.results.ix[:,:,'area'], axis=0).astype(np.float64))) < 5:
            print self.embryoLabels[embryo] + ' contains little embryo data (less than 5 time points).'
            print 'This is most likely to have occured due to it not being located succesfully.'
            print 'Consequently no results file will be produced, but a .log file will be to serve as a record.' 
            with open(savePath + self.embryoLabels[embryo] + 'dataset.HDF5_log', 'w') as f:
                f.write('Less than 5 time points containing data for ' + self.embryoLabels[embryo] + ".")
            savePath + self.embryoLabels[embryo] + 'dataset.HDF5_log'
            print self.embryoLabels[embryo] + "dataset.HDF5_log saved"
        else:
            ### Extract BLOCKWISE data and store as XArray
            blockWiseVals = self.results.ix[:,:,'blockWise'].values
            blockWiseArrayOne = np.zeros((blockWiseVals.shape[1],blockWiseVals.shape[0]))
            blockWiseArrayTwo = np.zeros((blockWiseVals.shape[1],blockWiseVals.shape[0],2,2))
            blockWiseArrayThree = np.zeros((blockWiseVals.shape[1],blockWiseVals.shape[0],4,4))
            blockWiseArrayFour = np.zeros((blockWiseVals.shape[1],blockWiseVals.shape[0],8,8))
            blockWiseArrayFive = np.zeros((blockWiseVals.shape[1],blockWiseVals.shape[0],16,16))
            blockWiseArrayFive = np.ascontiguousarray(blockWiseArrayFive)
            for t in np.arange(blockWiseVals.shape[1]):
                for f in range(blockWiseVals.shape[0]):  
                    if (type(blockWiseVals[f,t]) is list):
                        # Res 1x1
                        blockWiseArrayOne[t,f] = blockWiseVals[f,t][0]
                        # Res 2x2 : 16x16
                        blockWiseArrayTwo[t,f,:,:] = np.array(blockWiseVals[f,t][1])
                        blockWiseArrayThree[t,f,:,:] = np.array(blockWiseVals[f,t][2])
                        blockWiseArrayFour[t,f,:,:] = np.array(blockWiseVals[f,t][3])
                        blockWiseArrayFive[t,f,:,:] = np.array(blockWiseVals[f,t][4])
                    else:
                        blockWiseArrayOne[t,f] = np.NaN
                        blockWiseArrayTwo[t,f,:,:] = np.NaN
                        blockWiseArrayThree[t,f,:,:] = np.NaN
                        blockWiseArrayFour[t,f,:,:] = np.NaN
                        blockWiseArrayFive[t,f,:,:] = np.NaN
            self.blockWiseArrayOneXR = xr.DataArray(blockWiseArrayOne)
            self.blockWiseArrayTwoXR = xr.DataArray(blockWiseArrayTwo)
            self.blockWiseArrayThreeXR = xr.DataArray(blockWiseArrayThree)
            self.blockWiseArrayFourXR = xr.DataArray(blockWiseArrayFour)
            self.blockWiseArrayFiveXR = xr.DataArray(blockWiseArrayFive)  
    
            ### Generate FREQ data from Blockwise data and output as XArray.
            # Create output arrays and info necessary to make them
            frameRate = self.results.shape[1]/((np.float(self.results.ix[0,self.results.shape[1]-1,'elapsedTime'])-np.float(self.results.ix[0,0,'elapsedTime']))/1000)
            sampFreqs, powerSpect = signal.welch(self.blockWiseArrayOneXR.values[0,:],frameRate, scaling='spectrum', nfft=self.blockWiseArrayOneXR.shape[1])
            # Numpy arrays for other resolutions        
            self.freqOutputOne = np.zeros((self.results.shape[0],2,len(powerSpect)))
            self.freqOutputTwo = np.zeros((self.results.shape[0],2,2,2,len(powerSpect)))
            self.freqOutputThree = np.zeros((self.results.shape[0],4,4,2,len(powerSpect)))
            self.freqOutputFour = np.zeros((self.results.shape[0],8,8,2,len(powerSpect)))
            self.freqOutputFive = np.zeros((self.results.shape[0],16,16,2,len(powerSpect)))
            # Create a threading pool - limited to 12 as the default (on MacPro) cripples the computer (> 60GB RAM used).
            innerPool = pathos.multiprocessing.ThreadingPool(6)
            # Run freq analysis
            res = innerPool.map(self.parFreqAnalysis, range(self.freqOutputOne.shape[0]))
            # Unpack results
            for t in xrange(self.freqOutputOne.shape[0]):
                self.freqOutputOne[t,:,:] = res[t]['F1'][:,:]
                self.freqOutputTwo[t,:,:,:,:] = res[t]['F2'][:,:,:,:]
                self.freqOutputThree[t,:,:,:,:] = res[t]['F3'][:,:,:,:]
                self.freqOutputFour[t,:,:,:,:] = res[t]['F4'][:,:,:,:]
                self.freqOutputFive[t,:,:,:,:] = res[t]['F5'][:,:,:,:]
            # Store as XArrays
            self.freqOutputOneXR = xr.DataArray(self.freqOutputOne)
            self.freqOutputTwoXR = xr.DataArray(self.freqOutputTwo)
            self.freqOutputThreeXR = xr.DataArray(self.freqOutputThree)
            self.freqOutputFourXR = xr.DataArray(self.freqOutputFour)
            self.freqOutputFiveXR = xr.DataArray(self.freqOutputFive)
            ### Extract METADATA
            res = self.results
            meta = res.drop(['dateTime', 'elapsedTime', 'area', 'centroidX', 'centroidY',
                'solidity', 'aspect', 'extent', 'hullArea', 'bboxMincol',
                'bboxMinrow', 'bboxWidth', 'bboxHeight', 'embryoOutline',
                'eggRotBBox', 'eggBoxPoints', 'blockWise'], axis=2)
            #  Extract 'size & pos' data
            sizePos = self.results.ix[:,:,7:18].astype('float')
            sizePosArray = xr.DataArray(sizePos)
            # Extract ROTATED EGG BBOX  
            eggRotBB = self.results.ix[:,:,'eggRotBBox'].values
            eggRotBBArray = np.zeros((self.results.shape[1],self.results.shape[0],5))
            for t in range(self.results.shape[0]):
            	for f in range(self.results.shape[1]):
            		if type is list:
            			if (not math.isnan(float(eggRotBB[f,t][0]))):
            				eggRotBBArray[f,t,:] = eggRotBB[f,t]
            		else:
            			eggRotBBArray[f,t,:] = np.NaN
            eggRotBBArray = eggRotBBArray.swapaxes(0,1)
            eggRotBBXR = xr.DataArray(eggRotBBArray,{'Items':self.results.items,'Frame':self.results.major_axis,'EggRotBBox':['X','Y','W','H','A']})
            #  Extract SIZE AND POS data
            sizePos = self.results.ix[:,:,7:18].astype('float')
            sizePosArray = xr.DataArray(sizePos)    
             # Generate SUMMARY STATS
            timeSpStats, globalStats = self.generateSummaryStats()              
            #  NOW COMBINE and SAVE
            dateTime = self.results.ix[:,0,'dateTime']
            testDataSet = xr.Dataset({'Metadata':(['dateTime','frame','meta'],meta),
            'SizePos':(['dateTime','frame','sizePos'],sizePosArray),
            'EggBB':(['dateTime','frame','point','dim'],self.eggBBArray),
            'EmbryoOutline':(['dateTime','frame','dim','coord'],self.embryoOutlineXR),
            'EggRotBBox':(['dateTime','frame','eggRotBBox'],eggRotBBXR),
            'BlockWise_1x1':(['dateTime','frame'],self.blockWiseArrayOneXR),
            'BlockWise_2x2':(['dateTime','frame','X_2x2','Y_2x2'],self.blockWiseArrayTwoXR),
            'BlockWise_4x4':(['dateTime','frame','X_4x4','Y_4x4'],self.blockWiseArrayThreeXR),
            'BlockWise_8x8':(['dateTime','frame','X_8x8','Y_8x8'],self.blockWiseArrayFourXR),
            'BlockWise_16x16':(['dateTime','frame','X_16x16','Y_16x16'],self.blockWiseArrayFiveXR),
            'FreqOutput_1x1':(['dateTime','freq/Power','freqs'],self.freqOutputOneXR),
            'FreqOutput_2x2':(['dateTime','X_2x2','Y_2x2','freqPower','freqs'],self.freqOutputTwoXR),
            'FreqOutput_4x4':(['dateTime','X_4x4','Y_4x4','freqPower','freqs'],self.freqOutputThreeXR),
            'FreqOutput_8x8':(['dateTime','X_8x8','Y_8x8','freqPower','freqs'],self.freqOutputFourXR),
            'FreqOutput_16x16':(['dateTime','X_16x16','Y_16x16','freqPower','freqs'],self.freqOutputFiveXR),
            'TimeSpecificSummaryData':(['dateTime','timeSpecificMeasure'],timeSpStats.T),
            'GlobalSummaryData':(['globalMeasure', 'value'],xr.DataArray(pd.DataFrame(data = globalStats.values(), index= globalStats.keys())))},
            coords={'dateTime':dateTime,
            'frame':(list(np.arange(0,self.results.shape[1]))),
            'meta':(meta.minor_axis.values),
            'sizePos':(sizePosArray.coords[('dim_2')].values),
            'dim':['X','Y'],
            'point':['0','1','2','3'],
            'freqs':(list(np.arange(0,self.freqOutputFive.shape[4]))),
            'timeSpecificMeasure':['meanArea', 'minArea', 'maxArea','meanMinBB','meanMaxBB',
                                   'minMinBB','maxMaxBB','meanDistance','minDistance','maxDistance',
                                   'totalDistance','meanEggLength', 'meanEggWidth','Freq0-0.1Hz',
                                   'Freq0.1-0.3Hz','Freq0.3-0.5Hz','Freq0.5-0.7Hz',
                                   'Freq0.7-0.9Hz','Freq0.9-1.2Hz','Freq1.2-1.6Hz','Freq1.6-1.8Hz',
                                   'Freq1.8-2.2Hz','Freq2.2-3.0Hz','Freq3.0-4.0Hz','Freq4.0-5.0Hz', 
                                   'Freq5.0Hz-'],
            'globalMeasure':['MeanArGrowthRateSlope','MeanArGrowthRateIntercept','MeanArGrowthR',
                  'MeanArGrowthRatePval', 'MeanArGrowthSE','MinArGrowthRateSlope',
                  'MinArGrowthRateIntercept','MinArGrowthR','MinArGrowthRatePval',
                  'MinArGrowthSE','MaxArGrowthRateSlope','MaxArGrowthRateIntercept',
                  'MaxArGrowthR','MaxArGrowthRatePval','MaxArGrowthSE','MeanMinLenGrowthRateSlope', 
                  'MeanMinLenGrowthRateIntercept','MeanMinLenGrowthR','MeanMinLenGrowthRatePval',
                  'MeanMinLenGrowthSE','MinMinLenGrowthRateSlope','MinMinLenGrowthRateIntercept',
                  'MinMinLenGrowthR','MinMinLenGrowthRatePval','MinMinLenGrowthSE',
                  'MaxMaxLenGrowthRateSlope','MaxMaxLenGrowthRateIntercept',
                  'MaxMaxLenGrowthR','MaxMaxLenGrowthRatePval','MaxMaxLenGrowthSE']})                          
            testDataSet.to_netcdf(savePath + self.embryoLabels[embryo] + 'dataset.HDF5')
            print 'Saving ' + self.embryoLabels[embryo] + ' XArray Datset'
            print 'Phenome measures generated for ' + self.embryoLabels[embryo] + ' in {} s'.format(time.time()-ts)
            
#%%
#==============================================================================
#   Generate summary stats (called from calculatePhenomeMeasures) 
#==============================================================================   
    def generateSummaryStats(self):
        # Remove points containing no data..    
        minDistance = []
        maxDistance = []
        meanDistance = []
        totalDistance = []
        meanMinBB = []
        meanMaxBB = []
        minMinBB = []
        maxMaxBB = []
        mins = []
        maxs = []
        eggLength = []    
        eggWidth = []
        meanEggLength = []
        meanEggWidth = []
        #  EMBRYO OUTLINE       
        outlineVals = self.results.ix[:,:,'embryoOutline'].values
        arrLen = 10000
        lengths=[]
        outlineArray = np.zeros((self.results.shape[1],self.results.shape[0],2,arrLen))
        for t in range(self.results.shape[0]):
            for f in range(self.results.shape[1]):
                if (type(outlineVals[f,t]) is np.ndarray):
                    length = len(outlineVals[f,t][:,0,0])
                    lengths.append(length)
                    if length < arrLen:
                        outlineArray[f,t,0,0:length] = outlineVals[f,t][:,0,0]
                        outlineArray[f,t,1,0:length] = outlineVals[f,t][:,0,1]
                    else:
                        print 'Too many pixel coords for array'
                else:
                    outlineArray[f,t,:,:] = np.NaN
        # Reshape (time, frames..)
        outlineArray = outlineArray.swapaxes(0,1)
        # Catch instances when no embryo outlines are present
        if (len(lengths) > 0):
            # Crop to the maximum number of coordinates
            outlineArray = outlineArray[:,:,:,0:max(lengths)]
        else:
            outlineArray = outlineArray[:,:,:,0:1]
            outlineArray[:,:,:,:] = np.NaN
        # For some reason (prob simple mistake) cannot create a labelled XArray with these pixe data, however the default (non-labelled) creation below works fine.
        #embryoOutlineXR = xr.DataArray(outlineArray, {'Items':self.results.items,'Frame':self.results.major_axis,'embryoOutline':['X','Y'],'pixCoords':coord})
        self.embryoOutlineXR = xr.DataArray(outlineArray)
        eggBB = self.results.ix[:,:,'eggBoxPoints'].values
        #   EGG BOX POINTS
        eggArray = np.zeros((eggBB.shape[1],eggBB.shape[0],4,2))
        for t in range(eggBB.shape[1]):
            for f in range(eggBB.shape[0]):
                if (type(eggBB[f,t]) != list):
                    eggArray[t,f,:,:] = eggBB[f,t]
                else:
                    eggArray[t,f,:,:] = np.NaN   
        # Generate XArray
        self.eggBBArray = xr.DataArray(eggArray,{'Items':self.results.items,'Frame':self.results.major_axis,'Dim':['X','Y'],'Point':['0','1','2','3']})

        # Use a meshgrid to prevent the need for two loops to cycle over frames and timepoints.    
        frames, times = np.meshgrid(range(self.results.shape[1]),range(self.results.shape[0]))
        times = times.flatten()
        frames = frames.flatten()
        self.tmp = []              
        # Initiate tmpTime
        tmpTime = 0
        for t in xrange(len(times)):
        # Calculations not requiring looping over individual frames within this loop
          # Loop over times and if a new one is reached, proceed.
            if t is 0:
                # Calculate movement stats
                meanDisTmp, minDisTmp, maxDisTmp, totDisTmp = self.movementStats(times[t])
                meanDistance.append(meanDisTmp)
                minDistance.append(minDisTmp)
                maxDistance.append(maxDisTmp)
                totalDistance.append(totDisTmp)
                # Rotated bounding box dimensions
                mins = []
                maxs = []
                minDim, maxDim = self.calculateRotBBox(outlineArray,times[t], frames[t])
                mins.append(minDim)
                maxs.append(maxDim)
                # Egg length and width
                eggMin, eggMax = self.calculateEggRotBBox(eggArray,times[t], frames[t])
                eggWidth.append(eggMin)
                eggLength.append(eggMax)
            if times[t] != tmpTime:
                tmpTime = times[t]
                # Calculate movement stats
                meanDisTmp, minDisTmp, maxDisTmp, totDisTmp = self.movementStats(times[t])
                meanDistance.append(meanDisTmp)
                minDistance.append(minDisTmp)
                maxDistance.append(maxDisTmp)
                totalDistance.append(totDisTmp)
                self.tmp.append(times[t])
                # Rotated bounding box dimensions -  if moved onto a new time point, save the previous time data and start afresh.
                meanMinBB.append(np.nanmean(mins))
                meanMaxBB.append(np.nanmean(maxs))
                minMinBB.append(np.nanmin(mins))
                maxMaxBB.append(np.nanmax(maxs))            
                mins = []
                maxs = []
                minDim, maxDim = self.calculateRotBBox(outlineArray,times[t], frames[t])
                mins.append(minDim)
                maxs.append(maxDim)    
                # Egg length and width -  if moved onto a new time point, save the previous time data and start afresh.
                meanEggWidth.append(np.nanmean(eggWidth))
                meanEggLength.append(np.nanmean(eggLength))            
                eggWidth = []
                eggLength = []
                eggMin, eggMax = self.calculateEggRotBBox(eggArray,times[t], frames[t])
                eggWidth.append(eggMin)
                eggLength.append(eggMax)
            # If at the end of the last time point 
            if t == len(times)-1:
                # Add the last bounding box data to the BB stats.
                minDim, maxDim = self.calculateRotBBox(outlineArray,times[t], frames[t])
                mins.append(minDim)
                maxs.append(maxDim)
                meanMinBB.append(np.nanmean(mins))
                meanMaxBB.append(np.nanmean(maxs))
                minMinBB.append(np.nanmin(mins))
                maxMaxBB.append(np.nanmax(maxs))  
            else:
                # If not within a new timepoint continue calculating rotated bounding box dimensions and appending to list
                minDim, maxDim = self.calculateRotBBox(outlineArray,times[t], frames[t])
                mins.append(minDim)
                maxs.append(maxDim)
                eggMin, eggMax = self.calculateEggRotBBox(eggArray,times[t], frames[t])
                eggWidth.append(eggMin)
                eggLength.append(eggMax)
                
        # Calculate growth    
        # Area
        meanArGrowth, minArGrowth, maxArGrowth, meanArea, minArea, maxArea = self.calculateAreaGrowth()
        meanArGrowthSlope, meanArGrowthIntercept, meanArGrowthR, meanArGrowthPval,meanArGrowthSE = meanArGrowth
        minArGrowthSlope, minArGrowthIntercept, minArGrowthR, minArGrowthPval, minArGrowthSE = minArGrowth
        maxArGrowthSlope, maxArGrowthIntercept, maxArGrowthR, maxArGrowthPval, maxArGrowthSE = maxArGrowth      
        # Length
        meanMinLm, meanMaxLm, minMinLm, maxMaxLm = self.calculateLengthGrowth(meanMinBB,meanMaxBB,minMinBB,maxMaxBB)
        meanMinLenGrowthSlope, meanMinLenGrowthIntercept, meanMinLenGrowthR, meanMinLenGrowthPval,meanMinLenGrowthSE = meanMinLm
        meanMaxLenGrowthSlope, meanMaxLenGrowthIntercept, meanMaxLenGrowthR, meanMaxLenGrowthPval,meanMaxLenGrowthSE = meanMaxLm        
        minMinLenGrowthSlope, minMinLenGrowthIntercept, minMinLenGrowthR, minMinLenGrowthPval,minMinLenGrowthSE = minMinLm
        maxMaxLenGrowthSlope, maxMaxLenGrowthIntercept, maxMaxLenGrowthR, maxMaxLenGrowthPval,maxMaxLenGrowthSE = maxMaxLm
        # First find where is some data
        ind = np.argmax(~np.isnan(self.freqOutputFour[:,4,4,1,0]))
        bins = np.array([0,0.1,0.3,0.5,0.7,0.9,1.2,1.6,1.8,2.2,3.0,4.0,5.0,np.nanmax(self.freqOutputFour[:,:,:,0,:])])
        # Extract data for the first time point for which there is osme
        tmpData = self.freqOutputFour[ind,4,4,0,:]
        # Now use the bins specified above to create a list of labels applicable to each data point/frequency to 
        # inidcate to which bin it belongs
        binLabs = np.digitize(tmpData,bins)
        # Create some lists for saving output
        freqList_count = np.max(binLabs)
        freqList = [[] for i in range(0, freqList_count)]
        # Loop over the bins, calculate their energy and save to a list
        for i in range(len(bins)):
            # Create a filter to use to isolate the frequencies of relevance to the bin
            tmp = (binLabs==np.arange(1,len(binLabs))[i])
            if len(tmp) > 0:
                freqList[i].append(np.nansum(self.freqOutputFour[:,:,:,1,tmp],axis=(1,2,3)))
        # Debug
        self.freqList = freqList
        # Compile SummaryData output 
        timeSpecificData = pd.DataFrame(data = [meanArea, minArea, maxArea,
                                                meanMinBB,meanMaxBB,minMinBB,maxMaxBB,
                                                meanDistance,minDistance,maxDistance,
                                                totalDistance,meanEggLength, meanEggWidth,                                          
                                                list(np.squeeze(freqList[0])),list(np.squeeze(freqList[1])),
                                                list(np.squeeze(freqList[2])),list(np.squeeze(freqList[3])),
                                                list(np.squeeze(freqList[4])),list(np.squeeze(freqList[5])),
                                                list(np.squeeze(freqList[6])),list(np.squeeze(freqList[7])),
                                                list(np.squeeze(freqList[8])),list(np.squeeze(freqList[9])),
                                                list(np.squeeze(freqList[10])),list(np.squeeze(freqList[11])),
                                                list(np.squeeze(freqList[12])),list(np.squeeze(freqList[13]))],
                                                index = ['meanArea', 'minArea', 'maxArea',
                                                         'meanMinBB','meanMaxBB','minMinBB',
                                                         'maxMaxBB','meanDistance','minDistance',
                                                         'maxDistance','totalDistance',
                                                         'meanEggLength','meanEggWidth',
                                                         "FreqEnergy"+str(bins[0])+ ":"+ str(bins[1]),
                                                         "FreqEnergy"+str(bins[1])+ ":"+ str(bins[2]),
                                                         "FreqEnergy"+str(bins[2])+ ":"+ str(bins[3]),
                                                         "FreqEnergy"+str(bins[3])+ ":"+ str(bins[4]),
                                                         "FreqEnergy"+str(bins[4])+ ":"+ str(bins[5]),
                                                         "FreqEnergy"+str(bins[5])+ ":"+ str(bins[6]),
                                                         "FreqEnergy"+str(bins[6])+ ":"+ str(bins[7]),
                                                         "FreqEnergy"+str(bins[7])+ ":"+ str(bins[8]),
                                                         "FreqEnergy"+str(bins[8])+ ":"+ str(bins[9]),
                                                         "FreqEnergy"+str(bins[9])+ ":"+ str(bins[10]),
                                                         "FreqEnergy"+str(bins[10])+ ":"+ str(bins[11]),
                                                         "FreqEnergy"+str(bins[11])+ ":"+ str(bins[12]),
                                                         "FreqEnergy"+str(bins[12])+ ":"+ str(bins[13]),
                                                         "FreqEnergy"+str(bins[13])])
        # Debug
        self.timeSpData = timeSpecificData
        # Combine two last freq outputs to get relevant sums..
        timeSpecificData.ix[25,:] = timeSpecificData.ix[25,:] + timeSpecificData.ix[26,:]
        timeSpecificData = timeSpecificData.drop(timeSpecificData.index[len(timeSpecificData.index)-1])
        #timesp = timesp.drop('Freq5.0Hz-', axis=1)      
        # Add global data
        globalData = {'MeanArGrowthRateSlope': meanArGrowthSlope, 'MeanArGrowthRateIntercept': meanArGrowthIntercept, 'MeanArGrowthR': meanArGrowthR, 'MeanArGrowthRatePval': meanArGrowthPval, 'MeanArGrowthSE': meanArGrowthSE,
                      'MinArGrowthRateSlope': minArGrowthSlope, 'MinArGrowthRateIntercept': minArGrowthIntercept, 'MinArGrowthR': minArGrowthR, 'MinArGrowthRatePval': minArGrowthPval, 'MinArGrowthSE': minArGrowthSE,
                      'MaxArGrowthRateSlope': maxArGrowthSlope, 'MaxArGrowthRateIntercept': maxArGrowthIntercept, 'MaxArGrowthR': maxArGrowthR, 'MaxArGrowthRatePval': maxArGrowthPval, 'MaxArGrowthSE': maxArGrowthSE,
                      'MeanMinLenGrowthRateSlope': meanMinLenGrowthSlope, 'MeanMinLenGrowthRateIntercept': meanMinLenGrowthIntercept, 'MeanMinLenGrowthR': meanMinLenGrowthR, 'MeanMinLenGrowthRatePval': meanMinLenGrowthPval, 'MeanMinLenGrowthSE': meanMinLenGrowthSE,
                      'MinMinLenGrowthRateSlope': minMinLenGrowthSlope, 'MinMinLenGrowthRateIntercept': minMinLenGrowthIntercept, 'MinMinLenGrowthR': minMinLenGrowthR, 'MinMinLenGrowthRatePval': minMinLenGrowthPval, 'MinMinLenGrowthSE': minMinLenGrowthSE,
                      'MaxMaxLenGrowthRateSlope': maxMaxLenGrowthSlope, 'MaxMaxLenGrowthRateIntercept': maxMaxLenGrowthIntercept, 'MaxMaxLenGrowthR': maxMaxLenGrowthR, 'MaxMaxLenGrowthRatePval': maxMaxLenGrowthPval, 'MaxMaxLenGrowthSE': maxMaxLenGrowthSE}
                        
        return timeSpecificData, globalData
  
#%%
#==============================================================================
#   Quantify energy within different frequency bands in parallel 
#   (called from calculatePhenomeMeasures) generate freq data from blockwise data.
#==============================================================================
    def parFreqAnalysis(self,t):
        # Get necessary info..
        frameRate = self.results.shape[1]/((np.float(self.results.ix[0,self.results.shape[1]-1,'elapsedTime'])-np.float(self.results.ix[0,0,'elapsedTime']))/1000)
        sampFreqs, powerSpect = signal.welch(self.blockWiseArrayOneXR.values[t,:],frameRate, scaling='spectrum', nfft=self.blockWiseArrayOneXR.shape[1])  
        # Arrays for output
        ef1 = np.zeros((2,len(powerSpect)))
        ef2 = np.zeros((2,2,2,len(powerSpect)))
        ef3 = np.zeros((4,4,2,len(powerSpect)))
        ef4 = np.zeros((8,8,2,len(powerSpect)))
        ef5 = np.zeros((16,16,2,len(powerSpect)))
        if not np.isnan(self.blockWiseArrayOneXR.values[t,0]):
            # Create and run meshgrids to prevent need for nested loops..
            X2,Y2 = np.meshgrid(range(2),range(2))
            X2 = X2.flatten()
            Y2 = Y2.flatten()
            X3,Y3 = np.meshgrid(range(4),range(4))
            X3 = X3.flatten()
            Y3 = Y3.flatten()
            X4,Y4 = np.meshgrid(range(8),range(8))
            X4 = X4.flatten()
            Y4 = Y4.flatten()
            X5,Y5 = np.meshgrid(range(16),range(16))
            X5 = X5.flatten()
            Y5 = Y5.flatten()
            # Run analysis
            # 1x1
            ef1[0,:], ef1[1,:] = signal.welch(self.blockWiseArrayOneXR.values[t,:],frameRate, scaling='spectrum', nfft=self.results.shape[1])
            # 2x2
            for b in xrange(4):
                tmp = self.blockWiseArrayTwoXR.values
                ef2[X2[b],Y2[b],0,:], ef2[X2[b],Y2[b],1,:] = signal.welch(tmp[t,:,X2[b],Y2[b]],frameRate, scaling='spectrum', nfft=self.results.shape[1])
                # 4x4    
            for b in xrange(16):
                tmp = self.blockWiseArrayThreeXR.values
                ef3[X3[b],Y3[b],0,:], ef3[X3[b],Y3[b],1,:] = signal.welch(tmp[t,:,X3[b],Y3[b]],frameRate, scaling='spectrum', nfft=self.results.shape[1])
            # 8x8
            for b in xrange(64):
                tmp = self.blockWiseArrayFourXR.values
                ef4[X4[b],Y4[b],0,:], ef4[X4[b],Y4[b],1,:] = signal.welch(tmp[t,:,X4[b],Y4[b]],frameRate, scaling='spectrum', nfft=self.results.shape[1])
            # 16x16
            for b in xrange(256):
                tmp = self.blockWiseArrayFiveXR.values
                ef5[X5[b],Y5[b],0,:], ef5[X5[b],Y5[b],1,:] = signal.welch(tmp[t,:,X5[b],Y5[b]],frameRate, scaling='spectrum', nfft=self.results.shape[1])
        else:
            # If no blockwise signal assign NaN
            ef1[:,:] = np.NaN
            ef2[:,:,:,:] = np.NaN
            ef3[:,:,:,:] = np.NaN
            ef4[:,:,:,:] = np.NaN
            ef5[:,:,:,:] = np.NaN
        # Return numpy arrays in a dict for unpacking.
        return {'F1':ef1,'F2':ef2,'F3':ef3,'F4':ef4,'F5':ef5}

    #==============================================================================
    #   LENGTH GROWTH - Calculate mean, min and max slope stats (slope, intercept and pvalue)
    #==============================================================================
    def calculateLengthGrowth(self,meanMinBB,meanMaxBB,minMinBB,maxMaxBB):
        # Mask to remove nans Note: A rolling window is used here in developing the mask 
        # (and below) to remove outlying (temporally) data.
        mask = ~pd.Series(np.array(meanMinBB)).rolling(window=3).mean().isnull()
        if np.sum(~pd.Series(np.array(meanMinBB)).rolling(window=3).mean().isnull()) == 0:
            return float('nan'), float('nan'), float('nan'), float('nan') 
        # Date time is formatted appropriately
        filteredDateTimes = self.results.ix[:,0,'dateTime'].values[mask]
        DHMS = pd.to_datetime(filteredDateTimes)-pd.to_datetime(filteredDateTimes[0])
        fromExptStartMINS = np.array((DHMS/np.timedelta64(1,'m')).astype(np.float))
        # Generate growth data (using mask to remove nans, otherwise output for an embryo is NaN)       
        meanMinLenLm = scipy.stats.linregress(x=fromExptStartMINS, y = np.array(meanMinBB)[mask]*self.scale)
        meanMaxLenLm = scipy.stats.linregress(x=fromExptStartMINS, y = np.array(meanMaxBB)[mask]*self.scale)
        minMinLenLm = scipy.stats.linregress(x=fromExptStartMINS, y = np.array(minMinBB)[mask]*self.scale)        
        maxMaxLenLm = scipy.stats.linregress(x=fromExptStartMINS, y = np.array(maxMaxBB)[mask]*self.scale)        
        # Debug
        # plt.plot(fromExptStartMINS, np.array(self.meanMinBB)[mask]*self.scale)
        # plt.plot(fromExptStartMINS, meanMinLenLm.intercept + meanMinLenLm.slope*fromExptStartMINS, 'r', label='fitted line', color = 'black')
        return meanMinLenLm, meanMaxLenLm, minMinLenLm, maxMaxLenLm
        
    #==============================================================================
    #   AREA GROWTH - Calculate mean, min and max slope stats (slope, intercept and pvalue)
    #==============================================================================
    def calculateAreaGrowth(self):
        if self.mode is not 'xarray':
            # If too much missing data do not proceed
            if np.nansum(np.nanmean(self.results.ix[:,:,'area'], axis=0)) ==0:
            #np.sum(~pd.Series(np.nanmean(self.results.ix[:,:,'area'],axis=0)).rolling(window =3).mean().isnull()) ==0:
                return (float('nan'), float('nan'), float('nan'), float('nan'), float('nan')), (float('nan'), float('nan'),float('nan'), float('nan'), float('nan')), (float('nan'), float('nan'), float('nan'), float('nan'), float('nan')), [float('nan')]*self.results.shape[0], [float('nan')]*self.results.shape[0], [float('nan')]*self.results.shape[0]
            else:
                # Date time is formatted appropriately
                filteredDateTimes = self.results.ix[:,0,'dateTime'].values
                DHMS = pd.to_datetime(filteredDateTimes)-pd.to_datetime(filteredDateTimes[0])
                fromExptStartMINS = np.array((DHMS/np.timedelta64(1,'h')).astype(np.float))                           
                # Raw data - note conversion with metadata scale
                meanArea = np.nanmean(self.results.ix[:,:,'area'].values.astype(np.float), axis = 0)*(self.scale*self.scale)
                minArea = np.nanmin(self.results.ix[:,:,'area'].values.astype(np.float), axis = 0)*(self.scale*self.scale)
                maxArea = np.nanmax(self.results.ix[:,:,'area'].values.astype(np.float), axis = 0)*(self.scale*self.scale)
                
                # Generate growth data, including a scale conversion
                noDataMask = ~np.isnan(meanArea)
                meanLm = scipy.stats.linregress(x=fromExptStartMINS[noDataMask], y = np.log(meanArea)[noDataMask])
                minLm = scipy.stats.linregress(x=fromExptStartMINS[noDataMask], y = np.log(minArea)[noDataMask]) 
                maxLm = scipy.stats.linregress(x=fromExptStartMINS[noDataMask], y = np.log(maxArea)[noDataMask]) 
                
                # Debug
                # Plot data
                # plt.figure()
                # plt.plot(fromExptStartMINS, np.log(np.array(meanArea).astype(float)), color = 'black')
                # plt.plot(fromExptStartMINS,np.log(np.array(minArea).astype(float)), color = 'blue')
                # plt.plot(fromExptStartMINS,np.log(np.array(maxArea).astype(float)), color = 'orange')
                # Add model to plot
                # plt.plot(fromExptStartMINS, meanLm.intercept + meanLm.slope*fromExptStartMINS, 'r', label='fitted line', color = 'black')
                # plt.plot(fromExptStartMINS, minLm.intercept + minLm.slope*fromExptStartMINS, 'r', label='fitted line', color = 'blue')
                # plt.plot(fromExptStartMINS, maxLm.intercept + maxLm.slope*fromExptStartMINS, 'r', label='fitted line', color = 'orange')
                
        # Output growth model stats and raw data
        if self.mode is 'xarray':
            if np.nansum(np.nanmean(self.results['SizePos'].loc[:,:,'area'], axis=0)) ==0:
                #np.sum(~pd.Series(np.nanmean(self.results.ix[:,:,'area'],axis=0)).rolling(window =3).mean().isnull()) ==0:
                return (float('nan'), float('nan'), float('nan'), float('nan'), float('nan')), (float('nan'), float('nan'),float('nan'), float('nan'), float('nan')), (float('nan'), float('nan'), float('nan'), float('nan'), float('nan')), [float('nan')]*self.results.shape[0], [float('nan')]*self.results.shape[0], [float('nan')]*self.results.shape[0]
            else:
                #print 'here'
                ### Incorporate a mask to remove very high solidities, likely the entire egg - i.e. issues with segmentation?
                solidityMask = np.nanmean(self.results['SizePos'].loc[:,:,'solidity'].values,axis=1) < 0.99
                meanArea = np.nanmean(self.results['SizePos'].loc[:,:,'area'].values, axis=1)
                meanArea[~solidityMask] = np.NaN
                meanArea = meanArea*(self.scale*self.scale)
                minArea = np.nanmean(self.results['SizePos'].loc[:,:,'area'].values, axis=1)
                minArea[~solidityMask] = np.NaN
                minArea = minArea*(self.scale*self.scale)
                maxArea = np.nanmax(self.results['SizePos'].loc[:,:,'area'].values, axis=1)
                maxArea[~solidityMask] = np.NaN
                maxArea = maxArea*(self.scale*self.scale)
                # Date time is formatted appropriately
                filteredDateTimes = self.results['dateTime'].values
                DHMS = pd.to_datetime(filteredDateTimes)-pd.to_datetime(filteredDateTimes[0])
                fromExptStartMINS = np.array((DHMS/np.timedelta64(1,'h')).astype(np.float))             
                # Generate growth data, including a scale conversion
                noDataMask = ~np.isnan(meanArea)
                meanLm = scipy.stats.linregress(x=fromExptStartMINS[noDataMask], y = np.log(meanArea)[noDataMask])
                minLm = scipy.stats.linregress(x=fromExptStartMINS[noDataMask], y = np.log(minArea)[noDataMask]) 
                maxLm = scipy.stats.linregress(x=fromExptStartMINS[noDataMask], y = np.log(maxArea)[noDataMask]) 
                # Debug
                # Plot data
                plt.figure()
                plt.plot(fromExptStartMINS, np.log(np.array(meanArea).astype(float)), color = 'black')
                plt.plot(fromExptStartMINS,np.log(np.array(minArea).astype(float)), color = 'blue')
                plt.plot(fromExptStartMINS,np.log(np.array(maxArea).astype(float)), color = 'orange')
                # Add model to plot
                plt.plot(fromExptStartMINS, meanLm.intercept + meanLm.slope*fromExptStartMINS, 'r', label='fitted line', color = 'black')
                plt.plot(fromExptStartMINS, minLm.intercept + minLm.slope*fromExptStartMINS, 'r', label='fitted line', color = 'blue')
                plt.plot(fromExptStartMINS, maxLm.intercept + maxLm.slope*fromExptStartMINS, 'r', label='fitted line', color = 'orange')
                plt.title(str(self.embryo))
                # Output growth model stats and raw data
                return (meanLm.slope, meanLm.intercept, meanLm.rvalue, meanLm.pvalue, meanLm.stderr), (minLm.slope, minLm.intercept, minLm.rvalue, minLm.pvalue, minLm.stderr), (maxLm.slope, maxLm.intercept, maxLm.rvalue, maxLm.pvalue, maxLm.stderr), meanArea, minArea, maxArea
           
   #==============================================================================
   #   Determine rotated bounding box of embryo (to get min and max dimension)
   #==============================================================================
    def calculateRotBBox(self,coords,timesT,framesT):
        if (np.mean(coords[timesT,framesT,0,:]) > 0):
            # Get rotated bounding box
            (_,_),(w,h),_ = cv2.minAreaRect(np.array((coords[timesT,framesT,0,coords[timesT,framesT,0,:]>0],coords[timesT,framesT,1,coords[timesT,framesT,1,:]>0])).astype(np.int).T)
            minDim = np.min((w,h))
            maxDim = np.max((w,h))
        else:
            minDim = np.NaN
            maxDim = np.NaN
        return minDim, maxDim
       
   #==============================================================================
   #   Determine rotated bounding box of egg (to get min and max dimension)
   #==============================================================================
    def calculateEggRotBBox(self,eggBB,timesT,framesT):
        if (np.mean(eggBB[timesT,framesT,0,:]) > 0):
            # Get rotated bounding box
            (_,_),(w,h),_ = cv2.minAreaRect(np.array(eggBB[timesT,framesT,:,:].astype(np.int)))
            minDim = np.min((w,h))
            maxDim = np.max((w,h))
        else:
            minDim = np.NaN
            maxDim = np.NaN
        return minDim, maxDim

                   
    #==============================================================================
    #     MOVEMENT
    #==============================================================================
    def movementStats(self,t):
        # If not more than 10 data points missing for the time point then proceed with movement stat calculation.
        movementDataX = self.results.loc[:,:,'centroidX'].values.astype(np.float)
        movementDataY = self.results.loc[:,:,'centroidY'].values.astype(np.float)
        
        if (np.sum(np.isnan(movementDataX[:,t])) <10):
            distances = np.linalg.norm(np.array((movementDataX[:,t],movementDataY[:,t])).T[:-1] - np.array((movementDataX[:,t],movementDataY[:,t])).T[1:], axis=1)
            meanDist = np.nanmean((distances*self.scale))  
            minDist = np.nanmin((distances*self.scale))  
            maxDist = np.nanmax((distances*self.scale))  
            totalDist = np.nansum((distances*self.scale))  
        else:
            # If too much missing data then assign NaN
            meanDist = np.NaN
            minDist = np.NaN
            maxDist = np.NaN
            totalDist = np.NaN
        return meanDist, minDist, maxDist, totalDist    
    
