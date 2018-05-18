# Import dependencies
import cv2
import numpy as np
import pandas as pd
import glob
import os
import sys
import eggUI
import viewOutput
import time
import pathos
import json
import xarray as xr
from PyQt5.Qt import *
from PyQt5 import QtGui

class dataHandling(object):

#==============================================================================
#   Produce a list of embryos with results fileslarger than a particular size.
#==============================================================================
    def filterResultsList(self,fileSize):
        #self.resultsDir = os.path.dirname(self.parentPath + "phenomeData/")
        self.sizeFilteredEmbryos = []
        for f in range(len(self.embryoLabels)):
            self.embryo = self.embryoLabels[f]
            resultsStructure = self.resultsDir +'/' + self.embryo + '.pandas'
            size = os.path.getsize(resultsStructure)/1000000
            if size >= fileSize:
                self.sizeFilteredEmbryos.append(self.embryo)
                
#==============================================================================
#   Function to identify the egg 
#==============================================================================
    def seqImport(self, n):
        self.seq = np.zeros(shape=(self.results.shape[1],self.imImport(0,0).shape[0],self.imImport(0,0).shape[1]))
        for f in range(self.results.shape[1]):
            self.seq[f] = cv2.imread(self.results.iloc[n]['currentFolder'][f] + self.results.iloc[n]['file'][f],cv2.IMREAD_ANYDEPTH)
        ran = (self.seq.max() - self.seq.min()) / 255.
        self.seq = self.seq/ran
        self.seq = self.seq-self.seq.min()
        self.seq = np.ascontiguousarray(self.seq.astype(np.uint8))
        #np.ascontiguousarray
#==============================================================================
#   Import single image using results
#==============================================================================
    def imImport(self,n,m):
        im = cv2.imread(self.results.iloc[n]['currentFolder'][m] + self.results.iloc[n]['file'][m],cv2.IMREAD_ANYDEPTH)
        ran = (im.max()-im.min())/255.
        out = (im/ran)
        out = out-out.min()
        out = out.astype(np.uint8)
        return out

#==============================================================================
#   Get embryo labels from a parentPath        
#==============================================================================
    def getEmbryoLabels(self,parentPath):
        folders = glob.glob(parentPath + "*/*/")
        # Trouble with getctime .. doesn't seem well supported. Therefore switched.
        folders.sort(key=os.path.getmtime)  
        if len(folders)<1:
            print 'Data not found.'
        embryoLabels = []
        for p in range(len(folders)):
            embryoLabels.append(os.path.basename(os.path.normpath((folders[p]))))
        if self.species == 'rbalthica':
            # Reduce to only the unique labels
            embryoLabels = np.unique(embryoLabels)
            self.embryoLabels = embryoLabels[embryoLabels != 'BLANK']       
            # Use to remove 'problematic embryos'
            if self.exclude is not 'na':
                if type(self.exclude) is str:
                    self.embryoLabels = self.embryoLabels[self.embryoLabels != self.exclude]  
                else:
                    self.embryoLabels = self.embryoLabels
                    for p in range(len(self.exclude)):
                        self.embryoLabels = self.embryoLabels[self.embryoLabels != str(self.exclude[p])]  
        elif self.species == 'ogammarellus':
            # Add a filter to deal with 'copied data - such as Orchestia methods MS data'
            labs = []
            for p in range(len(folders)):
                labs.append(str(os.path.basename(os.path.normpath((folders[p])))).split(' ')[0])
            # Reduce to only the unique labels
            embryoLabels = np.unique(labs)
            self.embryoLabels = embryoLabels[embryoLabels != 'BLANK']       
#==============================================================================
#   Get embryo labels from a parentPath        
#==============================================================================
    def getEmbryoLabelsFromResultsFolder(self,resultsDir):
        if self.dataformat:
            folders = glob.glob(resultsDir + "/*.pandas")
            embryoLabels = []
            resultsDir +"/"
            for f in range(len(folders)):
                embryoLabels.append(folders[f].replace(".pandas", ""))
                #embryoLabels[f] = embryoLabels[f].replace(resultsDir +"/","")
                embryoLabels[f] = embryoLabels[f].replace(resultsDir,"")
            # Reduce to only the unique labels (shouldn't be a problem as from results).
            embryoLabels = np.unique(embryoLabels)
            self.embryoLabels = embryoLabels[embryoLabels != 'BLANK']    
             # Use to remove 'problematic embryos'
            if self.exclude is not 'na':
                if type(self.exclude) is str:
                    self.embryoLabels = self.embryoLabels[self.embryoLabels != self.exclude]  
                else:
                    self.embryoLabels = self.embryoLabels
                    for p in range(len(self.exclude)):
                        self.embryoLabels = self.embryoLabels[self.embryoLabels != str(self.exclude[p])]      
        else: 
            # If raw does not equal True (i.e. if data do not follow the typical raw MicroManager format
            # due perhaps to being copied from elesewhere).
            # For example Parent Folder/EmbryoA/Time1
            #                                  /Time2
            #                           /EmbryoB/Time1
            #                                  /Time2    
            
            folders = glob.glob(resultsDir + "/*.pandas")
            embryoLabels = []
            resultsDir +"/"
            for f in range(len(folders)):
                embryoLabels.append(folders[f].replace(".pandas", ""))
                embryoLabels[f] = embryoLabels[f].replace(resultsDir +"/","")
            # Reduce to only the unique labels (shouldn't be a problem as from results).
            embryoLabels = np.unique(embryoLabels)
            self.embryoLabels = embryoLabels[embryoLabels != 'BLANK']       
        # Debug
        print self.embryoLabels   
        # Debug
        # print self.embryoLabels      
        
#==============================================================================
#   Get embryo labels from a parentPath with XARRAY data 
#==============================================================================
    def getXREmbryoLabelsFromResultsFolder(self,resultsDir):
        # If data in normal MM format..
        folders = glob.glob(resultsDir + "/*.HDF5")
        embryoLabels = []
        resultsDir +"/"
        for f in range(len(folders)):
            embryoLabels.append(folders[f].replace(".HDF5", ""))
            embryoLabels.append(folders[f].replace("dataset", ""))
            embryoLabels[f] = embryoLabels[f].replace(resultsDir +"/","")
        # Reduce to only the unique labels (shouldn't be a problem as from results).
        embryoLabels = np.unique(embryoLabels)
        self.embryoLabels = embryoLabels[embryoLabels != 'BLANK']       
        # Debug
        # print self.embryoLabels          
        
    
#==============================================================================
#   Get folders for a particular embryo
#==============================================================================
    def getEmbryoFolders(self, parentPath, embryo):
        if self.mode == 'results':
            self.parentPath = parentPath
            self.embryo = embryo
            self.embryoFolders = self.results.ix[:,0,'currentFolder']
            self.getShortenedPaths()
        if (self.mode == 'resume'):
            self.parentPath = parentPath
            self.embryo = embryo
            if self.species == 'rbalthica':
                #self.embryoFolders = self.results.items
                # Oli changed (1211) to deal with corrupted HD (2017018)
                #self.embryoFolders = glob.glob(parentPath + "*/" + embryo +"/")
                self.embryoFolders = self.results.ix[:,0,'currentFolder']
                self.embryoFolders = list(self.embryoFolders)
                self.embryoFolders.sort(key=os.path.getmtime)
                #self.embryoFolders.sort(key=os.path.getmtime) 
                self.getShortenedPaths()
            elif self.species == 'ogammarellus':
                self.embryoFolders = self.results.items
                self.embryoFolders = glob.glob(parentPath + "*/" + embryo +" */")
                self.embryoFolders.sort(key=os.path.getmtime) 
                self.getShortenedPaths()
        if (self.mode == 'new'):
            if self.species == 'rbalthica':
                self.parentPath = parentPath
                self.embryo = embryo
                self.embryoFolders = glob.glob(parentPath + "*/" + embryo +"/")
                self.embryoFolders.sort(key=os.path.getmtime) 
                self.getShortenedPaths()
            elif self.species == 'ogammarellus':
                #self.embryoFolders = self.results.items
                self.embryoFolders = glob.glob(parentPath + "*/" + embryo +" */")
                self.embryoFolders.sort(key=os.path.getmtime) 
                self.getShortenedPaths()
#==============================================================================
#   Get list of shortenedPaths  
#==============================================================================
    def getShortenedPaths(self):
        if (self.mode == 'results')|(self.mode == 'resume'):
            self.shortenedPaths = self.results.items
        if self.mode == 'new': 
            self.shortenedPaths=[]
            for f in range(len(self.embryoFolders)):
                self.shortenedPaths.append(os.path.relpath(self.embryoFolders[f], self.parentPath))
                
#==============================================================================
#   Create a directory for saving results
#==============================================================================
    def createResultsFolder(self):
        self.resultsDir = os.path.dirname(self.parentPath + "phenomeData/")
        if not os.path.exists(self.resultsDir):
            os.makedirs(self.resultsDir)
            
#==============================================================================
#   Generate a results panel suitable for experiment AND then Run egg ID    
#==============================================================================
    def generateResultsAndFindEggs(self, parentPath, scale, eggInt=1234):
        self.eggInt = eggInt
        totts = time.time()
        self.compiledData = {}
        self.parentPath = parentPath
        self.scale = scale
        self.getEmbryoLabels(self.parentPath)
        self.createResultsFolder()
        for e in range(len(self.embryoLabels)):
            ts = time.time()
            print 'Initiation started for', self.embryoLabels[e]
            self.getEmbryoFolders(self.parentPath,self.embryoLabels[e])
            for f in range(len(self.embryoFolders)):
                print f, self.embryoFolders[f]
                self.currentFolder = self.embryoFolders[f]
                self.shortenedPath = self.shortenedPaths[f]
                self.parseMetadata()
                self.createResultsTable()
                self.compiledData[self.shortenedPath] = self.results
            print 'Initiation complete for', self.embryoLabels[e], ('in {} s'.format(time.time()-ts))    
            self.resultSheets = pd.Panel.from_dict(self.compiledData)
            self.resultSheets.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')
            self.results=[]
            self.resultSheets=[]
            self.compiledData={}
        self.saveMetadata()
        self.runEggID(self.eggInt)  
        print 'Egg identification and creation of results files is now complete. '
        print 'This took ', ('{} s'.format(time.time()-totts))


#==============================================================================
#   Just generate results - route for O. gammarellus analysis. Initiated upon creating
#   a 'new' experiment.
#
#   Generate a results panel suitable for experiment AND then Run egg ID    
#==============================================================================
    def generateResults(self, parentPath, scale):
        totts = time.time()
        self.compiledData = {}
        self.parentPath = parentPath
        self.scale = scale
        self.getEmbryoLabels(self.parentPath)
        self.createResultsFolder()
        for e in range(len(self.embryoLabels)):
            ts = time.time()
            print 'Initiation started for', self.embryoLabels[e]
            self.getEmbryoFolders(self.parentPath,self.embryoLabels[e])
            for f in range(len(self.embryoFolders)):
                print f, self.embryoFolders[f]
                self.currentFolder = self.embryoFolders[f]
                self.shortenedPath = self.shortenedPaths[f]
                self.parseMetadata()
                self.createResultsTable()
                self.compiledData[self.shortenedPath] = self.results
            print 'Results file created for ', self.embryoLabels[e], ('in {} s'.format(time.time()-ts))    
            self.resultSheets = pd.Panel.from_dict(self.compiledData)
            self.resultSheets.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')
            self.results=[]
            self.resultSheets=[]
            self.compiledData={}
        self.saveMetadata()
        #self.runEggID(self.eggInt)  
        print 'This took ', ('{} s'.format(time.time()-totts))
        print 'Egg identification and creation of results files is now complete. Use instance.analyseAllEmbryos function to continue'

#==============================================================================
#   Create results structure
#==============================================================================
    def createResultsTable(self):
        global da
        # Store embryo outline as numpy array, but results as panda.
        cols = ['embryo','currentFolder','parentPath','file','scale','UUID','dateTime','elapsedTime',
                'area','centroidX','centroidY','solidity','aspect','extent',
                'hullArea','bboxMincol','bboxMinrow','bboxWidth','bboxHeight', 'embryoOutline'
                ,'eggRotBBox','eggBoxPoints', 'blockWise']
        #, 'eggBBox'
        self.results = pd.DataFrame(index = range(int(len(self.filtFiles))), columns = cols)
        self.results.embryo = self.embryo
        self.results.currentFolder = self.currentFolder
        self.results.file = self.filtFiles
        self.results.parentPath = self.parentPath
        self.results.scale = self.scale
        self.results.dateTime = self.filtTimes
        self.results.elapsedTime = self.filtElapsedTimes
        self.results.UUID = self.filtUUID

#==============================================================================
# Load metadata      
#============================================================================== 
    def parseMetadata(self):   
        # Search currentFolder for metadata
        with open(glob.glob(self.currentFolder + "/*.txt")[0], 'r') as f: 
            jsonDict = json.loads(f.read())       
        times = []
        elapsedTimes = []
        files = []
        UUID = []
        for i in jsonDict.keys():
            if i!=u'Summary':
                times.append(jsonDict[i]['Time'])
                elapsedTimes.append(jsonDict[i]['ElapsedTime-ms'])
                files.append(jsonDict[i]['FileName'])
                UUID.append(jsonDict[i]['UUID'])
        meta = [(y,x,z,a) for (y,x,z,a) in sorted(zip(times,elapsedTimes,files,UUID))]
        self.filtTimes = []
        self.filtElapsedTimes = []
        self.filtFiles = []
        self.filtUUID = []
        for i in range(len(meta)):
            self.filtTimes.append(str(meta[i][0]))
            self.filtElapsedTimes.append(str(meta[i][1]))
            self.filtFiles.append(str(meta[i][2]))
            self.filtUUID.append(str(meta[i][3]))
         
#==============================================================================
#   Apply locateEgg to image sequence at specified interval
#==============================================================================
    def getIntervIndicesFromSequence(self, n):
        # n = Specifies how many frames to skip in between the locateEgg.
        #self.intN = self.results.shape[1]-1
        self.intN = n
        self.eggIDIms = []
        # If n = 1234 (default), only the first image from each sequence is taken.
        if self.intN == 1234:
            self.eggIDIms.append(0)
        # If an interval is specified...
        else:
            seqLength = self.results.shape[1]
            # Get number of IDs necessary, based on n and seqLength.
            eggIDNum = seqLength/n   
            # Get frame indices
            for i in range(eggIDNum):
                 self.eggIDIms.append(i*self.intN)        

#==============================================================================
#   Import a non continuous sequence             
#==============================================================================
    def nonContSeqImport(self, f):
        # Create np stack of appropriate size
        im = self.imImport(0,0)
        self.seq = np.zeros(shape=(int(len(self.eggIDIms)),im.shape[0],im.shape[1]))      
        # Loop over the eggIDIms
        for g in range(len(self.eggIDIms)):
            # Import image one by one, from correct folder (argument) and iterated eggIDIm.
            self.seq[g] = self.imImport(f,self.eggIDIms[g])
        self.seq = self.seq.astype(np.uint8)

#==============================================================================
# Load results for current embryo    
#==============================================================================
    def loadResults(self):
        # Iterate over embryo labels...
        resultsStructure = self.resultsDir +'/' + self.embryo + '.pandas'
        results = pd.read_pickle(resultsStructure)        
        # Reorder panels so they are chronological, using the first datetime from each panel.
        dates = pd.to_datetime(results.ix[:,0,6])
        # Sort them
        dates.sort_values(inplace=True)
        # Reindex axis
        self.results = results.reindex_axis(dates.index, copy='False')

#==============================================================================
#   Load results for current embryo - XARRAY
#==============================================================================
    def loadXRResults(self):
        # Iterate over embryo labels...
        if self.mode is 'xarray':
            if float(str(glob.glob(self.parentPath + self.embryo + "dataset*")[0]).find('log')) > 0:
                self.results = 'NoData'
            else:
                self.results = xr.open_mfdataset(self.parentPath + self.embryo + 'dataset.HDF5')

#==============================================================================
#   Load, but return results
#==============================================================================
    def returnResults(self, embryo):
        # Iterate over embryo labels...
        resultsStructure = self.resultsDir +'/' + embryo + '.pandas'
        results = pd.read_pickle(resultsStructure)        
        # Reorder panels so they are chronological, using the first datetime from each panel.
        dates = pd.to_datetime(results.ix[:,0,6])
        # Sort them
        dates.sort_values(inplace=True)
        # Reindex axis
        results = results.reindex_axis(dates.index, copy='False')
        return results
    
#==============================================================================
#   Apply egg ID
#==============================================================================
    # Run egg ID and add to the results table before saving
    # If no eggInt given then just run on the first image of each sequence
    def runEggID(self, eggInt = 1234):
        self.eggInt = eggInt
        for e in range(len(self.embryoLabels)):
            self.embryo = self.embryoLabels[e]
            self.getEmbryoFolders(self.parentPath, self.embryo)
            self.loadResults()
            #if self.eggInt ==1234:
                #self.eggIDIms = []
                #self.eggIDIms.append(0)
            self.getIntervIndicesFromSequence(int(self.eggInt))
            for f in range(len(self.embryoFolders)):
                print self.embryoFolders[f]
                self.nonContSeqImport(f)
                for g in range(len(self.eggIDIms)):
                    # If g = [0] i.e. eggInt = 1234, make g = 0.
                    if len(self.eggIDIms) == 1:
                        g =0
                    # NOW perform the eggID on each intervalled image
                    self.locateEgg(self.seq[g])
                    # AND store eggID output in results table
                    #self.results.ix[f,self.eggIDIms[g],'eggBBox'] = self.eggBBox
                    self.results.ix[f,self.eggIDIms[g],'eggRotBBox'] = self.eggRotBBox
                    self.results.ix[f,self.eggIDIms[g],'eggBoxPoints'] = self.boxPoints     
            self.results.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')
            
#==============================================================================
#   Fill missing egg measurements with ones from either earlier or later time points 
#   in the data panel series. Note: This is a suboptimal solution and should be 
#   avoided via manual checking using eggUI.
#==============================================================================
    def fillEggMeasurements(self):
        # Fill non calculated values frames from most nrecent previously calculated.
        for f in range(self.results.shape[0]):
             #self.results.ix[f,:,'eggBBox'] = self.results.ix[f,:,'eggBBox'].fillna(method = 'ffill')
             #self.results.ix[f,:,'eggBBox'] = self.results.ix[f,:,'eggBBox'].fillna(method = 'bfill')
             self.results.ix[f,:,'eggBoxPoints'] = self.results.ix[f,:,'eggBoxPoints'].fillna(method = 'ffill')
             self.results.ix[f,:,'eggBoxPoints'] = self.results.ix[f,:,'eggBoxPoints'].fillna(method = 'bfill')
             self.results.ix[f,:,'eggRotBBox'] = self.results.ix[f,:,'eggRotBBox'].fillna(method = 'ffill')
             self.results.ix[f,:,'eggRotBBox'] = self.results.ix[f,:,'eggRotBBox'].fillna(method = 'bfill')
             self.results.ix[f,:,'bboxMincol'] = self.results.ix[f,:,'bboxMincol'].fillna(method = 'ffill')
             self.results.ix[f,:,'bboxMinrow'] = self.results.ix[f,:,'bboxMinrow'].fillna(method = 'bfill')
        # If no value in first datapoint in first panel..
        for f in range(self.results.shape[0]):    
            null = self.results.ix[:,0,'eggBoxPoints'].isnull()
            if null[f]:
                for g in range(len(null)-f): 
                    if ~null[(f+g)]:
                        #self.results.ix[f,0,'eggBBox'] = self.results.ix[(f+g),0,'eggBBox']
                        #self.results.ix[f,:,'eggBBox'] = self.results.ix[f,:,'eggBBox'].fillna(method = 'ffill')
                        self.results.ix[f,0,'eggBoxPoints'] = self.results.ix[(f+g),0,'eggBoxPoints']
                        self.results.ix[f,:,'eggBoxPoints'] = self.results.ix[f,:,'eggBoxPoints'].fillna(method = 'ffill')
                        self.results.ix[f,0,'eggRotBBox'] = self.results.ix[(f+g),0,'eggRotBBox']
                        self.results.ix[f,:,'eggRotBBox'] = self.results.ix[f,:,'eggRotBBox'].fillna(method = 'ffill')
                        self.results.ix[f,0,'bboxMincol'] = self.results.ix[(f+g),0,'bboxMincol']
                        self.results.ix[f,:,'bboxMincol'] = self.results.ix[f,:,'bboxMincol'].fillna(method = 'ffill')        
                        self.results.ix[f,0,'bboxMinrow'] = self.results.ix[(f+g),0,'bboxMinrow']
                        self.results.ix[f,:,'bboxMinrow'] = self.results.ix[f,:,'bboxMinrow'].fillna(method = 'ffill')                               
                        break
        # Now repeat process but filling later panels in the dataframe with data from earlier panels.
        for f in range(self.results.shape[0]):   
            null = self.results.ix[:,0,'eggBoxPoints'].isnull()
            if null[f]:
                #print f
                for g in range(f,0,-1):
                    if ~null[g]:
                        #self.results.ix[f,0,'eggBBox'] = self.results.ix[g,0,'eggBBox']
                        #self.results.ix[f,:,'eggBBox'] = self.results.ix[f,:,'eggBBox'].fillna(method = 'ffill')
                        self.results.ix[f,0,'eggBoxPoints'] = self.results.ix[g,0,'eggBoxPoints']
                        self.results.ix[f,:,'eggBoxPoints'] = self.results.ix[f,:,'eggBoxPoints'].fillna(method = 'ffill')
                        self.results.ix[f,0,'eggRotBBox'] = self.results.ix[g,0,'eggRotBBox']
                        self.results.ix[f,:,'eggRotBBox'] = self.results.ix[f,:,'eggRotBBox'].fillna(method = 'ffill')
                        self.results.ix[f,0,'bboxMincol'] = self.results.ix[g,0,'bboxMincol']
                        self.results.ix[f,:,'bboxMincol'] = self.results.ix[f,:,'bboxMincol'].fillna(method = 'ffill')
                        self.results.ix[f,0,'bboxMinrow'] = self.results.ix[g,0,'bboxMinrow']
                        self.results.ix[f,:,'bboxMinrow'] = self.results.ix[f,:,'bboxMinrow'].fillna(method = 'ffill')                        
                        break

#==============================================================================
#   Apply embryo segmentation to all embryos
#==============================================================================
    def segmentAllEmbryos(self):
        for e in range(len(self.embryoLabels)): 
            print 'Starting segmentation of', self.embryoLabels[e]
            # Make embyo label current (for saving)
            self.embryo = self.embryoLabels[e]
            self.loadResults()
            # Interpolate egg measurements to whole dataframe
            self.fillEggMeasurements()
            self.getEmbryoFolders(self.parentPath,self.embryoLabels[e])
            self.runEmbryoSegmentation()
            print 'Finished segmentation of', self.embryoLabels[e]
            
#==============================================================================
#   Apply embryo segmentation to specific embryo
#==============================================================================
    def segmentSpecificEmbryos(self, embryo):
        # Make embyo label current (for saving)
        if isinstance(embryo, str):
            self.embryo = embryo
            self.loadResults()
            self.fillEggMeasurements()
            self.getEmbryoFolders(self.parentPath,self.embryo)
            self.runEmbryoSegmentation()
            print self.embryo, 'segmentation complete'
        if isinstance(embryo, list):
            for e in range(len(embryo)):
                self.embryo = embryo[e]
                self.loadResults()
                self.getEmbryoFolders(self.parentPath,self.embryo)
                self.runEmbryoSegmentation()
                print self.embryo, 'segmentation complete'

#==============================================================================
#   Worker function for parallel embryo segmentation                 
#==============================================================================
    def parallelSegmentAllEmbryos(self,e):     
        print 'Analysis started for:', self.embryoLabels[e]
        # Make embyo label current (for saving)
        ts = time.time()
        self.embryo = self.embryoLabels[e]
        self.loadResults()
        if self.species == 'rbalthica':
            # Interpolate egg measurements to whole dataframe
            self.fillEggMeasurements()
        self.getEmbryoFolders(self.parentPath,self.embryoLabels[e])        
        self.runParEmbryoSegmentation()
        print 'Analysis complete for', self.embryoLabels[e], ('in {} s'.format(time.time()-ts))        

#==============================================================================
#   Multiprocessing parallel embryo segmetnation function            
#==============================================================================
    def quantifyAllEmbryos(self,par=True,exclude='na'):
        # Exclude any embryos that user does not want to be analysed.
        if exclude is not 'na':
            if type(exclude) is str:
                self.embryoLabels = self.embryoLabels[self.embryoLabels != exclude]  
            else:
                for p in range(len(exclude)):
                    self.embryoLabels = self.embryoLabels[self.embryoLabels != str(exclude[p])]                  
        if par is True:
            # Four seems a good compromise for processing vs data IO on a 12 core
            # MacPro, with data on a 7200 RPM SATA, connected via USB3.0/eSATA or Thunderbolt.
            cpuCount = 4
            # Uncomment if you want to maximise cpu usage, note that data IO from
            # drives will likely become limiting and cause serious issues..
            #cpuCount = pathos.multiprocessing.cpu_count()
            self.getEmbryoLabels(self.parentPath)
            jobSize= len(self.embryoLabels)
            jobRange= range(len(self.embryoLabels))
            if cpuCount > len(self.embryoLabels):
                pool = pathos.multiprocessing.ProcessPool(jobSize)
                pool.map(self.parallelSegmentAllEmbryos, jobRange)        
            else:
                pool = pathos.multiprocessing.ProcessPool(cpuCount)
                pool.map(self.parallelSegmentAllEmbryos, jobRange)    
        else:
            # If par is not True then use non parallel version.
            self.segmentAllEmbryos()
#==============================================================================
#   Save metadata
#==============================================================================        
    def saveMetadata(self):
        if self.species == 'rbalthica':
            metadata = {'embryoLabels':self.embryoLabels, 'scale':self.scale, 'eggInt':self.eggInt}
        elif self.species == 'ogammarellus':
            metadata = {'embryoLabels':self.embryoLabels, 'scale':self.scale, 'eggInt':1234}
        np.save(self.resultsDir + "/phenomeMetadata", metadata)

#==============================================================================
#   Load metadata
#==============================================================================        
    def loadMetadata(self):
        self.metadata = np.load(self.resultsDir + "/phenomeMetadata.npy")
        self.scale = self.metadata[()]['scale']
        self.embryoLabels = self.metadata[()]['embryoLabels']
        # Use to remove 'problematic embryos'
        if self.exclude is not 'na':
            if type(self.exclude) is str:
                self.embryoLabels = self.embryoLabels[self.embryoLabels != self.exclude]  
            else:
                self.embryoLabels = self.embryoLabels
                for p in range(len(self.exclude)):
                    self.embryoLabels = self.embryoLabels[self.embryoLabels != str(self.exclude[p])]  
             
#==============================================================================
#   Multiprocessing parallel embryo segmetnation functions            
#==============================================================================
    # Performs analysis for each embryo, looping over embryoFolders.
    def initEmParallel(self,e):     
        print 'Initiation in progress for ', self.embryoLabels[e]
        # Make embyo label current (for saving)
        #ts = time.time()
        self.getEmbryoFolders(self.parentPath,self.embryoLabels[e])
        self.embryo = self.embryoLabels[e]
        for f in range(len(self.embryoFolders)):
            self.currentFolder = self.embryoFolders[f]
            self.shortenedPath = self.shortenedPaths[f]
            #self.shortenedPath = os.path.relpath(self.currentFolder, self.parentPath)
            self.parseMetadata()
            self.createResultsTable()
            self.compiledData[self.shortenedPath] = self.results
        self.resultSheets = pd.Panel.from_dict(self.compiledData)
        self.resultSheets.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')
        self.results=[]
        self.resultSheets=[]
        self.compiledData={}
        # Now the results file is created and saved open it and populate with egg measurements
        self.loadResults()
        self.getIntervIndicesFromSequence(int(self.eggInt))
        for f in range(len(self.embryoFolders)):
            self.nonContSeqImport(f)
            for g in range(len(self.eggIDIms)):
                if len(self.eggIDIms) == 1:
                    g =0
                # NOW perform the eggID on each intervalled image
                self.locateEgg(self.seq[g])
                # AND store eggID output in results table
                #self.results.ix[f,self.eggIDIms[g],'eggBBox'] = self.eggBBox
                self.results.ix[f,self.eggIDIms[g],'eggRotBBox'] = self.eggRotBBox
                self.results.ix[f,self.eggIDIms[g],'eggBoxPoints'] = self.boxPoints     
        self.results.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')                    

    # Creates worker pool and allocates jobs. 
    def parallelInitiation(self):
        self.getEmbryoLabels(self.parentPath)
        #cpuCount = pathos.multiprocessing.cpu_count()
        cpuCount = pathos.multiprocessing.cpu_count()
        jobSize= len(self.embryoLabels)
        jobRange= range(len(self.embryoLabels))
        if cpuCount > len(self.embryoLabels):
            pool = pathos.multiprocessing.ProcessPool(jobSize)
            pool.map(self.initEmParallel, jobRange)        
        else:
            pool = pathos.multiprocessing.ProcessPool(cpuCount)
            pool.map(self.initEmParallel, jobRange)                
            
    # Worker function ensuring appropriate functions are called before and after parallel processing.
    def parallelGenerateResultsAndFindEggs(self, parentPath, scale, eggInt=1234): 
        if self.species == 'rbalthica':
            totts = time.time()
            self.eggInt = eggInt
            self.compiledData = {}
            self.parentPath = parentPath
            self.scale = scale
            self.getEmbryoLabels(self.parentPath)
            # Debug
            print self.embryoLabels
            print len(self.embryoLabels)
            self.createResultsFolder()   
            self.parallelInitiation()
            self.saveMetadata()
            print 'Egg identification and creation of results files is now complete. '
            print 'This took ', ('{} s'.format(time.time()-totts))    
        
        
#==============================================================================
#   Return a pandas panel with blockwise and embryooutline removed - can be 
#   useful for some downstream, restricted analysis.
#==============================================================================
    def reduceData(self):
        res = self.results.drop('blockWise', axis=2)
        res = res.drop('embryoOutline', axis=2)
        return res

#==============================================================================
#     Save a restricted dataset (excluding blockwise data and embryo outlines
#     to a Numpy dictionary.
#==============================================================================
    def saveReducedResults(self, savePath):
    # Use function and add output to results dict.
        reducedResults = {}
        for f in range(len(self.embryoLabels)):
            self.embryo = self.embryoLabels[f]
            self.results = self.returnResults(self.embryo)
            reduced = self.reduceData()
            # Take embryo label and assign results
            #exec('self.%s = reduced' % self.embryo)
            print self.embryoLabels[f], 'loaded'
            reducedResults[self.embryoLabels[f]] = reduced
        # Get an appropriate name to save..
        out = self.parentPath.replace('/phenomeData/','')
        out = out.split('/')
        out = out[len(out)-1]
        # Finally save
        np.save(savePath + out + '_reducedPhenomeResults.npy', reducedResults)


#==============================================================================
#   Functionn to launch Egg ID UI 
#==============================================================================
    def validateEggs(self, eggInt = 1234):
        self.eggInt = eggInt
        app = 0
        #QtGui.QApplication.setGraphicsSystem('raster')
        app = QtGui.QApplication(sys.argv)
        self.UI = eggUI.eggUI()
        self.dataForUI(0)
        self.UI.embryoFolders = self.embryoFolders
        self.UI.showUI(self.UI.compSeq, self.results[:,self.eggIDIms,'eggRotBBox'].values, self.results[:,self.eggIDIms,'eggBoxPoints'].values,list(self.embryoLabels), self.eggInt)
        #instance1.UI.showUI(instance1.UI.compSeq, instance1.results[:,instance1.eggIDIms,'eggRotBBox'].values, instance1.results[:,instance1.eggIDIms,'eggBoxPoints'].values,list(instance1.embryoLabels), instance1.eggInt)
        self.UI.diag.imv.sigTimeChanged.connect(self.UI.updateOpenCVEggROICurrEmbryo)
        self.UI.diag.table.itemSelectionChanged.connect(self.supplyUINewEmbryoData)
        self.UI.approveROI_btn.clicked.connect(self.saveUpdatedROI)        
        # Update image when timeline slider is changed.
        # self.UI.diag.imv.timeLine.sigPositionChanged.connect(self.UI.updateImage)
        app.exec_()      
        
#==============================================================================
#  Functions for the validateEggs() user interface - to validate egg locations etc..
#==============================================================================
#==============================================================================
#   Load data for the eggUI for a particular embryo
#==============================================================================
    def dataForUI(self,e):
        # Make embyo label current (for saving)
        self.embryo = self.embryoLabels[e]
        # Load results
        self.loadResults()
        self.getEmbryoFolders(self.parentPath, self.embryo)
        # If self.eggInt = 1234, only check the first image of each image sequence
        # /time series. This is the default. However, users can check very image if 
        # desired by setting an appropriate eggInt.
        if self.eggInt ==1234:
            self.eggIDIms = []
            self.eggIDIms.append(0)
            self.intN = self.results.shape[1]-1
            self.getIntervIndicesFromSequence(int(self.eggInt))
        else:
            # Get indices for loading, based on interval (intN).
            self.getIntervIndicesFromSequence(int(self.eggInt))
        # Create an empty stack for the egg approval images. Use first im for dimensions
        im = self.imImport(0,0)
        self.UI.compSeq = np.zeros(shape=(int(len(self.eggIDIms)*len(self.embryoFolders)),im.shape[0],im.shape[1]))      
        self.UI.eggUIimPaths = []        
        for e in range(len(self.embryoFolders)):
            self.UI.eggUIimPaths[e*len(self.eggIDIms):(e*len(self.eggIDIms)+len(self.eggIDIms))] = self.results.iloc[e]['currentFolder'][self.eggIDIms] + self.results.iloc[e]['file'][self.eggIDIms]

#==============================================================================
#   Save changes to egg ROI to disk.
#==============================================================================
    def saveUpdatedROI(self):    
        # Save changes for each embryo's ROI
        for r in range(self.UI.eggRotBBox.shape[1]):    
            self.results.iloc[r]['eggRotBBox'][self.eggIDIms] = self.UI.eggRotBBox[:,r]
        for r in range(self.UI.eggBoxPoints.shape[1]):    
            self.results.iloc[r]['eggBoxPoints'][self.eggIDIms] = self.UI.eggBoxPoints[:,r]
        # Store changes to disk
        self.results.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')

#==============================================================================
#   When table row selection changes, load new embryo data
#==============================================================================
    def supplyUINewEmbryoData(self):
        # Debug
        # print 'supplyUINewEmbryoData', self.UI.diag.table.currentRow()
        currRow = self.UI.diag.table.currentRow()
        self.dataForUI(currRow)
        self.UI.updateUI(self.UI.compSeq, self.results[:,self.eggIDIms,'eggRotBBox'].values, self.results[:,self.eggIDIms,'eggBoxPoints'].values)
        #self.dataforViewOutput(currRow)


#==============================================================================
#   View Output functions    
#==============================================================================
    def viewOutput(self):
        app = 0
        app = QApplication(sys.argv)
        self.outputUI = viewOutput.viewOutput()
        self.outputUI.scale = self.scale
        self.dataforViewOutput(0)
        self.outputUI.showUI(list(self.embryoLabels))
        #self.dataForUI(0)
        #self.UI.outputUI(self.compSeq, self.results[:,self.eggIDIms,'eggRotBBox'].values, self.results[:,self.eggIDIms,'eggBoxPoints'].values,list(self.embryoLabels))
        #test.UI.diag.imv.sigTimeChanged.connect(updateAlert)
        #self.UI.diag.imv.sigTimeChanged.connect(self.UI.updateOpenCVEggROICurrEmbryo)
        self.outputUI.diag.table.itemSelectionChanged.connect(self.supplyoutputUINewEmbryoData)
        #self.UI.approveROI_btn.clicked.connect(self.saveUpdatedROI)
        #test.UI.roi.sigRegionChangeFinished.connect(roiChanged)
        app.exec_()     
        
        # Provide results data
    def dataforViewOutput(self,n=0):
        # Make embyo label current (for saving)
        self.embryo = self.embryoLabels[int(n)]
        # Load results
        self.loadResults()
        self.outputUI.results = self.results
        self.getEmbryoFolders(self.parentPath, self.embryo)
        # Send timeSeriesEmbryoBB data ...
        self.outputUI.embryoBBRange  = self.timeSeriesEmbryoBB()
        
    def supplyoutputUINewEmbryoData(self):
        # Debug
        # print 'supplyUINewEmbryoData', self.UI.diag.table.currentRow()
        currRow = self.outputUI.diag.table.currentRow()
        self.dataforViewOutput(currRow)
        self.outputUI.updateUI()

#==============================================================================
# Get Embryo bounding box max and min locations across entire data panel..
#==============================================================================
    def timeSeriesEmbryoBB(self):
        maxX = np.zeros(len(self.results))
        maxY = np.zeros(len(self.results))
        minX = np.zeros(len(self.results))
        minY = np.zeros(len(self.results))        
        for f in range(len(self.results)):
            self.currentFolder = self.embryoFolders[f]
            self.shortenedPath = self.shortenedPaths[f]
            #self.getShortenedPath()
            out = self.getEmbryoBB()
            maxX[f], maxY[f], minX[f], minY[f] = out['maxX'],out['maxY'],out['minX'],out['minY']
        maxX = max(maxX)
        maxY = max(maxY)
        minX = min(minX)
        minY = min(minY)
        
        return {'minX':minX, 'minY':minY, 'maxX':maxX, 'maxY':maxY}


# Get Egg bounding box max and min locations across entire data panel..
    def timeSeriesEggBB(self):
        maxX = np.zeros(len(self.results))
        maxY = np.zeros(len(self.results))
        minX = np.zeros(len(self.results))
        minY = np.zeros(len(self.results))        
        for f in range(len(self.results)):
            self.currentFolder = self.embryoFolders[f]
            self.shortenedPath = self.shortenedPaths[f]            
            out = self.getSeqEggBB()
            maxX[f], maxY[f], minX[f], minY[f] = out['eggMaxX'],out['eggMaxY'],out['eggMinX'],out['eggMinY']
        maxX = max(maxX)
        maxY = max(maxY)
        minX = min(minX)
        minY = min(minY)
        
        return {'minX':minX, 'minY':minY, 'maxX':maxX, 'maxY':maxY}

