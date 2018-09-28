# Import dependencies
import pandas as pd
import os
from imageAnalysis import imageAnalysis
from dataHandling import dataHandling
from dataAnalysis import dataAnalysis
from xarrayDataAnalysis import xarrayDataAnalysis
from dataIntegration import dataIntegration
import warnings
#from collapsedDataAnalysis import collapsedDataAnalysis

class EmbryoCV(imageAnalysis, dataHandling, dataAnalysis, xarrayDataAnalysis, dataIntegration):
    '''''
    - embryoCV class(parentPath, new = 'no')
        > Initialise with two argments parent path i.e. the folder containing the sequentially
        labelled MicroManager folders for each time point
        
                     /Volumes/DataDrive/Time_1/Position1
                    ................../Time_1/Position2
                    ................../Time_2/Position1
                    ................../Time_2/Position2
        
        parentPath = '/Volumes/DataDrive/'
        > Second argument specifies mode:
            - 'new' = a new analysis - build results files and locate eggs
            - 'resume' = load a previous analysis.
            - 'results' = load a previous analysis, but just from the results folder, 
            (i.e. the results folder is seperate from the original images)
           - 'xarray' = load an XArray version of the results. This is encouraged over
           summaryresults below as it is fast (data is accessed dynamically from disk)
           and it contains all of the results data in a format that can be easily added to
           with subsequent analysis.
           - 'summaryresults' = load a reduced version of the results (no embryo
            outlines or blockwise/freq data) from a numpy dict.
 
    - filterResultsList(fileSize)
        A function to generate a size filtered list of results files from the 
        results directory. Can be used if some of the results files are incomplete,
        owing to either a failed or prematurely aborted analysis. Provide the 
        file size to filter by in bytes.
        
    - seqImport(seqPath)
        Import a continuous sequence of images from specified dirctory. Images
        are converted from 16 bits to 8 bits, scaled and assigned to 'seq', a numpy 
        array.
        
    - imImport(n,m)
        Import a single image using 'currentFolder'and 'filename' taken from the 
        currently loaded results . n specifies the results panel and m the sequence
        number.
        
    - getEmbryoLabels(parentPath)
        Generate a list of embryo labels from the specified parent path. The embryo 
        list is assigned to embryoLabels.
        
    - getEmbryoLabelsFromResultsFolder(resultsDir)
        Use the results files in the results directory to generate a list of embryo 
        labels. These are assigned to embryoLabels.
    
    - getEmbryoFolders(parentPath, embryo)
        Generate a list of folders containing image sequences for a particular embryo.
        The embryo folder list is assigned to embryoFolders and is ordered by
        date created.
        
    - createResultsFolder()
        Generate a 'phenomeData' results folder in the parentPath.
    
    - generateResultsAndFindEggs(parentPath, scale, eggInt=1234):
        Generate results files for each embryo
    '''''    
#==============================================================================
#==============================================================================    
    def __init__(self,parentPath,mode='unspecified', scale='na',exclude = 'na',species='rbalthica', dataformat=True): 
        if mode == 'resume':
            self.mode = mode
            self.dataformat = dataformat
            self.parentPath = parentPath
            self.exclude = exclude
            if (species == 'rbalthica') or (species =='ogammarellus'):
                self.species = species
            else:
                print '** ' + str(species) + ' not currently supported' + ' **'
            self.getEmbryoLabels(self.parentPath)
            print len(self.embryoLabels), 'embryos identified'
            self.resultsDir = os.path.dirname(self.parentPath + "phenomeData/")
            self.loadMetadata()
            self.embryo = self.embryoLabels[0]
            #self.loadResults()
            #self.getEmbryoFolders(self.parentPath, self.embryo)
            print 'EmbryoCV launched in resume mode. Resume analysis at appropriate stage'
        elif mode == 'new':
            self.mode = mode
            self.dataformat = dataformat
            self.parentPath = parentPath
            self.resultsDir = os.path.dirname(self.parentPath + "phenomeData/")
            self.exclude = exclude
            if (species == 'rbalthica') or (species =='ogammarellus'):
                self.species = species
            else:
                print '** ' + str(species) + ' not currently supported' + ' **'
            self.getEmbryoLabels(self.parentPath)
            self.embryo = self.embryoLabels[0]
            if scale != 'na':
                self.scale = scale            
                print 'New instance of EmbryoCV launched.' 
                print len(self.embryoLabels), 'embryos identified'
                # If R.balthica - generate results table and locate eggs..
                if self.species == 'rbalthica':
                    print 'Identifying eggs and generating results tables'
                    # Non-parallel version - useful for debugging..
                    # self.generateResultsAndFindEggs(self.parentPath,self.scale)
                    self.parallelGenerateResultsAndFindEggs(self.parentPath,self.scale)
                # If O. gammarellus generate results table, but dont locate eggs
                # as egg and embryo are indistinguishable from each other, so diff. 
                # approach required.
                if self.species == 'ogammarellus':
                    self.generateResults(self.parentPath,self.scale)
            else:
                print 'Warning: scale not given!'
        elif mode == 'results':
            self.mode = mode
            self.dataformat = dataformat
            self.parentPath = parentPath
            self.resultsDir = self.parentPath
            self.exclude = exclude
            if (species == 'rbalthica') or (species =='ogammarellus'):
                self.species = species
            else:
                print '** ' + str(species) + ' not currently supported' + ' **'
            self.getEmbryoLabelsFromResultsFolder(self.parentPath)
            self.loadMetadata()
            if (species == 'rbalthica') or (species =='ogammarellus'):
                self.species = species
            else:
                print '** ' + str(species) + ' not currently supported' + ' **'
            print len(self.embryoLabels), 'embryos identified'
            self.resultsDir = os.path.dirname(self.parentPath)
            
        # Depreciated..
        #elif mode == 'summaryresults':
        #    self.mode = mode
        #    self.parentPath = parentPath
        #    self.collapsedAnalysis = collapsedDataAnalysis()
        #    self.collapsedAnalysis.parentPath = self.parentPath

        elif mode == 'xarray':
            self.mode = mode
            self.exclude = exclude
            self.parentPath = parentPath
            self.resultsDir = self.parentPath
            self.dataformat = dataformat
            if (species == 'rbalthica') or (species =='ogammarellus'):
                self.species = species
            else:
                print '** ' + str(species) + ' not currently supported' + ' **'
            self.getXREmbryoLabelsFromResultsFolder(self.parentPath)
            self.loadMetadata()
            print len(self.embryoLabels), 'embryos identified'
            self.xarray = xarrayDataAnalysis()
            #self.xarray.
            parentPath = self.parentPath
        else:
            print 'Invalid mode argument (resume, new or results).'
        
        # Surpress warning regarding ctypes.
        warnings.filterwarnings('ignore',category=pd.io.pytables.PerformanceWarning)


