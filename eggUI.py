from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
from scipy.spatial import distance as dist
import glob
import re
import os
from PyQt5 import QtGui
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import sys
import cv2
import pandas as pd
from PyQt5.Qt import *
import pyqtgraph as pg
#from PyQt4.Qt import *
#%%
class eggUI(QDialog):
    '''
    createOpenCVEggROI : take eggID defined ROIs and visualise
    
    '''
    sliderUpdate = QtCore.pyqtSignal()
    embryoUpdate = QtCore.pyqtSignal()
    keyPressed = QtCore.pyqtSignal()
    
    def __init__(self, parent=None):
        super(eggUI, self).__init__(parent)
        # Make QDialog
        self.diag = QtGui.QDialog()
        global parentPath, vidTime
        self.diag.setWindowTitle('Identify eggs')
        self.diag.imv = pg.ImageView()
        self.btn_save = QPushButton('Save', self)

#==============================================================================
# 
#==============================================================================
    def showUI(self,ims,eggRotBBox, eggBoxPoints, embryoLabels, eggInt):
        self.eggInt = eggInt
        self.embryoLabels = embryoLabels
        self.diag.setWindowTitle('Identify eggs')
        # Make ImageView
        self.diag.imv = pg.ImageView()
        self.diag.resize(1000,600)
        # Make ROI
        self.importOpenCVROIs(eggRotBBox, eggBoxPoints)
        if (eggRotBBox[0][0][0] != 'nan'): 
            self.createOpenCVEggROI()
            self.diag.imv.addItem(self.roi)
        # Remove buttons from ImageView widget
        self.diag.imv.ui.roiBtn.hide()
        self.diag.imv.ui.menuBtn.hide()
        # Make tableview
        self.diag.table = QtGui.QTableWidget()
        self.diag.table.setShowGrid(True)
        self.diag.table.setHorizontalHeaderLabels(['Embryo', 'Sorted'])
        # Sets different alignment data just on the first column
        self.diag.table.setRowCount(int(len(self.embryoLabels)))
        self.diag.table.setColumnCount(2)
        # Highlight first row
        self.diag.table.selectRow(0)
        # Make layout
        checkLayout = QGridLayout()
        # Deal with stretching for approrpraite formatting.
        checkLayout.setColumnStretch(0, 3)
        checkLayout.setColumnStretch(1, 1)
        checkLayout.setRowStretch(0, 1)
        checkLayout.setRowStretch(1, 3)
        # Add to layout
        checkLayout.addWidget(self.diag.imv,0,0,2,2)
        checkLayout.addWidget(self.diag.table,1,5)
        # Apply layout
        self.diag.setLayout(checkLayout)
        # Make buttons
        self.cpROI_btn = QtGui.QPushButton('&Copy ROI')
        self.cpROI_btn.setMinimumHeight(40);
        self.useCpROI_btn = QtGui.QPushButton('&Use Copied ROI')
        self.useCpROI_btn.setMinimumHeight(40);
        self.noEgg_btn = QtGui.QPushButton('&No Egg')
        self.noEgg_btn.setMinimumHeight(40);
        self.approveROI_btn = QtGui.QPushButton('&Approve ROIs')
        self.approveROI_btn.setMinimumHeight(40);
        self.exit_btn = QtGui.QPushButton('Exit')
        self.exit_btn.setMinimumHeight(40);
        # Make button layout
        self.btnLayout = QGridLayout()
        self.btnLayout.addWidget(self.cpROI_btn,0,0)
        self.btnLayout.addWidget(self.useCpROI_btn,0,1)
        self.btnLayout.addWidget(self.noEgg_btn,1,1)
        self.btnLayout.addWidget(self.approveROI_btn,1,0)
        # Exit button not implemented, just use window x (topRight).
        # self.btnLayout.addWidget(self.exit_btn,2,1)
        # Add button layout to GridLayout.
        checkLayout.addLayout(self.btnLayout,0,5)
        # Format images for pyqtgraph and put in ImageView
        # self.formatSequence(ims)
        self.imImport()
        self.diag.imv.setImage(self.compSeq)
        # Add the ROI to ImageItem        
        self.diag.show()
        # Call function to add data
        self.dataForTable()
        # Function for modifying the table when ROI is approved.
        self.approveROI_btn.clicked.connect(self.updateTable)
        # Copy current ROI
        self.cpROI_btn.clicked.connect(self.cpROI)
        # Apply copied ROI
        self.useCpROI_btn.clicked.connect(self.applyCopiedROI)
        # Assign nan to frames not containing egg
        self.noEgg_btn.clicked.connect(self.recordNoEgg)
        # Exit - prompt user to confirm
        #self.exit_btn.clicked.connect(self.closeEvent)
        # Connect changes in timeline so correct ROI is created and displayed.
        self.diag.imv.timeLine.sigPositionChanged.connect(self.updateOpenCVEggROICurrEmbryo)
        #self.diag.keyPressEvent(self.keyPressEvent)
#==============================================================================
#   Generate data for populating the embryo/approveROI table.
#==============================================================================
    def dataForTable(self):
        self.tableData = {'Embryo':list(self.embryoLabels),
        'ROI approved':['No'] * len(list(self.embryoLabels))}
        self.tableCols = [QtGui.QColor(0,0,100,120)]* len(list(self.embryoLabels))
        # Enter data onto Table
        horHeaders = []
        for n, key in enumerate(sorted(self.tableData.keys())):
            horHeaders.append(key)
            for m, item in enumerate(self.tableData[key]):
                newitem = QtGui.QTableWidgetItem(item)
                newitem.setBackground(QtGui.QColor(0,0,100,120))
                self.diag.table.setItem(m, n, newitem)
        # Add Header
        self.diag.table.setHorizontalHeaderLabels(horHeaders)  
        # Adjust size of Table
        self.diag.table.resizeRowsToContents()
        # self.diag.table.resizeColumnsToContents()

#==============================================================================
#   Update table when approve ROI button clicked.
#==============================================================================
    def updateTable(self):   
        self.tableData['ROI approved'][self.diag.table.currentRow()] = 'Approved'
        self.tableCols[self.diag.table.currentRow()] = QtGui.QColor(0,100,0,120)
        horHeaders = []
        for n, key in enumerate(sorted(self.tableData.keys())):
            horHeaders.append(key)
            for m, item in enumerate(self.tableData[key]):
                newitem = QtGui.QTableWidgetItem(item) 
                self.diag.table.setItem(m, n, newitem)
                newitem.setBackground(self.tableCols[m])
        #Add Header
        self.diag.table.setHorizontalHeaderLabels(horHeaders)  
        #Adjust size of Table
        self.diag.table.resizeRowsToContents()

#==============================================================================
#   Update the user interface        
#==============================================================================
    def updateUI(self,ims,eggRotBBox, eggBoxPoints):
        self.imImport()
        self.diag.imv.setImage(self.compSeq)
        self.importOpenCVROIs(eggRotBBox, eggBoxPoints)
        self.getSeqValsAndCurrROI()  
        self.updateOpenCVEggROINewEmbryo()
        # Add the ROI to ImageItem
        #self.diag.imv.addItem(self.roi)

#==============================================================================
#   Deal with data from the dataHandling class        
#==============================================================================
    def formatSequence(self,ims):
        # Format seq appropriately for pyqtgraph ROIs
        self.tSeqd = np.zeros_like(ims)
        for l in range(len(self.tSeqd)):
            self.tSeqd[l] = ims[l].T
 
#==============================================================================
#   Get folders for a particular embryo
#==============================================================================
    def getEmbryoFolders(self, parentPath, embryo):
        self.parentPath = parentPath
        self.embryo = embryo
        self.embryoFolders = glob.glob(parentPath + "*/" + embryo +"/")
        self.embryoFolders.sort(key=os.path.getctime) 
        
#==============================================================================
#   Get image
#==============================================================================
    def imImport(self):
        for f in range(len(self.eggUIimPaths)):
            im = cv2.imread(self.eggUIimPaths[f],cv2.IMREAD_ANYDEPTH)
            ran = (im.max()-im.min())/255.
            out = (im/ran)
            out = out-out.min()
            self.compSeq[int(f)] = out.astype(np.uint8)
            self.compSeq[f] = self.compSeq[f].T
#==============================================================================
#   Update image iteratively when slider moved
#==============================================================================
#==============================================================================
#     def updateImage(self):
#         self.getSeqValsAndCurrROI()
#         #self.UI.compSeq[e*len(self.eggIDIms):(e*len(self.eggIDIms)+len(self.eggIDIms))] = self.seq
#         #self.UI.comp(self.imImport(self.diag.imv.currentIndex()))
#         im = cv2.imread(self.eggUIimPaths[self.diag.imv.currentIndex],cv2.IMREAD_ANYDEPTH)
#         ran = (im.max()-im.min())/255.
#         out = (im/ran)
#         out = out-out.min()
#         self.compSeq[self.diag.imv.currentIndex] = out.astype(np.uint8)
#         self.diag.imv.setImage(self.compSeq.T)
#         self.diag.imv.show()
# #========
#==============================================================================
#==============================================================================
# ROI functions
#==============================================================================         
#==============================================================================
#   Import OpenCV determined ROIs from dataHandling instance. Called from showUI and updateUI.
#==============================================================================
    def importOpenCVROIs(self,eggRotBBox, eggBoxPoints):
        self.eggRotBBox = eggRotBBox
        self.eggBoxPoints = eggBoxPoints
        self.originalEggRotBBox = eggRotBBox.copy()
        self.originalEggBoxPoints = eggBoxPoints.copy()

#==============================================================================
#   Get index values for ROI data.     
#==============================================================================
    def getSeqValsAndCurrROI(self):
        # Calculate the indices for current frame
        if self.eggInt != 1234:
            self.divVal = self.diag.imv.currentIndex/float(len(self.eggRotBBox[1]))
            self.intDivVal = int(self.divVal)
            self.withinSeqVal = int((self.divVal - self.intDivVal)*len(self.eggRotBBox[self.intDivVal]))
            self.currROI_eggRotBBox = self.eggRotBBox[self.intDivVal,self.withinSeqVal]
            self.currROI_eggBoxPoints = self.eggBoxPoints[self.intDivVal,self.withinSeqVal]
        else:
            self.divVal = self.diag.imv.currentIndex
            self.intDivVal = int(self.divVal)
            self.currROI_eggRotBBox = self.eggRotBBox[0,self.intDivVal] 
            self.currROI_eggBoxPoints = self.eggBoxPoints[0,self.intDivVal] 
#==============================================================================
#   Generate a pyqtgraph ROI, using data from OpenCV.
#==============================================================================
    def createOpenCVEggROI(self):
        # Get relevant sequence position and ROI.
        self.getSeqValsAndCurrROI()   
        if (self.currROI_eggRotBBox[0] != 'nan'):
            # 0 or 90 degree angles seem very buggy. Shift to 1 and 89 as a bodge fix.
            if self.currROI_eggRotBBox[4] == -90:
                #self.currROI_eggRotBBox[4] = -89
                # Get rotated bounding box points
                ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
                # Get bottom most, and top most sorted corner points
                bottomMost = ySorted[:2, :]
                topMost = ySorted[2:, :]
                # Get bottom most
                bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
                (bl, br) = bottomMost
                # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
                # The point with the largest distance will be our bottom-right point
                D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
                (tl, tr) = topMost[np.argsort(D)[::-1], :]
                self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
            elif self.currROI_eggRotBBox[4] == -0:
                #self.currROI_eggRotBBox[4] = -1  
                ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
                # Get bottom most, and top most sorted corner points
                bottomMost = ySorted[:2, :]
                topMost = ySorted[2:, :]
                # Get bottom most
                bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
                (bl, br) = bottomMost
                # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
                # The point with the largest distance will be our bottom-right point
                D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
                (tl, tr) = topMost[np.argsort(D)[::-1], :]
                self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
            elif self.currROI_eggRotBBox[4] == -180:
                #self.currROI_eggRotBBox[4] = -179
                ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
                # Get bottom most, and top most sorted corner points
                bottomMost = ySorted[:2, :]
                topMost = ySorted[2:, :]
                # Get bottom most
                bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
                (bl, br) = bottomMost
                # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
                # The point with the largest distance will be our bottom-right point
                D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
                (tl, tr) = topMost[np.argsort(D)[::-1], :]
                self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
            
            else:
                # Get rotated bounding box points
                ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
                # Get bottom most, and top most sorted corner points
                bottomMost = ySorted[:2, :]
                topMost = ySorted[2:, :]
                # Get bottom most
                bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
                (bl, br) = bottomMost
                # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
                # The point with the largest distance will be our bottom-right point
                D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
                (tl, tr) = topMost[np.argsort(D)[::-1], :]
        
                # Make ROI - note non 0,or 90 degree angles, require different of the X size
                # Rectangular ROI used to enable more easy handling of corner handles for tracking user chagnges.
                if (self.currROI_eggRotBBox[4] == -90.0) | (self.currROI_eggRotBBox[4] == -0.0)| (self.currROI_eggRotBBox[4] == 0.0):
                    self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
                    # roi = pg.EllipseROI([bottomMost[0][0], bottomMost[0][1]], [eggRotBBox[vidTime][2], eggRotBBox[vidTime][3]])
                    # Debug
                    # print 'no angle'    
                else:
                    # Random angle ROIs
                    self.roi = pg.ROI([bottomMost[0][0], bottomMost[0][1]], [-self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
                    self.roi.setAngle(self.currROI_eggRotBBox[4], update=True)
                    # roi = pg.EllipseROI([bottomMost[0][0], bottomMost[0][1]], [-eggRotBBox[vidTime][2], eggRotBBox[vidTime][3]])
            # Add handles
            self.roi.addRotateHandle([1, 0],[0.5,0.5])
            self.roi.addRotateHandle([0, 1], [0.5,0.5])
            self.roi.addScaleHandle([1, 1], [0, 0])
            self.roi.addScaleHandle([0, 0], [1, 1])
            self.roi.setPen('y',width=3)
            self.roi.removable
            self.roi.invertible = 'True'
            # Make var for dealing with modifications to roi
            self.updatedEggROI=[]
            self.roi.sigRegionChangeFinished.connect(self.updateROI)
        #else:
#==============================================================================
#   Update the ROI for current embryo.
#==============================================================================    
    def updateOpenCVEggROICurrEmbryo(self):
        # Remove previous
        if (hasattr(self, 'roi')):
            self.diag.imv.removeItem(self.roi)
        # Get relevant video position and ROI.
        self.getSeqValsAndCurrROI()      
        # 0 or 90 degree angles seem very buggy. Shift to 1 and 89 as a bodge fix.
        if self.currROI_eggRotBBox[4] == -90:
            #self.currROI_eggRotBBox[4] = -89
            # Get rotated bounding box points
            ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
            # Get bottom most, and top most sorted corner points
            bottomMost = ySorted[:2, :]
            topMost = ySorted[2:, :]
            # Get bottom most
            bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
            (bl, br) = bottomMost
            # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
            # The point with the largest distance will be our bottom-right point
            D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
            (tl, tr) = topMost[np.argsort(D)[::-1], :]
            self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
        elif self.currROI_eggRotBBox[4] == -0:
            #self.currROI_eggRotBBox[4] = -1  
            ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
            # Get bottom most, and top most sorted corner points
            bottomMost = ySorted[:2, :]
            topMost = ySorted[2:, :]
            # Get bottom most
            bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
            (bl, br) = bottomMost
            # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
            # The point with the largest distance will be our bottom-right point
            D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
            (tl, tr) = topMost[np.argsort(D)[::-1], :]
            self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
        elif self.currROI_eggRotBBox[4] == -180:
            #self.currROI_eggRotBBox[4] = -179
            ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
            # Get bottom most, and top most sorted corner points
            bottomMost = ySorted[:2, :]
            topMost = ySorted[2:, :]
            # Get bottom most
            bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
            (bl, br) = bottomMost
            # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
            # The point with the largest distance will be our bottom-right point
            D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
            (tl, tr) = topMost[np.argsort(D)[::-1], :]
            self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
        
        else:
            # Get rotated bounding box points
            ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
            # Get bottom most, and top most sorted corner points
            bottomMost = ySorted[:2, :]
            topMost = ySorted[2:, :]
            # Get bottom most
            bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
            (bl, br) = bottomMost
            # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
            # The point with the largest distance will be our bottom-right point
            D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
            (tl, tr) = topMost[np.argsort(D)[::-1], :]
    
            # Make ROI - note non 0,or 90 degree angles, require different of the X size
            # Rectangular ROI used to enable more easy handling of corner handles for tracking user chagnges.
            if (self.currROI_eggRotBBox[4] == -90.0) | (self.currROI_eggRotBBox[4] == -0.0)| (self.currROI_eggRotBBox[4] == 0.0):
                self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
                # roi = pg.EllipseROI([bottomMost[0][0], bottomMost[0][1]], [eggRotBBox[vidTime][2], eggRotBBox[vidTime][3]])
                # Debug
                # print 'no angle'    
            else:
                # Random angle ROIs
                self.roi = pg.ROI([bottomMost[0][0], bottomMost[0][1]], [-self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
                self.roi.setAngle(self.currROI_eggRotBBox[4], update=True)
                # roi = pg.EllipseROI([bottomMost[0][0], bottomMost[0][1]], [-eggRotBBox[vidTime][2], eggRotBBox[vidTime][3]])
 
            # roi = pg.EllipseROI([bottomMost[0][0], bottomMost[0][1]], [-eggRotBBox[vidTime][2], eggRotBBox[vidTime][3]])
        # Add handles
        self.roi.addRotateHandle([1, 0],[0.5,0.5])
        self.roi.addRotateHandle([0, 1], [0.5,0.5])
        self.roi.addScaleHandle([1, 1], [0, 0])
        self.roi.addScaleHandle([0, 0], [1, 1])
        self.roi.setPen('y',width=3)
        self.roi.removable
        self.roi.invertible = 'True'
        # Make var for dealing with modifications to roi
        self.updatedEggROI=[]
        ### Still to do...
        self.diag.imv.addItem(self.roi)
        self.roi.sigRegionChangeFinished.connect(self.updateROI)

#==============================================================================
#   Update ROI for new embryo.        
#==============================================================================
    def updateOpenCVEggROINewEmbryo(self):
        # Remove old ROI
        if (hasattr(self, 'roi')):
            self.diag.imv.removeItem(self.roi)
        # Get relevant video position and ROI
        self.getSeqValsAndCurrROI()    
        # 0 or 90 degree angles seem very buggy. Shift to 1 and 89 as a bodge fix.        
        if self.currROI_eggRotBBox[4] == -90:
            #self.currROI_eggRotBBox[4] = -89
            # Get rotated bounding box points
            ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
            # Get bottom most, and top most sorted corner points
            bottomMost = ySorted[:2, :]
            topMost = ySorted[2:, :]
            # Get bottom most
            bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
            (bl, br) = bottomMost
            # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
            # The point with the largest distance will be our bottom-right point
            D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
            (tl, tr) = topMost[np.argsort(D)[::-1], :]
            self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
        elif self.currROI_eggRotBBox[4] == -0:
            #self.currROI_eggRotBBox[4] = -1  
            ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
            # Get bottom most, and top most sorted corner points
            bottomMost = ySorted[:2, :]
            topMost = ySorted[2:, :]
            # Get bottom most
            bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
            (bl, br) = bottomMost
            # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
            # The point with the largest distance will be our bottom-right point
            D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
            (tl, tr) = topMost[np.argsort(D)[::-1], :]
            self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
        elif self.currROI_eggRotBBox[4] == -180:
            #self.currROI_eggRotBBox[4] = -179
            ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
            # Get bottom most, and top most sorted corner points
            bottomMost = ySorted[:2, :]
            topMost = ySorted[2:, :]
            # Get bottom most
            bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
            (bl, br) = bottomMost
            # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
            # The point with the largest distance will be our bottom-right point
            D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
            (tl, tr) = topMost[np.argsort(D)[::-1], :]
            self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
        
        else:
            # Get rotated bounding box points
            ySorted = self.currROI_eggBoxPoints[np.argsort(self.currROI_eggBoxPoints[:, 1]), :]
            # Get bottom most, and top most sorted corner points
            bottomMost = ySorted[:2, :]
            topMost = ySorted[2:, :]
            # Get bottom most
            bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
            (bl, br) = bottomMost
            # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
            # The point with the largest distance will be our bottom-right point
            D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
            (tl, tr) = topMost[np.argsort(D)[::-1], :]
    
            # Make ROI - note non 0,or 90 degree angles, require different of the X size
            # Rectangular ROI used to enable more easy handling of corner handles for tracking user chagnges.
            if (self.currROI_eggRotBBox[4] == -90.0) | (self.currROI_eggRotBBox[4] == -0.0)| (self.currROI_eggRotBBox[4] == 0.0):
                self.roi = pg.ROI([bl[0], bl[1]], [self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
                # roi = pg.EllipseROI([bottomMost[0][0], bottomMost[0][1]], [eggRotBBox[vidTime][2], eggRotBBox[vidTime][3]])
                # Debug
                # print 'no angle'    
            else:
                # Random angle ROIs
                self.roi = pg.ROI([bottomMost[0][0], bottomMost[0][1]], [-self.currROI_eggRotBBox[2], self.currROI_eggRotBBox[3]])
                self.roi.setAngle(self.currROI_eggRotBBox[4], update=True)
        # Add handles
        self.roi.addRotateHandle([1, 0],[0.5,0.5])
        self.roi.addRotateHandle([0, 1], [0.5,0.5])
        self.roi.addScaleHandle([1, 1], [0, 0])
        self.roi.addScaleHandle([0, 0], [1, 1])
        self.roi.setPen('y',width=3)
        self.roi.removable
        self.roi.invertible = 'True'
        # Make var for dealing with modifications to roi
        self.updatedEggROI=[]
        ### Still to do...
        self.diag.imv.addItem(self.roi)
        self.roi.sigRegionChangeFinished.connect(self.updateROI)

#==============================================================================
#   Update ROI.        
#==============================================================================
    def updateROI(self):
        #global vidTime, xyPosHandles, ellipse, changeAngle, roiChanges,updatedEggROI, changeX, changeY, changeScaleX, changeScaleY, changeAngle
        # Get changes to ROI scale, angle and position
        roiChanges = self.roi.getGlobalTransform()    
        changeX = -roiChanges.getTranslation()[0]
        changeY = roiChanges.getTranslation()[1]
        changeScaleX = roiChanges.getScale()[0]
        changeScaleY = roiChanges.getScale()[1]
        changeAngle = roiChanges.getAngle()
        # Update ROI, either updating the previously updated or taking the unaltered ROI from OpenCV as a starting point.
        #if len(self.updatedEggROI) == 0:
        self.updatedEggROI = (((self.currROI_eggRotBBox[0]-changeX),(self.currROI_eggRotBBox[1]+changeY)),((max((self.currROI_eggRotBBox[3]*changeScaleX),(self.currROI_eggRotBBox[2]*changeScaleY))),(min((self.currROI_eggRotBBox[3]*changeScaleX),(self.currROI_eggRotBBox[2]*changeScaleY)))),self.currROI_eggRotBBox[4]+changeAngle)
        #else:
            #self.updatedEggROI = (((self.updatedEggROI[0][0]-changeX),(self.updatedEggROI[0][1]+changeY)),((max((self.updatedEggROI[1][0]*changeScaleX),(self.updatedEggROI[1][1]*changeScaleY))),(min((self.updatedEggROI[1][0]*changeScaleX),(self.updatedEggROI[1][1]*changeScaleY)))),self.updatedEggROI[2]+changeAngle)
        hh = self.roi.getHandles()
        hh = [self.roi.mapToItem(self.diag.imv.getImageItem(), h.pos()) for h in hh]
        # Handle on each corner. Get handle positions
        self.xyPosHandles =[]  
        for h in hh:
            self.xyPosHandles.append([h.x(),h.y()]) 
        (eggBBX, eggBBY), (eggBBW, eggBBH), eggBBAng = cv2.minAreaRect(np.array(self.xyPosHandles, dtype=np.int32) )
        if eggBBAng == -90:
            eggBBAng = -89
        elif eggBBAng == -180:
            eggBBAng = -179
        elif eggBBAng == -0:
            eggBBAng = -1  
            # Save updated 
        # If more than one frame eggID per sequence..
        if self.eggInt != 1234:
            self.eggRotBBox[self.intDivVal,self.withinSeqVal] = [eggBBX, eggBBY, eggBBW, eggBBH, eggBBAng]
            self.eggBoxPoints[self.intDivVal,self.withinSeqVal] = cv2.boxPoints(((eggBBX, eggBBY), (eggBBW, eggBBH), eggBBAng))
        # Otherwise just save simply
        else:
            self.eggRotBBox[0,self.intDivVal] = [eggBBX, eggBBY, eggBBW, eggBBH, eggBBAng]
            self.eggBoxPoints[0,self.intDivVal] = cv2.boxPoints(((eggBBX, eggBBY), (eggBBW, eggBBH), eggBBAng))

#==============================================================================
#   Copy ROI on button click.
#==============================================================================
    def cpROI(self):
        self.originalEggRotBBox = self.currROI_eggRotBBox
        self.originalEggBoxPoints = self.currROI_eggBoxPoints

#==============================================================================
#   Assign nan to current ROI if 'No Egg' button clicked
#==============================================================================
    def recordNoEgg(self):
        # Remove ROI
        self.diag.imv.removeItem(self.roi)
        # Store nans in place of ROI
        if self.eggInt != 1234:
            self.eggRotBBox[self.intDivVal,self.withinSeqVal] = [np.nan, np.nan, np.nan, np.nan, np.nan]
            self.eggBoxPoints[0,self.intDivVal] = [np.nan,np.nan,np.nan,np.nan]
        else:
            self.eggBoxPoints[0,self.intDivVal] = [np.nan,np.nan,np.nan,np.nan]
            self.eggRotBBox[0,self.intDivVal] = [np.nan, np.nan, np.nan, np.nan, np.nan]
        
#==============================================================================
#   Copy ROI on button click.
#==============================================================================
    def applyCopiedROI(self):
        self.getSeqValsAndCurrROI()
        # Store copied ROI to embryo sequence ROIs
        if self.eggInt != 1234:
            self.divVal = self.diag.imv.currentIndex/float(len(self.eggRotBBox[1]))
            self.intDivVal = int(self.divVal)
            self.withinSeqVal = int((self.divVal - self.intDivVal)*len(self.eggRotBBox[self.intDivVal]))
            self.eggRotBBox[self.intDivVal,self.withinSeqVal] = self.originalEggRotBBox
            self.eggBoxPoints[self.intDivVal,self.withinSeqVal] = self.originalEggBoxPoints
        else:
            self.divVal = self.diag.imv.currentIndex
            self.intDivVal = int(self.divVal)
            self.eggRotBBox[0,self.intDivVal] = self.originalEggRotBBox
            self.eggBoxPoints[0,self.intDivVal] = self.originalEggBoxPoints 
        self.updateOpenCVEggROICurrEmbryo()
        
        
 
#==============================================================================
#        
#==============================================================================

#==============================================================================
#   Close button - not implemented (hidden)
#==============================================================================
#==============================================================================
#     def closeEvent(self, event):
#     
#         quit_msg = "Are you sure you want to exit the program?"
#         reply = QtGui.QMessageBox.question(self, 'Message', 
#                          quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
#     
#         if reply == QtGui.QMessageBox.Yes:
#             #event.accept()
#             app.quit()
#         else:
#             event.ignore()
#         
#==============================================================================
#==============================================================================
#         #self.originalEggRotBBox = eggRotBBox.copy()
#         #self.originalEggBoxPoints = eggBoxPoints.copy()
#         #self.currROI_eggRotBBox = self.eggRotBBox[self.intDivVal,self.withinSeqVal]
#         #self.currROI_eggBoxPoints = self.eggBoxPoints[self.intDivVal,self.withinSeqVal]
# 
#         # Modified version of updateOpenCVEggROICurrEmbryo
#         # Remove previous
#         self.diag.imv.removeItem(self.roi)
#         # Get relevant video position and ROI.
#         self.getSeqValsAndCurrROI()      
#         # Get rotated bounding box points
#         ySorted = self.originalEggBoxPoints[np.argsort(self.originalEggBoxPoints[:, 1]), :]
#         # Get bottom most, and top most sorted corner points
#         bottomMost = ySorted[:2, :]
#         topMost = ySorted[2:, :]
#         # Get bottom most
#         bottomMost = bottomMost[np.argsort(bottomMost[:, 1]), :]
#         (bl, br) = bottomMost
#         # Use bottom-left coordinate as anchor to calculate the Euclidean distance between the
#         # The point with the largest distance will be our bottom-right point
#         D = dist.cdist(bl[np.newaxis], topMost, "euclidean")[0]
#         (tl, tr) = topMost[np.argsort(D)[::-1], :]
#         # Make ROI - note non 0,or 90 degree angles, require different of the X size
#         # Rectangular ROI used to enable more easy handling of corner handles for tracking user chagnges.
#         if (self.originalEggRotBBox[4] == -90.0) | (self.originalEggRotBBox[4] == -0.0)| (self.originalEggRotBBox[4] == 0.0):
#             self.roi = pg.ROI([bottomMost[0][0], bottomMost[0][1]], [self.originalEggRotBBox[2], self.originalEggRotBBox[3]])
#             # roi = pg.EllipseROI([bottomMost[0][0], bottomMost[0][1]], [eggRotBBox[vidTime][2], eggRotBBox[vidTime][3]])
#         else:
#             # Random angle ROIs
#             self.roi = pg.ROI([bottomMost[0][0], bottomMost[0][1]], [-self.originalEggRotBBox[2], self.originalEggRotBBox[3]])
#             self.roi.setAngle(self.originalEggRotBBox[4], update=True)
#             # roi = pg.EllipseROI([bottomMost[0][0], bottomMost[0][1]], [-eggRotBBox[vidTime][2], eggRotBBox[vidTime][3]])
#         # Add handles
#         self.roi.addRotateHandle([1, 0],[0.5,0.5])
#         self.roi.addRotateHandle([0, 1], [0.5,0.5])
#         self.roi.addScaleHandle([1, 1], [0, 0])
#         self.roi.addScaleHandle([0, 0], [1, 1])
#         self.roi.setPen('y',width=3)
#         self.roi.removable
#         self.roi.invertible = 'True'
#         # Make var for dealing with modifications to roi
#         self.updatedEggROI=[]
#         ### Still to do...
#         self.diag.imv.addItem(self.roi)
#         self.roi.sigRegionChangeFinished.connect(self.updateROI)
#==============================================================================
        
#===============