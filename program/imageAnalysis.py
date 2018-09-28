# Import dependencies
import cv2
import numpy as np
from skimage.segmentation import clear_border
import time
import skimage    
import math
import pathos

class imageAnalysis(object):

#==============================================================================
#   Identify the egg 
#==============================================================================
    def locateEgg(self,im):
        #self.eggBBox = []
        self.eggRotBBox = []
        self.boxPoints = []
        origIm = im.copy()
        #kernel = np.ones((4,4),np.uint8)   
        # Apply otsu thresholding to ID egg
        ret2,thresh = cv2.threshold(im,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
        # Open to remove all small/thin elements, including egg (usually leaving blobs and embryo.)
        s1 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (12,12))
        thresh2 = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, s1)
        # Debug
        # plt.imshow(thresh2)  
        # Now use inverted threshold of the opened image to remove embryo and other large items from original OTSU threshold.
        origIm[thresh == thresh2] =0
        # Apply adaptive threshold
        origIm = cv2.adaptiveThreshold(origIm, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 101,0) 
        # Use stats to remove non-egg candidates (rather than erode, which tends to remove too much egg)
        num, labels, stats, centroids = cv2.connectedComponentsWithStats(origIm)
        # 0 = n of blobs,  1 = labelled matrix, 2 = stats (leftMostX, topMostY, Width, Height, Area), 3 = centroids (X and Y
        # Calculate distance from centre of image        
        disFrCen = np.sqrt(((centroids[:,0] - (im.shape[0]/2))**2) + ((centroids[:,1] - (im.shape[1]/2))**2))
        # Get longest axis
        lon = np.amax(stats[:,2:3], axis=1)
        # Transpose stats array and add the extras
        stats = stats.transpose()
        stats = np.vstack([stats, centroids[:,0],centroids[:,1], disFrCen, lon, range(num)])
        # Now filter the blobs
        filtStats = stats[9,(stats[8,:]>20) & (stats[7,:]<(np.amax(im.shape[1:2])/0.8))& (stats[8,:]<(np.amax(origIm.shape[1:2])))]
        if filtStats.shape[0] ==0:
            #self.eggBBox=['nan','nan','nan','nan']
            self.eggRotBBox = ['nan','nan','nan','nan','nan']
            self.boxPoints = (('nan','nan'),('nan','nan'),('nan','nan'),('nan','nan'))
        else:
            #ind = int(ind)  
            mask = np.zeros_like(labels)
            for m in range(len(filtStats)):    
                mask[labels == int(filtStats[m])]=1
            # Debug
            # plt.imshow(mask)
            # Clear border
            origIm = clear_border(mask, buffer_size=20, bgval=0)
            # plt.imshow(testim)
            # Find contours
            contIm, contour,hier = cv2.findContours(origIm,cv2.RETR_CCOMP,cv2.CHAIN_APPROX_SIMPLE)
            
            # If no contours found then do not continue and store nans.
            if len(contour) ==0:
                #self.eggBBox=['nan','nan','nan','nan']
                self.eggRotBBox = ['nan','nan','nan','nan','nan']
                self.boxPoints = (('nan','nan'),('nan','nan'),('nan','nan'),('nan','nan'))
            else:
                # Collate contours
                self.eggCont = cv2.fitEllipse(np.concatenate(contour))
                # Draw ellipse on blank image
                eggEllipse = cv2.ellipse(np.zeros_like(im),self.eggCont,(1,1,1),-1)
                # Debug
                # np.zeros(im.shape,np.bool)
                # plt.imshow(cv2.ellipse(origIm,eggCont,(255,0,0),2))  
                # imcopy[j] = cv2.bitwise_and(thresh,imcopy[j])
                # pg.image(imcopy)
                # plt.imshow(cv2.drawContours(imcopy[j],np.concatenate(contour),-1,(255,0,0),20))
                # Use ellipse to get bounding box (both non rotated for im cropping and rotated for egg gui tweaking).
                _,tmp,_ = cv2.findContours(eggEllipse,cv2.RETR_CCOMP,cv2.CHAIN_APPROX_SIMPLE)
                # If nothing returned, then eggID unsuccesful so return NAs
                if len(tmp)!=1:
                    #self.eggBBox=['nan','nan','nan','nan']
                    self.eggRotBBox = ['nan','nan','nan','nan','nan']
                    self.boxPoints = (('nan','nan'),('nan','nan'),('nan','nan'),('nan','nan'))
                else:
                    # NonRotated BBox dimensions
                    #eggBBx,eggBBy,eggBBw,eggBBh = cv2.boundingRect(tmp[0])
                    #self.eggBBox = [eggBBx,eggBBy,eggBBw,eggBBh]
                    # Rotated BBox dimensions
                    (eggBBRotCentX, eggBBRotCentY), (eggBBRotWid, eggBBRotHei), eggBBRotAng = cv2.minAreaRect(tmp[0])
                    self.eggRotBBox = [eggBBRotCentX, eggBBRotCentY,eggBBRotWid, eggBBRotHei,eggBBRotAng]
                
                    # Get corners for rotated BB
                    self.boxPoints = cv2.boxPoints(cv2.minAreaRect(tmp[0]))                 
                #return {'successful':successful, 'eggCont':eggCont, 'eggBBox':eggBBox, 'eggRotBBox':eggRotBBox, 'eggCentXY':eggCentXY, 'boxPoints':boxPoints}

#===============================================================================
# Segment and measure embryo spatial characteristics
#==============================================================================         
    def segmentEmbryo(self,n):
        if self.species == 'rbalthica':
            # Oli - Note - 7 kernel size upto 3/11. Increased to 9 to try and overcome lack of
            # early stage segmentation.
            blobs = self.auto_thresh(cv2.medianBlur(self.eggSeq[n],9),0.01)
            cntMask = np.zeros((self.eggSeq[n].shape[0], self.eggSeq[n].shape[1]),dtype = np.uint8)
            im2, conts, hierarchy = cv2.findContours(blobs,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
            for cnt in conts:
                if cv2.arcLength(cnt, True) > 200:
                    cv2.drawContours(cntMask,[cnt],-1,(255,255,255),-1)
                    
            # Open to remove all small/thin elements, including egg (usually leaving blobs and embryo.)
            s1 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (12,12))
            cntMask = cv2.morphologyEx(cntMask, cv2.MORPH_OPEN, s1)
            # Get properties of thresholded blobs
            # 0 = n of blobs,  1 = labelled matrix, 2 = stats (leftMostX, topMostY, Width, Height, Area), 3 = centroids (X and Y)
            num, labels, stats, centroids = cv2.connectedComponentsWithStats(cntMask)
            # Calculate extent - region area/bounding box area
            extent = stats[:,4]/(stats[:,2].astype(float)*stats[:,3])
            # Calculate aspect ratio (L/W)
            asp = stats[:,2]/stats[:,3].astype(float)
            for aspRat in range(len(asp)):
                if asp[aspRat] < 1:
                    continue
                else:
                    asp[aspRat] = 1/asp[aspRat]
            # Format and combine extra descriptors
            filtStats = stats.transpose()
            filtStats = np.vstack([filtStats,extent,asp,range(len(stats))])
            # Filter blobs
            # Oli (16/11) changed min area from 5000 to 4000 to 3500 to 2000.
            ind = filtStats[7,(filtStats[4,]>2000) & (filtStats[5,]>0.5) & (filtStats[3,]<labels.shape[0])]
            if ind.shape[0] ==1:
                ind = int(ind)   
                # Debug - visualise output     
                # self.seq[n] = cv2.bitwise_and(self.seq[n],self.seq[n],mask= np.array(labels == ind,dtype='uint8'))
                im2, conts, hierarchy = cv2.findContours(np.array(labels == ind,dtype='uint8'),cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
                for cnt in conts:
                    # cv2.arcLength(cnt, True) > 200:
                    if len(conts) ==1:
                        cv2.drawContours(cntMask,[cnt],-1,(255,255,255),-1)
                        embryoOutline = cnt
                        # Moments
                        M = cv2.moments(embryoOutline)
                        # Centroid
                        cx = int(M['m10']/M['m00'])
                        cy = int(M['m01']/M['m00'])            
                        # Area
                        area = M['m00']
                        # Aspect ratio
                        x,y,w,h = cv2.boundingRect(embryoOutline)
                        aspect = float(w)/h
                        # Extent
                        rect_area = w*h
                        extent = float(area)/rect_area
                        # Solidity
                        hull = cv2.convexHull(embryoOutline)
                        hullArea = cv2.contourArea(hull)
                        solidity = float(area)/hullArea
                    else:
                        embryoOutline = np.NaN
                        cx = np.NaN
                        cy = np.NaN
                        area = np.NaN
                        x = np.NaN
                        y = np.NaN
                        w = np.NaN
                        h = np.NaN
                        solidity = np.NaN
                        hullArea = np.NaN
                        extent = np.NaN
                        aspect = np.NaN
            else:
                embryoOutline = np.NaN
                cx = np.NaN
                cy = np.NaN
                area = np.NaN
                x = np.NaN
                y = np.NaN
                w = np.NaN
                h = np.NaN
                solidity = np.NaN
                hullArea = np.NaN
                extent = np.NaN
                aspect = np.NaN
                #'Embryo':embryo,
        elif self.species == 'ogammarellus':
            # If not R.balthica, perform this segmentation approach... 
            # NOTE: This is a more crude approach - using the
            # bounding box of the egg, rather than a complete embryo segmentation
            # which for species such as O. gammarellus that completely fill the egg
            # capsule is not biologically relevant.)    
            s1 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (19,19))
            mask = cv2.morphologyEx(self.auto_thresh(self.seq[n]), cv2.MORPH_CLOSE, s1)
            s1 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3,3))
            mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, s1)
            cntMask = np.zeros((self.seq[n].shape[0], self.seq[n].shape[1]),dtype = np.uint8)
            im2, conts, hierarchy = cv2.findContours(mask,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
            for cnt in conts:
                if cv2.arcLength(cnt, True) > 200:
                    cv2.drawContours(cntMask,[cnt],-1,(255,255,255),-1)
                    # Get properties of thresholded blobs
                    # 0 = n of blobs,  1 = labelled matrix, 2 = stats (leftMostX, topMostY, Width, Height, Area), 3 = centroids (X and Y)
                    num, labels, stats, centroids = cv2.connectedComponentsWithStats(cntMask)
                    #plt.imshow(labels)
                    # Filter to big blob (gammarid)
                    # Filter blobs
                    #filtStats = stats[np.invert(stats[:,2] == self.seq[n].shape[1]),:]
                    if (stats.shape[0] is 2):
                        # Get index for filtering blob. Discounting the background i.e. size of im.
                        # ind = np.where(np.invert(stats[:,2] == self.seq[n].shape[1]))[0][0]
                        # maskArr = np.array(labels==ind)
                        # Uncomment to apply mask to original image
                        #im = cv2.bitwise_and(im,im, mask = maskArr.astype(np.uint8))
                        embryoOutline = cnt
                        # Get stats
                        # Moments
                        M = cv2.moments(embryoOutline)
                        # Centroid
                        cx = int(M['m10']/M['m00'])
                        cy = int(M['m01']/M['m00'])            
                        # Area
                        area = M['m00']
                        # Aspect ratio
                        x,y,w,h = cv2.boundingRect(embryoOutline)
                        aspect = float(w)/h
                        # Extent - region area/bounding box area
                        extent = area/(w*h)
                        # Solidity
                        hull = cv2.convexHull(cnt)
                        hullArea = cv2.contourArea(hull)
                        solidity = float(area)/hullArea
                    else:
                        embryoOutline = np.NaN
                        cx = np.NaN
                        cy = np.NaN
                        area = np.NaN
                        x = np.NaN
                        y = np.NaN
                        w = np.NaN
                        h = np.NaN
                        solidity = np.NaN
                        hullArea = np.NaN
                        extent = np.NaN
                        aspect = np.NaN
        return {'n':n,'outline':embryoOutline,'centX': cx,'centY': cy, 'ar': area,'aspect': aspect, 'extent': extent,'hullArea': hullArea,'solidity': solidity,'bboxMinCol': x,'bboxMinRow': y,'bboxWidth': w,'bboxHeight': h}    

#==============================================================================#==============================================================================
#   Run embryo segmentation
#==============================================================================
    def runEmbryoSegmentation(self):            
        # Loop over the folders for an embryo and apply segmentation.
        for f in range(len(self.embryoFolders)):
            # Debug
            #print f
            totts = time.time()
            self.currentFolder = self.embryoFolders[f]
            self.shortenedPath = self.shortenedPaths[f]
            print 'Starting: ', self.embryoFolders[f]
            # If no egg dims in first row, do not continue with attempting segmentation.       
            if (self.results.loc[self.shortenedPath]['eggRotBBox'][0] == ['nan','nan','nan','nan','nan'])|(math.isnan(float(self.results.loc[self.shortenedPath]['eggRotBBox'][0][0])))| (self.results.loc[self.shortenedPath]['eggRotBBox'][0] == 'nan')| (self.results.loc[self.shortenedPath]['eggRotBBox'][0] == np.NaN):
                self.results.loc[self.shortenedPath]['embryoOutline'] = np.NaN
                self.results.loc[self.shortenedPath]['centroidX'] = np.NaN
                self.results.loc[self.shortenedPath]['centroidY'] = np.NaN
                self.results.loc[self.shortenedPath]['area'] = np.NaN
                self.results.loc[self.shortenedPath]['aspect'] = np.NaN
                self.results.loc[self.shortenedPath]['extent'] = np.NaN
                self.results.loc[self.shortenedPath]['hullArea'] = np.NaN
                self.results.loc[self.shortenedPath]['solidity'] = np.NaN
                self.results.loc[self.shortenedPath]['bboxMincol'] = np.NaN
                self.results.loc[self.shortenedPath]['bboxMinrow'] = np.NaN
                self.results.loc[self.shortenedPath]['bboxWidth'] = np.NaN
                self.results.loc[self.shortenedPath]['bboxHeight'] = np.NaN
                self.results.loc[self.shortenedPath]['blockWise'] = np.NaN
                print 'No egg for', self.results.loc[self.shortenedPath]['currentFolder'][0]
            else:
                self.seqImport(f)
                # Mask seq using egg masks
                #self.getSeqEggBB()
                self.applyEggMaskAndCrop()
                # Create a threading pool
                innerPool = pathos.multiprocessing.ThreadingPool()
                # Run segmentation
                res = innerPool.map(self.segmentEmbryo, range(self.eggSeq.shape[0]))
                # Save to dataframe
                for g in range(len(res)):
                    self.results.loc[self.shortenedPath]['embryoOutline'][g] = res[g]['outline']
                    self.results.loc[self.shortenedPath]['centroidX'][g] = res[g]['centX']
                    self.results.loc[self.shortenedPath]['centroidY'][g] = res[g]['centY']
                    self.results.loc[self.shortenedPath]['area'][g] = res[g]['ar']
                    self.results.loc[self.shortenedPath]['aspect'][g] = res[g]['aspect']
                    self.results.loc[self.shortenedPath]['extent'][g] = res[g]['extent']
                    self.results.loc[self.shortenedPath]['hullArea'][g] = res[g]['hullArea']
                    self.results.loc[self.shortenedPath]['solidity'][g] = res[g]['solidity']
                    self.results.loc[self.shortenedPath]['bboxMincol'][g] = res[g]['bboxMinCol']
                    self.results.loc[self.shortenedPath]['bboxMinrow'][g] = res[g]['bboxMinRow']
                    self.results.loc[self.shortenedPath]['bboxWidth'][g] = res[g]['bboxWidth']
                    self.results.loc[self.shortenedPath]['bboxHeight'][g] = res[g]['bboxHeight']                      
                # And create a cropped seq
                #self.emSeq = self.seq[:,self.minY:self.maxY, self.minX:self.maxX]
                res = innerPool.map(self.runNestedMeanWindowCalc, range(self.seq.shape[0]))
                for g in range(len(res)):
                    self.results.loc[self.shortenedPath]['blockWise'][g] = res[g]
                print self.embryo, ':', self.shortenedPath, 'analysed in {} s'.format(time.time()-totts) 
        # Save
        self.results.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')

#==============================================================================
#   Segment/analyse particular time points and add these to the results.
#==============================================================================
    def segmentSpecificTimePoints(self):
        print 'Loading results for: ' + self.embryo
        # Load embryo
        self.loadResults()
        self.getEmbryoFolders(self.parentPath, self.embryo)
        # Identify missing data
        indicies = np.isnan(np.nanmean(self.results.ix[:,:,'area'].values.astype(np.float),axis=0))
        inds = np.nonzero(indicies)[0]
        # Loop over missing time points and perform analysis.
        print str(len(inds)) + ' time points with missing data identified'
        for t in range(len(inds)):
            f = inds[t]
            totts = time.time()
            self.currentFolder = self.embryoFolders[f]
            self.shortenedPath = self.shortenedPaths[f]
            #self.getShortenedPath()
            #self.shortenedPath = os.path.relpath(self.currentFolder, self.parentPath)
            print 'Starting: ', self.embryoFolders[f]
            # If no egg dims in first row, do not continue with attempting segmentation.       
            if (self.results.loc[self.shortenedPath]['eggRotBBox'][0] == ['nan','nan','nan','nan','nan'])|(math.isnan(float(self.results.loc[self.shortenedPath]['eggRotBBox'][0][0])))| (self.results.loc[self.shortenedPath]['eggRotBBox'][0] == 'nan')| (self.results.loc[self.shortenedPath]['eggRotBBox'][0] == np.NaN):
                self.results.loc[self.shortenedPath]['embryoOutline'] = np.NaN
                self.results.loc[self.shortenedPath]['centroidX'] = np.NaN
                self.results.loc[self.shortenedPath]['centroidY'] = np.NaN
                self.results.loc[self.shortenedPath]['area'] = np.NaN
                self.results.loc[self.shortenedPath]['aspect'] = np.NaN
                self.results.loc[self.shortenedPath]['extent'] = np.NaN
                self.results.loc[self.shortenedPath]['hullArea'] = np.NaN
                self.results.loc[self.shortenedPath]['solidity'] = np.NaN
                self.results.loc[self.shortenedPath]['bboxMincol'] = np.NaN
                self.results.loc[self.shortenedPath]['bboxMinrow'] = np.NaN
                self.results.loc[self.shortenedPath]['bboxWidth'] = np.NaN
                self.results.loc[self.shortenedPath]['bboxHeight'] = np.NaN
                self.results.loc[self.shortenedPath]['blockWise'] = np.NaN
                print 'No egg for', self.results.loc[self.shortenedPath]['currentFolder'][0]
            else:
                self.seqImport(f)
                # Mask seq using egg masks
                #self.getSeqEggBB()
                self.applyEggMaskAndCrop()
                # Create a threading pool
                innerPool = pathos.multiprocessing.ThreadingPool()
                # Run segmentation
                res = innerPool.map(self.segmentEmbryo, range(self.eggSeq.shape[0]))
                # Save to dataframe
                for g in range(len(res)):
                    self.results.loc[self.shortenedPath]['embryoOutline'][g] = res[g]['outline']
                    self.results.loc[self.shortenedPath]['centroidX'][g] = res[g]['centX']
                    self.results.loc[self.shortenedPath]['centroidY'][g] = res[g]['centY']
                    self.results.loc[self.shortenedPath]['area'][g] = res[g]['ar']
                    self.results.loc[self.shortenedPath]['aspect'][g] = res[g]['aspect']
                    self.results.loc[self.shortenedPath]['extent'][g] = res[g]['extent']
                    self.results.loc[self.shortenedPath]['hullArea'][g] = res[g]['hullArea']
                    self.results.loc[self.shortenedPath]['solidity'][g] = res[g]['solidity']
                    self.results.loc[self.shortenedPath]['bboxMincol'][g] = res[g]['bboxMinCol']
                    self.results.loc[self.shortenedPath]['bboxMinrow'][g] = res[g]['bboxMinRow']
                    self.results.loc[self.shortenedPath]['bboxWidth'][g] = res[g]['bboxWidth']
                    self.results.loc[self.shortenedPath]['bboxHeight'][g] = res[g]['bboxHeight']                      
                # And create a cropped seq
                #self.emSeq = self.seq[:,self.minY:self.maxY, self.minX:self.maxX]
                res = innerPool.map(self.runNestedMeanWindowCalc, range(self.seq.shape[0]))
                for g in range(len(res)):
                    self.results.loc[self.shortenedPath]['blockWise'][g] = res[g]
                print self.embryo, ':', self.shortenedPath, 'analysed in {} s'.format(time.time()-totts) 
                
        tmp = len(np.nonzero(np.isnan(np.nanmean(self.results.ix[:,:,'area'].values.astype(np.float),axis=0)))[0])
        print str(tmp) + ' time points with missing data identified AFTER re-analysis.'
        # Save
        self.results.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')
        print self.embryo + ' analysis complete'
#==============================================================================
#   Parallel segmentation
#==============================================================================
    def runParEmbryoSegmentation(self):
        if self.species == 'rbalthica':
            # Loop over the folders for an embryo and apply segmentation.
            totts = time.time()
            for f in range(len(self.embryoFolders)):
                # Debug
                print 'Starting: ', self.embryoFolders[f]
                ts = time.time()
                self.currentFolder = self.embryoFolders[f]
                self.shortenedPath = self.shortenedPaths[f]
                #self.getShortenedPath()
                #self.shortenedPath = os.path.relpath(self.currentFolder, self.parentPath)
                # If no egg dims in first row, do not continue with attempting segmentation.       
                if (self.results.loc[self.shortenedPath]['eggRotBBox'][0] == ['nan','nan','nan','nan','nan'])|(math.isnan(float(self.results.loc[self.shortenedPath]['eggRotBBox'][0][0])))|(self.results.loc[self.shortenedPath]['eggRotBBox'][0] == 'nan')| (self.results.loc[self.shortenedPath]['eggRotBBox'][0] == np.NaN):
                    self.results.loc[self.shortenedPath]['embryoOutline'] = np.NaN
                    self.results.loc[self.shortenedPath]['centroidX'] = np.NaN
                    self.results.loc[self.shortenedPath]['centroidY'] = np.NaN
                    self.results.loc[self.shortenedPath]['area'] = np.NaN
                    self.results.loc[self.shortenedPath]['aspect'] = np.NaN
                    self.results.loc[self.shortenedPath]['extent'] = np.NaN
                    self.results.loc[self.shortenedPath]['hullArea'] = np.NaN
                    self.results.loc[self.shortenedPath]['solidity'] = np.NaN
                    self.results.loc[self.shortenedPath]['bboxMincol'] = np.NaN
                    self.results.loc[self.shortenedPath]['bboxMinrow'] = np.NaN
                    self.results.loc[self.shortenedPath]['bboxWidth'] = np.NaN
                    self.results.loc[self.shortenedPath]['bboxHeight'] = np.NaN
                    self.results.loc[self.shortenedPath]['blockWise'] = np.NaN
                    print 'No egg for', self.results.loc[self.shortenedPath]['currentFolder'][0]
                else:
                    self.seqImport(f)
                    self.applyEggMaskAndCrop()
                    print 'debug starting segmentation'
                    # Create a threading pool
                    innerPool = pathos.multiprocessing.ThreadingPool()
                    # Run segmentation
                    res = innerPool.map(self.segmentEmbryo, range(self.eggSeq.shape[0]))
                    # Save to dataframe
                    for g in range(len(res)):
                        self.results.loc[self.shortenedPath]['embryoOutline'][g] = res[g]['outline']
                        self.results.loc[self.shortenedPath]['centroidX'][g] = res[g]['centX']
                        self.results.loc[self.shortenedPath]['centroidY'][g] = res[g]['centY']
                        self.results.loc[self.shortenedPath]['area'][g] = res[g]['ar']
                        self.results.loc[self.shortenedPath]['aspect'][g] = res[g]['aspect']
                        self.results.loc[self.shortenedPath]['extent'][g] = res[g]['extent']
                        self.results.loc[self.shortenedPath]['hullArea'][g] = res[g]['hullArea']
                        self.results.loc[self.shortenedPath]['solidity'][g] = res[g]['solidity']
                        self.results.loc[self.shortenedPath]['bboxMincol'][g] = res[g]['bboxMinCol']
                        self.results.loc[self.shortenedPath]['bboxMinrow'][g] = res[g]['bboxMinRow']
                        self.results.loc[self.shortenedPath]['bboxWidth'][g] = res[g]['bboxWidth']
                        self.results.loc[self.shortenedPath]['bboxHeight'][g] = res[g]['bboxHeight']                      
                    # And create a cropped seq
                    #self.emSeq = self.seq[:,self.minY:self.maxY, self.minX:self.maxX]
                    res = innerPool.map(self.runNestedMeanWindowCalc, range(self.seq.shape[0]))                 
                    # Put data into Pandas DP
                    for g in range(len(res)):
                        self.results.loc[self.shortenedPath]['blockWise'][g] = res[g]
                    print self.embryo, ':', self.shortenedPath, 'analysed in {} s'.format(time.time()-ts) 
            print self.embryo, ':', ' results currently being saved.'    
            # Save
            self.results.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')
            print self.embryo,':',' results saved!' 
            print 'Analysis complete for', self.embryo, ' in {} s'.format(time.time()-totts) 
        # If Orchestia ..
        elif self.species == 'ogammarellus':
            # Loop over the folders for an embryo and apply segmentation.
            for f in range(len(self.embryoFolders)):
                print 'Starting: ', self.embryoFolders[f]
                ts = time.time()
                self.currentFolder = self.embryoFolders[f]
                self.shortenedPath = self.shortenedPaths[f]
                # Debug
                print 'Analysing: ', self.embryoFolders[f]
                ts = time.time()
                self.currentFolder = self.embryoFolders[f]
                self.shortenedPath = self.shortenedPaths[f]
                self.seqImport(f)
                #self.applyEggMaskAndCrop()
                # Create a threading pool
                innerPool = pathos.multiprocessing.ThreadingPool()
                # Run segmentation
                res = innerPool.map(self.segmentEmbryo, range(self.seq.shape[0]))
                # Save to dataframe
                for g in range(len(res)):
                    self.results.loc[self.shortenedPath]['embryoOutline'][g] = res[g]['outline']
                    self.results.loc[self.shortenedPath]['centroidX'][g] = res[g]['centX']
                    self.results.loc[self.shortenedPath]['centroidY'][g] = res[g]['centY']
                    self.results.loc[self.shortenedPath]['area'][g] = res[g]['ar']
                    self.results.loc[self.shortenedPath]['aspect'][g] = res[g]['aspect']
                    self.results.loc[self.shortenedPath]['extent'][g] = res[g]['extent']
                    self.results.loc[self.shortenedPath]['hullArea'][g] = res[g]['hullArea']
                    self.results.loc[self.shortenedPath]['solidity'][g] = res[g]['solidity']
                    self.results.loc[self.shortenedPath]['bboxMincol'][g] = res[g]['bboxMinCol']
                    self.results.loc[self.shortenedPath]['bboxMinrow'][g] = res[g]['bboxMinRow']
                    self.results.loc[self.shortenedPath]['bboxWidth'][g] = res[g]['bboxWidth']
                    self.results.loc[self.shortenedPath]['bboxHeight'][g] = res[g]['bboxHeight']                      
                # And create a cropped seq
                #self.emSeq = self.seq[:,self.minY:self.maxY, self.minX:self.maxX]
                res = innerPool.map(self.runNestedMeanWindowCalc, range(self.seq.shape[0]))
                # Put data into Pandas DP
                for g in range(len(res)):
                    self.results.loc[self.shortenedPath]['blockWise'][g] = res[g]
                print self.embryo, ':', self.shortenedPath, 'analysed in {} s'.format(time.time()-ts) 
        print self.embryo, ':', ' results currently being saved.'    
        # Save
        self.results.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')
        print self.embryo,':',' results saved!' 
        print 'Analysis complete for', self.embryo, ' in {} s'.format(time.time()-totts) 
     
                    
#==============================================================================
#     def parSeg(self, n):
#         # Segment embryo
#         self.segmentEmbryo(n)
#         return {'n':n, 'Embryo':self.embryo,'outline':self.embryoOutline,'centX': self.cx,'centY': self.cy, 'ar': self.area,'aspect': self.aspect, 'extent': self.extent,'hullArea': self.hullArea,'solidity': self.solidity,'bboxMinCol': self.x,'bboxMinRow': self.y,'bboxWidth': self.w,'bboxHeight': self.h}    
# 
#==============================================================================
        # Save
        # self.results.to_pickle(self.resultsDir + "/" + self.embryo + '.pandas')

#==============================================================================
#   Mask outside of egg from seq stack         
#==============================================================================
    def applyEggMask(self):
        # Take eggRotBBox from results, format appropriately and apply mask to self.seq
        for g in range(self.seq.shape[0]):
            # Apply egg mask.
            formattedRotBBox = ((self.results[self.shortenedPath]['eggRotBBox'][g][0],self.results[self.shortenedPath]['eggRotBBox'][g][1]),
            (self.results[self.shortenedPath]['eggRotBBox'][g][2],self.results[self.shortenedPath]['eggRotBBox'][g][3]),
            (self.results[self.shortenedPath]['eggRotBBox'][g][4]))
            out = cv2.ellipse(np.zeros_like(self.seq[g]),formattedRotBBox,(255,255,255),-1)
            self.seq[g] = cv2.bitwise_and(self.seq[g], out)
 
#==============================================================================
#   Mask outside of egg from seq stack         
#==============================================================================
    def applyEggMaskAndCrop(self):
        # Reduce seq to egg BB and apply egg mask.
        out = self.getSeqEggBB()
        eggMinX, eggMinY, eggMaxX, eggMaxY = out['eggMinX'], out['eggMinY'], out['eggMaxX'], out['eggMaxY']
        formattedRotBBox = ((self.results[self.shortenedPath]['eggRotBBox'][0][0]-eggMinX,self.results[self.shortenedPath]['eggRotBBox'][0][1]-eggMinY),
        (self.results[self.shortenedPath]['eggRotBBox'][0][2],self.results[self.shortenedPath]['eggRotBBox'][0][3]),
        (self.results[self.shortenedPath]['eggRotBBox'][0][4]))
        self.eggSeq = self.seq[:,eggMinY:eggMaxY,eggMinX:eggMaxX]
        for g in range(self.eggSeq.shape[0]):
            out = cv2.ellipse(np.zeros_like(self.eggSeq[g]),formattedRotBBox,(255,255,255),-1)
            self.eggSeq[g] = cv2.bitwise_and(self.eggSeq[g], out)
            
#==============================================================================
#   Get egg ROI
#==============================================================================
    def getSeqEggBB(self):
        #self.results.ix[self.shortenedPath,:,'eggBoxPoints']
        # Currently assumes egg ROI fixed during sequence..
        eggMinX = int(np.array(self.results.ix[self.shortenedPath,:,'eggBoxPoints'][0])[:,0].min())
        eggMaxX = int(np.array(self.results.ix[self.shortenedPath,:,'eggBoxPoints'][0])[:,0].max())
        eggMinY = int(np.array(self.results.ix[self.shortenedPath,:,'eggBoxPoints'][0])[:,1].min())
        eggMaxY = int(np.array(self.results.ix[self.shortenedPath,:,'eggBoxPoints'][0])[:,1].max())
        # Check if min or max are outside im region, in which restrict..
        if (eggMinX < 0) | (eggMinY < 0) | (eggMinX < 0) | (eggMinX < 0):
            if (eggMinX < 0):
                eggMinX = 0
            if (eggMinY < 0):
                eggMinY = 0
            if (eggMaxX > self.seq.shape[2]):
                eggMaxX = self.seq.shape[2]
            if (eggMaxY > self.seq.shape[1]):
                eggMaxY = self.seq.shape[1]                
        return {'eggMinX':eggMinX, 'eggMinY':eggMinY, 'eggMaxX':eggMaxX, 'eggMaxY':eggMaxY}
        #cv2.boundingRect()
#==============================================================================
#   Crop seq to embryo BB
#==============================================================================
    def getEmbryoBB(self):   
        # If Radix balthica...
        if self.species == 'rbalthica':
            # Check for python nan
            tmp=[]
            for f in self.results.ix[self.shortenedPath,:,'bboxWidth'].values:
                tmp.append(math.isnan(f))
                
            if np.sum(tmp)>0|np.sum(self.results.ix[self.shortenedPath,:,'bboxWidth'] =='nan')|np.sum(self.results.ix[self.shortenedPath,:,'bboxWidth'] ==np.NaN) >0|np.sum(self.results.ix[self.shortenedPath,:,'eggBoxPoints']==np.NaN) >0:
                maxX = np.NaN
                maxY = np.NaN
                minX = np.NaN
                minY = np.NaN
            
            else:
                # Get egg BB to add
                out = self.getSeqEggBB()
                eggMinX, eggMinY, eggMaxX, eggMaxY = out['eggMinX'], out['eggMinY'], out['eggMaxX'], out['eggMaxY']
                maxX = eggMinX + np.max(self.results.ix[self.shortenedPath,:,'bboxWidth'] + self.results.ix[self.shortenedPath,:,'bboxMincol'])
                maxY = eggMinY + np.max(self.results.ix[self.shortenedPath,:,'bboxHeight'] + self.results.ix[self.shortenedPath,:,'bboxMinrow'])
                minX = eggMinX + self.results.ix[self.shortenedPath,:,'bboxMincol'].min()
                minY = eggMinY + self.results.ix[self.shortenedPath,:,'bboxMinrow'].min()
        
        # If Orchestia gammarellus..
        elif self.species == 'ogammarellus':
            # Check for python nan
            maxX = np.max(self.results.ix[self.shortenedPath,:,'bboxWidth'] + self.results.ix[self.shortenedPath,:,'bboxMincol'])
            maxY = np.max(self.results.ix[self.shortenedPath,:,'bboxHeight'] + self.results.ix[self.shortenedPath,:,'bboxMinrow'])
            minX = self.results.ix[self.shortenedPath,:,'bboxMincol'].min()
            minY = self.results.ix[self.shortenedPath,:,'bboxMinrow'].min()
            
        return {'minX':minX, 'minY':minY, 'maxX':maxX, 'maxY':maxY}

#==============================================================================
#   Get mean window calculation at multiple resolutions.. nested.
#==============================================================================
    def runNestedMeanWindowCalc(self,f):   
        out = self.getEmbryoBB()
        minX, minY, maxX, maxY = out['minX'], out['minY'], out['maxX'], out['maxY']
        
        if (minX != np.NaN) and (not math.isnan(minX)):
            #im = self.seq[f,self.minY:self.maxY, self.minX:self.maxX]
            # Get shape dims
            imy = maxY - minY
            imx = maxX - minX
    
            # To be able to divide by 2,4,8,16 the embryoBB may need to be expanded slightly.
            # Use modulo 16 to find out to what value the range must be expanded to make this possible.
            expandY = (imy%16)
            expandX = (imx%16)
            skipY = False
            if expandY !=0:
                # Test whether subtraction is required and if so modify expandY
                if ((imy-imy%16)%16)==0:
                    expandY = 16 - (imy%16)
                # If either minX or minY ==1 or max im size
                if ((minY-(expandY/2)) <= 0) | ((maxY+(expandY/2))>self.seq.shape[1]):
                    if ((minY-(expandY/2)) <= 0):
                        maxY = (maxY + expandY)
                        imy = maxY - minY
                        skipY = True
                    elif ((maxY+(expandY/2))>self.seq.shape[1]):
                        minY = (minY - expandY)  
                        imy = maxY - minY
                        skipY = True
                # If expansion/contraction does not reach img limits...
                if not skipY:
                    if (not((minY-(expandY/2)) <= 0)) | (not((maxY+(expandY/2))>self.seq.shape[1])):
                        # Even
                        if expandY%2 ==0:
                            minY = int(minY - (expandY/2))
                            maxY = int(maxY + (expandY/2))
                            imy = maxY - minY
                        # Odd
                        else:
                            minY = int(minY - (math.ceil(expandY/2.)))
                            maxY = int(maxY + (math.floor(expandY/2.)))
                            imy = maxY - minY              
                
            skipX = False
            if expandX !=0:
                # Test whether subtraction is required and if so modify expandY
                if ((imx-imx%16)%16)==0:
                    expandX = 16 - (imx%16)
                # If either minX or minY ==1 or max im size
                if ((minX-(expandX/2)) <= 0) | ((maxX+(expandX/2))>self.seq.shape[2]):
                    if ((minX-(expandX/2)) <= 0):
                        maxX = (maxX + expandX)
                        imx = maxX - minX
                        skipX = True
                    elif ((maxX+(expandX/2))>self.seq.shape[2]):
                        minX = (minX - expandX)  
                        imx = maxX - minX
                        skipX = True
                # If expansion/contraction does not reach img limits...
                if not skipX:
                    if (not((minX-(expandX/2)) <= 0)) | (not((maxX+(expandX/2))>self.seq.shape[2])):
                    # Even
                        if expandX%2 ==0:
                            minX = int(minX - (expandX/2))
                            maxX = int(maxX + (expandX/2))
                            imx = maxX - minX
                        # Odd
                        else:
                            minX = int(minX - (math.ceil(expandX/2.)))
                            maxX = int(maxX + (math.floor(expandX/2.)))
                            imx = maxX - minX         
            # Debug
            imx = imx
            imy = imy
            
            # Extract frame relevant to f
            im = np.ascontiguousarray(self.seq[f,minY:maxY, minX:maxX])
            minWinSize =60
            # 60 Pixel window calculations
            rangeY = maxY - minY
            rangeX = maxX - minX
            winNoY = rangeY/minWinSize
            winNoX = rangeX/minWinSize
            
            # List to collect results
            nestMeanWindOut = []
            # First get overall embryoBB mean
            nestMeanWindOut.append(im.mean())
            
            # Now get 2 x 2 means   
            out = skimage.util.view_as_blocks(im,block_shape = (imy/2,imx/2))
            nestMeanWindOut.append(((out[0,0].mean(),out[0,1].mean()),(out[1,0].mean(),out[1,1].mean())))
    
            # Now get 4 x 4 means   
            out = skimage.util.view_as_blocks(im,block_shape = (imy/4,imx/4))
            nestMeanWindOut.append(((out[0,0].mean(),out[0,1].mean(), out[0,2].mean(),out[0,3].mean()),
                                    (out[1,0].mean(),out[1,1].mean(),out[1,2].mean(),out[1,3].mean()),
                                    (out[2,0].mean(),out[2,1].mean(),out[2,2].mean(),out[2,3].mean()),
                                    (out[3,0].mean(),out[3,1].mean(),out[3,2].mean(),out[3,3].mean())))
    
            # Now get 8 x 8 means   
            out = skimage.util.view_as_blocks(im,block_shape = (imy/8,imx/8))
            nestMeanWindOut.append(((out[0,0].mean(),out[0,1].mean(), out[0,2].mean(),out[0,3].mean(),
                                    out[0,4].mean(),out[0,5].mean(),out[0,6].mean(),out[0,7].mean()),
                                    (out[1,0].mean(),out[1,1].mean(), out[1,2].mean(),out[1,3].mean(),
                                    out[1,4].mean(),out[1,5].mean(),out[1,6].mean(),out[1,7].mean()),
                                    (out[2,0].mean(),out[2,1].mean(), out[2,2].mean(),out[2,3].mean(),
                                    out[2,4].mean(),out[2,5].mean(),out[2,6].mean(),out[2,7].mean()),
                                    (out[3,0].mean(),out[3,1].mean(), out[3,2].mean(),out[3,3].mean(),
                                    out[3,4].mean(),out[3,5].mean(),out[3,6].mean(),out[3,7].mean()),
                                    (out[4,0].mean(),out[4,1].mean(), out[4,2].mean(),out[4,3].mean(),
                                    out[4,4].mean(),out[4,5].mean(),out[4,6].mean(),out[4,7].mean()),
                                    (out[5,0].mean(),out[5,1].mean(), out[5,2].mean(),out[5,3].mean(),
                                    out[5,4].mean(),out[5,5].mean(),out[5,6].mean(),out[5,7].mean()),
                                    (out[6,0].mean(),out[6,1].mean(), out[6,2].mean(),out[6,3].mean(),
                                    out[6,4].mean(),out[6,5].mean(),out[6,6].mean(),out[6,7].mean()),
                                    (out[7,0].mean(),out[7,1].mean(), out[7,2].mean(),out[7,3].mean(),
                                    out[7,4].mean(),out[7,5].mean(),out[7,6].mean(),out[7,7].mean())))
    
            # Now get 16 x 16 means   
            out = skimage.util.view_as_blocks(im,block_shape = (imy/16,imx/16))
            nestMeanWindOut.append(((out[0,0].mean(),out[0,1].mean(), out[0,2].mean(),out[0,3].mean(),
                                    out[0,4].mean(),out[0,5].mean(),out[0,6].mean(),out[0,7].mean(),
                                    out[0,8].mean(),out[0,9].mean(), out[0,10].mean(),out[0,11].mean(),
                                    out[0,12].mean(),out[0,13].mean(),out[0,14].mean(),out[0,15].mean()),
                                    (out[1,0].mean(),out[1,1].mean(), out[1,2].mean(),out[1,3].mean(),
                                    out[1,4].mean(),out[1,5].mean(),out[1,6].mean(),out[1,7].mean(),
                                    out[1,8].mean(),out[1,9].mean(), out[1,10].mean(),out[1,11].mean(),
                                    out[1,12].mean(),out[1,13].mean(),out[1,14].mean(),out[1,15].mean()),
                                    (out[2,0].mean(),out[2,1].mean(), out[2,2].mean(),out[2,3].mean(),
                                    out[2,4].mean(),out[2,5].mean(),out[2,6].mean(),out[2,7].mean(),
                                    out[2,8].mean(),out[2,9].mean(), out[2,10].mean(),out[2,11].mean(),
                                    out[2,12].mean(),out[2,13].mean(),out[2,14].mean(),out[2,15].mean()),
                                    (out[3,0].mean(),out[3,1].mean(), out[3,2].mean(),out[3,3].mean(),
                                    out[3,4].mean(),out[3,5].mean(),out[3,6].mean(),out[3,7].mean(),
                                    out[3,8].mean(),out[3,9].mean(), out[3,10].mean(),out[3,11].mean(),
                                    out[3,12].mean(),out[3,13].mean(),out[3,14].mean(),out[3,15].mean()),
                                    (out[4,0].mean(),out[4,1].mean(), out[4,2].mean(),out[4,3].mean(),
                                    out[4,4].mean(),out[4,5].mean(),out[4,6].mean(),out[4,7].mean(),
                                    out[4,8].mean(),out[4,9].mean(), out[4,10].mean(),out[4,11].mean(),
                                    out[4,12].mean(),out[4,13].mean(),out[4,14].mean(),out[4,15].mean()),
                                    (out[5,0].mean(),out[5,1].mean(), out[5,2].mean(),out[5,3].mean(),
                                    out[5,4].mean(),out[5,5].mean(),out[5,6].mean(),out[5,7].mean(),
                                    out[5,8].mean(),out[5,9].mean(), out[5,10].mean(),out[5,11].mean(),
                                    out[5,12].mean(),out[5,13].mean(),out[5,14].mean(),out[5,15].mean()),
                                    (out[6,0].mean(),out[6,1].mean(), out[6,2].mean(),out[6,3].mean(),
                                    out[6,4].mean(),out[6,5].mean(),out[6,6].mean(),out[6,7].mean(),
                                    out[6,8].mean(),out[6,9].mean(), out[6,10].mean(),out[6,11].mean(),
                                    out[6,12].mean(),out[6,13].mean(),out[6,14].mean(),out[6,15].mean()),
                                    (out[7,0].mean(),out[7,1].mean(), out[7,2].mean(),out[7,3].mean(),
                                    out[7,4].mean(),out[7,5].mean(),out[7,6].mean(),out[7,7].mean(),
                                    out[7,8].mean(),out[7,9].mean(), out[7,10].mean(),out[7,11].mean(),
                                    out[7,12].mean(),out[7,13].mean(),out[7,14].mean(),out[7,15].mean()),
                                    (out[8,0].mean(),out[8,1].mean(), out[8,2].mean(),out[8,3].mean(),
                                    out[8,4].mean(),out[8,5].mean(),out[8,6].mean(),out[8,7].mean(),
                                    out[8,8].mean(),out[8,9].mean(), out[8,10].mean(),out[8,11].mean(),
                                    out[8,12].mean(),out[8,13].mean(),out[8,14].mean(),out[8,15].mean()),
                                    (out[9,0].mean(),out[9,1].mean(), out[9,2].mean(),out[9,3].mean(),
                                    out[9,4].mean(),out[9,5].mean(),out[9,6].mean(),out[9,7].mean(),
                                    out[9,8].mean(),out[9,9].mean(), out[9,10].mean(),out[9,11].mean(),
                                    out[9,12].mean(),out[9,13].mean(),out[9,14].mean(),out[9,15].mean()),
                                    (out[10,0].mean(),out[10,1].mean(), out[10,2].mean(),out[10,3].mean(),
                                    out[10,4].mean(),out[10,5].mean(),out[10,6].mean(),out[10,7].mean(),
                                    out[10,8].mean(),out[10,9].mean(), out[10,10].mean(),out[10,11].mean(),
                                    out[10,12].mean(),out[10,13].mean(),out[10,14].mean(),out[10,15].mean()),
                                    (out[11,0].mean(),out[11,1].mean(), out[11,2].mean(),out[11,3].mean(),
                                    out[11,4].mean(),out[11,5].mean(),out[11,6].mean(),out[11,7].mean(),
                                    out[11,8].mean(),out[11,9].mean(), out[11,10].mean(),out[11,11].mean(),
                                    out[11,12].mean(),out[11,13].mean(),out[11,14].mean(),out[11,15].mean()),
                                    (out[12,0].mean(),out[12,1].mean(), out[12,2].mean(),out[12,3].mean(),
                                    out[12,4].mean(),out[12,5].mean(),out[12,6].mean(),out[12,7].mean(),
                                    out[12,8].mean(),out[12,9].mean(), out[12,10].mean(),out[12,11].mean(),
                                    out[12,12].mean(),out[12,13].mean(),out[12,14].mean(),out[12,15].mean()),
                                    (out[13,0].mean(),out[13,1].mean(), out[13,2].mean(),out[13,3].mean(),
                                    out[13,4].mean(),out[13,5].mean(),out[13,6].mean(),out[13,7].mean(),
                                    out[13,8].mean(),out[13,9].mean(), out[13,10].mean(),out[13,11].mean(),
                                    out[13,12].mean(),out[13,13].mean(),out[13,14].mean(),out[13,15].mean()),
                                    (out[14,0].mean(),out[14,1].mean(), out[14,2].mean(),out[14,3].mean(),
                                    out[14,4].mean(),out[14,5].mean(),out[14,6].mean(),out[14,7].mean(),
                                    out[14,8].mean(),out[14,9].mean(), out[14,10].mean(),out[14,11].mean(),
                                    out[14,12].mean(),out[14,13].mean(),out[14,14].mean(),out[14,15].mean()),
                                    (out[15,0].mean(),out[15,1].mean(), out[15,2].mean(),out[15,3].mean(),
                                    out[15,4].mean(),out[15,5].mean(),out[15,6].mean(),out[15,7].mean(),
                                    out[15,8].mean(),out[15,9].mean(), out[15,10].mean(),out[15,11].mean(),
                                    out[15,12].mean(),out[15,13].mean(),out[15,14].mean(),out[15,15].mean())))
        else:
            nestMeanWindOut = np.NaN
            
        return nestMeanWindOut
    
#%%
#==============================================================================
#   Image segmentation methods
#==============================================================================
#%%
# Adaptive form of thresholding.
    @staticmethod
    def auto_thresh(image,sigma = 0.1): 
        v = np.mean(image[image!=0])
        lower = int(max(0, (1.0 - sigma) * v))
        upper = int(min(255, (1.0 + sigma) * v))   
        ret2,thresh = cv2.threshold(image,lower, upper,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
        # ret2,thresh = cv2.threshold(blur,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
        # Debug
        # plt.imshow(thresh)
        return thresh
    
    #%% 
    @staticmethod     
    def auto_canny(image, sigma=0.33):
        # compute the median of the single channel pixel intensities
        v = np.median(image[image!=0])
        # apply automatic Canny edge detection using the computed median
        lower = int(max(0, (1.0 - sigma) * v))
        upper = int(min(255, (1.0 + sigma) * v))
        edged = cv2.Canny(image, lower, upper)
        # return the edged image
        return edged
        
    #%% 
    @staticmethod     
    def auto_varSeg(image, sigma=0.33):
        # compute the median of the single channel pixel intensities
        v = np.mean(image[image!=0])
        # apply automatic Canny edge detection using the computed median
        lower = int(max(0, (1.0 - sigma) * v))
        upper = int(min(255, (1.0 + sigma) * v))
        edged = cv2.Canny(image, lower, upper)
        # return the edged image
        return edged    
        
