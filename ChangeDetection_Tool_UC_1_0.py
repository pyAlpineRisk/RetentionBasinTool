# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   ChangeDetection Tool                        *
*   Nicole Kamp                                         *
*   Jänner 2021                                                            *
*                                                                         *
***************************************************************************
"""

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProject,
                       QgsVectorLayer,
                       QgsTextFormat,
                       QgsExpression,
                       QgsFeatureRequest,
                       QgsFeature,
                       QgsGeometry,
                       QgsPoint,
                       QgsPointXY,
                       QgsVectorFileWriter,
                       QgsRasterBandStats,
                       QgsColorRampShader,
                       QgsRasterTransparency,
                       QgsFillSymbol,
                       QgsRasterShader,
                       QgsSingleBandPseudoColorRenderer,
                       QgsWkbTypes,
                       QgsVectorLayerSimpleLabeling,
                       QgsPalLayerSettings,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterFile,
                       QgsRasterLayer,
                       QgsCoordinateReferenceSystem,
                       QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterExtent,
                       QgsProcessingMultiStepFeedback,
                       QgsMessageLog)
from qgis import processing
from qgis.analysis import QgsRasterCalculatorEntry, QgsRasterCalculator
from qgis.PyQt.QtWidgets import QApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtGui import QColor
from qgis.utils import iface

from PIL import Image
Image.MAX_IMAGE_PIXELS = None
import osgeo.gdal as gdal
import numpy as np
import csv

from osgeo import ogr, osr
from shapely.geometry import Polygon
from matplotlib import rcParams
import matplotlib.pyplot as plt 
import pandas as pd
import subprocess
from subprocess import call

import string, os, sys, copy, shutil, math, numpy, time, datetime
from time import *
from sys import *


#class ProjectUTM33N

#class AddSurfaceInfo

#class AddFields

#class LoopHeight


class ChangeDetectionProcessingAlgorithm(QgsProcessingAlgorithm):
    """
Python Script zur Ermittlung von Steinschlag-Impact Punkten auf Basis des Pauschalgefälles.
Eingabe der Input- und Output-Parameter
    """
    INPUT_Shape = 'INPUT_SHP'
    INPUT_ALS_new = 'DGM1'
    INPUT_ALS_old = 'DGM2'
    INPUT_UC = 'UC'
    #INPUT_extent = 'AUSSCHNITT'
    TEMP = 'TEMP'
    
    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return ChangeDetectionProcessingAlgorithm()

    def name(self):
        return 'changedetectiontool'

    def displayName(self):
        return self.tr('QCD_ChangeDetectionUC')

    def group(self):
        return self.tr('pyAlpineRisk')

    def groupId(self):
        return 'pyAlpineRiskScripts'

    def shortHelpString(self):
        return self.tr("GIS-Tool zur Ermittlung der Veränderung von zwei Höhenmodellen (Kamp, 2021).")
        
        
    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_Shape,
                self.tr('StudyArea'),
                [QgsProcessing.TypeVectorAnyGeometry],
                defaultValue='G:/01_dissertation_neu/06_tools/01_lidar/21_QGIS_ChangeDetection/00_data_neu/vectordata/mask_SA1_4_teilA.shp'
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_ALS_new,
                self.tr('Geländemodell_neu'),
                #extension='tiff',
                #fileFilter="tiff (*.tif)",
                defaultValue='G:/01_dissertation_neu/06_tools/01_lidar/21_QGIS_ChangeDetection/00_data_neu/rasterdata/SA1/PC2_dtm50cm_rds_SA1_ULS.tif'
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_ALS_old,
                self.tr('Geländemodell_alt'),
                #extension='tiff',
                #fileFilter="tiff (*.tif)",
                defaultValue='G:/01_dissertation_neu/06_tools/01_lidar/21_QGIS_ChangeDetection/00_data_neu/rasterdata/SA1/PC1_dtm50cm_rds_SA1_ULS.tif'
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_UC,
                self.tr('UC_Threshold'),
                #extension='tiff',
                #fileFilter="tiff (*.tif)",
                defaultValue='G:/01_dissertation_neu/06_tools/01_lidar/21_QGIS_ChangeDetection/00_data_neu/rasterdata/SA1/UC_50cm_SA1_ULS.tif'
            )
        )

              
        self.addParameter(
            QgsProcessingParameterFolderDestination(
                self.TEMP, 
                self.tr('Output-Ordner'), 
                defaultValue='C:/temp/CD'
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        feedback = QgsProcessingMultiStepFeedback(1, feedback)
        results = {}
        outputs = {}
        
        
        # Pfade definieren + Timestamp
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        temp_path = str(parameters[self.TEMP])+'/temp_'+str(timestamp)
        final_path = str(parameters[self.TEMP])+'/final_'+str(timestamp)
        
        
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
        if not os.path.exists(final_path):
            os.makedirs(final_path)
        
        #feedback.pushInfo(str(parameters[self.INPUT_extent]))


        ## Buffer
        buffer = out = temp_path + '/' + 'extent.shp'
        alg_params = {
            'DISSOLVE': False,
            'DISTANCE': 10,
            'END_CAP_STYLE': 2,
            'INPUT': str(parameters[self.INPUT_Shape]),
            'JOIN_STYLE': 1,
            'MITER_LIMIT': 2,
            'SEGMENTS': 1,
            'OUTPUT': buffer
        }
        processing.run('native:buffer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        ## Clip Rasters
        out_new = temp_path + '/' + 'dhm_new.tif'
        out_old = temp_path + '/' + 'dhm_old.tif'
        out_uc = temp_path + '/' + 'dhm_uc.tif'
        alg_params = {
            'ALPHA_BAND': False,
            'CROP_TO_CUTLINE': True,
            'DATA_TYPE': 0,
            'EXTRA': '',
            'INPUT': str(parameters[self.INPUT_ALS_new]),
            'KEEP_RESOLUTION': False,
            'MASK': str(buffer),
            'MULTITHREADING': False,
            'NODATA': None,
            'OPTIONS': '',
            'SET_RESOLUTION': False,
            'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:32633'),
            'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:32633'),
            'X_RESOLUTION': None,
            'Y_RESOLUTION': None,
            'OUTPUT': out_new
        }
        processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        alg_params = {
            'ALPHA_BAND': False,
            'CROP_TO_CUTLINE': True,
            'DATA_TYPE': 0,
            'EXTRA': '',
            'INPUT': str(parameters[self.INPUT_ALS_old]),
            'KEEP_RESOLUTION': False,
            'MASK': str(buffer),
            'MULTITHREADING': False,
            'NODATA': None,
            'OPTIONS': '',
            'SET_RESOLUTION': False,
            'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:32633'),
            'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:32633'),
            'X_RESOLUTION': None,
            'Y_RESOLUTION': None,
            'OUTPUT': out_old
        }
        processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        alg_params = {
            'ALPHA_BAND': False,
            'CROP_TO_CUTLINE': True,
            'DATA_TYPE': 0,
            'EXTRA': '',
            'INPUT': str(parameters[self.INPUT_UC]),
            'KEEP_RESOLUTION': False,
            'MASK': str(buffer),
            'MULTITHREADING': False,
            'NODATA': None,
            'OPTIONS': '',
            'SET_RESOLUTION': False,
            'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:32633'),
            'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:32633'),
            'X_RESOLUTION': None,
            'Y_RESOLUTION': None,
            'OUTPUT': out_uc
        }
        processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        del(out_uc)
        
        
        ## Cell Size/ Extent
        raster_cs = gdal.Open(out_new)
        gt_cs =raster_cs.GetGeoTransform() 
        cs = gt_cs[1]
        
        del raster_cs
        del gt_cs
        
        dtm_clip_gdal = gdal.Open(out_new)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
        band1 = dtm_clip_gdal.GetRasterBand(1)
        im = Image.open(out_new)
        pix = im.load()
        gt = dtm_clip_gdal.GetGeoTransform()
        feedback.pushInfo(str(gt))
        width = dtm_clip_gdal.RasterXSize
        height = dtm_clip_gdal.RasterYSize
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5] 
        maxx = gt[0] + width*gt[1] + height*gt[2]
        maxy = gt[3]
        extent=str(minx)+","+str(maxx)+","+str(miny)+","+str(maxy)
        extent2=(minx,maxx,miny,maxy)
        feedback.pushInfo(str(extent))
        feedback.pushInfo(str(width))
        feedback.pushInfo(str(height))
        
        # Create empty grids
        grid_diff = np.zeros(shape=(width,height), dtype=np.float32)
        grid_diff_uc = np.zeros(shape=(width,height), dtype=np.float32)
               

        ## DOD
        out2_old = temp_path + '/' + 'diff.tif'
 
        ##Rasters to Arrays
        new_rds = gdal.Open(out_new)
        format = "GTiff"
        driver1 = gdal.GetDriverByName( format )
        band_new = new_rds.GetRasterBand(1)
        im_new = Image.open(out_new)
        pix_new = im_new.load()
        
        old_rds = gdal.Open(out_old)
        format = "GTiff"
        driver2 = gdal.GetDriverByName( format )
        band_old = old_rds.GetRasterBand(1)
        im_old = Image.open(out_old)
        pix_old = im_old.load()

                
        for row in range(0, width):
            for col in range(0, height):
                val_new = pix_new[row, col]
                val_old = pix_old[row, col]
                if val_new > -9999.0:
                    val_diff = val_new-val_old
                    #feedback.pushInfo(str(val_new))
                    grid_diff[row, col] = val_diff
        
        grid_diff=np.flip(grid_diff,1)
        grid_diff=np.rot90(grid_diff)
        imsave = Image.fromarray(grid_diff, mode='F')
        imsave.save(out2_old, "TIFF")
        
        ## Raster georeferenzieren
        out2 = final_path + '/' + 'diff.tif'
        
        src_ds = gdal.Open(out2_old)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
                
        dst_ds = driver.CreateCopy(out2, src_ds, 0)
        dst_ds.SetGeoTransform(gt)
        epsg = 32633
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
        dest_wkt = srs.ExportToWkt()
        dst_ds.SetProjection(dest_wkt)

        # Close files
        dst_ds = None
        src_ds = None
        
        ## ChangeDetection
        analyse = open(final_path + '/' + 'change_detection.txt', "w")
        analyse_csv = open(final_path + '/' + 'change_detection.csv', "w")
        analyse_csv.write(" ; AOI [m2]; Detectable Change [m2]; Detectable Change [%]; Surface Lowering [m2]; Surface Raising [m2]; Surface Lowering [m3]; Surface Raising [m3]; Volume of Difference [m3]; Net Volume of Difference [m3]" + "\n")
 
 
        ##Rasters to Arrays
        out_uc = temp_path + '/' + 'dhm_uc.tif'
                
        diff_rds = gdal.Open(out2)
        format = "GTiff"
        driver3 = gdal.GetDriverByName( format )
        band_diff = diff_rds.GetRasterBand(1)
        im_diff = Image.open(out2)
        pix_diff = im_diff.load()
        
        uc_rds = gdal.Open(out_uc)
        format = "GTiff"
        driver4 = gdal.GetDriverByName( format )
        band_uc = uc_rds.GetRasterBand(1)
        im_uc = Image.open(out_uc)
        pix_uc = im_uc.load()
        
        ## Classification1 - 0 % 0.1
        countAOI=0
        countPosGes = 0
        countNegGes = 0
        countPos = 0
        countNeg = 0
        countPosArea = 0
        countNegArea = 0
        countPosAreaGes = 0
        countNegAreaGes = 0
        countAcc = 0
        countCell = 0
        countCellVol = 0
        
        for row in range(0, width):
            for col in range(0, height):
                diff = pix_diff[row, col]
                if diff > -9999.0 and diff < 100:
                    countAOI = countAOI+(cs*cs)
                if diff < 0 and diff > -100:
                    volNegGes = cs * cs * (abs(diff))
                    countNegGes = countNegGes + volNegGes
                    countNegAreaGes = countNegAreaGes + (cs * cs)
                if diff > 0 and diff < 100:
                    volPosGes = cs * cs * (abs(diff))
                    countPosGes = countPosGes + volPosGes
                    countPosAreaGes = countPosAreaGes + (cs * cs)
                if diff < -0.1 and diff > -100: # 0.1 m ist der Standardfehler zwischen den 2 Modellen
                    volNeg = cs * cs * (abs(diff)-0.1)# -0.1 Standardfehler wird bei GCD-Tool nicht abgezogen
                    countNeg = countNeg + volNeg
                    countNegArea = countNegArea + (cs * cs)  
                    countAcc = countAcc + (cs * cs * 0.1)
                if diff > 0.1 and diff < 100:
                    volPos = cs * cs * (diff-0.1)
                    countPos = countPos + (volPos)
                    countPosArea = countPosArea + (cs * cs) 
                    countAcc = countAcc + (cs * cs * 0.1)
                if diff < 0 and diff > -0.1:
                    countCell = countCell + 1
                    countCellVol = countCellVol + (cs * cs * abs(diff))
                    #print countCell
                if diff > 0 and diff < 0.1:
                    countCell = countCell + 1
                    countCellVol = countCellVol + (cs * cs * abs(diff))
            
        analyse.write("Total Area [m2] of Interest: " + str(countAOI)+"\n")
        analyse.write("Total Area [km2] of Interest: " + str(countAOI/1000000)+"\n")
        analyse.write("\n") 
        analyse.write("Total Area [m2] of Detectable Change: " + str(countPosAreaGes+countNegAreaGes)+"\n")
        analyse.write("Total Area of Interest [%] with Detectable Change: " + str(((countPosAreaGes+countNegAreaGes)/countAOI)*100)+"\n")
        analyse.write("Total Area [m2] of Surface Lowering: " + str(countNegAreaGes)+"\n")
        analyse.write("Total Area [m2] of Surface Raising: " + str(countPosAreaGes)+"\n")
        analyse.write("Total Volume [m3] of Surface Lowering: " + str(countNegGes)+"\n")
        analyse.write("Total Volume [m3] of Surface Raising: " + str(countPosGes)+"\n")
        analyse.write("Total Volume [m3] of Difference: " + str(countPosGes+countNegGes)+"\n")
        analyse.write("Total Net Volume [m3] of Difference: " + str(countPosGes-countNegGes)+"\n")
        analyse.write("\n")
        analyse.write("\n")
        
        analyse_csv.write("0.0; "+str(countAOI)+"; "+str(countPosAreaGes+countNegAreaGes)+"; "+str(((countPosAreaGes+countNegAreaGes)/countAOI)*100)+"; "+str(countNegAreaGes)+"; "+str(countPosAreaGes)+"; "+str(countNegGes)+"; "+str(countPosGes)+"; "+str(countPosGes+countNegGes)+"; "+str(countPosGes-countNegGes)+"\n")
        net_v_0 = countPosGes-countNegGes
        sr_0 = countPosGes
        sl_0 = countNegGes
        
        analyse.write("Analysis with Threshold of +/- 0.1"+"\n")
        analyse.write("Thresholded (0.1) Area [m2] of Detectable Change: " + str(countPosArea+countNegArea)+"\n")
        analyse.write("Thresholded (0.1) Area [km2] of Detectable Change: " + str((countPosArea+countNegArea)/1000000)+"\n")
        analyse.write("Thresholded (0.1) Area of Interest [%] with Detectable Change: " + str(((countPosArea+countNegArea)/countAOI)*100)+"\n")
        analyse.write("Thresholded (0.1) Area [m2] of Surface Lowering: " + str(countNegArea)+"\n")
        analyse.write("Thresholded (0.1) Area [km2] of Surface Lowering: " + str(countNegArea/1000000)+"\n")
        analyse.write("Thresholded (0.1) Area [m2] of Surface Raising: " + str(countPosArea)+"\n")
        analyse.write("Thresholded (0.1) Area [km2] of Surface Raising: " + str(countPosArea/1000000)+"\n")
        analyse.write("Thresholded (0.1) Volume [m3] of Surface Lowering: " + str(countNeg)+"\n")
        analyse.write("Thresholded (0.1) Volume [m3] of Surface Raising: " + str(countPos)+"\n")
        analyse.write("Thresholded (0.1) Volume [m3] of Difference: " + str(countPos+countNeg)+"\n")
        analyse.write("Thresholded (0.1) Net Volume [m3] of Difference: " + str(countPos-countNeg)+"\n")
        analyse.write("\n")
        analyse.write("\n") 
        analyse.write("Volume [m3] of Error within Threshold of -0.1 and 0.1: " + str(countAcc)+"\n")
        analyse.write("Count of Cells within -0.1 and 0.1: " + str(countCell)+"\n")
        analyse.write("Volume [m3] of Cells between -0.1 and 0.1: " + str(countCellVol)+"\n")
        analyse.write("Percent [%] of Error of Cells between -0.1 and 0.1: " + str((countCellVol/(countPosGes+countNegGes))*100)+"\n")
        analyse.write("\n")
        analyse.write("\n") 
        
        analyse_csv.write("0.1; "+str(countAOI)+"; "+str(countPosArea+countNegArea)+"; "+str(((countPosArea+countNegArea)/countAOI)*100)+"; "+str(countNegArea)+"; "+str(countPosArea)+"; "+str(countNeg)+"; "+str(countPos)+"; "+str(countPos+countNeg)+"; "+str(countPos-countNeg)+"\n")
        net_v_1 = countPos-countNeg
        sr_1 = countPos
        sl_1 = countNeg

        del countPos
        del countNeg
        del countPosArea
        del countNegArea
        del countAcc
        del countCell
        del countCellVol
 

        ## Classification2 - 0.3
        countPos = 0
        countNeg = 0
        countPosArea = 0
        countNegArea = 0
        countAcc = 0
        countCell = 0
        countCellVol = 0
    
        for row in range(0, width):
            for col in range(0, height):
                diff = pix_diff[row, col]
                ES = pix_uc[row, col]
                ES2 = abs(ES)
                if diff < -0.3 and diff > -100: # 0.1 ist der Standardfehler zwischen den 2 Modellen
                    volNeg = cs * cs * (abs(diff)-0.3)
                    countNeg = countNeg + volNeg
                    countNegArea = countNegArea + (cs * cs)  
                    countAcc = countAcc + (cs * cs * 0.3)
                if diff > 0.3 and diff < 100:
                    volPos = cs * cs * (diff-0.3)
                    countPos = countPos + (volPos)
                    countPosArea = countPosArea + (cs * cs) 
                    countAcc = countAcc + (cs * cs * 0.3)
                if diff < 0 and diff > -0.3:
                    countCell = countCell + 1
                    countCellVol = countCellVol + (cs * cs * abs(diff))
                if diff > 0 and diff < 0.3:
                    countCell = countCell + 1
                    countCellVol = countCellVol + (cs * cs * diff)

        analyse.write("Analysis with Threshold of +/- 0.3"+"\n")
        analyse.write("Thresholded (0.3) Area [m2] of Detectable Change: " + str(countPosArea+countNegArea)+"\n")
        analyse.write("Thresholded (0.3) Area [km2] of Detectable Change: " + str((countPosArea+countNegArea)/1000000)+"\n")
        analyse.write("Thresholded (0.3) of Interest [%] with Detectable Change: " + str(((countPosArea+countNegArea)/countAOI)*100)+"\n")
        analyse.write("Thresholded (0.3) Area [m2] of Surface Lowering: " + str(countNegArea)+"\n")
        analyse.write("Thresholded (0.3) Area [km2] of Surface Lowering: " + str(countNegArea/1000000)+"\n")
        analyse.write("Thresholded (0.3) Area [m2] of Surface Raising: " + str(countPosArea)+"\n")
        analyse.write("Thresholded (0.3) Area [km2] of Surface Raising: " + str(countPosArea/1000000)+"\n")
        analyse.write("Thresholded (0.3) Volume [m3] of Surface Lowering: " + str(countNeg)+"\n")
        analyse.write("Thresholded (0.3) Volume [m3] of Surface Raising: " + str(countPos)+"\n")
        analyse.write("Thresholded (0.3) Volume [m3] of Difference: " + str(countPos+countNeg)+"\n")
        analyse.write("Thresholded (0.3) Net Volume [m3] of Difference: " + str(countPos-countNeg)+"\n")
        analyse.write("\n") 
        analyse.write("Volume [m3] of Error within Threshold of -0.3 and 0.3: " + str(countAcc)+"\n")
        analyse.write("Count of Cells within -0.3 and 0.3: " + str(countCell)+"\n")
        analyse.write("Volume [m3] of Cells between -0.3 and 0.3: " + str(countCellVol)+"\n")
        analyse.write("Percent [%] of Error of Cells between -0.3 and 0.3: " + str((countCellVol/(countPosGes+countNegGes))*100)+"\n")
        analyse.write("\n")
        analyse.write("\n") 
        
        analyse_csv.write("0.3; "+str(countAOI)+"; "+str(countPosArea+countNegArea)+"; "+str(((countPosArea+countNegArea)/countAOI)*100)+"; "+str(countNegArea)+"; "+str(countPosArea)+"; "+str(countNeg)+"; "+str(countPos)+"; "+str(countPos+countNeg)+"; "+str(countPos-countNeg)+"\n")
        net_v_2 = countPos-countNeg
        sr_2 = countPos
        sl_2 = countNeg

        del countPos
        del countNeg
        del countPosArea
        del countNegArea
        del countAcc
        del countCell
        del countCellVol


        ## Classification3 - UC
        countAOI=0
        countPosGes = 0
        countNegGes = 0
        countPos = 0
        countNeg = 0
        countPosArea = 0
        countNegArea = 0
        countES = 0
        countESneg = 0
        countESpos = 0
        countAcc = 0
        countCell = 0
        countCellVol = 0
        for row in range(0, width):
            for col in range(0, height):
                diff = pix_diff[row, col]
                ES = pix_uc[row, col]
                ES2 = abs(ES)
                if diff > -9999.0 and diff < 100:
                    countAOI = countAOI+(cs*cs)
                if diff < -ES2 and diff > -100.0:
                    grid_diff_uc[row,col]=diff+ES2
                    volESneg = cs * cs * (abs(diff)-ES2)
                    countESneg = countESneg + volESneg
                    countNegArea = countNegArea + (cs * cs)
                    countAcc = countAcc + (cs * cs * ES2)
                if diff > ES2 and diff < 100:
                    grid_diff_uc[row,col]=diff-ES2
                    volESpos = cs * cs * (diff-ES2)
                    countESpos = countESpos + volESpos
                    countPosArea = countPosArea + (cs * cs) 
                    countAcc = countAcc + (cs * cs * 0.1)
                if diff < 0 and diff > -ES2:
                    grid_diff_uc[row,col]=0.0
                    countCell = countCell + 1
                    countCellVol = countCellVol + (cs * cs * abs(diff))
                    #print countCell
                if diff > 0 and diff < ES2:
                    grid_diff_uc[row,col]=0.0
                    countCell = countCell + 1
                    countCellVol = countCellVol + (cs * cs * abs(diff))


        analyse.write("Analysis including Uncertainty Analysis"+"\n")
        analyse.write("Thresholded Area [m2] of Detectable Change: " + str(countPosArea+countNegArea)+"\n")
        analyse.write("Thresholded Area of Interest [%] with Detectable Change: " + str(((countPosArea+countNegArea)/countAOI)*100)+"\n")
        analyse.write("Thresholded Volume [m3] of Surface Lowering: " + str(countESneg)+"\n")
        analyse.write("Thresholded Volume [m3] of Surface Raising: " + str(countESpos)+"\n")
        analyse.write("Total Volume [m3] of Difference: " + str(countESpos+countESneg)+"\n")
        analyse.write("Total Net Volume [m3] of Difference: " + str(countESpos-countESneg)+"\n")
        analyse.write("\n")
        analyse.write("\n") 
        analyse.write("Volume [m3] of Error within Threshold of UC: " + str(countAcc)+"\n")
        analyse.write("Count of Cells within UC: " + str(countCell)+"\n")
        analyse.write("Volume [m3] of Cells between UC: " + str(countCellVol)+"\n")
        analyse.write("Percent [%] of Error of Cells between UC: " + str((countCellVol/(countESpos+countESneg))*100)+"\n")
        analyse.write("\n")
        analyse.write("\n")
       
        analyse_csv.write("UC; "+str(countAOI)+"; "+str(countPosArea+countNegArea)+"; "+str(((countPosArea+countNegArea)/countAOI)*100)+"; ; ; "+str(countESneg)+"; "+str(countESpos)+"; "+str(countESpos+countESneg)+"; "+str(countESpos-countESneg)+"\n")
        net_v_3 = countESpos-countESneg
        sr_3 = countESpos
        sl_3 = countESneg
        
        ## BarChart
        plotdata = pd.DataFrame({
            'surface raising': [sr_0, sr_1, sr_2, sr_3],
            'surface lowering': [-sl_0, -sl_1, -sl_2, -sl_3]
            }, 
            index=['raw','0.1 m','0.3 m','UC model']) 
        New_Colors = ['silver','dimgrey']
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Verdana']
        plotdata.plot(kind="bar", stacked=True, width = 0.5, color=New_Colors,  align='center', widtH=0.5)
        #plt.xticks(np.arange(0, 4, step=1))
        plt.xticks(rotation=5, horizontalalignment="center")
        #plt.yticks(rotation=0, horizontalalignment="center")
        plt.title('volume changes [m3]', fontsize=10)
        plt.xlabel('threshold', fontsize=10)
        #plt.ylabel('volume changes [m3]', fontsize=10)
        #plt.ylabel('net volume [m3]', fontsize=10)
        plt.grid(True)

        plt.savefig(final_path + '/' + 'change_detection_plot.png')

        ## Grid UC
        out3 = temp_path + '/' + 'diff_uc_old.tif'
        out4 = final_path + '/' + 'diff_uc.tif'
        grid_diff_uc=np.flip(grid_diff_uc,1)
        grid_diff_uc=np.rot90(grid_diff_uc)
        imsave = Image.fromarray(grid_diff_uc, mode='F')
        imsave.save(out3, "TIFF")
        
        ## Raster georeferenzieren       
        src_ds = gdal.Open(out3)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
                
        dst_ds = driver.CreateCopy(out4, src_ds, 0)
        dst_ds.SetGeoTransform(gt)
        epsg = 32633
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
        dest_wkt = srs.ExportToWkt()
        dst_ds.SetProjection(dest_wkt)

        # Close files
        dst_ds = None
        src_ds = None

        ## Add Layer
        root = QgsProject.instance().layerTreeRoot()
        mygroup = root.findGroup("CD_Tool")
        rlayer = QgsRasterLayer(out4, 'Raster')
        rprovider = rlayer.dataProvider()
        colDic = {'f':'#d7191c','f1':'#eb6640','f2': '#feb165', 'f3': '#ffdc96', 'f4':'#ffffff' , 'f5':'#ffffff', 'f6':'#9cd3a7', 'f7':'#5ea7b1','f8':'#2b83ba', 'f9':'#1e5c83' }
        valueList = [-3, -2, -1, -0.5, -0.05, 0.05, 0.5, 1, 2, 3]
        lst = [QgsColorRampShader.ColorRampItem(valueList[0], QColor(colDic['f'])),
            QgsColorRampShader.ColorRampItem(valueList[1], QColor(colDic['f1'])),
            QgsColorRampShader.ColorRampItem(valueList[2], QColor(colDic['f2'])),
            QgsColorRampShader.ColorRampItem(valueList[3], QColor(colDic['f3'])),
            QgsColorRampShader.ColorRampItem(valueList[4], QColor(colDic['f4'])),
            QgsColorRampShader.ColorRampItem(valueList[5], QColor(colDic['f5'])),
            QgsColorRampShader.ColorRampItem(valueList[6], QColor(colDic['f6'])),
            QgsColorRampShader.ColorRampItem(valueList[7], QColor(colDic['f7'])),
            QgsColorRampShader.ColorRampItem(valueList[8], QColor(colDic['f8'])),
            QgsColorRampShader.ColorRampItem(valueList[9], QColor(colDic['f9']))]
        my_shader = QgsRasterShader()
        my_colorramp = QgsColorRampShader()
        #fcn = QgsColorRampShader()
        #fcn.setColorRampType(QgsColorRampShader.Interpolated)
        #lst = [ QgsColorRampShader.ColorRampItem(0, QColor(0,255,0)),QgsColorRampShader.ColorRampItem(255, QColor(255,255,0)) ]
        my_colorramp.setColorRampItemList(lst)
        my_colorramp.setColorRampType(QgsColorRampShader.Interpolated)
        my_shader.setRasterShaderFunction(my_colorramp)
        renderer = QgsSingleBandPseudoColorRenderer(rlayer.dataProvider(), 1, my_shader)
        rasterTransparency = QgsRasterTransparency()
        myTransparentSingleValuePixelList = []
        myTransparentPixel1 = QgsRasterTransparency.TransparentSingleValuePixel()
        myTransparentPixel1.min = -0.05
        myTransparentPixel1.max = 0.05
        myTransparentPixel1.percentTransparent = 100
        myTransparentSingleValuePixelList.append(myTransparentPixel1)
        rasterTransparency.setTransparentSingleValuePixelList(myTransparentSingleValuePixelList)
        renderer.setRasterTransparency(rasterTransparency)
        rlayer.setRenderer(renderer)
        rlayer.triggerRepaint()
        QgsProject.instance().addMapLayer(rlayer)
        mygroup.addLayer(rlayer)

        outputs['LastStep'] = out4
        results['Tool abgeschlossen'] = outputs['LastStep']
        return results



