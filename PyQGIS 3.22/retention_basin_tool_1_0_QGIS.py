# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   Ermittlung der Kubatur von Rückhaltebecken                            *
*   pyAlpineRisk                                                          *
*   QGIS 3.22                                                             *
*   Nicole Kamp & Franz Langegger                                         *
*   niki.kamp@gmail.com                                                   *
*   Mai 2020                                                              *
*                                                                         *
***************************************************************************
"""

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProject,
                       QgsVectorLayer,
                       QgsTextFormat,
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

from PIL import Image
Image.MAX_IMAGE_PIXELS = None
import osgeo.gdal as gdal
import numpy as np

from osgeo import ogr, osr
from shapely.geometry import Polygon
import matplotlib
import subprocess
from subprocess import call

import string, os, sys, copy, shutil, math, numpy, time, datetime
from time import *
from sys import *


#class ProjectUTM33N

#class AddSurfaceInfo

#class AddFields

#class LoopHeight


class KubaturProcessingAlgorithm(QgsProcessingAlgorithm):
    """
Python Script zur Berechnung des Volumens des Geschieberückhaltebeckens.
Eingabe der Input- und Output-Parameter
    """
    INPUT_Shape = 'INPUT_SHP'
    INPUT_ALS = 'DGM'
    INPUT_hoehe = 'HOEHE'
    INPUT_prozent = 'PROZENT'
    #INPUT_extent = 'AUSSCHNITT'
    TEMP = 'TEMP'
    
    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return KubaturProcessingAlgorithm()

    def name(self):
        return 'kubaturtool'

    def displayName(self):
        return self.tr('KubaturTool_1_0')

    def group(self):
        return self.tr('pyAlpineRisk')

    def groupId(self):
        return 'pyAlpineRiskScripts'

    def shortHelpString(self):
        return self.tr("Ermittlung des Rückhaltevolumens\nDas Kubatur Tool ermöglicht, ausgehend von einer Abfragelinie, die Ermittlung der Kubatur innerhalb einer generierten Verschneidungsfläche. Aufgrund der Abfragelinie wird eine künstliche Oberfläche erzeugt, die sich aus der Lage und Höhe sowie dem Verlandungsgefälle ergibt und mit dem Geländemodell (ALS 1m) verschnitten wird. Durch eine Rasteranalyse wird das Differenzvolumen ermittelt und die Verschneidung als Flächenpolygon bereitgestellt. Die temporären Arbeitsdaten sowie die Ergebnisse werden getrennt nach Datum und Uhrzeit im vordefinierten Ausgabeverzeichnis „C:/temp/KT/“ gespeichert.\nArbeitsschritte:\n1. Abfragelinie erstellen (Neues Linienobjekt - keine BWK Linien verwenden!) \n2. Kubatur Tool öffnen (Bedienfeld – Verarbeitungswerkzeuge – Skripte)\n3. Sperrenhöhe und Verlandungsgefälle festlegen\n4. Kubatur Tool starten\n5. Das Ergebnis (Verschneidungspolygon und Differenzvolumen) wird automatisch zur aktuellen Ansicht hinzugefügt")

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_Shape,
                self.tr('Abfragelinie'),
                [QgsProcessing.TypeVectorAnyGeometry]
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_ALS,
                self.tr('Geländemodell'),
                #extension='tiff',
                #fileFilter="tiff (*.tif)",
                defaultValue='/DGM/ALS_DGM_1m.tif'
            )
        )

        self.addParameter(QgsProcessingParameterNumber(
            self.INPUT_hoehe, 
            self.tr('Sperrenhoehe'),
            QgsProcessingParameterNumber.Double,
            3.0
            )
        )
        
        self.addParameter(QgsProcessingParameterNumber(
            self.INPUT_prozent, 
            self.tr('Verlandungsgefaelle'),
            QgsProcessingParameterNumber.Double,
            3.0
            )
        )
              
        self.addParameter(
            QgsProcessingParameterFolderDestination(
                self.TEMP, 
                self.tr('Output-Ordner'), 
                defaultValue='C:/temp/KT'
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        feedback = QgsProcessingMultiStepFeedback(1, feedback)
        results = {}
        outputs = {}
        
        laenge = 1000
        laenge_k = 10
        
        # Pfade definieren + Timestamp
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        temp_path = str(parameters[self.TEMP])+'/temp_'+str(timestamp)
        final_path = str(parameters[self.TEMP])+'/final_'+str(timestamp)
        
        # Alle Daten im Folder löschen
        #for root, dirs, files in os.walk(final_path):
        #    for f in files:
        #        os.unlink(os.path.join(root, f))
        #    for d in dirs:
        #        shutil.rmtree(os.path.join(root, d))
        
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
        if not os.path.exists(final_path):
            os.makedirs(final_path)
        
        #feedback.pushInfo(str(parameters[self.INPUT_extent]))
        
        ## Buffer
        buffer = out = temp_path + '/' + 'b1000.shp'
        alg_params = {
            'DISSOLVE': False,
            'DISTANCE': 700,
            'END_CAP_STYLE': 2,
            'INPUT': str(parameters[self.INPUT_Shape]),
            'JOIN_STYLE': 1,
            'MITER_LIMIT': 2,
            'SEGMENTS': 1,
            'OUTPUT': buffer
        }
        processing.run('native:buffer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        # Raster auf Layermaske zuschneiden
        out = temp_path + '/' + 'clip_vrt.tif'
        alg_params = {
            'ALPHA_BAND': False,
            'CROP_TO_CUTLINE': True,
            'DATA_TYPE': 0,
            'EXTRA': '',
            'INPUT': str(parameters[self.INPUT_ALS]),
            'KEEP_RESOLUTION': False,
            'MASK': str(buffer),
            'MULTITHREADING': False,
            'NODATA': None,
            'OPTIONS': '',
            'SET_RESOLUTION': True,
            'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:2056'),
            'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:2056'),
            'X_RESOLUTION': 1,
            'Y_RESOLUTION': 1,
            'OUTPUT': out
        }
        processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        
        ## Cell Size
        raster_cs = gdal.Open(out)
        gt_cs =raster_cs.GetGeoTransform() 
        cs = gt_cs[1]
        
        del raster_cs
        del gt_cs
        
        ## 1. Profil-Punkte aus einer Linie generieren
        out1 = temp_path + '/' + 'out1.shp'
        
        alg_params = {
            'DEM': str(out),
            'LINES': str(parameters[self.INPUT_Shape]),
            'NAME': 'profil',
            'SPLIT         ': True,
            'VALUES': [],
            'PROFILE': out1,
            'PROFILES': out1
        }
        processing.run('saga:profilesfromlines', alg_params, context=context, feedback=feedback, is_child_algorithm=True)


        ## 2. Liste mit Hoehenwerten aus den Profil-Punkten
        listeZ =[]
        ds=ogr.Open(out1)
        lyr=ds.GetLayer()
        for feat in lyr:
            zVal = feat.GetField("Z")
            #feedback.pushInfo(str(zVal))
            listeZ.append(zVal)
        
        del out1


        ## 3. Hoehe der Sperre, Abfragehoehe, Hoehe100 und ZMin berechnen
        zmin = min(listeZ)
        objhoehe = zmin + parameters[self.INPUT_hoehe]
        abfragehoehe = objhoehe + (laenge * parameters[self.INPUT_prozent]/100)
        hoehe100 = zmin + 100
        feedback.pushInfo(str(objhoehe))
        feedback.pushInfo(str(abfragehoehe))
        feedback.pushInfo(str(hoehe100))
        feedback.pushInfo(str(zmin))
  
         # Felder berechnen
        outputs = {}
        out2 = temp_path + '/' + 'out2.shp'
        out3 = temp_path + '/' + 'out3.shp'
        out4 = temp_path + '/' + 'out4.shp'
        out5 = final_path + '/' + 'abfragelinie.shp'
        
        # SperreH
        alg_params = {
            'FIELD_LENGTH': 10,
            'FIELD_NAME': 'SperreH',
            'FIELD_PRECISION': 3,
            'FIELD_TYPE': 1,
            'FORMULA': 'value = ' + str(objhoehe),
            'GLOBAL': '',
            'INPUT': str(parameters[self.INPUT_Shape]),
            'OUTPUT': out2
        }
        processing.run('qgis:advancedpythonfieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        # AbfrageH
        alg_params = {
            'FIELD_LENGTH': 10,
            'FIELD_NAME': 'AbfrageH',
            'FIELD_PRECISION': 3,
            'FIELD_TYPE': 1,
            'FORMULA': 'value = ' + str(abfragehoehe),
            'GLOBAL': '',
            'INPUT': str(out2),
            'OUTPUT': out3
        }
        processing.run('qgis:advancedpythonfieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        # Hoehe100
        alg_params = {
            'FIELD_LENGTH': 10,
            'FIELD_NAME': 'H100',
            'FIELD_PRECISION': 3,
            'FIELD_TYPE': 1,
            'FORMULA': 'value = ' + str(hoehe100),
            'GLOBAL': '',
            'INPUT': str(out3),
            'OUTPUT': out4
        }
        processing.run('qgis:advancedpythonfieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
               
        # ZMIN
        alg_params = {
            'FIELD_LENGTH': 10,
            'FIELD_NAME': 'ZMIN',
            'FIELD_PRECISION': 3,
            'FIELD_TYPE': 1,
            'FORMULA': 'value = ' + str(zmin),
            'GLOBAL': '',
            'INPUT': str(out4),
            'OUTPUT': out5
        }
        processing.run('qgis:advancedpythonfieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        feedback.pushInfo(str("si"))

        del out2, out3, out5
        feedback.pushInfo(str("ok"))
      
      
        ## 4. Rechteck + Linien berechnen
        #outputs = {}
        
        # Koordinaten berechnen
        ds=ogr.Open(out4)
        lyr=ds.GetLayer()
        pts_list=[]
        for i in range(lyr.GetFeatureCount()):
            feat=lyr.GetFeature(i)
            geom=feat.GetGeometryRef()
            ptcount = geom.GetPointCount()
            #feedback.pushInfo(str(geom.GetPointCount()))
            firstpoint=geom.GetPoint(0)
            lastpoint=geom.GetPoint(geom.GetPointCount()-1)
            x1 = firstpoint[0]
            y1 = firstpoint[1]
            x2 = lastpoint[0]
            y2 = lastpoint[1]
            
            for n in range(0, ptcount):
                point=geom.GetPoint(n)
                pts_list.append(point)
        feedback.pushInfo(str(pts_list))
            
        y_neu = y2 - y1
        x_neu = x2 - x1
        phi = math.radians(90)
        r1 = math.atan2 (y_neu, x_neu)
        r2 = r1 + phi
        r3 = r1 - phi
        x3 = x2 + laenge * math.cos(r2) 
        y3 = y2 + laenge * math.sin(r2) 
        x4 = x1 + laenge * math.cos(r2)
        y4 = y1 + laenge * math.sin(r2)
        x5 = x2 + laenge * math.cos(r3)
        y5 = y2 + laenge * math.sin(r3)
        x6 = x1 + laenge * math.cos(r3)
        y6 = y1 + laenge * math.sin(r3)
        
        # Vergleichslinie
        x3k = x2 + laenge_k * math.cos(r2) 
        y3k = y2 + laenge_k * math.sin(r2) 
        x4k = x1 + laenge_k * math.cos(r2)
        y4k = y1 + laenge_k * math.sin(r2)
        x5k = x2 + laenge_k * math.cos(r3)
        y5k = y2 + laenge_k * math.sin(r3)
        x6k = x1 + laenge_k * math.cos(r3)
        y6k = y1 + laenge_k * math.sin(r3)
        
        z12 = float(objhoehe)
        z3456 = float(abfragehoehe)
        extent=str(x3)+","+str(x6)+","+str(y4)+","+str(y5) 
        feedback.pushInfo(str(extent))
               
        #del out4
 
        # Output-Pfade definieren
        out6 = temp_path + '/' + 'out6.shp'
        out7 = temp_path + '/' + 'out7.shp'
        outL1 = temp_path + '/' + 'line1.shp'
        outL2 = temp_path + '/' + 'line2.shp'
        outL1_neu = temp_path + '/' + 'line1_pkt.shp'
        outL2_neu = temp_path + '/' + 'line2_pkt.shp'
        
        # Line1 anlegen
        driver = ogr.GetDriverByName('Esri Shapefile')
        ds = driver.CreateDataSource(outL1)
        layer = ds.CreateLayer('', None, ogr.wkbLineString)
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        defn = layer.GetLayerDefn()
        feat = ogr.Feature(defn)
        feat.SetField('id', 1)
        geom = ogr.Geometry(ogr.wkbLineString)
        geom.AddPoint(x4k, y4k)
        geom.AddPoint(x3k, y3k)
        feat.SetGeometry(geom)
        layer.CreateFeature(feat)
        feat = geom = None 
        ds = layer = feat = geom = None 
               
        # Profil-Punkte aus Linie1 generieren
        alg_params = {
            'DEM': str(out),
            'LINES': str(outL1),
            'NAME': 'profil',
            'SPLIT         ': True,
            'VALUES': [],
            'PROFILE': outL1_neu,
            'PROFILES': outL1_neu
        }
        processing.run('saga:profilesfromlines', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        # Liste mit Hoehenwerten aus den Profil-Punkten
        listeZ_L1 =[]
        ds=ogr.Open(outL1_neu)
        lyr=ds.GetLayer()
        for feat in lyr:
            zVal = feat.GetField("Z")
            listeZ_L1.append(zVal) 

        # Line2 anlegen
        driver = ogr.GetDriverByName('Esri Shapefile')
        ds = driver.CreateDataSource(outL2)
        layer = ds.CreateLayer('', None, ogr.wkbLineString)
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        defn = layer.GetLayerDefn()
        feat = ogr.Feature(defn)
        feat.SetField('id', 1)
        geom = ogr.Geometry(ogr.wkbLineString)
        geom.AddPoint(x5k, y5k)
        geom.AddPoint(x6k, y6k)
        feat.SetGeometry(geom)
        layer.CreateFeature(feat)
        feat = geom = None 
        ds = layer = feat = geom = None
        
        # Profil-Punkte aus Linie2 generieren
        alg_params = {
            'DEM': str(out),
            'LINES': str(outL2),
            'NAME': 'profil',
            'SPLIT         ': True,
            'VALUES': [],
            'PROFILE': outL2_neu,
            'PROFILES': outL2_neu
        }
        processing.run('saga:profilesfromlines', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        # Liste mit Hoehenwerten aus den Profil-Punkten
        listeZ_L2 =[]
        ds=ogr.Open(outL2_neu)
        lyr=ds.GetLayer()
        for feat in lyr:
            zVal = feat.GetField("Z")
            #feedback.pushInfo(str(zVal))
            listeZ_L2.append(zVal) 
        
        
        # Polygon anlegen - Check, welche Linie höher liegt + Koordinatenliste
        driver = ogr.GetDriverByName('Esri Shapefile')
        ds = driver.CreateDataSource(out6)
        layer = ds.CreateLayer('', None, ogr.wkbPolygon25D) #wkbPolygon
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        defn = layer.GetLayerDefn()
        feat = ogr.Feature(defn)
        feat.SetField('id', 1)
        
        pts_length = len(pts_list)
        feedback.pushInfo(str(pts_length))
        zmin1 = min(listeZ_L1)
        zmin2 = min(listeZ_L2)
        if pts_length == 2:
            if zmin1 > zmin2:
                poly1 = Polygon([(x1,y1,z12), (x4,y4,z3456), (x3,y3,z3456), (x2,y2,z12), (x1,y1,z12)])
            else:
                poly1 = Polygon([(x1,y1,z12), (x2,y2,z12), (x5,y5,z3456), (x6,y6,z3456), (x1,y1,z12)])
            
                    # Rechteck1 anlegen
                    #driver = ogr.GetDriverByName('Esri Shapefile')
                    #ds = driver.CreateDataSource(out6)
                    #layer = ds.CreateLayer('', None, ogr.wkbPolygon25D) #wkbPolygon
                    #layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
                    #defn = layer.GetLayerDefn()
                    #feat = ogr.Feature(defn)
                    #feat.SetField('id', 1)
            geom = ogr.CreateGeometryFromWkb(poly1.wkb)
            feat.SetGeometry(geom)
            layer.CreateFeature(feat)
            feat = geom = None 
            ds = layer = feat = geom = None 
        
        if pts_length > 2:
            coords0 = pts_list[0]
            x0 = coords0[0]
            y0 = coords0[1]
            z0=z12
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for n in range (0,pts_length-1):
                coords = pts_list[n]
                xn = coords[0]
                yn = coords[1]
                zn=z12
                ring.AddPoint(xn, yn, zn)
            if zmin1 > zmin2:
                ring.AddPoint(x3,y3,z3456)
                ring.AddPoint(x4,y4,z3456)
                ring.AddPoint(x0, y0, z0)
            else:
                ring.AddPoint(x5,y5,z3456)
                ring.AddPoint(x6,y6,z3456)
                ring.AddPoint(x0, y0, z0)
            poly = ogr.Geometry(ogr.wkbPolygon25D)
            poly.AddGeometry(ring)
            poly1 = poly.ExportToWkt()
            
            # Rechteck1 anlegen
            geom = ogr.CreateGeometryFromWkt(poly1)
            #geom = ogr.CreateGeometryFromWkb(poly1.wkb)
            feat.SetGeometry(geom)
            layer.CreateFeature(feat)
            feat = geom = None 
            ds = layer = feat = geom = None 
        
        feedback.pushInfo("TEST")
        feedback.pushInfo(str(poly1))
                

        # Rechteck1 anlegen
        #driver = ogr.GetDriverByName('Esri Shapefile')
        #ds = driver.CreateDataSource(out6)
        #layer = ds.CreateLayer('', None, ogr.wkbPolygon25D) #wkbPolygon
        #layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        #defn = layer.GetLayerDefn()
        #feat = ogr.Feature(defn)
        #feat.SetField('id', 1)
        #geom = ogr.CreateGeometryFromWkb(poly1.wkb)
        #feat.SetGeometry(geom)
        #layer.CreateFeature(feat)
        #feat = geom = None 
        #ds = layer = feat = geom = None 
        
        # AbfrageHoehe als Attribut anhaengen
        alg_params= {
            'FIELD_LENGTH': 10,
            'FIELD_NAME': 'AbfrageH',
            'FIELD_PRECISION': 3,
            'FIELD_TYPE': 1,
            'FORMULA': 'value = ' + str(abfragehoehe),
            'GLOBAL': '',
            'INPUT': str(out6),
            'OUTPUT': out7
        }
        processing.run('qgis:advancedpythonfieldcalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        del out6, outL1, outL2, outL1_neu, outL2_neu
 
 
        # 5. Aus 3D Rechteck Raster berechnen
        out8 = temp_path + '/' + 'out8.tif'
 
        alg_params = {
            'EXTENT': str(extent),
            'INTERPOLATION_DATA': str(out7)+'::~::1::~::-1::~::2',
            'METHOD': 0,
            'PIXEL_SIZE': int(cs),
            'OUTPUT': out8
        }
        processing.run('qgis:tininterpolation', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        
        ## 6. Differenzmodell und Stauraum berechnen       
        # Raster mit Rechteck zuschneiden und Auflösung ändern
        out9 = temp_path + '/' + 'out9.tif'
        out10 = temp_path + '/' + 'out10.tif'
        out11 = temp_path + '/' + 'out11.tif'
               
        # DGM klippen
        OutTile = gdal.Warp(str(out10), 
            str(out), 
            cutlineDSName=str(out7),
            cropToCutline=True,
            dstNodata = 0)
        OutTile = None
                        
        # Geklipptes DGM öffnen
        dtm_clip_gdal = gdal.Open(out10)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
        band1 = dtm_clip_gdal.GetRasterBand(1)
        im = Image.open(out10)
        pix = im.load()
        gt = dtm_clip_gdal.GetGeoTransform()
        feedback.pushInfo(str(gt))
        width = dtm_clip_gdal.RasterXSize
        height = dtm_clip_gdal.RasterYSize
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5] 
        maxx = gt[0] + width*gt[1] + height*gt[2]
        maxy = gt[3]
        extent1=str(minx)+","+str(maxx)+","+str(miny)+","+str(maxy)
        extent2=(minx,maxx,miny,maxy)
        feedback.pushInfo(str(extent1))

        # Snap Raster
        flaeche3d = gdal.Open(out8)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
        OutTile = gdal.Warp(out9, flaeche3d,
            format=format, outputBounds=[minx, miny, maxx, maxy], 
            xRes=int(cs), yRes=int(cs))
        OutTile = None
       
        alg_params = {'INPUT_A' : str(out10),
            'BAND_A' : 1,
            'INPUT_B' : str(out9),
            'BAND_B' : 1,
            'FORMULA' : '(A - B)',  
            'OUTPUT' : out11
        }
        processing.run('gdal:rastercalculator', alg_params, context=context, feedback=feedback, is_child_algorithm=True) 
        
        del out7, out8, out9, out10
       
        # Raster to Array
        kubatur = gdal.Open(out11)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
        band1 = kubatur.GetRasterBand(1)
        im = Image.open(out11)
        pix = im.load()
        gt = kubatur.GetGeoTransform()
        width = kubatur.RasterXSize
        height = kubatur.RasterYSize
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5] 
        maxx = gt[0] + width*gt[1] + height*gt[2]
        maxy = gt[3]
        extent=str(minx)+","+str(maxx)+","+str(miny)+","+str(maxy) 
        
        # Leeres Grid anlegen
        grid1 = np.zeros(shape=(width,height), dtype=np.float32)
        grid2 = np.zeros(shape=(width,height), dtype=np.float32)
        
        kubatur_gesamt = 0
        feedback.pushInfo('ok')
        
        for row in range(0, width):
            for col in range(0, height):
                val1 =pix[row, col]
                if val1 < 0:
                    nval = (val1*-1)*int(cs)*int(cs)
                    nval2 = 1.0
                    kubatur_gesamt = round(kubatur_gesamt+nval, 1)
                    grid1[row,col]=nval
                    grid2[row,col]=nval2
                else:
                    nval = np.nan
                    nval2 = 2.0 #np.nan
                    grid1[row,col]=nval
                    grid2[row,col]=nval2

        #feedback.pushInfo(str(grid1.shape))
        #feedback.pushInfo(str(grid1))
        
        # Array to Raster
        rasterinfo = gdal.Open(out11)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
        band1 = rasterinfo.GetRasterBand(1)
        gt = rasterinfo.GetGeoTransform()
        feedback.pushInfo(str(gt))
        width = rasterinfo.RasterXSize
        height = rasterinfo.RasterYSize
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5] 
        maxx = gt[0] + width*gt[1] + height*gt[2]
        maxy = gt[3] 
        
        # Beide Raster abspeichern
        grid1=np.flip(grid1,1)
        grid1=np.rot90(grid1)
        out12 = temp_path + '/' + 'out12.tif'
        imsave = Image.fromarray(grid1, mode='F')
        imsave.save(out12, "TIFF")
        
        grid2=np.flip(grid2,1)
        grid2=np.rot90(grid2)
        out13 = temp_path + '/' + 'out13.tif'
        imsave = Image.fromarray(grid2, mode='F')
        imsave.save(out13, "TIFF")
        
        del imsave, out11
        
        # Raster1 georeferenzieren
        out12 = temp_path + '/' + 'out12.tif'
        out14 = final_path + '/' + 'kubatur.tif'
        src_ds = gdal.Open(out12)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
                
        dst_ds = driver.CreateCopy(out14, src_ds, 0)
        dst_ds.SetGeoTransform(gt)
        epsg = 2056
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
        dest_wkt = srs.ExportToWkt()
        dst_ds.SetProjection(dest_wkt)

        # Close files
        dst_ds = None
        src_ds = None
        
         # Raster2 georeferenzieren
        out13 = temp_path + '/' + 'out13.tif'
        out15 = temp_path + '/' + 'mask.tif'
        src_ds = gdal.Open(out13)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
                
        dst_ds = driver.CreateCopy(out15, src_ds, 0)
        dst_ds.SetGeoTransform(gt)
        epsg = 2056
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
        dest_wkt = srs.ExportToWkt()
        dst_ds.SetProjection(dest_wkt)

        # Close files
        dst_ds = None
        src_ds = None       
        
        feedback.pushInfo('JUHU')
        
        ## Raster in Polygon umwandeln
        out15 = temp_path + '/' + 'mask.tif'
        out16 = temp_path + '/' + 'stauraum_alt.shp'
        src_ds = gdal.Open(out15)
        srcband = src_ds.GetRasterBand(1)
        driver = ogr.GetDriverByName('Esri Shapefile')
        ds = driver.CreateDataSource(out16)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(2056)
        dst_ds = ds.CreateLayer('', srs=srs)
        newField = ogr.FieldDefn('id', ogr.OFTInteger)
        #newField.SetPrecision(1)
        dst_ds.CreateField(newField)
        newField2 = ogr.FieldDefn('kubatur', ogr.OFTReal)
        dst_ds.CreateField(newField2)
        gdal.Polygonize(srcband, None, dst_ds, 0, [], callback=None)
        feat = geom = None 
        ds = dst_ds = feat = geom = None
               
        # Dissolve
        out17 = temp_path + '/' + 'stauraum.shp'
        alg_params = {
            'COMPUTE_AREA': True,
            'COMPUTE_STATISTICS': False,
            'COUNT_FEATURES': False,
            'EXPLODE_COLLECTIONS': False,
            'FIELD': 'id',
            'GEOMETRY': 'geometry',
            'INPUT': str(out16),
            'KEEP_ATTRIBUTES': True,
            'OPTIONS': '',
            'STATISTICS_ATTRIBUTE': None,
            'OUTPUT': out17
        }
        processing.run('gdal:dissolve', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
      
        ds=ogr.Open(out17, 1)
        layer=ds.GetLayer()
        for i in range(layer.GetFeatureCount()):
            feat=layer.GetFeature(i)
            geom=feat.GetGeometryRef()
            area = geom.GetArea()
            volume = kubatur_gesamt
            feat.SetField('kubatur', volume)
            layer.SetFeature(feat)
        
        ds = layer = feat = geom = None
        del i
            
        ds=ogr.Open(out17, 1)
        layer=ds.GetLayer()
        for feat in layer:
            fid = feat.GetFID()
            if feat.GetField('id') == 2:
                #feedback.pushInfo(str(fid))
                layer.DeleteFeature(feat.GetFID())
        ds.ExecuteSQL('REPACK lowercase')
 
        ds = layer = feat = geom = None
        del out12, out13, out14, out15, out16
        
        ## Temporäre Daten löschen
        #for root, dirs, files in os.walk(temp_path):
        #    for f in files:
        #        os.unlink(os.path.join(root, f))
        #    for d in dirs:
        #        shutil.rmtree(os.path.join(root, d))

        # Smoothing
        out18 =temp_path + '/' + 'stauraum_geglaettet.shp'
        alg_params = {
            'INPUT': str(out17),
            'ITERATIONS': 10,
            'MAX_ANGLE': 100,
            'OFFSET': 0.4,
            'OUTPUT': out18
        }
        processing.run('native:smoothgeometry', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        feedback.pushInfo('Stauraum: ' + str(kubatur_gesamt) + ' m³')
        
        # Vereinfachen
        out19 = final_path + '/' + 'stauraum.shp'
        alg_params = {
            'INPUT': str(out18),
            'METHOD': 1,
            'TOLERANCE': 1,
            'OUTPUT': out19
        }
        outputs['Vereinfachen'] = processing.run('native:simplifygeometries', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        
        ## Layer hinzufügen
        root = QgsProject.instance().layerTreeRoot()
        mygroup = root.findGroup("Kubatur-Tool - Ergebnisse")
        layer = QgsVectorLayer(out19, "Stauraum", "ogr")
        text_format = QgsTextFormat()
        label = QgsPalLayerSettings()
        label.fieldName = 'kubatur'
        label.enabled = True
        label.setFormat(text_format)
        labeler = QgsVectorLayerSimpleLabeling(label)
        layer.setLabelsEnabled(True)
        layer.setLabeling(labeler)
        #layer.triggerRepaint()
        layer.renderer().symbol().setColor(QColor("blue"))
        layer.triggerRepaint()
        #QgsProject.instance().removeMapLayers(layer)
        QgsProject.instance().addMapLayer(layer, False)
        mygroup.addLayer(layer)

        
        results['KubaturTool'] = outputs['Vereinfachen']['OUTPUT']
        return results

##       ## --------------------------------------------------------------------------------------------------------------## 
##        ## --------------------------------------------------------------------------------------------------------------## 
##        ## 13. Add Layers to Map
##        root = QgsProject.instance().layerTreeRoot()
##        mygroup = root.insertGroup(0,"Kubatur_Tool")        
##        
##        ## --------------------------------------------------------------------------------------------------------------## 
##        ## 13.1. Add Impact-Points
##        layer1 = QgsVectorLayer(pts, "Impact-Points", "ogr")
##
##        # Symbology
##        symbol = QgsSymbol.defaultSymbol(layer1.geometryType())    
##        
##        cat0_wert = 0
##        cat0_symbol = QgsMarkerSymbol.createSimple({'color': 'black', 'size': '1', 'outline_color': 'black'})
##        cat0 = QgsRendererCategory(cat0_wert, cat0_symbol, 'Start')
##        
##        categories =(cat0)
##        renderer = QgsCategorizedSymbolRenderer('grad_class', categories)
##        if renderer is not None:
##            layer1.setRenderer(renderer)
##        layer1.triggerRepaint()
##        
##        QgsProject.instance().addMapLayer(layer1, False)
##        mygroup.addLayer(layer1)
##
##        ## --------------------------------------------------------------------------------------------------------------## 
##        ## 13.2. Fliesspfade
##        layer2 = QgsVectorLayer(line_path, "Fliesspfade", "ogr")
##        label = QgsPalLayerSettings()
##        labeler2 = QgsVectorLayerSimpleLabeling(label)
##        layer2.setLabelsEnabled(True)
##        layer2.setLabeling(labeler2)
##        layer2.renderer().symbol().setColor(QColor("#33b5f5"))
##        layer2.triggerRepaint()
##        QgsProject.instance().addMapLayer(layer2, False)
##        mygroup.addLayer(layer2)



