# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   QGIS RetentionBasin Tool  v.1.1.0-beta                                *
#   QQGIS-Version: 3.28.15-Firenze                                        * 
*   pyAlpineRisk                                                          *
*   Nicole Kamp & Franz Langegger                                         *
*   niki.kamp@gmail.com                                                   *
*   www.nicolekamp.com                                                    *
*   May 2020                                                              *
*   revised 03/2024                                                       *
*                                                                         *
*   Bugs:                                                                 *
*      1. Labeling with m³ not working yet                                *
*      2. Layers in group cannot be switched on and off                   *
***************************************************************************
"""


from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProject,
                       QgsVectorLayer,
                       QgsTextFormat,
                       QgsVectorLayerSimpleLabeling,
                       QgsPalLayerSettings,
                       QgsTextBufferSettings,
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
from qgis.PyQt.QtGui import QColor, QFont

from PIL import Image
Image.MAX_IMAGE_PIXELS = None
import osgeo.gdal as gdal
import numpy as np

from osgeo import ogr, osr
from osgeo import gdal
from shapely.geometry import Polygon
import matplotlib
import subprocess
from subprocess import call

import string, os, sys, copy, shutil, math, numpy, time, datetime
from time import *
from sys import *

import geopandas as gpd
from shapely.geometry import Point
from shapely.ops import unary_union


##---------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------##
###FUNCTIONS
# Get Elevation Value
def get_elevation(x, y, band, gt):
    px = int((x - gt[0]) / gt[1])
    py = int((y - gt[3]) / gt[5])

    # Are coordinates within raster
    if px < 0 or py < 0 or px >= band.XSize or py >= band.YSize:
        return float('nan')  # Outside the raster
    
    elevation = band.ReadAsArray(px, py, 1, 1)
    if elevation is None:
        return float('nan')

    return elevation[0][0]

##---------------------------------------------------------------------------------------------------------------------------##
#Save Raster
def save_grid_as_raster(filename, grid, geotransform, projection, nodata=None):
    # Assuming grid is a 2D numpy array of the same dimensions as the source raster
    rows, cols = grid.shape
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform(geotransform)
    outRaster.SetProjection(projection)
    outband = outRaster.GetRasterBand(1)
    if nodata is not None:
        outband.SetNoDataValue(nodata)
    outband.WriteArray(grid)
    outband.FlushCache()

##---------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------------------------------------------##
class KubaturProcessingAlgorithm(QgsProcessingAlgorithm):
    """
Python script to calculate the volume of the retention basin.
Input of the input and output parameters
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
        return 'retentionbasintool'

    def displayName(self):
        return self.tr('RetentionBasinTool_v_1_1_0-beta')

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

    ##---------------------------------------------------------------------------------------------------------------------------##
    ##---------------------------------------------------------------------------------------------------------------------------##
    def processAlgorithm(self, parameters, context, feedback):
        feedback = QgsProcessingMultiStepFeedback(1, feedback)
        results = {}
        outputs = {}
        
        laenge = 1000
        laenge_k = 10
        
        
        ## PREWORK
        ## Define paths + timestamp
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        temp_path = str(parameters[self.TEMP])+'/temp_'+str(timestamp)
        final_path = str(parameters[self.TEMP])+'/final_'+str(timestamp)
        
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
        if not os.path.exists(final_path):
            os.makedirs(final_path)
        
        ## Cell Size and EPSG-Code
        raster_cs = gdal.Open(str(parameters[self.INPUT_ALS]))
        gt_cs =raster_cs.GetGeoTransform()
        proj = osr.SpatialReference(wkt=raster_cs.GetProjection()) 
        cs = gt_cs[1]
        epsg = proj.GetAttrValue('AUTHORITY',1)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg))
               
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
        
        ## Crop raster to layer mask
        out = temp_path + '/' + 'clip_vrt.tif'
        alg_params = {
            'INPUT': str(parameters[self.INPUT_ALS]),
            'MASK':buffer,
            'SOURCE_CRS':QgsCoordinateReferenceSystem(epsg),
            'TARGET_CRS':QgsCoordinateReferenceSystem(epsg),
            'TARGET_EXTENT':None,
            'NODATA':None,
            'ALPHA_BAND':False,
            'CROP_TO_CUTLINE':True,
            'KEEP_RESOLUTION':True,
            'SET_RESOLUTION':False,
            'X_RESOLUTION':None,
            'Y_RESOLUTION':None,
            'MULTITHREADING':False,
            'OPTIONS':'',
            'DATA_TYPE':0,
            'EXTRA':'',
            'OUTPUT':out
        }
        processing.run("gdal:cliprasterbymasklayer", alg_params, context=context, feedback=feedback, is_child_algorithm=True)
               

        ##---------------------------------------------------------------------------------------------------------------------------##
        ## 1. generate profile points from a line
        out_pts = temp_path + '/' + 'out_pts.shp'
        out1 = temp_path + '/' + 'out_pts_Z.shp'
        alg_params = {
            'INPUT': str(parameters[self.INPUT_Shape]),
            'DISTANCE':cs,
            'START_OFFSET':0,
            'END_OFFSET':0,
            'OUTPUT':out_pts
        }
        processing.run('native:pointsalonglines', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        
        # 2. Get Height 
        punkte_gdf = gpd.read_file(out_pts)
        ds = gdal.Open(out)
        band = ds.GetRasterBand(1)
        gt = ds.GetGeoTransform()
        
        # Make sure that the CRS of the shapefile and the raster match
        raster_crs = osr.SpatialReference(wkt=ds.GetProjection())
        shapefile_crs = osr.SpatialReference(wkt=punkte_gdf.crs.to_wkt())
        if not raster_crs.IsSame(shapefile_crs):
            print("CRS stimmt nicht überein. Eine Transformation ist erforderlich.")
        else:
            print("CRS stimmt überein.")

        # Query elevation values for each point and save them in the GeoDataFrame
        punkte_gdf['Z'] = punkte_gdf.apply(
            lambda row: get_elevation(row.geometry.x, row.geometry.y, band, gt), axis=1)

        punkte_gdf.to_file(out1)
        

        ## 3. list with height values from the profile points
        listeZ =[]
        ds=ogr.Open(out1)
        lyr=ds.GetLayer()
        for feat in lyr:
            zVal = feat.GetField("Z")
            listeZ.append(zVal)
        del out1


        ## 4. calculate height of protective structure, query height, height100 and ZMin
        zmin = min(listeZ)
        objhoehe = zmin + parameters[self.INPUT_hoehe]
        abfragehoehe = objhoehe + (laenge * parameters[self.INPUT_prozent]/100)
        hoehe100 = zmin + 100
        #feedback.pushInfo(str(objhoehe))
        #feedback.pushInfo(str(abfragehoehe))
        #feedback.pushInfo(str(hoehe100))
        #feedback.pushInfo(str(zmin))
  
         # Calculate fields
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

        del out2, out3, out5
      
      
        ## 5. calculate rectangle + lines
        # Calculate coordinates
        ds=ogr.Open(out4)
        lyr=ds.GetLayer()
        pts_list=[]
        for i in range(lyr.GetFeatureCount()):
            feat=lyr.GetFeature(i)
            geom=feat.GetGeometryRef()
            ptcount = geom.GetPointCount()
            firstpoint=geom.GetPoint(0)
            lastpoint=geom.GetPoint(geom.GetPointCount()-1)
            x1 = firstpoint[0]
            y1 = firstpoint[1]
            x2 = lastpoint[0]
            y2 = lastpoint[1]
            
            for n in range(0, ptcount):
                point=geom.GetPoint(n)
                pts_list.append(point)
        #feedback.pushInfo(str(pts_list))
            
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
        
        # Reference line
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
                
        # Define output paths
        out6 = temp_path + '/' + 'out6.shp'
        out7 = temp_path + '/' + 'out7.shp'
        outL1 = temp_path + '/' + 'line1.shp'
        outL2 = temp_path + '/' + 'line2.shp'
        outL1_neu = temp_path + '/' + 'line1_pkt.shp'
        outL2_neu = temp_path + '/' + 'line2_pkt.shp'
        
        # Create Line1
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
               
        # Generate profile points from line1 #adapt
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

        # List with height values from the profile points
        listeZ_L1 =[]
        ds=ogr.Open(outL1_neu)
        lyr=ds.GetLayer()
        for feat in lyr:
            zVal = feat.GetField("Z")
            listeZ_L1.append(zVal) 

        # Create Line2
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
        
        # Generate profile points from line2 #adapt
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

        # List with height values from the profile points
        listeZ_L2 =[]
        ds=ogr.Open(outL2_neu)
        lyr=ds.GetLayer()
        for feat in lyr:
            zVal = feat.GetField("Z")
            listeZ_L2.append(zVal) 
        
        # Create polygon - check which line is higher + list of coordinates
        driver = ogr.GetDriverByName('Esri Shapefile')
        ds = driver.CreateDataSource(out6)
        layer = ds.CreateLayer('', None, ogr.wkbPolygon25D) #wkbPolygon
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        defn = layer.GetLayerDefn()
        feat = ogr.Feature(defn)
        feat.SetField('id', 1)
        
        pts_length = len(pts_list)
        #feedback.pushInfo(str(pts_length))
        zmin1 = min(listeZ_L1)
        zmin2 = min(listeZ_L2)
        if pts_length == 2:
            if zmin1 > zmin2:
                poly1 = Polygon([(x1,y1,z12), (x4,y4,z3456), (x3,y3,z3456), (x2,y2,z12), (x1,y1,z12)])
            else:
                poly1 = Polygon([(x1,y1,z12), (x2,y2,z12), (x5,y5,z3456), (x6,y6,z3456), (x1,y1,z12)])

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
            
            # Create rectangle1
            geom = ogr.CreateGeometryFromWkt(poly1)
            feat.SetGeometry(geom)
            layer.CreateFeature(feat)
            feat = geom = None 
            ds = layer = feat = geom = None 
                
        # Append AbfrageHoehe as attribute
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
 
 
        # 6. Calculate raster from 3D rectangle
        out8 = temp_path + '/' + 'out8.tif'
        alg_params = {
            'INTERPOLATION_DATA':str(out7)+'::~::0::~::1::~::0', #str(out7)+'::~::1::~::-1::~::2',
            'METHOD':0,
            'EXTENT':str(extent),
            'PIXEL_SIZE':float(cs),
            'OUTPUT':out8
        }
        processing.run('qgis:tininterpolation', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        
        ## 7. Calculate difference model and retention capacity  
        # Crop raster with rectangle and change resolution
        out9 = temp_path + '/' + 'out9.tif'
        out10 = temp_path + '/' + 'out10.tif'
        out11 = temp_path + '/' + 'out11.tif'
               
        # Clip DTM
        OutTile = gdal.Warp(str(out10), 
            str(out), 
            cutlineDSName=str(out7),
            cropToCutline=True,
            dstNodata = 0)
        OutTile = None
                        
        # Open clipped DTM
        dtm_clip_gdal = gdal.Open(out10)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
        band1 = dtm_clip_gdal.GetRasterBand(1)
        gt = dtm_clip_gdal.GetGeoTransform()
        width = dtm_clip_gdal.RasterXSize
        height = dtm_clip_gdal.RasterYSize
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5] 
        maxx = gt[0] + width*gt[1] + height*gt[2]
        maxy = gt[3]
        extent1=str(minx)+","+str(maxx)+","+str(miny)+","+str(maxy)
        extent2=(minx,maxx,miny,maxy)
        feedback.pushInfo(str(extent1))

        # Snap raster
        flaeche3d = gdal.Open(out8)
        format = "GTiff"
        options = ["-of", "GTiff", "-te", str(minx), str(miny), str(maxx), str(maxy), "-tr", str(cs), str(cs)]
        OutTile = gdal.Warp(out9, flaeche3d, options=options)
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
       
        # Raster to array
        kubatur = gdal.Open(out11)
        projection = kubatur.GetProjection()
        band1 = kubatur.GetRasterBand(1)
        array = band1.ReadAsArray()
        gt = kubatur.GetGeoTransform()
        width = kubatur.RasterXSize
        height = kubatur.RasterYSize
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5] 
        maxx = gt[0] + width*gt[1] + height*gt[2]
        maxy = gt[3]
        extent=str(minx)+","+str(maxx)+","+str(miny)+","+str(maxy) 
        
        # Create empty raster
        grid1 = np.zeros(shape=(height, width), dtype=np.float32)
        grid2 = np.zeros(shape=(height, width), dtype=np.float32)
        
        kubatur_gesamt = 0.0
        
        for row in range(height):
            for col in range(width):
                val1 = array[row, col]
                if val1 < 0:
                    nval = (val1 * -1) * float(cs) * float(cs)
                    kubatur_gesamt += nval
                    grid1[row, col] = nval
                    grid2[row, col] = 1.0  # Mark as type 1
                else:
                    #grid1[row, col] = np.nan  # Use NaN for non-negative values
                    grid2[row, col] = 2.0  # Mark as type 2

        kubatur_gesamt = round(kubatur_gesamt, 1)
        
        #Save raster
        out12 = temp_path + '/' + 'out12.tif'
        out13 = temp_path + '/' + 'out13.tif'
        save_grid_as_raster(out12, grid1, gt, projection, nodata=np.nan)
        save_grid_as_raster(out13, grid2, gt, projection, nodata=np.nan)
                
        # Georeference raster 1
        out12 = temp_path + '/' + 'out12.tif'
        out14 = final_path + '/' + 'kubatur.tif'
        src_ds = gdal.Open(out12)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
                
        dst_ds = driver.CreateCopy(out14, src_ds, 0)
        dst_ds.SetGeoTransform(gt)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg))
        dest_wkt = srs.ExportToWkt()
        dst_ds.SetProjection(dest_wkt)

        # Close files
        dst_ds = None
        src_ds = None
        
         # Georeference raster 2
        out13 = temp_path + '/' + 'out13.tif'
        out15 = temp_path + '/' + 'mask.tif'
        src_ds = gdal.Open(out13)
        format = "GTiff"
        driver = gdal.GetDriverByName( format )
                
        dst_ds = driver.CreateCopy(out15, src_ds, 0)
        dst_ds.SetGeoTransform(gt)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg))
        dest_wkt = srs.ExportToWkt()
        dst_ds.SetProjection(dest_wkt)

        # Close files
        dst_ds = None
        src_ds = None       

        
        ## Raster to polygon
        out15 = temp_path + '/' + 'mask.tif'
        out16 = temp_path + '/' + 'stauraum_alt.shp'
        src_ds = gdal.Open(out15)
        srcband = src_ds.GetRasterBand(1)
        driver = ogr.GetDriverByName('Esri Shapefile')
        ds = driver.CreateDataSource(out16)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg))
        dst_ds = ds.CreateLayer('', srs=srs)
        newField = ogr.FieldDefn('id', ogr.OFTInteger)
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
        
        # Simplify
        out19 = final_path + '/' + 'stauraum.shp'
        alg_params = {
            'INPUT': str(out18),
            'METHOD': 1,
            'TOLERANCE': 1,
            'OUTPUT': out19
        }
        outputs['Vereinfachen'] = processing.run('native:simplifygeometries', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        
        ##---------------------------------------------------------------------------------------------------------------------------##
        ## Add Layers
        root = QgsProject.instance().layerTreeRoot()
        mygroup = root.findGroup("Kubatur-Tool - Ergebnisse")
        
        # If the group doesn't exist, create it
        if mygroup is None:
            mygroup = root.insertGroup(0, "Kubatur-Tool - Ergebnisse")
        else:
            clone = mygroup.clone()
            parent = mygroup.parent() if mygroup.parent() else root
            parent.removeChildNode(mygroup)
            newGroup = parent.insertChildNode(0, clone)
            mygroup = newGroup or clone

        layer = QgsVectorLayer(out19, "Stauraum", "ogr")
        if not layer.isValid():
            feedback.pushInfo("Layer failed to load!")

        # Label Settings       
        # Text-Format
        text_format = QgsTextFormat()
        text_format.setFont(QFont("Arial", 10))
        text_format.setColor(QColor(0, 0, 255))  # Blau
        
        # Halo
        buffer_settings = QgsTextBufferSettings()
        buffer_settings.setEnabled(True)
        buffer_settings.setSize(1)  # Größe des Halos
        buffer_settings.setColor(QColor("white"))  # Farbe des Halos
        text_format.setBuffer(buffer_settings)
        
        #Label Setting
        label = QgsPalLayerSettings()
        label.enabled = True
        #label.fieldName = 'concat("kubatur", \' m³\')' #error
        label.fieldName = 'kubatur'
        
        label.setFormat(text_format)
        
        labeler = QgsVectorLayerSimpleLabeling(label)
        layer.setLabelsEnabled(True)
        layer.setLabeling(labeler)
        
        # Set the colour with 50% transparency
        symbol = layer.renderer().symbol()
        symbol.setColor(QColor(0, 0, 255, 127))  # RGBA für Blau mit 50% Transparenz

        layer.triggerRepaint()
        QgsProject.instance().addMapLayer(layer, False)
        mygroup.addLayer(layer)

        
        results['KubaturTool'] = outputs['Vereinfachen']['OUTPUT']
        return results



