# -*- coding: cp1252 -*-
## --------------------------------------------------------------------------- ##
##
## 
## Python Script zur Berechnung des Volumens des Geschieberückhaltebeckens
##
## Software: ArcGIS 10.4.1
##
##
## Description: 
##
## Copyright © Franz Langegger & Nicole Kamp, Graz, August 2016, Version 1.0.
## niki.kamp@gmail.com
## Revision: 

##
## --------------------------------------------------------------------------- ##

# Import modules
import arcpy
from arcpy import env
from arcpy.sa import *
import string, os, sys, copy, shutil, arcgisscripting
from time import *
from sys import *

# Create the Geoprocessor object
gp = arcgisscripting.create()

#Meldungen auch am Bildschirm (True)
Bildschirm = False

arcpy.ResetEnvironments()
arcpy.env.overwriteOutput = True

# ---------------------------------------------------------------------------#
# ---------------------------------------------------------------------------#
def gueltig(pfad):
    if arcpy.Exists (pfad):
        print "Vorgabe gültig:"+pfad
        protokoll1.write("Vorgabe gültig:"+pfad+"\n")
    else:
        print "**"
        print "..Fehlerhafte Variable:"+pfad
        print "**"
        protokoll1.write("\n**\n")
        protokoll1.write("..Fehlerhafte Variable:"+pfad+"\n")
        protokoll1.write("**\n")
    return


# ---------------------------------------------------------------------------#
# ---------------------------------------------------------------------------#
#### INPUT
## test or not
test=0
if test==0:
    shape = str(sys.argv[1])  
    hoehe = float(sys.argv[2])
    prozent = float(sys.argv[3])
    input_als = str(sys.argv[4])

elif test==1:
    shape = "J:/01_dissertation_neu/06_tools/01_lidar/03_volumen_rückhaltebecken/Abfragelinie.shp"
    hoehe = 4
    prozent = 0


ordner =  "C:/ArcGISTemp/Rueckhaltevolumen/temp"


ordner_neu = ordner + "/work"
if os.path.exists(ordner_neu): arcpy.Delete_management(ordner_neu)

if not os.path.exists(ordner_neu): os.makedirs(ordner_neu)

output = ordner + "/output"
if os.path.exists(output): arcpy.Delete_management(output)

if not os.path.exists(output): os.makedirs(output)


protokoll = output + "/Volumen_Geschieberueckhaltebeckens.txt"
protokoll1 = open(protokoll, "w")

info = output + "/Volumen_Geschieberueckhaltebeckens_INFO.txt"
info1 = open(info, "w")



# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
##Starting Times
time = strftime("Starting Time: %x %X %Z") + "\n-----------------------------------------------\n"
protokoll1.write(time)

status = gp.CheckOutExtension("spatial")
if not (status == 'CheckedOut'):
    print "Derzeit steht keine Spatial Analyst Lizenz bereit"
    protokoll1.write("\nDerzeit steht keine Spatial Analyst Lizenz bereit\n")
    sys.exit(1)

status = gp.CheckOutExtension("3D")
if not (status == 'CheckedOut'):
    print "Derzeit steht keine 3D Analyst Lizenz bereit"
    protokoll1.write("\nDerzeit steht keine Spatial Analyst Lizenz bereit\n")
    sys.exit(1)

# Gueltige Variablen
gueltig(input_als)


# ---------------------------------------------------------------------------#
## Create CreatePersonalGDB
arcpy.AddMessage("Create GDB...")

if arcpy.Exists (ordner_neu + "/data.gdb"):
	protokoll1.write("alte Geodatenbank vorhanden --> loeschen eingeleitet" + "\n")
	arcpy.Delete_management(ordner_neu + "/data.gdb")
arcpy.CreateFileGDB_management(out_folder_path= ordner_neu, out_name="data", out_version="CURRENT")
	
report = "Process: Create GDB sucessful" + "\n"
protokoll1.write(report)


## 1. SHAPE in UTM33N umprojezieren
arcpy.AddMessage("1. SHAPE in UTM33N umprojezieren...")

shape_select = ordner_neu + "/data.gdb/shape_lambert_select"

arcpy.Select_analysis(in_features=shape, out_feature_class=shape_select, where_clause="")
utm33n = "PROJCS['WGS_1984_UTM_Zone_33N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',00],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',15.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]"
transform = "MGI_To_WGS_1984_8"
lambert = "PROJCS['MGI_Austria_Lambert',GEOGCS['GCS_MGI',DATUM['D_MGI',SPHEROID['Bessel_1841',6377397.155,299.1528128]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',400000.0],PARAMETER['False_Northing',400000.0],PARAMETER['Central_Meridian',13.33333333333333],PARAMETER['Standard_Parallel_1',46.0],PARAMETER['Standard_Parallel_2',49.0],PARAMETER['Latitude_Of_Origin',47.5],UNIT['Meter',1.0]]"
arcpy.Project_management(in_dataset=shape_select, out_dataset=ordner_neu + "/data.gdb/shape_utm33n", out_coor_system=utm33n, transform_method=transform, in_coor_system=lambert, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="")

report = "Shape erfolgreich in UTM33N umprojeziert" + "\n"
protokoll1.write(report)


## 2. ADD SURFACE INFORMATION
arcpy.AddMessage("2. ADD SURFACE INFORMATION...")

shape_neu = ordner_neu + "/data.gdb/shape_utm33n"
arcpy.AddSurfaceInformation_3d(in_feature_class=shape_neu, in_surface=input_als, out_property="Z_MIN", method="BILINEAR", sample_distance="", z_factor="1", pyramid_level_resolution="0", noise_filtering="NO_FILTER")

report = "Prozess: Add Surface Information erfolgreich" + "\n"
protokoll1.write(report)


## 3. ADD FIELDS (SperreH, AbfrageH, H100)
arcpy.AddMessage("3. ADD FIELDS (SperreH, AbfrageH, H100)...")

arcpy.AddField_management(in_table=shape_neu, field_name="SperreH", field_type="DOUBLE", field_precision="", field_scale="", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table=shape_neu, field_name="AbfrageH", field_type="DOUBLE", field_precision="", field_scale="", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table=shape_neu, field_name="H100", field_type="DOUBLE", field_precision="", field_scale="", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")

hoehe_sperre = "!Z_MIN!+"+ str(hoehe)
arcpy.CalculateField_management(in_table=shape_neu, field="SperreH", expression=hoehe_sperre, expression_type="PYTHON_9.3", code_block="")

wert = 300 * prozent/100
hoehe_abfrage = "!SperreH!+"+str(wert)
arcpy.CalculateField_management(in_table=shape_neu, field="AbfrageH", expression=hoehe_abfrage, expression_type="PYTHON_9.3", code_block="")

hoehe100 = "!Z_MIN!+100"
arcpy.CalculateField_management(in_table=shape_neu, field="H100", expression=hoehe100, expression_type="PYTHON_9.3", code_block="")
    

## Loop - Hoehe
i = 1
liste_tin=[]
while i < hoehe:
    liste_tin.append(i)
    arcpy.AddField_management(in_table=shape_neu, field_name="SperreH"+str(i), field_type="DOUBLE", field_precision="", field_scale="", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")

    hoehe_sperre = "!Z_MIN!+"+ str(i)
    arcpy.CalculateField_management(in_table=shape_neu, field="SperreH"+str(i), expression=hoehe_sperre, expression_type="PYTHON_9.3", code_block="")

    i = i + 1

report = "Prozess: Add Fields erfolgreich" + "\n"
protokoll1.write(report)


## 4. UPDATE CURSOR - Erstellen eines Rechtecks (Berechnungsfläche)
## 300M fuer Berechnungsflaeche

arcpy.AddMessage("4. Erzeugen der Berechnungsflaeche ...")

fields = ['xstart','xend','ystart','yend']

# Add fields to your FC
for field in fields:
    arcpy.AddField_management(shape_neu,str(field),"DOUBLE")

with arcpy.da.UpdateCursor(shape_neu, ('xstart','xend','ystart','yend', "SHAPE@")) as cursor:
    for row in cursor:
        row[0] = row[4].firstPoint.X
        x1 = row[0]
        row[1] = row[4].lastPoint.X
        x2 = row[1]
        row[2] = row[4].firstPoint.Y
        y1 = row[2]
        row[3] = row[4].lastPoint.Y
        y2 = row[3]
        cursor.updateRow(row)

# Berechnung der neuen Koordinaten für die Fläche
y_neu = y2 - y1
x_neu = x2 - x1
phi = math.radians(90)
r1 = math.atan2 (y_neu, x_neu)
r2 = r1 + phi
r3 = r1 - phi
x3 = x2 + 300 * math.cos(r2)
y3 = y2 + 300 * math.sin(r2)
x4 = x1 + 300 * math.cos(r2)
y4 = y1 + 300 * math.sin(r2)
x5 = x2 + 300 * math.cos(r3)
y5 = y2 + 300 * math.sin(r3)
x6 = x1 + 300 * math.cos(r3)
y6 = y1 + 300 * math.sin(r3)        


# Pfad und Name des Geodatensatzes, wo das Rechteck verspeichert wird
rechteck = ordner_neu + "/data.gdb/flaeche_utm33n"

# Der Shape Rechteck wird geloescht, falls er bereits existiert
pruefObGibt = gp.exists(rechteck)
if pruefObGibt :
    gp.delete(rechteck)

# Anlegen eines leeren Shape rechteck im Koordinatensystem UTMN33N
arcpy.CreateFeatureclass_management(out_path=ordner_neu + "/data.gdb", out_name="flaeche_utm33n", geometry_type="POLYGON", template="", has_m="DISABLED", has_z="DISABLED", spatial_reference=utm33n, config_keyword="", spatial_grid_1="0", spatial_grid_2="0", spatial_grid_3="0")

# Rechteck in Shape schreiben
aGPPnt = gp.CreateObject("Point")
RectangleAsGPArray = gp.CreateObject("Array")
cur = gp.InsertCursor(rechteck)
futter = cur.NewRow()
RectangleAsList = [[x6,y6],[x4,y4],[x3,y3],[x5,y5],[x6,y6]]
for idx in range(len(RectangleAsList)):
    aGPPnt.id = idx
    aGPPnt.x = RectangleAsList[idx][0]
    aGPPnt.y = RectangleAsList[idx][1]
    RectangleAsGPArray.add(aGPPnt)
futter.shape = RectangleAsGPArray
RectangleAsGPArray.RemoveAll()
cur.InsertRow(futter)
cur = None

arcpy.AddField_management(in_table=rechteck, field_name="AbfrageH", field_type="DOUBLE", field_precision="", field_scale="", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")

field = ['AbfrageH']
with arcpy.da.UpdateCursor(shape_neu, field) as cursor:
    for row in cursor:
        hoehe_abfrage = row[0]

with arcpy.da.UpdateCursor(rechteck, field) as cursor:
    for row in cursor:
        row[0] = hoehe_abfrage
        cursor.updateRow(row)

report = "Prozess: Erzeugen des Rechtecks inkl. Anlegen einer neuen Attibutspalte abgeschlossen. " + "\n"
protokoll1.write(report)


## 5. Create TIN
arcpy.AddMessage("5. Create TIN...")

tin = ordner_neu + "/tin"
arcpy.CreateTin_3d(out_tin=tin, spatial_reference=utm33n, in_features=rechteck+" AbfrageH Soft_Clip <None>;" + shape_neu + " SperreH Hard_Line <None>", constrained_delaunay="DELAUNAY")

for j in liste_tin:
    tin_h = ordner_neu + "/tin"+str(j)
    arcpy.CreateTin_3d(out_tin=tin_h, spatial_reference=utm33n, in_features=rechteck+" AbfrageH Soft_Clip <None>;" + shape_neu + " SperreH"+str(j)+" Hard_Line <None>", constrained_delaunay="DELAUNAY")

report = "Prozess: Create TIN erfolgreich" + "\n"
protokoll1.write(report)



## 6. Extract ALS Raster by MASK (Berechnungsflaeche)
# Set Snap Raster environment
arcpy.env.snapRaster = input_als

arcpy.AddMessage("6. Extract ALS Raster by MASK (Berechnungsflaeche)...")

als_neu = ordner_neu + "/data.gdb/als"
arcpy.gp.ExtractByMask_sa(input_als, rechteck, als_neu)

report = "Prozess: Extract by Mask erfolgreich" + "\n"
protokoll1.write(report)


## 7. ALS Raster to TIN
arcpy.AddMessage("7. ALS Raster to TIN...")

als_pt = ordner_neu + "/data.gdb/als_pt"
als_tin = ordner_neu + "/als_tin"

arcpy.RasterToPoint_conversion(in_raster=als_neu, out_point_features=als_pt, raster_field="Value")

arcpy.CreateTin_3d(out_tin=als_tin, spatial_reference=utm33n, in_features=als_pt +" grid_code Mass_Points <None>;"+ shape_neu +" H100 Hard_Line <None>", constrained_delaunay="DELAUNAY")

report = "Prozess: Raster to TIN erfolgreich" + "\n"
protokoll1.write(report)


## 8. Surface Difference
arcpy.AddMessage("8. Surface Difference...")

volume_t1 = ordner_neu + "/data.gdb/geschieberueckhaltevolumen_t1"

arcpy.SurfaceDifference_3d(in_surface=als_tin, in_reference_surface=tin, out_feature_class=volume_t1, pyramid_level_resolution="0", reference_pyramid_level_resolution="0", out_raster="", raster_cell_size="10", out_tin_folder="", out_tin_basename="")

for k in liste_tin:
    volume_t1_h = ordner_neu + "/data.gdb/geschieberueckhaltevolumen_t1_"+str(k)
    tin_h = ordner_neu + "/tin"+str(k)
    arcpy.SurfaceDifference_3d(in_surface=als_tin, in_reference_surface=tin_h, out_feature_class=volume_t1_h, pyramid_level_resolution="0", reference_pyramid_level_resolution="0", out_raster="", raster_cell_size="10", out_tin_folder="", out_tin_basename="")



report = "Prozess: Surface Difference erfolgreich" + "\n"
protokoll1.write(report)



## 9. Update Shape
arcpy.AddMessage("9. Update Shape...")

volume_t2 = ordner_neu + "/data.gdb/geschieberueckhaltevolumen_t2"

if arcpy.Exists (volume_t2):
	arcpy.Delete_management(volume_t2)


arcpy.MakeFeatureLayer_management(in_features=volume_t1, out_layer="geschieberueckhaltevolumen_layer", where_clause="", workspace="", field_info="OID OID VISIBLE NONE;Shape Shape VISIBLE NONE;Volume Volume VISIBLE NONE;SArea SArea VISIBLE NONE;Code Code VISIBLE NONE;Shape_Length Shape_Length VISIBLE NONE;Shape_Area Shape_Area VISIBLE NONE")
arcpy.SelectLayerByLocation_management(in_layer="geschieberueckhaltevolumen_layer", overlap_type="WITHIN_A_DISTANCE", select_features=shape_neu, search_distance="4 Meters", selection_type="NEW_SELECTION", invert_spatial_relationship="NOT_INVERT")
arcpy.Select_analysis(in_features="geschieberueckhaltevolumen_layer", out_feature_class=volume_t2, where_clause="")


field1 = ['CODE']
with arcpy.da.UpdateCursor(volume_t2, field1) as cursor:
    for row in cursor:
        if row[0] == 1:
            cursor.deleteRow()

list=[]
field2 = ['Volume']
with arcpy.da.UpdateCursor(volume_t2, field2) as cursor:
    for row in cursor:
        list.append(row[0])

max_value = max(list)

field3 = ['Volume']
with arcpy.da.UpdateCursor(volume_t2, field3) as cursor:
    for row in cursor:
        if row[0] == max_value:
            cursor.deleteRow()


for l in liste_tin:
    volume_t2_h = ordner_neu + "/data.gdb/geschieberueckhaltevolumen_t2_"+str(l)
    volume_t1_h = ordner_neu + "/data.gdb/geschieberueckhaltevolumen_t1_"+str(l)
    arcpy.MakeFeatureLayer_management(in_features=volume_t1_h, out_layer="geschieberueckhaltevolumen_layer"+str(l), where_clause="", workspace="", field_info="OID OID VISIBLE NONE;Shape Shape VISIBLE NONE;Volume Volume VISIBLE NONE;SArea SArea VISIBLE NONE;Code Code VISIBLE NONE;Shape_Length Shape_Length VISIBLE NONE;Shape_Area Shape_Area VISIBLE NONE")
    arcpy.SelectLayerByLocation_management(in_layer="geschieberueckhaltevolumen_layer"+str(l), overlap_type="WITHIN_A_DISTANCE", select_features=shape_neu, search_distance="4 Meters", selection_type="NEW_SELECTION", invert_spatial_relationship="NOT_INVERT")
    arcpy.Select_analysis(in_features="geschieberueckhaltevolumen_layer"+str(l), out_feature_class=volume_t2_h, where_clause="")


    field1 = ['CODE']
    with arcpy.da.UpdateCursor(volume_t2_h, field1) as cursor:
        for row in cursor:
            if row[0] == 1:
                cursor.deleteRow()

    list=[]
    field2 = ['Volume']
    with arcpy.da.UpdateCursor(volume_t2_h, field2) as cursor:
        for row in cursor:
            list.append(row[0])

    max_value = max(list)

    field3 = ['Volume']
    with arcpy.da.UpdateCursor(volume_t2_h, field3) as cursor:
        for row in cursor:
            if row[0] == max_value:
                cursor.deleteRow()

report = "Prozess: Update Shape erfolgreich" + "\n"
protokoll1.write(report)



## 10. Dissolve
arcpy.AddMessage("10. Dissolve...")
volume = output + "/geschieberueckhaltevolumen.shp"
arcpy.Dissolve_management(in_features=volume_t2, out_feature_class=volume, dissolve_field="", statistics_fields="Volume SUM", multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES")

for m in liste_tin:
    volume_t2_h = ordner_neu + "/data.gdb/geschieberueckhaltevolumen_t2_"+str(m)
    volume_t3_h = ordner_neu + "/data.gdb/geschieberueckhaltevolumen_t3_"+str(m)
    arcpy.Dissolve_management(in_features=volume_t2_h, out_feature_class=volume_t3_h, dissolve_field="", statistics_fields="Volume SUM", multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES")

report = "Prozess: Dissolve erfolgreich" + "\n"
protokoll1.write(report)



## 11. TXT erstellen und Volumen aus Shape lesen
arcpy.AddMessage("11. TXT erstellen und Volumen aus Shape lesen...")
header = "h\tVolumen\n"
info1.write(header)

field = ['SUM_Volume']

for l in liste_tin:
    volume_t3_h = ordner_neu + "/data.gdb/geschieberueckhaltevolumen_t3_"+str(l)
    with arcpy.da.UpdateCursor(volume_t3_h, field) as cursor:
        for row in cursor:
            value = row[0]/1000
            txt = str(value)
            txt_neu = txt.replace(',','.')
            info_txt = str(l)+"\t"+txt_neu+"\n"
            info1.write(info_txt)


with arcpy.da.UpdateCursor(volume, field) as cursor:
    for row in cursor:
        value = row[0]/1000
        txt = str(value)
        txt_neu = txt.replace(',','.')
        info_txt = str(hoehe)+"\t"+txt_neu+"\n"
        info1.write(info_txt)

info1.close()

report = "Prozess: Erstellen einer INFO TXT war erfolgreich" + "\n"
protokoll1.write(report)



## 12. Merge
arcpy.AddMessage("12. Merge...")

volume_gesamt = output + "/geschieberueckhaltevolumen_gesamt.shp"

gesamt = str(volume)
for n in liste_tin[::-1]:
    txt = ";"+str(ordner_neu + "/data.gdb/geschieberueckhaltevolumen_t3_"+str(n))
    gesamt = gesamt+txt

arcpy.Merge_management(inputs=gesamt, output=volume_gesamt, field_mappings="")

report = "Prozess: Merge erfolgreich" + "\n"
protokoll1.write(report)



        
# 13. SHAPE in AUSTRIA Lambert umprojezieren
arcpy.AddMessage("SHAPE in AUSTRIA Lambert umprojezieren...")

volume_lambert = output + "/geschieberueckhaltevolumen_lambert.shp"
volume_gesamt_lambert = output + "/geschieberueckhaltevolumen_gesamt_lambert.shp"

utm33n = "PROJCS['WGS_1984_UTM_Zone_33N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',15.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]"
transform ="MGI_To_WGS_1984_8"
lambert = "PROJCS['MGI_Austria_Lambert',GEOGCS['GCS_MGI',DATUM['D_MGI',SPHEROID['Bessel_1841',6377397.155,299.1528128]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',400000.0],PARAMETER['False_Northing',400000.0],PARAMETER['Central_Meridian',13.33333333333333],PARAMETER['Standard_Parallel_1',46.0],PARAMETER['Standard_Parallel_2',49.0],PARAMETER['Latitude_Of_Origin',47.5],UNIT['Meter',1.0]]"
arcpy.Project_management(in_dataset=volume, out_dataset=volume_lambert, out_coor_system=lambert, transform_method=transform, in_coor_system=utm33n, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="")
arcpy.Project_management(in_dataset=volume_gesamt, out_dataset=volume_gesamt_lambert, out_coor_system=lambert, transform_method=transform, in_coor_system=utm33n, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="")


report = "Shape erfolgreich in AUSTRIA Lambert umprojeziert" + "\n"
protokoll1.write(report)



print "Volumen des Geschieberueckhaltebeckens erfolgreich berechnet!"
report = "Volumen des Geschieberueckhaltebeckens erfolgreich berechnet!" + "\n"
protokoll1.write(report)
protokoll1.flush()
