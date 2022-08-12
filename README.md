# pyAlpineRisk
Natural hazard management of alpine torrents.

In a time of increasing meteorological extreme events, high-resolution terrain models are an important tool for managing natural hazards. Combining simplified models of the terrain with different spatial analysis methods, various information and parameters can be obtained, which make it easier to answer specific hydrological questions. pyAlpineRisk is a Python-based toolset that can be used in natural hazard management to support a variety of workflows. The collection of tools ranges from the first assessment of torrent catchments to the classification of torrents and the utilization of sediment material.

1. Change Detection Tool

2. Retention Basin Tool
GIS tool for determining the retention basin of flood retention structures

This tool allows an automated calculation of the volume within a generated intersection area, based on a query line. Based on the query line, detailed LiDAR terrain models are used to determine location and elevation formations and to generate an artificial surface. Using this generated surface and the LiDAR data, a difference model is created and the volume is calculated from the positive value areas. As a result, the tool delivers a corresponding intersection area and the sums from the volume calculation. The field of application ranges from reservoir surveys for water retention structures to the detection and analysis of terrain depressions to large-scale valley and basin floor analyses. This tool was implemented using the Python programming language and is available for the QGIS and ArcMap programs.

# case examples

input parameters:

![B1](https://user-images.githubusercontent.com/52344347/184308968-6d689638-a457-4606-97c2-0eb275bc241e.jpg)

results:

![B2](https://user-images.githubusercontent.com/52344347/184309283-4f7a7b2e-472d-4c50-8a58-77cc8f0bde50.jpg)

![B3](https://user-images.githubusercontent.com/52344347/184309293-4278199e-f25b-48db-8a19-c8435dd84722.jpg)


# Installation/Application
The scripts are written for PyQGIS 3.16/3.22 and can be used by installing QGIS 3.16/3.22 or above.

To install QGIS tools developed for QGIS 3.x, copy them into
~/AppData/Roaming/QGIS/QGIS3/profiles/default/processing/scripts or in the upper part of the toolbox dialog you can add the scripts with ![mIconPythonFile](https://user-images.githubusercontent.com/52344347/136413201-b4a1f7d3-4053-4aa6-b11c-9433ae617057.png) Scripts - Add Script to Toolbox ...

After that the tools can be found in the QGIS "Processing Toolbox" - Scripts

# License
The tools in this repository are free software; you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation.

![QCD_Tool](https://user-images.githubusercontent.com/52344347/136409495-74d62525-a26a-4dbd-8a8f-b2b18425acca.png)
