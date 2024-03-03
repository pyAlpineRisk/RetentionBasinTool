<h1>pyAlpineRisk</h1>
<p>
<u>Natural hazard management of alpine torrents</u>
</p>

<p>In a time of increasing meteorological extreme events, high-resolution terrain models are an important tool for managing natural hazards. Combining simplified models of the terrain with different spatial analysis methods, various information and parameters can be obtained, which make it easier to answer specific hydrological questions. pyAlpineRisk is a Python-based toolset that can be used in natural hazard management to support a variety of workflows. The collection of tools ranges from the first assessment of torrent catchments to the classification of torrents and the utilization of sediment material.</p>
 
<h2>Retention Basin Tool</h2>
<p>GIS tool for determining the retention basin of flood retention structures. This tool allows an automated calculation of the volume within a generated intersection area, based on a query line. Based on the query line, detailed LiDAR terrain models are used to determine location and elevation formations and to generate an artificial surface. Using this generated surface and the LiDAR data, a difference model is created and the volume is calculated from the positive value areas. As a result, the tool delivers a corresponding intersection area and the sums from the volume calculation. The field of application ranges from reservoir surveys for water retention structures to the detection and analysis of terrain depressions to large-scale valley and basin floor analyses. This tool was implemented using the Python programming language and is available for the QGIS and ArcMap programs.</p>
<p><strong>QGIS-Note:</strong> Add Layerfile (pyAlpineRisk/layerfiles/layerfile.qlr) to your QGIS-Project! </p>

<h2>Case Examples</h2>

<i lang="id">Input Parameters:</i>

![B1](https://user-images.githubusercontent.com/52344347/184308968-6d689638-a457-4606-97c2-0eb275bc241e.jpg)

<i lang="id">Results:</i>

![B2](https://user-images.githubusercontent.com/52344347/184309283-4f7a7b2e-472d-4c50-8a58-77cc8f0bde50.jpg)

![B3](https://user-images.githubusercontent.com/52344347/184309293-4278199e-f25b-48db-8a19-c8435dd84722.jpg)


<h2>Installation/Application</h2>
<p>The latest release is written for PyQGIS 3.28 and can be used by installing QGIS 3.28 or above.

To install QGIS tools developed for QGIS 3.x, copy them into
~/AppData/Roaming/QGIS/QGIS3/profiles/default/processing/scripts or in the upper part of the toolbox dialog you can add the scripts with ![mIconPythonFile](https://user-images.githubusercontent.com/52344347/136413201-b4a1f7d3-4053-4aa6-b11c-9433ae617057.png) Scripts - Add Script to Toolbox ...

After that the tools can be found in the QGIS "Processing Toolbox" - Scripts</p>
