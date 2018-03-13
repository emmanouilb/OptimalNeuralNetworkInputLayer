@echo off

c:
cd c:
cd C:\Users\eurico\SunspotAnalysis\SpatialTemporalFeatureSelectionOptimalApproach

"C:\Program Files (x86)\IrfanView\i_view32.exe" sunspots_original.ppm   /gamma=5 /resize=(1888,250) /convert=sunspots_original.jpg /jpgq=100
"C:\Program Files (x86)\IrfanView\i_view32.exe" sunspots_forecast.ppm   /gamma=5 /resize=(1888,250) /convert=sunspots_forecast.jpg /jpgq=100
"C:\Program Files (x86)\IrfanView\i_view32.exe" /panorama=(2,sunspots_original.jpg,sunspots_forecast.jpg)  /jpgq=100  /convert=sunspots_both.jpg 

