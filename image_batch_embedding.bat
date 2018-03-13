@echo off

c:
cd c:
cd C:\Users\eurico\SunspotAnalysis\SpatialTemporalFeatureSelectionOptimalApproach

"C:\Program Files (x86)\IrfanView\i_view32.exe" sunspots_original_embedding.ppm   /gamma=5 /resize=(1888,250)  /convert=sunspots_original_embedding.jpg /jpgq=100
"C:\Program Files (x86)\IrfanView\i_view32.exe" sunspots_forecast_embedding.ppm   /gamma=5 /resize=(1888,250) /convert=sunspots_forecast_embedding.jpg /jpgq=100
"C:\Program Files (x86)\IrfanView\i_view32.exe" /panorama=(2,sunspots_original_embedding.jpg,sunspots_forecast_embedding.jpg)  /jpgq=100  /convert=sunspots_both_embedding.jpg 

