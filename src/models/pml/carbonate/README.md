### Using the Mauna Loa CO2 observations as GOTM input

The [data](https://scripps.ucsd.edu/programs/keelingcurve/wp-content/plugins/sio-bluemoon/graphs/mlo_full_record.png) are avilable as download from 
[Mauna Loa CO2 monthly mean data](ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt)

After download use the following command to convert to GOTM format:

./mouna_loa2gotm.py co2_mm_mlo.txt > CO2.dat

