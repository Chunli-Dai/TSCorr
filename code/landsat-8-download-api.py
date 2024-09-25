''' Once you have anaconda/miniconda python (version 3.8 preferably) installed, install the following required packages:
    satsearch: pip install sat-search : required
    geopandas: conda install geopandas: option to work with shapefile, else just pass bbox to query

'''

import os
from satsearch import Search
import geopandas as gpd
import json

stac_url = 'https://earth-search.aws.element84.com/v0'

# I. Get json representation of shapefile polygon/polyline (but closed)
# Should have only one polygon
shpfile = '/fs/project/howat.4/landsat8_chunli/Elliot.shp'
aoi = gpd.read_file(shpfile)
shp_json = aoi.to_json()
geom = json.loads(shp_json)['features'][0]['geometry']

# II. Get the number of remote images
search = Search.search(url=stac_url,
                       datetime='2020-10-01/2021-11-01',
                       query=["eo:cloud_cover<50"],
                       collections=['landsat-8-l1-c2'],
                       intersects=geom
                      )
#You can replace intersects = geom with bbox=[-110, 39.5, -105, 40.5]
# That way you need not worry about shapefile/geopandas processing to get json reprentation (ie Part I)
print(f'intersects search: {search.found()} items')

# III. Get the items on remote server (AWS)
items = search.items()

# Download data, one band at a time
items.download('B8', filename_template='./data/${id}') ## Download Band 08 data
#items.download('B8', filename_template='/fs/project/howat.4/landsat8_chunli/test/${id}') ## Download Band 08 data
items.download('MTL', filename_template='./data/${id}') # Download metadata file
# And so on ...
# filename_template is just the full path to where you want to save file
# $(id) is required to match this original name format : LC08_L1TP_139045_20170304_20170316_01_T1 


