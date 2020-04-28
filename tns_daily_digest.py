# Script to identify new discoveries associated with nearby galaxies and subsequent TNS classifications.

import json
import requests
import datetime
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.ned import Ned
from astroquery.ned.core import RemoteServiceError
from astropy.cosmology import FlatLambdaCDM

from time import sleep
from time import time as pytime

# TNS API key and URLs
TNS_API_KEY = "ab4e069b026bd426d1ca3e37cb46c64a21c5a455"
TNS_API_URL = "https://wis-tns.weizmann.ac.il/api/get"

# Define cosmology used to calculate distances to galaxy
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Method to use TNS API to get sources
def tns_api_call(url, data):
    
    send_data = [
        ('api_key', (None, TNS_API_KEY)),
        ('data', (None, json.dumps(data)))
    ]
    
    response = requests.post(url, files=send_data)
    
    if response.status_code != 200:
        
        logger.warning("TNS search return error code {}".format(
            response.status_code))
        
        return
    
    return json.loads(response.text)['data']['reply']

# Method to query NED for associated galaxies
# Offset is distance from galaxy in arcminutes

def query_region_NED(ra, dec, radius):
    
    pos = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
    
    try:
        # Query NED
        results = Ned.query_region(pos, radius=radius*u.arcmin)
        
        # Convert results to pandas dataframe
        results_df = results.to_pandas()
    
    except:
        
        # Return blank entry if no match found
        results_df = pd.DataFrame([], columns=['Object Name',
                                               'Type',
                                               'Magnitude and Filter',
                                               'RA',
                                               'Dec',
                                               'Redshift',
                                               'Separation'])
    
    # Format strings to UTF-8
    for col in results_df:
        if col in ['Object Name', 'Type', 'Magnitude and Filter']:
            results_df[col] = results_df[col].str.decode('utf-8')
    
    # Define columns
    cols = ['ned_name', 'ned_host_ra', 'ned_host_dec', 'ned_type',
            'ned_mag_filter', 'ned_host_redshift', 'ned_offset (arcmin)']
    
    # Rename columns
    results_df = results_df.rename(columns={'Object Name':'ned_name',
                                            'Type':'ned_type',
                                            'Magnitude and Filter':'ned_mag_filter',
                                            'RA':'ned_host_ra',
                                            'DEC':'ned_host_dec',
                                            'Redshift':'ned_host_redshift',
                                            'Separation':'ned_offset (arcmin)'})
    
    results_df = results_df[cols]
    
    # Calculate luminosity distance from host redshift, Mpc
    results_df['luminosity_distance (Mpc)'] = cosmo.luminosity_distance(results_df['ned_host_redshift']).value
    
    # Keep nearest galaxy
    results_df = results_df[
        results_df['ned_offset (arcmin)']==min(results_df['ned_offset (arcmin)'])
    ].reset_index(drop=True)
    
    return results_df.to_dict(orient='records')[0]

# Create search parameters to get objects reported previous night
today = datetime.datetime.today()
timestamp = '{0}-{1}-{2}'.format(today.year, today.month, today.day)
search_data = dict(
    public_timestamp=timestamp)

# Execute query to search objects realeased from timestamp
search_url = '{}/search'.format(TNS_API_URL)
search_reply = tns_api_call(search_url, search_data)

# Get object details from TNS
object_search_url = '{}/object'.format(TNS_API_URL)

tns_objects = []
print('Getting object details ...')

for i, obj in enumerate(search_reply):
    
    print('Progress = {0:.1f} % ..'.format((i+1)*100/len(search_reply)),
         end='\r')
    
    objname = obj['objname']
    search_data = dict(objname=objname)
    
    obj_search_reply = tns_api_call(object_search_url, search_data)
    
    # Format some columns
    discmagfilter_family = obj_search_reply['discmagfilter']['family']
    discmagfilter_name = obj_search_reply['discmagfilter']['name']
    obj_search_reply['discmagfilter'] = '{0}-{1}'.format(
        discmagfilter_family, discmagfilter_name)
    
    object_type = obj_search_reply['object_type']['name']
    obj_search_reply['object_type'] = object_type
    
    reporting_group = obj_search_reply['reporting_group']['group_name']
    obj_search_reply['reporting_group'] = reporting_group
    
    disc_data_src = obj_search_reply['discovery_data_source']['group_name']
    obj_search_reply['discovery_data_source'] = disc_data_src
    
    tns_objects.append(obj_search_reply)

results = list()

# Search NED for nearby galaxies for TNS objects
print('Searching NED for nearby galaxies ...')
for i, obj in enumerate(tns_objects):
    
    print('Progress = {0:.1f} % ..'.format((i+1)*100/len(tns_objects)),
         end='\r')
    ra = obj['radeg']
    dec = obj['decdeg']
    
    ned_results = query_region_NED(ra, dec, radius=1.0)
    
    # Add host galaxy information to TNS object information
    obj.update(ned_results)
    results.append(obj)

# Get objects with a luminosity distance
obj_df = pd.DataFrame(results)
obj_df = obj_df[obj_df['luminosity_distance (Mpc)' ]> 0].reset_index(drop=True)

# Save to .csv file
obj_df.to_csv('tns_results-{}.csv'.format(timestamp), index=None)