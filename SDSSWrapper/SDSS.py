import requests
import io
import pandas as pd
from astropy import units as u
import astropy.coordinates as c
from astropy.coordinates import SkyCoord
import numpy as np

def multi_assign_floats(keys, values, dictionary):
    """
    Appends values to a dictionary consisting of only lists
    """
    for i, key in enumerate(keys):
        for row in values:
            dictionary[key].append(row[i])
    return dictionary


def calc_luminosity(redshift, flux_density, wl):
    """
    Takes a redshift and flux density and returns an approximate luminosity.
    """
    distance = c.Distance(z=redshift)
    distance = distance.value
    distance = distance * 3.086 * 10**19
    flux = flux_density * 3.631 * (10 ** -32) * wl
    lum = 4 * np.pi * flux * (distance ** 2)
    return lum

class SDSS_API:
    """
    Basic class for wrapping the rectangular search SDSS
    """

    def __init__(self):
        self.rec_url = 'http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/RectangularSearch'
        self.sql_url = 'http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/SqlSearch'

    def rectangular_search(self, ra_range, dec_range, searchtype='equitorial', limit='10', whichquery='irspectra'):
        """
        ra_range and dec_range should be in tuple form, (lower bound, upper bound)
        """
        prams = {'min_ra': str(ra_range[0]), 'max_ra': str(ra_range[1]), 'min_dec': str(dec_range[0]), 'max_dec': str(
            dec_range[1]), 'searchtype': searchtype, 'limit': limit, 'format': 'csv', 'whichquery': whichquery}
        return(requests.get(self.rec_url, prams).text)

    def sql_search(self, cmd, fmat='csv', whichquery='irspectra', lim='10'):
        prams = {'cmd': cmd, 'format': fmat, 'searchtype': 'equitorial',
                 'whichquery': whichquery, 'limit': lim}
        r = requests.get(self.sql_url, prams)
        search_results = {}
        text = r.text.replace('#Table1\n','')
        split = text.split('\n')
        split = split[:-1]
        headers = split[0].split(',')
        data = [row.split(',') for row in split[1:]]
        for header in headers:
            search_results[header] = []
        search_results = multi_assign_floats(headers, data, search_results)
        return search_results


    def get_luminosity(self, galaxy):
        """
        Returns the luminosity of a specific galaxy

        """
        freq = 3.282159601489 * 10 ** 14
        command = 'select z, sky_z \n from photoObj \n where ObjID = ' + str(galaxy)
        r = self.sql_search(command)
        z, flux = r['z'][0], r['sky_z'][0]
        lum = calc_luminosity(float(z), float(flux), freq)
        return lum
        

    def close_galaxy_search(self, ra, dec, z, interval, z_interval):
        command = 'select \n z, ra, dec, bestObjID \n from \n specObj \n where \n class = \'galaxy\'  \n and zWarning = 0 \n and ra<' + \
            str(ra + interval) + '\n and ra>' + str(ra - interval) + '\n and dec<' + \
            str(dec + interval) + '\n and dec>' + str(dec - interval) + '\n and z >' + str(z-z_interval) + '\n and z <' + str(z+z_interval)
        search_results = self.sql_search(command)
        galaxies = pd.DataFrame.from_dict(search_results)
        return galaxies

    def find_closest_galaxy(self, ra, dec, z):
        if ra < 0:
            ra += 360
        c_sne = SkyCoord(ra, dec, unit='deg')
        near_galaxies = self.close_galaxy_search(ra, dec, z, interval=0.5, z_interval=0.05)
        redshifts = near_galaxies['z'].tolist()
        gal_ra = near_galaxies['ra'].tolist()
        gal_dec = near_galaxies['dec'].tolist()
        gal_ids = near_galaxies['bestObjID'].tolist()
        if not len(redshifts) != 0:
            raise ValueError("No galaxies found within the interval!")
        distance = 100.0
        closest_galaxy = ''
        for z, r, d, objid in zip(redshifts, gal_ra, gal_dec, gal_ids):
            c = SkyCoord(r, d, unit='deg')
            sep = c.separation(c_sne).radian
            if sep < distance:
                distance = sep
                closest_galaxy = objid
        if closest_galaxy == 0:
            raise ValueError('ObjID is 0 so returns multiple galaxies')
        return near_galaxies.loc[near_galaxies['bestObjID'] == closest_galaxy], distance



ss = SDSS_API()
print(ss.get_luminosity(1237663784201748556))