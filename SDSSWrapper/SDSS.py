import requests
import pandas as pd
import astropy.coordinates as c
from astropy.coordinates import SkyCoord
import numpy as np


test_galaxy = 1237663784201748556 #the objID of a random galaxy for testing stuff :):):)

def multi_assign_floats(keys, values, dictionary):
    """
    Appends values to a dictionary consisting of only lists
    Returns the appended dict

    keys: list
        Keys for the values to be appended

    values: list
        Values to be appended

    dictionary: dict
        Dictionary to be appended to
    """
    for i, key in enumerate(keys):
        for row in values:
            dictionary[key].append(row[i])
    return dictionary


def calc_luminosity(redshift, flux_density, wl):
    """
    Takes a redshift and flux density and returns an approximate luminosity.
    Supperceded by HELMS_Wrapper stuff
    """
    distance = c.Distance(z=redshift)
    distance = distance.value
    distance = distance * 3.086 * 10**19
    flux = flux_density * 3.631 * (10 ** -32) * wl
    lum = 4 * np.pi * flux * (distance ** 2)
    return lum

class SDSS_API:
    """
    Basic class for wrapping the SDSS API searches.
    Specific functions are used depending on the aim of the search.
    """

    def __init__(self):
        self.rec_url = 'http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/RectangularSearch'
        self.sql_url = 'http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/SqlSearch'

    def rectangular_search(self, ra_range, dec_range, searchtype='equitorial', limit='10', whichquery='irspectra'):
        """
        ra_range and dec_range should be in tuple form, (lower bound, upper bound)
        Largely a useless function, SQL search is vastly superior in flexibility.
        """
        prams = {'min_ra': str(ra_range[0]), 'max_ra': str(ra_range[1]), 'min_dec': str(dec_range[0]), 'max_dec': str(
            dec_range[1]), 'searchtype': searchtype, 'limit': limit, 'format': 'csv', 'whichquery': whichquery}
        return(requests.get(self.rec_url, prams).text)

    def sql_search(self, cmd, whichquery='irspectra', lim='10'):
        """
        Searches the SDSS SQL database and returns a dictionary of the results under the searched
        headers.

        cmd: str
            The SQL command to be executed. Should be correctly formatted for direct application.

        whichquery: str 
            Some random crap required by the SQL searcher. Included here because it might be important
            but for now basically just leave on default.

        lim: str
            The limit on number of displayed results from the SQL search. Default 10.
        """
        prams = {'cmd': cmd, 'format': 'csv', 'searchtype': 'equitorial',
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
        Returns the luminosity of a specific galaxy in the z band.
        For far IR measurements, this function is superceded by the 
        HELMS_wrapper one.
        """
        freq = 3.282159601489 * 10 ** 14
        command = 'select z, sky_z \n from photoObj \n where ObjID = ' + str(galaxy)
        r = self.sql_search(command)
        z, flux = r['z'][0], r['sky_z'][0]
        lum = calc_luminosity(float(z), float(flux), freq)
        return lum       

    def galaxy_radius_search(self, galaxy, degrees=True):
        """
        Finds the petrosian radius of a given galaxy. Returns
        a float.

        galaxy: int, str
            The galaxy in questions' object ID

        degrees: bool
            If true, returns the results in degrees instead of arcesconds
            (Default true)
        """
        command = 'select s.z, p.petroRad_z \n from photoObj p \n JOIN specObj s ON s.bestObjID = p.objID \n where s.bestObjID = ' + str(galaxy)
        r = self.sql_search(command)
        z, rad = float(r['z'][0]), float(r['petroRad_z'][0])
        if degrees:
            rad = rad * 0.000277778
        return z, rad

    def close_galaxy_search(self, ra, dec, z, interval, z_interval):
        """
        SQL searches for galaxies within a certain distance of a given location.
        Returns a data frame of the closest galaxies including their z, ra, dec 
        and object ID

        ra, dec: float
            Right ascension and declination (deg)

        z: float
            Redshift

        interval: float
            The angular interval (deg) to be searched. ra and dec will be searched
            based on their values ± this interval

        z_interval: float
            The z interval to be searched (z±z_interval)
        """
        command = 'select \n z, ra, dec, bestObjID \n from \n specObj \n where \n class = \'galaxy\'  \n and zWarning = 0 \n and ra<' + \
            str(ra + interval) + '\n and ra>' + str(ra - interval) + '\n and dec<' + \
            str(dec + interval) + '\n and dec>' + str(dec - interval) + '\n and z >' + str(z-z_interval) + '\n and z <' + str(z+z_interval)
        search_results = self.sql_search(command)
        galaxies = pd.DataFrame.from_dict(search_results)
        return galaxies

    def find_closest_galaxy(self, ra, dec, z):
        """
        Finds the closest galaxy for a given position.
        Returns the object ID and angular separation of the closest galaxy.

        ra, dec: float
            Right ascension and declination (deg)

        z: float
            Redshift
        """
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

