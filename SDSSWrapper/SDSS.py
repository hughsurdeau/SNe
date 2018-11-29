import requests
import io
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord


class SkySearch:
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
        # theres some random crap at the start I have to cut off for use as a data frame
        return r.content[8:]

    def close_galaxy_search(self, ra, dec, interval):
        command = 'select \n z, ra, dec, bestObjID \n from \n specObj \n where \n class = \'galaxy\'  \n and zWarning = 0 \n and ra<' + \
            str(ra + interval) + '\n and ra>' + str(ra - interval) + '\n and dec<' + \
            str(dec + interval) + '\n and dec>' + str(dec - interval)
        csv = self.sql_search(command)
        galaxies = pd.read_csv(io.StringIO(csv.decode('utf-8')))
        return galaxies

    def find_closest_galaxy(self, ra, dec):
        if ra < 0:
            ra += 360
        c_sne = SkyCoord(ra, dec, unit='deg')
        near_galaxies = self.close_galaxy_search(ra, dec, interval=0.5)
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
        return near_galaxies.loc[near_galaxies['bestObjID'] == closest_galaxy], distance


