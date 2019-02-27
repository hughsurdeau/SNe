import pandas as pd
import pickle
import math
from basic_funcs import *
from astropy.table import Table
from SDSSWrapper import SDSS as ss
from HELMS_wrapper import *


sne_cat_file = os.getcwd() + '/CSV/test_data.csv'
cleaned_sne_cat_file = os.getcwd() + '/CSV/filtered_data.csv'
cwd = os.getcwd()



class SNeCatalogue:
    """
    A general class for working with the SNe Catalogs.

    self.sn_type: str
        The type of SNe considered in the catalog 

    self.locs: list of tuples
        A list of the locations of the galaxies in (ra, dec)

    self.cid: list
        A list of the CID (SNe ID number) of the SNe

    self.data_file: str
        Location of the SNe catalog
    """
    def __init__(self, sn_type, locs=[], cid=[], data_file=sne_cat_file):
        if not all(isinstance(x, i) for x, i in zip([sn_type, locs, cid, data_file], [str, list, list, str])):
            raise TypeError(
                "sn_type must be a string, locs and cid must be lists!")
        if not (data_file[-4:] == 'fits' or data_file[-3:] == 'csv'):
            raise ValueError("Data file must be in the fits or csv format")
        self.sn_type = sn_type
        self.locs = locs
        self.searcher = ss.SDSS_API()
        self.cid = cid
        self.data_file = data_file

    def convert_units(self, HMS=True, rounding=True, decs=6):
        """
        Converts the coordinates fo the self.locs between HMS and degrees.

        HMS: bool
            True = start coords are in HMS format.

        rounding: bool
            True = outputs rounded

        decs: int
            Number of decimals to round to
        """
        updated_coords = []
        if HMS:
            for coord in self.locs:
                updated_coords.append(
                    HMS2deg(ra=coord[0], dec=coord[1], delim=':', rounding=rounding, decimals=decs))
        else:
            for coord in self.locs:
                updated_coords.append(
                    deg2HMS(ra=coord[0], dec=coord[1], rounding=rounding))
        self.locs = updated_coords

    def generate_dataframe(self):
        """
        Generates a dataframe from the datafile and returns the dataframe.
        File must be in fits/csv format
        """
        if self.data_file[-4:] == 'fits':
            dat = Table.read(self.data_file, format='fits')
            source_df = dat.to_pandas()
            self.dataframe = source_df
            return source_df
        else:
            self.dataframe = pd.read_csv(self.data_file)

    def search_dataframe(self, sne_type=None, set_locs=True):
        """
        Searches the inherent dataframe for a certain type of 
        SNe. 
        Returns a data frame with: CID, locations, type, z and z_err

        sne_type: string
            Determines the type of SNe the function will search for.
            If none, will use the instance's SNe type (self.sn_type)

        set_locs: bool
            If true, the instance's locc and CID (self.locs/self.cid)
            will be set to the locs and CID of the search.
        """
        if sne_type == None:
            sne_type = self.sn_type
        data_keys = ['cid', 'locs', 'type_dr', 'z', 'z_best_err']
        obs_sne = {}
        for key in data_keys:
            obs_sne[key] = []
        for index, row in self.dataframe.iterrows():
            if sne_type in row['type_dr'] and row['z_best_with_BOSS_HELIO'] != -9.0:
                obs_sne['cid'].append(row['cid'])
                obs_sne['type_dr'].append(row['type_dr'])
                obs_sne['z'].append(row['z_best_with_BOSS_HELIO'])
                obs_sne['z_best_err'].append(row['z_best_err_with_BOSS_HELIO'])
                obs_sne['locs'].append((row['peak_ra'], row['peak_decl']))
        observable_sne = pd.DataFrame.from_dict(obs_sne)
        if set_locs:
            self.locs = obs_sne['locs']
            self.cid = obs_sne['cid']
        return observable_sne

    def get_cid(self, ra_key='peak_ra'):
        """
        Returns the CID of the supernovae found in self.locs

        ra_key: str
            The key for right ascension
        """
        zipped_values = list(zip(list(self.dataframe['cid'].values), list(self.dataframe[ra_key].values)))
        locations = [loc[0] for loc in self.locs]
        ra = []
        for loc in locations:
            if loc > 180.0:
                    nn = loc - 360
                    nnn = round(nn, 6)
                    ra.append(nnn)
            else:
                ra.append(loc)
        cids = y_search(ra, zipped_values)
        self.cid = cids
        return cids

    def get_redshifts(self, z_key='z_best_with_BOSS_HELIO', z_err_key='z_best_err_with_BOSS_HELIO'):
        """
        Returns a list of tuples including all key data for SNe
        Columns are: cid, z, z_err, ra, decl, type
        """
        if len(self.cid) == 0:
            self.get_cid()
        cid_df = self.dataframe.loc[self.dataframe['cid'].isin(c for c in self.cid)]
        redshifts = zip(self.cid, cid_df[z_key].tolist(), cid_df[z_err_key].tolist(), cid_df['peak_ra'], cid_df['peak_decl'], cid_df['type_dr'])
        self.redshifts = list(redshifts)

    def make_regions(self, filename):
        """
        Creats a region file of the locations included.
        """
        if not filename[-3:] == 'reg':
            raise ValueError('Filename must be a region file!')
        with open(filename, 'w') as f:
            f.write('# Region file format: DS9 version 4.1 \n')
            f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
            f.write('fk5 \n')
            for loc in self.locs:
                f.write('circle(%s,%s,50")' % (loc[0], loc[1]))
                f.write('\n')

    def import_regions(self, region_file):
        """
        Gets locations from a region file. Skips first 3 lines of file!
        """
        coords = []
        with open(region_file, 'r') as file:
            for i, line in enumerate(file):
                if i >2:
                    splitted = line.split('(')[1]
                    split_2 = splitted.split(')')[0]
                    coord = split_2.split(',')
                    coords.append((coord[0], coord[1]))
        self.locs = coords

    def find_host_galaxy(self, ra, dec, z):
        """
        Returns the ra, dec, radius, objID and angular distance of the SNe's host galaxy
        """
        df, distance = self.searcher.find_closest_galaxy(ra, dec, z)
        gal_z, gal_rad = self.searcher.galaxy_radius_search(int(df['bestObjID']))
        gal_data = [float(df['ra']), float(df['dec']), int(df['bestObjID']), float(df['z']), gal_rad, distance]
        return gal_data

    def glume(self, gal_id, loc, z, rad=None):
        if rad is None:
            redshift, rad = self.searcher.galaxy_radius_search(gal_id)
        flux_d = get_flux_density(loc[0], loc[1], rad)
        lum = get_lum(flux_d, z)
        return lum

    def get_core_metrics(self):
        """
        Returns a summary of all the key metrics needed for analysis
        """
        headers = ['cid', 'type', 'loc', 'gal_lum', 'z', 'rad_dist']
        metrics = {}
        for header in headers:
            metrics[header] = []
        for row in self.redshifts:
            loc = (row[3], row[4])
            print("d")
            try:
                gldat = self.find_host_galaxy(loc[0], loc[1], row[1])
            except(TypeError, ValueError, IndexError) as e:
                continue
            lum = self.glume(gldat[2], (gldat[0], gldat[1]), gldat[3], gldat[4])
            rad_pc = gldat[5]/gldat[4]
            data = [row[0], row[5], loc, lum, (row[1], row[2]), rad_pc]
            metrics = multi_assign(headers, data, metrics)
        self.metrics = pd.DataFrame.from_dict(metrics)
        return metrics



class FullTypeCatalogue(SNeCatalogue):
    """
    A subclass of SNeCatalogue which aims to more easily
    manage the processing of the entire list of a given type
    of SNe.

    Initializes with self.locs and self.cid set to the locations 
    of every SNe candidate of the given type. 
    """

    def __init__(self, sn_type, locs=[], cid=[], data_file= os.getcwd() + '/CSV/test_data.csv'):
        SNeCatalogue.__init__(self, sn_type, locs=[], cid=[], data_file= os.getcwd() + '/CSV/test_data.csv')
        self.generate_dataframe()
        self.type_catalogue = self.search_dataframe()
        self.host_galaxies = ''

    def generate_galaxies(self):
        """
        Returns a dictionary containing SNe host galaxy data.
        Headers: cid, gal_locs, gal_z, sne_z, sne_locs, sep, objID
        """
        data_keys = ['cid', 'gal_locs', 'gal_z', 'sne_z', 'sne_locs', 'sep', 'objID']
        host_galaxies = {}
        for key in data_keys:
            host_galaxies[key] = []
        for index, row in self.type_catalogue.iterrows():
            try:
                print(index/len(self.type_catalogue))
                r, distance = self.searcher.find_closest_galaxy(row['locs'][0], row['locs'][1], row['z'])
                vals = [row['cid'], (float(r['ra']), float(r['dec'])), float(r['z']), row['z'], row['locs'], distance, int(r['bestObjID'])]
                host_galaxies = multi_assign(data_keys, vals, host_galaxies)
            except (ValueError, TypeError) as e:
                print(e)
                pass
        self.host_galaxies = host_galaxies

    def lum_bin(self):
        if self.host_galaxies == '':
            self.generate_galaxies()
        hg = self.host_galaxies
        z_lum = []
        for objid, z, locs in zip(list(hg['objID']), list(hg['gal_z']), list(hg['gal_locs'])):
            try:
                lum = self.glume(objid, locs, z)
                if not math.isnan(lum):
                    z_lum.append((z, lum))
            except (TypeError, IndexError) as e:
                pass
        return z_lum



class Universe_Cataloger(SNeCatalogue):
    """
    A class for collecting data from the universe_map_supernovae catalogue
    """
    def __init__(self, tp='', df='rough_cleaned_ud.csv'):
        SNeCatalogue.__init__(self, sn_type=tp, locs=[], cid=[], data_file=df)
        self.new_database = pd.read_csv(self.data_file)
        self.old_database = pd.read_csv(cleaned_sne_cat_file)

    def clean_data(self):
        """
        Removes any SNe outside of the Helms map
        """
        to_drop = []
        for index, row in self.new_database.iterrows():
            print(index)
            try:
                fd = check_if_flux(row['ra'], row['dec'], 0.01)
                if fd == 'Outside':
                    to_drop.append(index)
            except ValueError:
                print('Fucked')
                to_drop.append(index)
        self.new_database = self.new_database.drop(to_drop)

    def check_duplicates(self, rounding=0.000001):
        """
        Ensures that no duplicates from the main catalogue remain
        Rounds to 5 decimal places to check 
        """
        to_drop = []
        old_ra = sorted([float(r) for r in self.old_database['peak_ra'].tolist()])
        old_dec = sorted([float(r) for r in self.old_database['peak_decl'].tolist()])
        for index, row in self.new_database.iterrows():
            ra = row['ra']
            dec = row['dec']
            for r in old_ra:
                if ra+rounding > r:
                    dif = abs(ra-r)
                    if dif < rounding:
                        for d in old_dec:
                            if dec+rounding > d:
                                dif = abs(ra-r)
                                if dif < rounding:
                                    to_drop.append(index)
                                else:
                                    break
                    else:
                        break
        self.new_database = self.new_database.drop(to_drop)

