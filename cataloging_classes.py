from basic_funcs import *
from astropy.table import Table
from SDSSWrapper import SDSS as ss
import pandas as pd
import os
import pickle


class SNeCatalogue:

    def __init__(self, sn_type, locs=[], cid=[], data_file='/Users/hughsurdeau/Desktop/Imperial/Year 4/Project/code/CSV/test_data.csv'):
        if not all(isinstance(x, i) for x, i in zip([sn_type, locs, cid, data_file], [str, list, list, str])):
            raise TypeError(
                "sn_type must be a string, locs and cid must be lists!")
        if not (data_file[-4:] == 'fits' or data_file[-3:] == 'csv'):
            raise ValueError("Data file must be in the fits or csv format")
        self.sn_type = sn_type
        self.locs = locs
        self.cid = cid
        self.data_file = data_file

    def convert_units(self, HMS=True, rounding=True, decs=6):
        """
        Converts the coordinates fo the SNe between HMS and degrees.
        HMS=True means that the start coords are in HMS format.
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
        File must be in fits/csv format!
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
        SNe. If set_locs=True, the instance's locs and CID will 
        be set as those of the result of the search.
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
        Returns a list of tuples, 1st item is CID, 2nd is z and 3rd is z_err. 
        """
        cid_df = self.dataframe.loc[self.dataframe['cid'].isin(
            c for c in self.cid)]
        redshifts = zip(
            self.cid, cid_df[z_key].tolist(), cid_df[z_err_key].tolist())
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


class FullTypeCatalogue(SNeCatalogue):

    def __init__(self, sn_type, locs=[], cid=[], data_file='/Users/hughsurdeau/Desktop/Imperial/Year 4/Project/code/CSV/test_data.csv'):
        SNeCatalogue.__init__(self, sn_type, locs=[], cid=[], data_file='/Users/hughsurdeau/Desktop/Imperial/Year 4/Project/code/CSV/test_data.csv')
        self.generate_dataframe()
        self.type_catalogue = self.search_dataframe()
        self.searcher = ss.SkySearch()
        self.host_galaxies = ''

    def generate_galaxies(self):
        data_keys = ['cid', 'gal_locs', 'gal_z', 'sne_z', 'sne_locs', 'sep', 'objID']
        host_galaxies = {}
        for key in data_keys:
            host_galaxies[key] = []
        for index, row in self.type_catalogue.iterrows():
            try:
                r, distance = self.searcher.find_closest_galaxy(row['locs'][0], row['locs'][1], row['z'])
                vals = [row['cid'], (float(r['ra']), float(r['dec'])), float(r['z']), row['z'], row['locs'], distance, int(r['bestObjID'])]
                host_galaxies = multi_assign(data_keys, vals, host_galaxies)
            except ValueError:
                pass
        self.host_galaxies = host_galaxies

    def lum_bin(self):
        if self.host_galaxies == '':
            self.generate_galaxies()
        hg = self.host_galaxies
        z_lum = []
        for objid, z in zip(list(hg['objID']), list(hg['gal_z'])):
            lum = self.searcher.get_luminosity(objid)
            z_lum.append((z, lum))
        return z_lum










s = FullTypeCatalogue('SNIc')
print(s.searcher)
print(s.lum_bin())
