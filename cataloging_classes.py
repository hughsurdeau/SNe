from basic_funcs import *
from astropy.table import Table
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

    def get_cid_old(self, ra_key='peak_ra'):
        """
        Returns the list of SNe to get a list of CID if not already known.
        """
        ra = []
        for loc in self.locs:
            if loc[0] > 180:
                neg = str(round(loc[0]-360, 6))
                if len(neg.split('.')[1]) < 6:
                    neg = neg + '0'*(6-len(neg.split('.')[1]))
                ra.append(neg)
            else:
                neg = str(loc[0])
                if len(neg.split('.')[1]) < 6:
                    neg = neg + '0'
                ra.append(neg)
        cid_df = self.dataframe.loc[self.dataframe[ra_key].isin(ra)]
        cid_list = cid_df['cid'].tolist()
        self.cid = cid_list

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

