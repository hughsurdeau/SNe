def deg2HMS(ra='', dec='', round=False):
    """
    Converts a degree value to an HMS value. Can return right ascension and declenation,
    or just one.
    """
    RA, DEC, rs, ds = '', '', '', ''
    if dec:
        if str(dec)[0] == '-':
            ds, dec = '-', abs(dec)
        deg = int(dec)
        decM = abs(int((dec - deg) * 60))
        if round:
            decS = int((abs((dec - deg) * 60) - decM) * 60)
        else:
            decS = (abs((dec - deg) * 60) - decM) * 60
        DEC = '{0}{1} {2} {3}'.format(ds, deg, decM, decS)

    if ra:
        if str(ra)[0] == '-':
            rs, ra = '-', abs(ra)
        raH = int(ra / 15)
        raM = int(((ra / 15) - raH) * 60)
        if round:
            raS = int(((((ra / 15) - raH) * 60) - raM) * 60)
        else:
            raS = ((((ra / 15) - raH) * 60) - raM) * 60
        RA = '{0}{1} {2} {3}'.format(rs, raH, raM, raS)

    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC

def HMS2deg(ra='', dec='', delim=':', rounding=True, decimals=6):
    """
    Converts an HMS value to degrees. Can return right ascension and declenation,
    or just one.
    """
    if dec:
        DEC = 0.0
        dec_components = dec.split(delim)
        DEC += abs(float(dec_components[0])) * 15.0
        DEC += float(dec_components[1]) * 0.25
        DEC += float(dec_components[2]) * (0.25 / 60)
        if dec_components[0][0] == '-':
            DEC *= -1
        if rounding:
            DEC = round(DEC, decimals)
    if ra:
        RA = 0.0
        ra_components = ra.split(delim)
        RA += abs(float(ra_components[0])) * 15.0
        RA += float(ra_components[1]) * 0.25
        RA += float(ra_components[2]) * (0.25 / 60)
        if ra_components[0][0] == '-':
            RA *= -1
        if rounding:
            RA = round(RA, decimals)

    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC

def y_search(search_list, tuple_list):
    y = []
    for item in tuple_list:
        if round(item[1], 6) in search_list:
            y.append(item[0])
    return y

def multi_assign(keys, values, dictionary):
    """
    Appends values to a dictionary consisting of only lists
    """
    for key, value in zip(keys,values):
        dictionary[key].append(value)
    return dictionary

def z_splitter(bins, data, z_key=0):
    """
    Splits data into bins based on their redshift
    Returns a dictionary of the binned data

    bins: int or list
        Determines the bins the data will be split into. 
        If int - splits the redshifts evenly into the number
        of bins
        If list - splits the redshifts into bins based on the
        given list

    data: list
        A list of tuples which contains the data that will be 
        split

    z_key: int
        Containt the tuple key of the redshift values within
        data
        Default = 0

    """
    bans = {}
    if type(bins) == int:
        z_max = max([row[z_key] for row in data])
        interval = z_max/bins
        bons = []
        for i in range(bins):
            bons.append((i+1)*interval)
            bans[(i+1)*interval] = []
        bins = bons
    else:
        for b in bins:
            bans[b] = []
    for row in data:
        for b in bons:
            if row[z_key] < b:
                bans[b].append(row)
                break
    return bans

def max_list_of_lists(lists):
    maxs = []
    for l in lists:
        maxs.append(max(l))
    grand_max = max(maxs)
    return grand_max
