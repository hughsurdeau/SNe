from astropy.io import fits
from regions.core import PixCoord
from regions import CirclePixelRegion, CircleSkyRegion
from astropy.units import Quantity
import numpy as np
from astropy.coordinates import Angle, SkyCoord
import astropy.wcs as wcs
from astropy import constants as const
from basic_funcs import *
import astropy.coordinates as c
import os


cwd = os.getcwd()
cat_dir = cwd + '/Catalogs/'
cat_250um = cat_dir + 'HeLMS_PSW_SANEPIC_v0.3.fits'
cat_350um = cat_dir + 'HeLMS_PMW_SANEPIC_v0.3.fits'
cat_500um = cat_dir + 'HeLMS_PLW_SANEPIC_v0.3.fits'
pixel_size = 6.0000012




hdul = fits.open(cat_250um)
img_dat = hdul[1].data
pixel_convert = wcs.WCS(hdul[1])


def get_flux_density(ra, dec, radius, catalog=cat_250um):
    """
    Returns flux density in the unit used for the catalog
    (Jy in the case of the HELMS maps)

    ra: float
        Right ascension of the object (degrees).

    dec: float
        Declination of the object (degrees).

    radius: float
        The radius of the galaxy (degrees).

    catalog: cat_250um, cat_350um, cat_500um
        The catalog used in the search. All
        variables are defined above.

    """
    hdul = fits.open(catalog)
    img_dat = hdul[1].data
    pixel_convert = wcs.WCS(hdul[1])
    center = SkyCoord(ra, dec, unit='deg')
    radius = Angle(radius, 'deg')
    region = CircleSkyRegion(center, radius)
    pixel_region = region.to_pixel(pixel_convert)
    if pixel_region.radius > 500:
        pixel_region.radius = 50
    mask = pixel_region.to_mask(mode='exact')
    data = mask.cutout(img_dat)
    flux_density = max_list_of_lists(data)
    return flux_density

def get_lum(flux_density, redshift, band=250):
    """
    Returns the wavelength specific luminosity given its flux density.
    Returns luminosity in watts.

    flux_density: float
        The flux density of the object at the given wavelength. Assumed
        to be in Jy (HELMS survey is in Jy)

    redshift: float
        The redshift of the object

    band: int
        The wavelength of the flux density. In micrometres.
    
    """
    if not band in [250, 350, 500]:
        raise ValueError("Band must be 250, 350 of 500 microns!")
    freq = const.c / (band * 10 ** -6)
    d = c.Distance(z=redshift)
    distance = d.value
    flux = flux_density * freq
    lum = 4 * np.pi * flux * (distance ** 2)
    return lum.value

def check_if_flux(ra, dec, radius, catalog=cat_250um):
    """
    Returns flux density in the unit used for the catalog
    (Jy in the case of the HELMS maps)
    Returns 'Outside' if the area is outside the image.

    ra: float
        Right ascension of the object (degrees).

    dec: float
        Declination of the object (degrees).

    radius: float
        The radius of the galaxy (degrees).

    catalog: cat_250um, cat_350um, cat_500um
        The catalog used in the search. All
        variables are defined above.

    """
    center = SkyCoord(ra, dec, unit='deg')
    radius = Quantity(radius, 'deg')
    region = CircleSkyRegion(center, radius)
    pixel_region = region.to_pixel(pixel_convert)
    if pixel_region.radius > 500:
        pixel_region.radius = 20
    mask = pixel_region.to_mask(mode='exact')
    data = mask.cutout(img_dat)
    if data is not None:
        flux_density = max_list_of_lists(data)
    else:
        flux_density = 'Outside'
    return flux_density






