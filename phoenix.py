from astropy.utils.data import download_file
from astropy.io import fits
import numpy as np

phoenix_model_temps = np.array(
    [2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100,
     3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
     4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900,
     5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800,
     5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700,
     6800, 6900, 7000, 7200, 7400, 7600, 7800, 8000, 8200,
     8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000,
     10200, 10400, 10600, 10800, 11000, 11200, 11400, 11600, 11800,
     12000, 12500, 13000, 13500, 14000, 14500, 15000])


def get_any_metallicity_url(T_eff, log_g, z):
    closest_grid_temperature = phoenix_model_temps[np.argmin(np.abs(phoenix_model_temps - T_eff))]

    if z > 0:
        z = "+{0:1.1f}".format(z)
    elif z == 0:
        z = "-{0:1.1f}".format(z)
    else:
        z = "{0:1.1f}".format(z)

    url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
           'PHOENIX-ACES-AGSS-COND-2011/Z{z}/lte{T_eff:05d}-'
           '{log_g:1.2f}{z}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
           ).format(T_eff=closest_grid_temperature, log_g=log_g, z=z)
    return url


def get_phoenix_model_wavelengths(cache=True):
    """
    Return the wavelength grid that the PHOENIX models were computed on,
    transformed into wavelength units in air (not vacuum).
    """
    wavelength_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/'
                      'HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
    wavelength_path = download_file(wavelength_url, cache=cache, timeout=30)
    wavelengths_vacuum = fits.getdata(wavelength_path)

    # Wavelengths are provided at vacuum wavelengths. For ground-based
    # observations convert this to wavelengths in air, as described in
    # Husser 2013, Eqns. 8-10:
    sigma_2 = (10**4 / wavelengths_vacuum)**2
    f = (1.0 + 0.05792105/(238.0185 - sigma_2) + 0.00167917 /
         (57.362 - sigma_2))
    wavelengths_air = wavelengths_vacuum / f
    return wavelengths_air