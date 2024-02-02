# DeRedden Toolset
This toolset can be used either as a standalone to be applied to SED data stored in a file, or as a module which can be loaded into a python program. It uses the [Schlafly and Finkbeiner (2011)](https://ui.adsabs.harvard.edu/abs/2011ApJ...737..103S/abstract) estimates of the dust extinction.

## Important info
You must use accurate naming convention to use this toolset correctly.

* For sources: 
    * Names must be resolvable with the [CDS name resolver](https://cds.unistra.fr/cgi-bin/Sesame)
* For filters: 
    * Facility names and filter IDs must match those found on the [SVO filter profile service](http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse). 
        * You need to find the fields titled "Obs. Facility" and "Filter ID" respectively. Note that some relevant filters may be listed under Misc, but Misc will not be the Obs. Facility.
        * Note that sometimes the observatory used to take the data will use another observatory's established filters. ASAS-SN is a [prime example](https://www.zooniverse.org/projects/tharinduj/citizen-asas-sn/about/research) of this.

## To use standalone
Work in an environment that has all the dependencies installed (astroquery, astropy, numpy)

Edit your `config_dereddn.ini` file to use the correct source name, and point it to the data you want to read in (`data_file`: an example file is at `SED_data/B2_1811+31.dat`), and give it the name of the output file (`output_file`: should always be of the ecsv format).

Adapt your new `data_file` to match the example at `SED_data/B2_1811+31.dat`. F(mJy) is the spectral flux density in mJansky. Do not use any other unit, or any other sed type. Note that you do not have to provide any information about the wavelength / frequency of the measurement - this is by design, as that information is redundant with the filter. Error columns should be strictly positive and are the size of the error bars. For facility name and filter ID columns, see **Important info** section above for how to identify the correct names.

Once you have done all this, you can run the code using:
```
python deredden.py
```
which will output the corrected fluxes to a file in the `results` directory, reporting all fluxes in [erg / (cm2 s)].

If this file needs to be read in later, it can be done so quickly and easily with:
```
from astropy.table import QTable

QTable.read('/path/to/file.ecsv', format='ascii.ecsv')
```

## To use as a module
An easy way in Python 3 to load a module (from wherever) into your python file (without fiddling with your $PATH) is to do the following:
```
import importlib.util

spec = importlib.util.spec_from_file_location("deredden", "/path/to/deredden.py")
deredden = importlib.util.module_from_spec(spec)
spec.loader.exec_module(deredden)
```

Now it can easily be used! Simply write down all the information you need and then run the dereddening on it. As an example, the following is for ASAS-SN g-filter data; note, however, that ASAS-SN uses [Sloan g-band filters](https://www.zooniverse.org/projects/tharinduj/citizen-asas-sn/about/research), so that is what I must use to identify the filter system.
```
import astropy.units as u

src_name = 'B2 1811+31'
facility_name = 'Sloan'
filter_ID = 'SLOAN/SDSS.g'

# Produce the S&F A_lambda for the source and filter
A_lambda = deredden.deredden(src_name, facility_name, filter_ID)

# Correct a flux measured by ASAS-SN using the g filter.
## It is essential that you provide the units for this correctly.
fluxdensity = 1.6369274 * u.mJy

## You can correct it and get it back in the same format as before, in this case, spectral flux density (in mJy)
flux_extcorr_Jy = deredden.deredden(src_name, facility_name, filter_ID, val=fluxdensity, convertJy2Eflux=False)

## Or you can correct it and get it back in energy flux [erg / (cm2 s)]
flux_extcorr_eflux = deredden.deredden(src_name, facility_name, filter_ID, val=fluxdensity, convertJy2Eflux=True)
```

One can also access the generic filters, by only specifying the filter_ID and not providing a facility_name:
```
deredden.deredden(src_name, filter_ID='Generic/Cousins.R', val=fluxdensity, convertJy2Eflux=True)
```
