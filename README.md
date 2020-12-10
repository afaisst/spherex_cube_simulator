# SPHEREx GalSim Cube Simulator

## Summary

This script simulates images in the SPHEREx filters using realistic SEDs and galaxy shapes from COSMOS. These cubes can be interpolated at the SPHEREx filter wavelengths and used by the SPHEREx simulator to create realistic LVFs.

The fluxes for each of the GalSim COSMOS galaxies are computed based on the EL-COSMOS SEDs (Saito et al. 2020) by convolving them with 100 SPHEREx filters. 

The script works as a wrapper for GalSim and is very simple: It creates a GalSim image for each of the SPHEREx filters by choosing random galaxies from the GalSim COSMOS shape catalog, matching them to the SPHEREx flux file (containing the same galaxies), adjusting the fluxes, and painting them on an image. All the simulated images for each of the input filters are collected in a FITS file (one HDU for each filter). The HDU names correspond to the name of the filter and also include the central wavelength ("LAMBDA" header keyword). 

## External Data and Requirements

### Galsim

GalSim has to be installed to run this wrapper. It can be downloaded from the GalSim GitHub repository:
https://github.com/GalSim-developers/GalSim

Furthermore, the GalSim COSMOS shape catalog (the "25.2AB" version) has to be downloaded (this can be done in GalSim).

### External flux and filter catalog

Two external catalogs are needed to run the wrapper, which define the SPHEREx filters and contain the photometry for the COSMOS galaxies in the GalSim COSMOS shape catalog:

* __SPHEREx filter file__: A file containing the SPHEREx filters and central wavelengths, which were used to calculate the fluxes for each of the galaxies in the GalSim COSMOS shape catalog. This catalog can be downloaded here: [LINK]

* __SPHEREX SED file__: A catalog containing the fluxes of the GalSim COSMOS shape catalog. This catalog can be downloaded here: [LINK]


## Usage

The cubes are created by calling the function main.py:

```
main(spherex_filter_name, output_name , params)
```

where

* spherex_filter_name: the SPHEREx filter name as listed in the SPHEREx filter file (e.g., "spherex_lvf1_m1")
* output_name: name of the output FITS image and truth catalogs (e.g., "test.fits")
* params: a Python dictionary containing the parameters to create the simulated image (see below)

The simulation parameter dictionary contains information on how to create the simualted image as well as the path to the external catalogs.
Here is an example for the "params" dictionary

```
params = {"pixel_scale":0.62,               # pixel scale [arcsec/px]
          "image_size_arcmin":1,            # image size [arcminutes]
          "noise_sigma":1e-4,               # sigma, in the units of the image
          "obj_density_per_arcminsq":10,    # number density of galaxies [gal/arcmin2] 
          "center_ra":10,                   # The RA of the center of the image on the sky [degrees]
          "center_dec":2,                   # The Dec of the center of the image on the sky [degrees]
          "random_seed":1,                  # Random seed
          "psf_fwhm":0.8,                    # PSF FWHM [arcsec]
          "output_path":"../output/",        # Output path to save images and truth catalogs
          "spherex_filter_file":"../../external_catalogs/spherex_lvf_filters_with_centwave_2020Dec3.txt", # Location of SPHEREx filter file
          "spherex_sed_file":"../../external_catalogs/SPHEREx_fluxes_matched_galsim_2020Dec8.fits", # Location of SPHEREx flux file
          "galsim_shape_directory":'/Users/afaisst/Work/Tools/GalSim/cosmos_data/COSMOS_25.2_training_sample/' # Path to GalSim shape directory containing the shape catalog
}
```

The script called run.py shows how the images for different filters can be created. If the "output_name" is the same, the script adds the new image to the existing images as subsequent HDUs. Each of the HDU is named after the filter (e.g., "spherex_lvf1_m1") and also has the central wavelength of the filter in the header keyword "LAMBDA". In addition, a truth catalog is saved with the namig convention "[output_name]_[SPHEREx filter].csv". The file contains the IDENT (ID used by the GalSim COSMOS shape catalog) as well as the input coordinates (RA/DEC) and fluxes (in units of the image as well as AB magnitude).

Here is the run.py example script:

```
from main import *


## Simulations parameters
params = {"pixel_scale":0.62,               # pixel scale [arcsec/px]
          "image_size_arcmin":1,            # image size [arcminutes]
          "noise_sigma":1e-4,               # sigma, in the units of the image
          "obj_density_per_arcminsq":10,    # number density of galaxies [gal/arcmin2] 
          "center_ra":10,                   # The RA of the center of the image on the sky [degrees]
          "center_dec":2,                   # The Dec of the center of the image on the sky [degrees]
          "random_seed":1,                  # Random seed
          "psf_fwhm":0.8,                    # PSF FWHM [arcsec]
          "output_path":"../output/",        # Output path to save images and truth catalogs
          "spherex_filter_file":"../../external_catalogs/spherex_lvf_filters_with_centwave_2020Dec3.txt", # Location of SPHEREx filter file
          "spherex_sed_file":"../../external_catalogs/SPHEREx_fluxes_matched_galsim_2020Dec8.fits", # Location of SPHEREx flux file
          "galsim_shape_directory":'/Users/afaisst/Work/Tools/GalSim/cosmos_data/COSMOS_25.2_training_sample/' # Path to GalSim shape directory containing the shape catalog
}


## Which filters to run it
spherex_filters = ascii.read("../../external_catalogs/spherex_lvf_filters_with_centwave_2020Dec3.txt" , names=["name","lam"])
filters = spherex_filters["name"][0:2]

for filt in filters:
    main(spherex_filter_name=filt , output_name = "test.fits" , params=params)

```

Running this script would create an image "test.fits" with 2 HDU extensions, each containing one filter. In addition, truth catalogs are saved with the naming convension "[output_name]_[SPHEREx filter].csv". The file contains the IDENT (ID used by the GalSim COSMOS shape catalog) as well as the input coordinates (RA/DEC) and fluxes (in units of the image as well as AB magnitude).




## References

Saito et al., 2020, MNRAS, 464, 199 [https://ui.adsabs.harvard.edu/abs/2020MNRAS.494..199S/abstract]
