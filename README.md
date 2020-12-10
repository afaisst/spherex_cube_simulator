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
main(spherex_filter_name=filt , output_name = "test.fits" , params=params)
```

where

* spherex_filter_name: the SPHEREx filter name as listed in the SPHEREx filter file (e.g., "spherex_lvf1_m1")
* output_name: name of the output FITS image and truth catalogs
* params: a dictionary containing the parameters to create the simulated image (see below)

adfadf

## References

Saito et al., 2020, MNRAS, 464, 199 [https://ui.adsabs.harvard.edu/abs/2020MNRAS.494..199S/abstract]
