from main import *


## Simulations parameters
params = {"pixel_scale":0.62,               # pixel scale [arcsec/px]
          "image_size_arcmin":1,            # image size [arcminutes]
          "noise_sigma":1e-4,               # sigma, in the units of the image
          "obj_density_per_arcminsq":10,    # number density of galaxies [gal/arcmin2]
          "star_density_per_arcminsq":0,    # number density of stars [stars/arcmins]
          "star_mags_range":[15,19.5],      # magnitude limits for stars [bright , faint] in AB mags 
          "center_ra":10,                   # The RA of the center of the image on the sky [degrees]
          "center_dec":2,                   # The Dec of the center of the image on the sky [degrees]
          "random_seed":1,                  # Random seed
          "psf_fwhm":0.8,                    # PSF FWHM [arcsec]
          "output_path":"../output/",        # Output path to save images and truth catalogs
          "spherex_filter_file":"../../external_catalogs/spherex_lvf_filters_2020Dec2.txt", # Location of SPHEREx filter file
          "spherex_sed_file":"../../external_catalogs/SPHEREx_fluxes_matched_galsim_2020Dec8.fits", # Location of SPHEREx flux file
          "galsim_shape_directory":"/Users/afaisst/Work/Tools/GalSim/cosmos_data/COSMOS_25.2_training_sample/", # Path to GalSim shape directory containing the shape catalog
          "galsim_input_catalog_name":"real_galaxy_catalog_25.2.fits", # Catalog name in galsim_shape_directory from where the galaxies are drawn
          "psf_file_name":"none" # file name for PSF. If that is set, then psf_fwhm is not used. Else, set to "none"
}


## Which filters to run it
spherex_filters = ascii.read("../../external_catalogs/spherex_lvf_filters_2020Dec2.txt" , names=["name","lam"])
#filters = ["spherex_lvf1_m1" , "spherex_lvf1_m3"]
filters = spherex_filters["name"][0:2]

print(filters)

for filt in filters:
    print(filt)
    main(spherex_filter_name=filt , output_name = "test.fits" , params=params)

