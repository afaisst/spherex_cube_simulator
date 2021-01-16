# Copyright (c) 2012-2019 by the GalSim developers team on GitHub
# https://github.com/GalSim-developers
#
# This file is part of GalSim: The modular galaxy image simulation toolkit.
# https://github.com/GalSim-developers/GalSim
#
# GalSim is free software: redistribution and use in source and binary forms,
# with or without modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions, and the disclaimer given in the accompanying LICENSE
#    file.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions, and the disclaimer given in the documentation
#    and/or other materials provided with the distribution.
#
""""
This YAML file creates an HSC image.
"""

import sys
import os
import math
import numpy as np
import logging
import time
import galsim

from astropy.io import fits, ascii
from astropy.table import Table


def main(spherex_filter_name, output_name, params):
    '''
    Main script to run the SPHEREx GalSim simulation.
    USAGE: main(spherex_filter_name , image_output_filename , params)
    where
        - spherex_filter_name: SPHEREx filter name as in the SPHEREx filter file (e.g., spherex_lvf1_m1)
        - image_output_filename: file name of FITS file (note that if the file name is the same, new images are added to the HDU list)
        - params: a dictionary with parameters to set up the simulation (see below)

    OUTPUT: Outputs a FITS file with the simulated image in the HDU list with the name of the SPHEREx filter. It also creates
            "truth" catalogs in the flux units of the image and magnitude. They are saved as [NAME]_[SPHEREx-FILTER].csv
    '''
    
    ## DEFINITIONS ########
    # Define some parameters here
    pixel_scale = params["pixel_scale"]                             # arcsec/pixel
    image_size_arcsec = params["image_size_arcmin"]*60              # image size in arcseconds
    noise_sigma = params["noise_sigma"]                             # ADU  (Just use simple Gaussian noise here. Variance would be the square)
    obj_density_per_arcminsq = params["obj_density_per_arcminsq"]   # number density of galaxies per arcmin squared 
    center_ra = params["center_ra"]*galsim.hours                    # The RA of the center of the image on the sky
    center_dec = params["center_dec"]*galsim.degrees                # The Dec of the center of the image on the sky
    random_seed = params["random_seed"]                             # random seed
    psf_fwhm = params["psf_fwhm"]                                   # PSF FWHM (gaussian) in arcseconds
    image_output_filename = os.path.join( params["output_path"] , output_name )
    #####################

    ## Set up Logger
    #logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logging.basicConfig(format="%(message)s", level=logging.INFO, filename=image_output_filename.replace(".fits",".log"))
    logger = logging.getLogger("SPHEREx_sim for filter %s" % spherex_filter_name)
    logger.info("SPHEREx_sim for filter %s" % spherex_filter_name )

    ## LOAD EXTERNAL CATALOGS ###

    ## Load the central wavelength of the SPHEREx filters
    # name - name of the filter
    # lam - wavelength in Angstroms
    spherex_filters = ascii.read(params["spherex_filter_file"] , names=["name","lam"])
    print("%g filters found in the filter file" % len(spherex_filters))

    ## Load SPHEREx SED catalog matched to GalSim shape catalog
    # col10 - col109 are the 100 SPHEREx filters
    sed_catalog = Table.read(params["spherex_sed_file"])
    print("Number of galaxies in the SED catalog: %g" % len(sed_catalog))

    ###########################

    # size of image in pixels
    image_size = int(image_size_arcsec/pixel_scale)

    # Number of objects
    nobj_tot = int(obj_density_per_arcminsq * image_size_arcsec**2/3600)
    logger.info("Simulating %g objects" % (nobj_tot) )

    # set numpy random seed
    np.random.seed = random_seed

    logger.info('Starting')

    ## Read in galaxy catalog
    # The COSMOSCatalog contains stamps and parametric fits to the COSMOS galaxies. Here
    # down to a magnitude of 25.2AB. (Fill in rest with point sources for example)
    cosmos_cat_name = 'real_galaxy_catalog_25.2.fits'
    cosmos_cat_path = params["galsim_shape_directory"]
    cosmos_cat = galsim.COSMOSCatalog(cosmos_cat_name, dir=cosmos_cat_path , exclusion_level="marginal")
    logger.info('Read in %d galaxies from catalog', cosmos_cat.nobjects)

    ## Load PSF
    # Note that one can define the pixel scale of the PSF. Usually it's the same
    # as the image
    psf = galsim.Gaussian(fwhm=psf_fwhm,flux=1)
    logger.info('Create Gaussian PSF')

    ## Setup the image:
    full_image = galsim.ImageF(image_size, image_size)

    # Set the lower-left corner pixel
    full_image.setOrigin(1,1)

    # Initiate the random number
    rng = galsim.BaseDeviate(random_seed)

    # We keep track of how much noise is already in the image from the RealGalaxies.
    # The default initial value is all pixels = 0.
    noise_image = galsim.ImageF(image_size, image_size)
    noise_image.setOrigin(1,1)

    ## Set up the WCS system
    theta = 0.0 * galsim.degrees
    # ( dudx  dudy ) = ( cos(theta)  -sin(theta) ) * pixel_scale
    # ( dvdx  dvdy )   ( sin(theta)   cos(theta) )
    # Aside: You can call numpy trig functions on Angle objects directly, rather than getting
    #        their values in radians first.  Or, if you prefer, you can write things like
    #        theta.sin() or theta.cos(), which are equivalent.
    dudx = np.cos(theta) * pixel_scale
    dudy = -np.sin(theta) * pixel_scale
    dvdx = np.sin(theta) * pixel_scale
    dvdy = np.cos(theta) * pixel_scale
    affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=full_image.true_center)

    # We can also put it on the celestial sphere to give it a bit more realism.
    # The TAN projection takes a (u,v) coordinate system on a tangent plane and projects
    # that plane onto the sky using a given point as the tangent point.  The tangent
    # point should be given as a CelestialCoord.
    sky_center = galsim.CelestialCoord(ra=center_ra, dec=center_dec)

    # The third parameter, units, defaults to arcsec, but we make it explicit here.
    # It sets the angular units of the (u,v) intermediate coordinate system.
    wcs = galsim.TanWCS(affine, sky_center, units=galsim.arcsec)
    full_image.wcs = wcs

    # prepare truth catalog
    truth_catalog = Table( names=["IDENT","NUMBER","ra","dec","fluxtot","magtot","theta"] ,
                            dtype=[np.int , np.int , np.float64 , np.float64, np.float64, np.float64 , np.float64])

    ## Now we need to loop over our GALAXIES:
    logger.info('Creating galaxies')
    time1 = time.time()
    for k in range(nobj_tot):

        # The usual random number generator using a different seed for each galaxy.
        ud = galsim.UniformDeviate(random_seed+k+1)

        # Choose a random RA, Dec around the sky_center.
        # Note that for this to come out close to a square shape, we need to account for the
        # cos(dec) part of the metric: ds^2 = dr^2 + r^2 d(dec)^2 + r^2 cos^2(dec) d(ra)^2
        # So need to calculate dec first.
        dec = center_dec + (ud()-0.5) * image_size_arcsec * galsim.arcsec
        ra = center_ra + (ud()-0.5) * image_size_arcsec / np.cos(dec) * galsim.arcsec
        world_pos = galsim.CelestialCoord(ra,dec)

        # We will need the image position as well, so use the wcs to get that
        image_pos = wcs.toImage(world_pos)
        
        # pick a galaxy from the COSMOS catalog
        gal = cosmos_cat.makeGalaxy(gal_type='parametric', rng=ud)
        gal_index = gal.index.copy() # this is the galaxy index linked to the cosmos_cat (created by GalSim)
        
        # NOTE: in the following, use gal.index and cosmos_cat.real_cat.KEY to access the parameters in the real catalog. 
        #print( "Flux comparison (catalog / from mag / used): (%g , %g , %g)" % (cosmos_cat.real_cat.stamp_flux[gal_index] , 10**(-0.4*(cosmos_cat.real_cat.mag[gal_index]  - 25.94734)) , gal.flux) )
        
        # Apply a random rotation
        gal_theta = ud()*2.0*np.pi*galsim.radians
        gal = gal.rotate(gal_theta)

        # Rescale the flux of the galaxy ----
        # NOTE: we cannot assign a flux via gal.flux = A, but we can do gal *= B for a factor B.
        # Therefore we have to compute the flux relative to the ACS flux to get it in ACS counts.
        # After that, we have to apply the zeropoint correction.

        # 0) get the flux of this galaxy in this SPHEREx filter
        spherex_filter_id = np.where( spherex_filters["name"] == spherex_filter_name )[0] # from 0 to 99
        if len(spherex_filter_id) == 0:
            logger.info("Filter %s was not found. Abort." % spherex_filter_name)
            quit()
        idx = np.where( sed_catalog["IDENT"] == int(cosmos_cat.real_cat.ident[gal_index]) )[0] # row number of galaxy in SED catalog
        if len(idx) == 0: # if a galaxy does not have an SED, just skip it.
            continue
        gal_SPHEREX_fnu = float(sed_catalog["col%g" % (10 + int(spherex_filter_id)) ][idx]) # TODO: variable with filter
        

        # 1) scale the flux according to current filter.
        gal_ACS_counts_new = gal_SPHEREX_fnu * 10**(-0.4*(-48.6-25.94734)) # converts the SPHEREx f_nu [erg/s/cm2/A] to the ACS counts.
        flux_scaling_sed = gal_ACS_counts_new / gal.flux # flux scaling for SED

        # 2) flux scaling to change zero point from counts to MJy/sr
        flux_scaling_zp = 10**(-0.4*(25.94734 - 8.9) ) / (pixel_scale**2) / 2.350443e-5 # from HST counts to MJy/sr

        # 3) apply the flux scaling
        gal *= flux_scaling_sed # scale flux to reflect SED
        gal *= flux_scaling_zp # scale flux to reflect SPHEREx zeropoint (flux is not in SPHEREx flux units)
        
        # -------

        # Convolve with the PSF.
        final = galsim.Convolve(psf, gal)

        # Account for the fractional part of the position
        # cf. demo9.py for an explanation of this nominal position stuff.
        x_nominal = image_pos.x + 0.5
        y_nominal = image_pos.y + 0.5
        ix_nominal = int(math.floor(x_nominal+0.5))
        iy_nominal = int(math.floor(y_nominal+0.5))
        dx = x_nominal - ix_nominal
        dy = y_nominal - iy_nominal
        offset = galsim.PositionD(dx,dy)

        # We use method='no_pixel' here because the PSF image that we are using includes the
        # pixel response already.
        stamp = final.drawImage(wcs=wcs.local(image_pos), offset=offset, method='no_pixel')

        # Recenter the stamp at the desired position:
        stamp.setCenter(ix_nominal,iy_nominal)

        # Find the overlapping bounds:
        bounds = stamp.bounds & full_image.bounds

        # Finally, add the stamp to the full image.
        full_image[bounds] += stamp[bounds]

        # update truth catalog
        truth_catalog.add_row( [int(cosmos_cat.real_cat.ident[gal_index]) , int(sed_catalog["NUMBER"][idx]), ra.deg ,
                                dec.deg ,
                                gal.flux , 
                                -2.5 * np.log10(gal.flux/flux_scaling_zp) + 25.94734 , 
                                float(gal_theta.deg)] )


    time2 = time.time()
    tot_time = time2-time1
    logger.info('Galaxies created in t=%f s', tot_time)

    ## Now add Gaussian noise with this variance to the final image.  We have to do this step
    # at the end, rather than adding to individual postage stamps, in order to get the noise
    # level right in the overlap regions between postage stamps.
    
    ## If we want to add correlated noise of an image, we can compute this here.
    # For this we would need a noise image. For example a residual image or something
    # similar. For now, we skip correlated noise and just use plain gaussian noise.
    #cn = galsim.CorrelatedNoise(image, rng=rng)
    #full_image.addNoise(cn)
    noise = galsim.GaussianNoise(rng, sigma=noise_sigma)
    full_image.addNoise(noise)
    logger.info('Added noise to final large image')

    ## Write image to disk. If an image is already available, add this image to the HDU.
    # Else create a new image first.
    if os.path.exists(image_output_filename): # load existing image if exists, else create a new image
        with fits.open(image_output_filename) as hdul: 
            galsim.fits.write(full_image , hdu_list=hdul) # add to existing HDUL
            hdul[-1].header["EXTNAME"] = spherex_filter_name # add extension name
            hdul[-1].header["LAMBDA"] = str(float(spherex_filters["lam"][spherex_filter_id])) # add filter lambda
            hdul.writeto(image_output_filename , overwrite=True) # write
    else:
        full_image.write(image_output_filename) # create a new image
        with fits.open(image_output_filename) as hdul:
            hdul[0].header["EXTNAME"] = spherex_filter_name # add extension name
            hdul[0].header["LAMBDA"] = str(float(spherex_filters["lam"][spherex_filter_id])) # add filter lambda
            hdul.writeto(image_output_filename , overwrite=True) # write

    logger.info('Wrote image to %r',image_output_filename)

    ## write truth catalog
    truth_catalog_filename = image_output_filename.replace(".fits" , "_%s.csv" % spherex_filter_name)
    truth_catalog.write(truth_catalog_filename , overwrite=True , format="csv" )
    logger.info('Wrote truth catalog to %r',truth_catalog_filename)

    print("--- ALL DONE! ----")

#####
