import os
import numpy as np
from astropy.io import fits, ascii

def makeCube(filename , save):
    '''
    This function converts a FITS MEF into a FITS cube.
    USAGE: output = makeCube(filename , save)
    where
        - filename: file name of FITS MEF [string]
        - save: If set to "True", then the cube is saved to disk with filename_cube = filename.replace(".fits" , "_cube.fits"). If set to "False", then the cube is directly returned as a numpy 3D array.  [True/False]
    
    OUTPUT: Either "True" if save=True or a numpy 3D array with shape (nbr extensions, x_dim , y_dim).
    '''

    ## Open file
    with fits.open(filename) as hdul:

        for ii,hdu in enumerate(hdul):

            print("%g/%g " % (ii , len(hdul))  , end=" ")

            if ii == 0:
                cube = np.zeros( (1,hdu.data.shape[0],hdu.data.shape[1]) )
                cube[0,:,:] = hdu.data.copy()
                hdr = hdu.header.copy()
            else:
                tmp = np.zeros( (1,hdu.data.shape[0],hdu.data.shape[1]) )
                tmp[0,:,:] = hdu.data.copy()
                cube = np.concatenate( (cube,tmp) )

    ## save cube or return it
    if save:
        hdu_new = fits.PrimaryHDU()
        hdu_new.data = cube.copy()
        hdu_new.header = hdr.copy()
        hdu_new.header["EXTNAME"] = "CUBE"
        del hdu_new.header["LAMBDA"]

        hdul_new = fits.HDUList([hdu_new])
        hdul_new.writeto( filename.replace(".fits","_cube.fits"),  overwrite=True )
        return(True)
    else:
        return(cube)


## EXAMPLE
#makeCube("/Users/afaisst/Work/SphereX/modules/B_Perform_Forced_Photometry/GalSim/spherex_cube_simulator/output/run_all_6arcmin/cube_6arcmin.fits", save=True)