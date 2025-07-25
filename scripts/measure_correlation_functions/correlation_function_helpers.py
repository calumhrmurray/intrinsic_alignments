import astropy.io.fits as fits
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import treecorr
import matplotlib.pyplot as plt
import healpy as hp
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM( 70 , 0.3 )


def process_ng_rpar_bin( shape_catalogue , 
                         position_catalogue ,
                         random_position_catalogue ,
                         min_rpar, 
                         max_rpar ,
                         min_sep, 
                         max_sep, 
                         nbins , 
                         bin_type = 'Linear', 
                         metric  = 'Rperp'):
    
    print('Running between rpar =', min_rpar, 'and rpar =', max_rpar)

    # Create the NNCorrelation objects
    ng = treecorr.NGCorrelation( min_sep = min_sep, 
                                 max_sep = max_sep, 
                                 nbins = nbins, 
                                 min_rpar = min_rpar, 
                                 max_rpar = max_rpar, 
                                 bin_type = bin_type) 
    rg = treecorr.NGCorrelation( min_sep = min_sep, 
                                 max_sep = max_sep, 
                                 nbins = nbins, 
                                 min_rpar = min_rpar, 
                                 max_rpar = max_rpar, 
                                 bin_type = bin_type) 

    # Process the position and random catalogues
    ng.process( position_catalogue, shape_catalogue , metric = metric )
    rg.process( random_position_catalogue, shape_catalogue, metric = metric )

    # Calculate the Landy-Szalay estimator
    xi_p , xi_x , var_xi = ng.calculateXi( rg = rg )
    r = np.exp(ng.meanlogr)

    return r , xi_p , xi_x, var_xi

def calculate_xi_2d( shape_catalogue , 
                     position_catalogue, 
                     random_position_catalogue,
                     min_rpar,
                     max_rpar,
                     min_sep,
                     max_sep,
                     nbins,
                     bin_type ):
    
    # Initialize lists to store results
    xi_gn_p_results = []
    xi_gn_x_results = []
    xi_var_results = []
    r_results = []

    # Iterate over rpar bins
    for i in range(len(rpar_bins) - 1):

        min_rpar = rpar_bins[i]
        max_rpar = rpar_bins[i + 1]

        r , w_xi_p , w_xi_x, w_var = process_ng_rpar_bin( shape_catalogue , 
                                                          position_catalogue, 
                                                          random_position_catalogue, 
                                                          min_rpar,
                                                          max_rpar,
                                                          min_sep,
                                                          max_sep,
                                                          nbins,
                                                          bin_type )

        # Store the results
        xi_gn_p_results.append( w_xi_p )
        xi_gn_x_results.append( w_xi_x )
        xi_var_results.append( w_var )
        r_results.append(r)

    return np.array( xi_gn_p_results ), np.array( xi_gn_x_results ), np.array( xi_var_results ), np.array( r_results )




