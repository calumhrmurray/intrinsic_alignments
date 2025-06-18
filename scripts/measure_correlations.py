import numpy as np
import treecorr     
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import pandas as pd

cosmo = FlatLambdaCDM( 70 , 0.3 )

desi_catalogue_folder = '/n17data/murray/desi_data/DESI/catalogs/'
desi_recon_catalogue_folder = '/n17data/murray/desi_data/DESI/results/catalogs_rec/'
desi_rsd_removed_catalogue_folder = '/n17data/murray/desi_data/DESI/results_rsd_removal_only/catalogs_rec/'

def match_shapes_to_positions( config ):
    """
    Match the shapes to the positions by repeating the shapes for each position.

    Returns the various reconstructed matched catalogues as well.
    """
    shapes = fits.open( config['general']['shape_catalogue_file'] )[1].data

    position_ngc = fits.open( desi_catalogue_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data
    position_sgc = fits.open( desi_catalogue_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    galaxies = np.concatenate(( position_ngc , position_sgc ))

    position_ngc = fits.open( desi_recon_catalogue_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data
    position_sgc = fits.open( desi_recon_catalogue_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    recon_galaxies = np.concatenate(( position_ngc , position_sgc ))

    position_ngc = fits.open( desi_rsd_removed_catalogue_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data
    position_sgc = fits.open( desi_rsd_removed_catalogue_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    rsd_removed_galaxies = np.concatenate(( position_ngc , position_sgc ))

    # Check the length of the combined catalogue
    print('Number of galaxies:' , len(galaxies) )

    # Create SkyCoord objects for matching
    galaxy_coords = coord.SkyCoord(ra=galaxies['RA']*u.degree, dec=galaxies['DEC']*u.degree)
    shape_coords = coord.SkyCoord(ra=shapes['RA']*u.degree, dec=shapes['Dec']*u.degree)

    print('Performing matching of galaxies to shapes...')
    # Perform the matching
    idx, d2d, d3d = galaxy_coords.match_to_catalog_sky(shape_coords)

    # Apply a maximum separation criterion (e.g., 1 arcsecond)
    max_sep = 1 * u.arcsec
    matches = d2d < max_sep

    # Filter the matched catalogues
    matched_galaxy_ra = galaxies['RA'][matches]
    matched_galaxy_dec = galaxies['DEC'][matches]
    matched_galaxy_redshift = galaxies['Z'][matches]
    matched_shape_ra = shapes['RA'][idx[matches]]
    matched_shape_dec = shapes['Dec'][idx[matches]]
    matched_shape_e1 = shapes['e1'][idx[matches]]
    matched_shape_e2 = shapes['e2'][idx[matches]]
    matched_shape_w = shapes['w_iv'][idx[matches]]

    # Create a new FITS table with the matched columns 
    cols = [
        fits.Column(name='RA', format='E', array=matched_galaxy_ra),
        fits.Column(name='Dec', format='E', array=matched_galaxy_dec),
        fits.Column(name='e1', format='E', array=matched_shape_e1),
        fits.Column(name='e2', format='E', array=matched_shape_e2),
        fits.Column(name='w_iv', format='E', array=matched_shape_w),
        fits.Column(name='galaxy_RA', format='E', array=matched_galaxy_ra),
        fits.Column(name='galaxy_Dec', format='E', array=matched_galaxy_dec),
        fits.Column(name='redshift', format='E', array=matched_galaxy_redshift)
    ]

    print('Writing matched catalogues to files...')

    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.writeto('/n17data/murray/desi_data/DESI/shape_catalogs/' + config['general']['position_tracer'] + '_unions_desi_matched_observed.fits', overwrite=True)

    # Filter the matched catalogues
    recon_matched_galaxy_ra = recon_galaxies['RA'][matches]
    recon_matched_galaxy_dec = recon_galaxies['DEC'][matches]
    recon_matched_galaxy_redshift = recon_galaxies['Z'][matches]
    recon_matched_shape_ra = shapes['RA'][idx[matches]]
    recon_matched_shape_dec = shapes['Dec'][idx[matches]]
    recon_matched_shape_e1 = shapes['e1'][idx[matches]]
    recon_matched_shape_e2 = shapes['e2'][idx[matches]]
    recon_matched_shape_w = shapes['w_iv'][idx[matches]]

    # Create a new FITS table with the matched columns
    cols = [
        fits.Column(name='RA', format='E', array=recon_matched_galaxy_ra),
        fits.Column(name='Dec', format='E', array=recon_matched_galaxy_dec),
        fits.Column(name='e1', format='E', array=recon_matched_shape_e1),
        fits.Column(name='e2', format='E', array=recon_matched_shape_e2),
        fits.Column(name='w_iv', format='E', array=recon_matched_shape_w),
        fits.Column(name='galaxy_RA', format='E', array=recon_matched_galaxy_ra),
        fits.Column(name='galaxy_Dec', format='E', array=recon_matched_galaxy_dec),
        fits.Column(name='redshift', format='E', array=recon_matched_galaxy_redshift)
    ]

    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.writeto('/n17data/murray/desi_data/DESI/shape_catalogs/' + config['general']['position_tracer'] + '_unions_desi_matched_reconstructed.fits', overwrite=True)

        # Filter the matched catalogues
    rsd_removed_matched_galaxy_ra = rsd_removed_galaxies['RA'][matches]
    rsd_removed_matched_galaxy_dec = rsd_removed_galaxies['DEC'][matches]
    rsd_removed_matched_galaxy_redshift = rsd_removed_galaxies['Z'][matches]
    rsd_removed_matched_shape_ra = shapes['RA'][idx[matches]]
    rsd_removed_matched_shape_dec = shapes['Dec'][idx[matches]]
    rsd_removed_matched_shape_e1 = shapes['e1'][idx[matches]]
    rsd_removed_matched_shape_e2 = shapes['e2'][idx[matches]]
    rsd_removed_matched_shape_w = shapes['w_iv'][idx[matches]]

    # Create a new FITS table with the matched columns
    cols = [
        fits.Column(name='RA', format='E', array=rsd_removed_matched_galaxy_ra),
        fits.Column(name='Dec', format='E', array=rsd_removed_matched_galaxy_dec),
        fits.Column(name='e1', format='E', array=rsd_removed_matched_shape_e1),
        fits.Column(name='e2', format='E', array=rsd_removed_matched_shape_e2),
        fits.Column(name='w_iv', format='E', array=rsd_removed_matched_shape_w),
        fits.Column(name='galaxy_RA', format='E', array=rsd_removed_matched_galaxy_ra),
        fits.Column(name='galaxy_Dec', format='E', array=rsd_removed_matched_galaxy_dec),
        fits.Column(name='redshift', format='E', array=rsd_removed_matched_galaxy_redshift)
    ]

    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.writeto('/n17data/murray/desi_data/DESI/shape_catalogs/' + config['general']['position_tracer'] + '_unions_desi_matched_rsd_removed.fits', overwrite=True)

    return 

def process_ng_rpar_bin( shape_catalogue , position_catalogue , random_position_catalogue , config , min_rpar , max_rpar ):
    print('Running between rpar =', min_rpar, 'and rpar =', max_rpar)

    # Create the NNCorrelation objects
    ng = treecorr.NGCorrelation(min_sep=float( config['treecorr']['min_rperp']), 
                                max_sep=float(config['treecorr']['max_rperp']), 
                                nbins=int(config['treecorr']['n_rperp_bins']), 
                                min_rpar=min_rpar, 
                                max_rpar=max_rpar, 
                                bin_type=config['treecorr']['bin_type'],
                                bin_slop=float(config['treecorr']['bin_slop']),
                                var_method=config['treecorr']['var_method'])#,
                                #angle_slop = float(config['treecorr']['angle_slop']))

    rg = treecorr.NGCorrelation(min_sep=float(config['treecorr']['min_rperp']), 
                                max_sep=float(config['treecorr']['max_rperp']), 
                                nbins=int(config['treecorr']['n_rperp_bins']), 
                                min_rpar=min_rpar, 
                                max_rpar=max_rpar, 
                                bin_type=config['treecorr']['bin_type'],
                                bin_slop=float(config['treecorr']['bin_slop']),
                                var_method=config['treecorr']['var_method'] )#,
                                #angle_slop = float(config['treecorr']['angle_slop']))

    # Process the position and random catalogues
    ng.process( position_catalogue, shape_catalogue , metric='Rperp')
    rg.process( random_position_catalogue, shape_catalogue, metric='Rperp')

    # Calculate the Landy-Szalay estimator
    xi_p , xi_x , var_xi = ng.calculateXi( rg = rg )
    r = np.exp( ng.meanlogr )

    return r , xi_p , xi_x , var_xi

def calculate_cartesian_coordinates( catalogue , ra_col , dec_col , z_col ):
    """
    Calculate the cartesian coordinates for the given catalogue.
    """
    # Convert RA, Dec, and redshift to Cartesian coordinates
    ra = catalogue[ ra_col ] * np.pi / 180.0  # Convert degrees to radians
    dec = catalogue[ dec_col ] * np.pi / 180.0  # Convert degrees to radians
    redshift = catalogue[ z_col ]

    d = cosmo.comoving_distance( redshift ).value

    # Calculate Cartesian coordinates
    x = d * np.cos(dec) * np.cos(ra)
    y = d * np.cos(dec) * np.sin(ra)
    z = d * np.sin(dec)

    catalogue['x'] = x
    catalogue['y'] = y
    catalogue['z'] = z

    return catalogue

def calculate_correlations( config  ):
    """
    Calculate the correlation between shapes and positions.
    """

    print('Calculating the correlation functions...')

    # load positions, randoms, shapes
    if config['general']['position_type'] == 'observed':
        cat_folder = '/n17data/murray/desi_data/DESI/catalogs/'
    elif config['general']['position_type'] == 'reconstructed':
        cat_folder = '/n17data/murray/desi_data/DESI/results/catalogs_rec/'
    elif config['general']['position_type'] == 'rsd_removed':
        cat_folder = '/n17data/murray/desi_data/DESI/results_rsd_removal_only/catalogs_rec/'

    position_ngc = fits.open( cat_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data
    position_sgc = fits.open( cat_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    positions = np.concatenate(( position_ngc , position_sgc ))

    randoms_ngc = fits.open( cat_folder + config['general']['position_tracer'] + '_NGC_' + config['general']['random_index'] +'_clustering.ran.fits' )[1].data
    randoms_sgc = fits.open( cat_folder + config['general']['position_tracer'] + '_SGC_' + config['general']['random_index'] +'_clustering.ran.fits' )[1].data

    # Combine the DESI north and south random catalogues
    randoms = np.concatenate(( randoms_ngc , randoms_sgc ))

    shapes = fits.open('/n17data/murray/desi_data/DESI/shape_catalogs/' + config['general']['shape_tracer'] + '_unions_desi_matched_' + config['general']['shape_type'] + '.fits')[1].data

    # convert from records -> pd.DataFrame
    shapes = pd.DataFrame.from_records( shapes )
    positions = pd.DataFrame.from_records( positions )
    randoms = pd.DataFrame.from_records( randoms )

    shapes = calculate_cartesian_coordinates( shapes , 'RA' , 'Dec' , 'redshift' )
    positions = calculate_cartesian_coordinates( positions , 'RA' , 'DEC' , 'Z' )
    randoms = calculate_cartesian_coordinates( randoms , 'RA' , 'DEC' , 'Z' )

    print('Calculate the cartesian coordinates for the positions, shapes and randoms...')
    position_catalogue = treecorr.Catalog( x=positions['x'], 
                                           y=positions['y'], 
                                           z=positions['z'],
                                           w=positions['WEIGHT'],
                                           npatch =  float( config['treecorr']['npatch'] ) )

    random_position_catalogue = treecorr.Catalog( x=randoms['x'], 
                                                  y=randoms['y'], 
                                                  z=randoms['z'],
                                                  w=randoms['WEIGHT'],
                                                  patch_centers = position_catalogue.patch_centers )
    
    shape_catalogue = treecorr.Catalog( x=shapes['x'], 
                                        y=shapes['y'], 
                                        z=shapes['z'], 
                                        g1 = shapes['e1'],
                                        g2 = shapes['e2'], 
                                        w=shapes['w_iv'] ,
                                        patch_centers = position_catalogue.patch_centers )
        
    # Initialize lists to store results
    xi_gn_p_results = []
    xi_gn_x_results = []
    xi_gn_var_results = []
    r_results = []

    rpar_bins = np.linspace(
        float(config['treecorr']['min_rpar']),
        float(config['treecorr']['max_rpar']),
        int( config['treecorr']['n_rpar_bins'] )
    )
    
    # Iterate over rpar bins
    for i in range(len(rpar_bins) - 1):

        min_rpar = rpar_bins[i]
        max_rpar = rpar_bins[i + 1]

        r , xi_p , xi_x, var_xi = process_ng_rpar_bin( shape_catalogue , 
                                                       position_catalogue , 
                                                       random_position_catalogue , 
                                                       config ,
                                                       min_rpar,
                                                       max_rpar )

        # Store the results
        xi_gn_p_results.append( xi_p )
        xi_gn_x_results.append( xi_x )
        xi_gn_var_results.append( var_xi )
        r_results.append( r )

    r_results = np.array( r_results )
    xi_gn_p_results = np.array( xi_gn_p_results )
    xi_gn_x_results = np.array( xi_gn_x_results )
    xi_gn_var_results = np.array( xi_gn_var_results )

    # save the correlation functions
    np.save( config['general']['correlation_function_folder'] + config['general']['position_tracer'] + '_xi_gn_p_' + config['general']['shape_type'] + '_' + config['general']['position_type'] + '.npy', xi_gn_p_results)
    np.save( config['general']['correlation_function_folder'] + config['general']['position_tracer'] + '_xi_gn_x_' + config['general']['shape_type'] + '_' + config['general']['position_type'] + '.npy', xi_gn_x_results)
    np.save( config['general']['correlation_function_folder'] + config['general']['position_tracer'] + '_xi_gn_var_' + config['general']['shape_type'] + '_' + config['general']['position_type'] + '.npy', xi_gn_var_results)
    np.save( config['general']['correlation_function_folder'] + config['general']['position_tracer'] + '_r_' + config['general']['shape_type'] + '_' + config['general']['position_type'] + '.npy', r_results)

    return 


def calculate_shifted_position_correlations( rsd_removed_shape_catalogue, 
                                             recon_shape_catalogue,
                                             rsd_removed_position_catalogue , 
                                             recon_position_catalogue , 
                                             random_catalogue,
                                             tau ):
    return