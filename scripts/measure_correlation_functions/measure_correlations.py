import numpy as np
import treecorr     
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import pandas as pd
from scipy.spatial import cKDTree

cosmo = FlatLambdaCDM( 70 , 0.3 )

desi_catalogue_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/catalogs/'
desi_recon_catalogue_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results/catalogs_rec/'
desi_rsd_removed_catalogue_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results_rsd_removal_only/catalogs_rec/'

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
    # Add only the requested columns: mag, flux_radius, snr, w_des, e1_uncal, e2_uncal
    matched_shape_mag = shapes['mag'][idx[matches]]
    matched_shape_flux_radius = shapes['FLUX_RADIUS'][idx[matches]]
    matched_shape_snr = shapes['snr'][idx[matches]]
    matched_shape_w_des = shapes['w_des'][idx[matches]]
    matched_shape_e1_uncal = shapes['e1_uncal'][idx[matches]]
    matched_shape_e2_uncal = shapes['e2_uncal'][idx[matches]]

    # Create a new FITS table with the matched columns, including extra shape properties
    cols = [
        fits.Column(name='RA', format='E', array=matched_galaxy_ra),
        fits.Column(name='Dec', format='E', array=matched_galaxy_dec),
        fits.Column(name='e1', format='E', array=matched_shape_e1),
        fits.Column(name='e2', format='E', array=matched_shape_e2),
        fits.Column(name='w_iv', format='E', array=matched_shape_w),
        fits.Column(name='galaxy_RA', format='E', array=matched_galaxy_ra),
        fits.Column(name='galaxy_Dec', format='E', array=matched_galaxy_dec),
        fits.Column(name='redshift', format='E', array=matched_galaxy_redshift),
        fits.Column(name='mag', format='E', array=matched_shape_mag),
        fits.Column(name='flux_radius', format='E', array=matched_shape_flux_radius),
        fits.Column(name='snr', format='E', array=matched_shape_snr),
        fits.Column(name='w_des', format='E', array=matched_shape_w_des),
        fits.Column(name='e1_uncal', format='E', array=matched_shape_e1_uncal),
        fits.Column(name='e2_uncal', format='E', array=matched_shape_e2_uncal)
    ]

    print('Writing matched catalogues to files...')

    # Write observed matched catalogue with extra shape columns
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.writeto('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['position_tracer'] + '_unions_desi_matched_observed.fits', overwrite=True)

    # Filter the matched catalogues for reconstructed
    recon_matched_galaxy_ra = recon_galaxies['RA'][matches]
    recon_matched_galaxy_dec = recon_galaxies['DEC'][matches]
    recon_matched_galaxy_redshift = recon_galaxies['Z'][matches]
    recon_matched_shape_ra = shapes['RA'][idx[matches]]
    recon_matched_shape_dec = shapes['Dec'][idx[matches]]
    recon_matched_shape_e1 = shapes['e1'][idx[matches]]
    recon_matched_shape_e2 = shapes['e2'][idx[matches]]
    recon_matched_shape_w = shapes['w_iv'][idx[matches]]
    recon_matched_shape_mag = shapes['mag'][idx[matches]]
    recon_matched_shape_flux_radius = shapes['FLUX_RADIUS'][idx[matches]]
    recon_matched_shape_snr = shapes['snr'][idx[matches]]
    recon_matched_shape_w_des = shapes['w_des'][idx[matches]]
    recon_matched_shape_e1_uncal = shapes['e1_uncal'][idx[matches]]
    recon_matched_shape_e2_uncal = shapes['e2_uncal'][idx[matches]]

    # Create a new FITS table with the matched columns and extra shape properties for reconstructed
    cols_recon = [
        fits.Column(name='RA', format='E', array=recon_matched_galaxy_ra),
        fits.Column(name='Dec', format='E', array=recon_matched_galaxy_dec),
        fits.Column(name='e1', format='E', array=recon_matched_shape_e1),
        fits.Column(name='e2', format='E', array=recon_matched_shape_e2),
        fits.Column(name='w_iv', format='E', array=recon_matched_shape_w),
        fits.Column(name='galaxy_RA', format='E', array=recon_matched_galaxy_ra),
        fits.Column(name='galaxy_Dec', format='E', array=recon_matched_galaxy_dec),
        fits.Column(name='redshift', format='E', array=recon_matched_galaxy_redshift),
        fits.Column(name='mag', format='E', array=recon_matched_shape_mag),
        fits.Column(name='flux_radius', format='E', array=recon_matched_shape_flux_radius),
        fits.Column(name='snr', format='E', array=recon_matched_shape_snr),
        fits.Column(name='w_des', format='E', array=recon_matched_shape_w_des),
        fits.Column(name='e1_uncal', format='E', array=recon_matched_shape_e1_uncal),
        fits.Column(name='e2_uncal', format='E', array=recon_matched_shape_e2_uncal)
    ]

    hdu = fits.BinTableHDU.from_columns(cols_recon)
    hdu.writeto('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['position_tracer'] + '_unions_desi_matched_reconstructed.fits', overwrite=True)

    # Filter the matched catalogues for rsd_removed
    rsd_removed_matched_galaxy_ra = rsd_removed_galaxies['RA'][matches]
    rsd_removed_matched_galaxy_dec = rsd_removed_galaxies['DEC'][matches]
    rsd_removed_matched_galaxy_redshift = rsd_removed_galaxies['Z'][matches]
    rsd_removed_matched_shape_ra = shapes['RA'][idx[matches]]
    rsd_removed_matched_shape_dec = shapes['Dec'][idx[matches]]
    rsd_removed_matched_shape_e1 = shapes['e1'][idx[matches]]
    rsd_removed_matched_shape_e2 = shapes['e2'][idx[matches]]
    rsd_removed_matched_shape_w = shapes['w_iv'][idx[matches]]
    rsd_removed_matched_shape_mag = shapes['mag'][idx[matches]]
    rsd_removed_matched_shape_flux_radius = shapes['FLUX_RADIUS'][idx[matches]]
    rsd_removed_matched_shape_snr = shapes['snr'][idx[matches]]
    rsd_removed_matched_shape_w_des = shapes['w_des'][idx[matches]]
    rsd_removed_matched_shape_e1_uncal = shapes['e1_uncal'][idx[matches]]
    rsd_removed_matched_shape_e2_uncal = shapes['e2_uncal'][idx[matches]]

    # Create a new FITS table with the matched columns and extra shape properties for rsd_removed
    cols_rsd = [
        fits.Column(name='RA', format='E', array=rsd_removed_matched_galaxy_ra),
        fits.Column(name='Dec', format='E', array=rsd_removed_matched_galaxy_dec),
        fits.Column(name='e1', format='E', array=rsd_removed_matched_shape_e1),
        fits.Column(name='e2', format='E', array=rsd_removed_matched_shape_e2),
        fits.Column(name='w_iv', format='E', array=rsd_removed_matched_shape_w),
        fits.Column(name='galaxy_RA', format='E', array=rsd_removed_matched_galaxy_ra),
        fits.Column(name='galaxy_Dec', format='E', array=rsd_removed_matched_galaxy_dec),
        fits.Column(name='redshift', format='E', array=rsd_removed_matched_galaxy_redshift),
        fits.Column(name='mag', format='E', array=rsd_removed_matched_shape_mag),
        fits.Column(name='flux_radius', format='E', array=rsd_removed_matched_shape_flux_radius),
        fits.Column(name='snr', format='E', array=rsd_removed_matched_shape_snr),
        fits.Column(name='w_des', format='E', array=rsd_removed_matched_shape_w_des),
        fits.Column(name='e1_uncal', format='E', array=rsd_removed_matched_shape_e1_uncal),
        fits.Column(name='e2_uncal', format='E', array=rsd_removed_matched_shape_e2_uncal)
    ]

    hdu = fits.BinTableHDU.from_columns(cols_rsd)
    hdu.writeto('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['position_tracer'] + '_unions_desi_matched_rsd_removed.fits', overwrite=True)

    return 

desi_catalogue_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/catalogs/'
desi_recon_catalogue_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results/catalogs_rec/'
desi_rsd_removed_catalogue_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results_rsd_removal_only/catalogs_rec/'

def calculate_displacement_vectors( config ):
    """
    Calculate the displacement vectors between the RSD-removed and reconstructed catalogues
    for both positions and shapes. Saves the displacement arrays to disk.
    """
    recon_position_ngc = fits.open( desi_recon_catalogue_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data
    recon_position_sgc = fits.open( desi_recon_catalogue_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    recon_positions = np.concatenate(( recon_position_ngc , recon_position_sgc ))

    rsd_removed_position_ngc = fits.open( desi_rsd_removed_catalogue_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data
    rsd_removed_position_sgc = fits.open( desi_rsd_removed_catalogue_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    rsd_removed_positions = np.concatenate(( rsd_removed_position_ngc , rsd_removed_position_sgc ))

    # calculate cartesian coordinates
    recon_positions = pd.DataFrame.from_records( recon_positions )
    rsd_removed_positions = pd.DataFrame.from_records( rsd_removed_positions )

    recon_positions = calculate_cartesian_coordinates( recon_positions , 'RA' , 'DEC' , 'Z' )
    rsd_removed_positions = calculate_cartesian_coordinates( rsd_removed_positions , 'RA' , 'DEC' , 'Z' )

    r_recon_positions = np.array( [ recon_positions['x'] , recon_positions['y'] , recon_positions['z'] ])
    r_rsd_removed_positions = np.array( [ rsd_removed_positions['x'] , rsd_removed_positions['y'] , rsd_removed_positions['z'] ])

    # Calculate the displacement vector in Cartesian coordinates
    dr = r_recon_positions - r_rsd_removed_positions  # shape: (3, N)

    # Get RA and Dec of reconstructed positions in radians
    ra_rad = np.deg2rad(recon_positions['RA'])
    dec_rad = np.deg2rad(recon_positions['DEC'])

    # Define basis vectors at each reconstructed position
    los = r_recon_positions / np.linalg.norm(r_recon_positions, axis=0)  # shape: (3, N)
    east = np.array([-np.sin(ra_rad), np.cos(ra_rad), np.zeros_like(ra_rad)])  # shape: (3, N)
    north = np.array([
        -np.sin(dec_rad) * np.cos(ra_rad),
        -np.sin(dec_rad) * np.sin(ra_rad),
        np.cos(dec_rad)
    ])  # shape: (3, N)

    # Project dr onto each basis vector
    drp = np.sum(dr * los, axis=0)    # shape: (N,)
    dr1 = np.sum(dr * east, axis=0)  # shape: (N,)
    dr2 = np.sum(dr * north, axis=0)  # shape: (N,)

    # Create a new FITS table with RA, Dec, redshift, dRA, dDec, dZ, and weights
    displacement_cols = [
        fits.Column(name='RA', format='E', array=recon_positions['RA']),
        fits.Column(name='Dec', format='E', array=recon_positions['DEC']),
        fits.Column(name='redshift', format='E', array=recon_positions['Z']),
        fits.Column(name='rsd_removed_RA', format='E', array=rsd_removed_positions['RA']),
        fits.Column(name='rsd_removed_Dec', format='E', array=rsd_removed_positions['DEC']),
        fits.Column(name='rsd_removed_redshift', format='E', array=rsd_removed_positions['Z']),
        fits.Column(name='dr1', format='E', array=dr1),
        fits.Column(name='dr2', format='E', array=dr2),
        fits.Column(name='drp', format='E', array=drp),
        fits.Column(name='WEIGHT', format='E', array=recon_positions['WEIGHT']),
        fits.Column(name='rsd_removed_WEIGHT', format='E', array=rsd_removed_positions['WEIGHT'])
    ]
    displacement_hdu = fits.BinTableHDU.from_columns(displacement_cols)
    displacement_hdu.writeto( config['general']['velocity_catalogue_folder'] + config['general']['position_tracer'] + '_displacements.fits', overwrite=True)

    recon_shapes = fits.open('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['shape_tracer'] + '_unions_desi_matched_reconstructed.fits')[1].data
    rsd_removed_shapes = fits.open('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['shape_tracer'] + '_unions_desi_matched_rsd_removed.fits')[1].data

    # calculate cartesian coordinates for shapes
    recon_shapes = pd.DataFrame.from_records( recon_shapes )
    rsd_removed_shapes = pd.DataFrame.from_records( rsd_removed_shapes )

    recon_shapes = calculate_cartesian_coordinates( recon_shapes , 'RA' , 'Dec' , 'redshift' )
    rsd_removed_shapes = calculate_cartesian_coordinates( rsd_removed_shapes , 'RA' , 'Dec' , 'redshift' )

    r_recon_shapes = np.array( [ recon_shapes['x'] , recon_shapes['y'] , recon_shapes['z'] ])
    r_rsd_removed_shapes = np.array( [ rsd_removed_shapes['x'] , rsd_removed_shapes['y'] , rsd_removed_shapes['z'] ])

    # Calculate the displacement vector in Cartesian coordinates for shapes 
    dr_shapes = r_recon_shapes - r_rsd_removed_shapes  # shape: (3, N)

    # Define basis vectors at each reconstructed shape position
    ra_rad_shapes = np.deg2rad(recon_shapes['RA'])  
    dec_rad_shapes = np.deg2rad(recon_shapes['Dec'])

    los_shapes = r_recon_shapes / np.linalg.norm(r_recon_shapes, axis=0)  # shape: (3, N)
    east_shapes = np.array([-np.sin(ra_rad_shapes), np.cos(ra_rad_shapes), np.zeros_like(ra_rad_shapes)])  # shape: (3, N)
    north_shapes = np.array([
        -np.sin(dec_rad_shapes) * np.cos(ra_rad_shapes),
        -np.sin(dec_rad_shapes) * np.sin(ra_rad_shapes),
        np.cos(dec_rad_shapes)
    ])  # shape: (3, N)

    # Project dr_shapes onto each basis vector
    drp_shapes = np.sum(dr_shapes * los_shapes, axis=0)    # shape: (N,)
    dr1_shapes = np.sum(dr_shapes * east_shapes, axis=0)    # shape: (N,)
    dr2_shapes = np.sum(dr_shapes * north_shapes, axis=0)   # shape: (N,)

    # Create a new FITS table for shape displacements
    shape_displacement_cols = [
        fits.Column(name='RA', format='E', array=recon_shapes['RA']),
        fits.Column(name='Dec', format='E', array=recon_shapes['Dec']),
        fits.Column(name='redshift', format='E', array=recon_shapes['redshift']),
        fits.Column(name='rsd_removed_RA', format='E', array=rsd_removed_shapes['RA']),
        fits.Column(name='rsd_removed_Dec', format='E', array=rsd_removed_shapes['Dec']),
        fits.Column(name='rsd_removed_redshift', format='E', array=rsd_removed_shapes['redshift']),
        fits.Column(name='e1', format='E', array=recon_shapes['e1']),
        fits.Column(name='e2', format='E', array=recon_shapes['e2']),
        fits.Column(name='w_iv', format='E', array=recon_shapes['w_iv']),
        fits.Column(name='dr1', format='E', array=dr1_shapes),
        fits.Column(name='dr2', format='E', array=dr2_shapes),
        fits.Column(name='drp', format='E', array=drp_shapes),
    ]
    shape_displacement_hdu = fits.BinTableHDU.from_columns(shape_displacement_cols)
    shape_displacement_hdu.writeto(
        config['general']['velocity_catalogue_folder'] + config['general']['shape_tracer'] + '_shape_displacements.fits',
        overwrite=True
    )

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

def process_vg_rpar_bin( displacement_catalogue , shape_catalogue, config, min_rpar, max_rpar):
    print('Running between rpar =', min_rpar, 'and rpar =', max_rpar)

    # Use the VGCorrelation for spin-1 (vector) and spin-2 (shear) cross-correlation
    vg = treecorr.VGCorrelation(
        min_sep=float(config['treecorr']['min_rperp']),
        max_sep=float(config['treecorr']['max_rperp']),
        nbins=int(config['treecorr']['n_rperp_bins']),
        min_rpar=min_rpar,
        max_rpar=max_rpar,
        bin_type=config['treecorr']['bin_type'],
        bin_slop=float(config['treecorr']['bin_slop']),
        var_method=config['treecorr']['var_method']
    )   

    # Process the spin-1 and shape (spin-2) catalogues
    vg.process( displacement_catalogue , shape_catalogue, metric='Rperp')

    xi_p, xi_x, var_xi = vg.calculateXi()
    r = np.exp(vg.meanlogr)

    return r, xi_p, xi_x, var_xi

def process_nv_rpar_bin( position_catalogue , displacement_catalogue, config, min_rpar, max_rpar):
    print('Running between rpar =', min_rpar, 'and rpar =', max_rpar)

    # Use the VGCorrelation for spin-1 (vector) and spin-2 (shear) cross-correlation
    nv = treecorr.NVCorrelation(
        min_sep=float(config['treecorr']['min_rperp']),
        max_sep=float(config['treecorr']['max_rperp']),
        nbins=int(config['treecorr']['n_rperp_bins']),
        min_rpar=min_rpar,
        max_rpar=max_rpar,
        bin_type=config['treecorr']['bin_type'],
        bin_slop=float(config['treecorr']['bin_slop']),
        var_method=config['treecorr']['var_method']
    )   

    # Process the spin-1 and shape (spin-2) catalogues
    nv.process( position_catalogue , displacement_catalogue, metric='Rperp')

    return nv.xi , nv.xi_im , nv.varxi

def process_nr_rpar_bin( position_catalogue , size_catalogue, config, min_rpar, max_rpar):
    print('Running between rpar =', min_rpar, 'and rpar =', max_rpar)

    # Use the VGCorrelation for spin-1 (vector) and spin-2 (shear) cross-correlation
    nr = treecorr.NKCorrelation(
        min_sep=float(config['treecorr']['min_rperp']),
        max_sep=float(config['treecorr']['max_rperp']),
        nbins=int(config['treecorr']['n_rperp_bins']),
        min_rpar=min_rpar,
        max_rpar=max_rpar,
        bin_type=config['treecorr']['bin_type'],
        bin_slop=float(config['treecorr']['bin_slop']),
        var_method=config['treecorr']['var_method']
    )   

    nr.process( position_catalogue , size_catalogue, metric='Rperp')

    return nr.xi , nr.varxi

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

def calculate_correlations(config):
    """
    Calculate the correlation between shapes and positions, with optional ellipticity filtering.
    """

    print('Calculating the correlation functions...')

    # load positions, randoms, shapes
    if config['general']['position_type'] == 'observed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/catalogs/'
    elif config['general']['position_type'] == 'reconstructed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results/catalogs_rec/'
    elif config['general']['position_type'] == 'rsd_removed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results_rsd_removal_only/catalogs_rec/'

    position_ngc = fits.open(cat_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits')[1].data
    position_sgc = fits.open(cat_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits')[1].data

    # Combine the DESI north and south galaxy catalogues
    positions = np.concatenate((position_ngc, position_sgc))

    randoms_ngc = fits.open(cat_folder + config['general']['position_tracer'] + '_NGC_' + config['general']['random_index'] + '_clustering.ran.fits')[1].data
    randoms_sgc = fits.open(cat_folder + config['general']['position_tracer'] + '_SGC_' + config['general']['random_index'] + '_clustering.ran.fits')[1].data

    # Combine the DESI north and south random catalogues
    randoms = np.concatenate((randoms_ngc, randoms_sgc))

    shapes = fits.open('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['shape_tracer'] + '_unions_desi_matched_' + config['general']['shape_type'] + '.fits')[1].data

    # convert from records -> pd.DataFrame
    shapes = pd.DataFrame.from_records(shapes)
    positions = pd.DataFrame.from_records(positions)
    randoms = pd.DataFrame.from_records(randoms)

    # Optionally filter shapes by ellipticity
    e_filter_on = config.has_section('ellipticity_filter') and config.getboolean('ellipticity_filter', 'enabled')
    e_filter_str = ""
    if e_filter_on:
        e_min = config.getfloat('ellipticity_filter', 'e_min')
        e_max = config.getfloat('ellipticity_filter', 'e_max')
        print( e_min , e_max )
        e = np.sqrt(shapes['e1']**2 + shapes['e2']**2)
        mask = (e >= e_min) & (e <= e_max)
        print(f'Applying ellipticity filter: {e_min:.2f} <= e <= {e_max:.2f}')
        print(f'Using {np.sum(mask)} shapes out of {len(shapes)} total shapes after filtering.')
        #shapes = shapes[mask].reset_index(drop=True)
        e_filter_str = f"_e{e_min:.2f}-{e_max:.2f}"

    shapes = calculate_cartesian_coordinates(shapes, 'RA', 'Dec', 'redshift')
    positions = calculate_cartesian_coordinates(positions, 'RA', 'DEC', 'Z')
    randoms = calculate_cartesian_coordinates(randoms, 'RA', 'DEC', 'Z')

    print('Calculate the cartesian coordinates for the positions, shapes and randoms...')
    position_catalogue = treecorr.Catalog(
        x=positions['x'],
        y=positions['y'],
        z=positions['z'],
        w=positions['WEIGHT'],
        npatch=float(config['treecorr']['npatch'])
    )

    random_position_catalogue = treecorr.Catalog(
        x=randoms['x'],
        y=randoms['y'],
        z=randoms['z'],
        w=randoms['WEIGHT'],
        patch_centers=position_catalogue.patch_centers
    )

    shape_catalogue = treecorr.Catalog(
        x=shapes['x'][mask] if e_filter_on else shapes['x'],
        y=shapes['y'][mask] if e_filter_on else shapes['y'],
        z=shapes['z'][mask] if e_filter_on else shapes['z'],
        g1=shapes['e1'][mask] if e_filter_on else shapes['e1'],
        g2=shapes['e2'][mask] if e_filter_on else shapes['e2'],
        w=shapes['w_iv'][mask] if e_filter_on else shapes['w_iv'],
        patch_centers=position_catalogue.patch_centers
    )

    # Initialize lists to store results
    xi_gn_p_results = []
    xi_gn_x_results = []
    xi_gn_var_results = []
    r_results = []

    rpar_bins = np.linspace(
        float(config['treecorr']['min_rpar']),
        float(config['treecorr']['max_rpar']),
        int(config['treecorr']['n_rpar_bins'])
    )

    # Iterate over rpar bins
    for i in range(len(rpar_bins) - 1):
        min_rpar = rpar_bins[i]
        max_rpar = rpar_bins[i + 1]

        r, xi_p, xi_x, var_xi = process_ng_rpar_bin(
            shape_catalogue,
            position_catalogue,
            random_position_catalogue,
            config,
            min_rpar,
            max_rpar
        )

        # Store the results
        xi_gn_p_results.append(xi_p)
        xi_gn_x_results.append(xi_x)
        xi_gn_var_results.append(var_xi)
        r_results.append(r)

    r_results = np.array(r_results)
    xi_gn_p_results = np.array(xi_gn_p_results)
    xi_gn_x_results = np.array(xi_gn_x_results)
    xi_gn_var_results = np.array(xi_gn_var_results)

    # save the correlation functions, include e-filter in filename if used
    base = (
        config['general']['correlation_function_folder']
        + config['general']['position_tracer']
        + '_xi_gn'
    )

    np.save(base + '_p_' + config['general']['shape_type'] + '_' + config['general']['position_type'] + e_filter_str + '.npy', xi_gn_p_results)
    np.save(base + '_x_' + config['general']['shape_type'] + '_' + config['general']['position_type'] + e_filter_str + '.npy', xi_gn_x_results)
    np.save(base + '_var_' + config['general']['shape_type'] + '_' + config['general']['position_type'] + e_filter_str + '.npy', xi_gn_var_results)
    np.save(base + '_r_' + config['general']['shape_type'] + '_' + config['general']['position_type'] + e_filter_str + '.npy', r_results)

    return

def calculate_shear_displacement_correlations( config  ):
    """
    Calculate the correlation between shapes and positions.
    """

    print('Calculating the correlation functions...')

    # load displacements and shapes
    displacements = fits.open( config['general']['velocity_catalogue_folder'] + config['general']['position_tracer'] + '_displacements.fits' )[1].data

    shapes = fits.open('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['shape_tracer'] + '_unions_desi_matched_' + config['general']['shape_type'] + '.fits')[1].data

    # convert from records -> pd.DataFrame
    shapes = pd.DataFrame.from_records( shapes )
    displacements = pd.DataFrame.from_records( displacements )

    shapes = calculate_cartesian_coordinates( shapes , 'RA' , 'Dec' , 'redshift' )
    displacements = calculate_cartesian_coordinates( displacements , 'RA' , 'Dec' , 'redshift' )

    print('Calculate the cartesian coordinates for the positions, shapes and randoms...')
    displacement_catalogue = treecorr.Catalog( x=displacements['x'], 
                                                y=displacements['y'], 
                                                z=displacements['z'],
                                                w=displacements['WEIGHT'],
                                                v1 = displacements['dRA'],
                                                v2 = displacements['dDec'], 
                                                npatch =  float( config['treecorr']['npatch'] ) )

    shape_catalogue = treecorr.Catalog( x=shapes['x'], 
                                        y=shapes['y'], 
                                        z=shapes['z'], 
                                        g1 = shapes['e1'],
                                        g2 = shapes['e2'], 
                                        w=shapes['w_iv'] ,
                                        patch_centers = displacement_catalogue.patch_centers )
        
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

        r , xi_p , xi_x, var_xi = process_vg_rpar_bin( displacement_catalogue , 
                                                       shape_catalogue , 
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
    np.save( config['general']['correlation_function_folder'] + 'velocity_'+ config['general']['position_tracer'] + '_xi_gn_p.npy', xi_gn_p_results)
    np.save( config['general']['correlation_function_folder'] + 'velocity_'+ config['general']['position_tracer'] + '_xi_gn_x.npy', xi_gn_x_results)
    np.save( config['general']['correlation_function_folder'] + 'velocity_'+ config['general']['position_tracer'] + '_xi_gn_var.npy', xi_gn_var_results)
    np.save( config['general']['correlation_function_folder'] + 'velocity_'+ config['general']['position_tracer'] + '_r.npy', r_results)

    return 

def calculate_count_displacement_correlations( config  ):
    """
    Calculate the correlation between shapes and positions.
    """

    print('Calculating the correlation functions...')

    # load positions, displacements
    displacements = fits.open( config['general']['velocity_catalogue_folder'] + config['general']['position_tracer'] + '_displacements.fits' )[1].data

    if config['general']['position_type'] == 'observed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/catalogs/'
    elif config['general']['position_type'] == 'reconstructed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results/catalogs_rec/'
    elif config['general']['position_type'] == 'rsd_removed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results_rsd_removal_only/catalogs_rec/'

    position_ngc = fits.open( cat_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data
    position_sgc = fits.open( cat_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    positions = np.concatenate(( position_ngc , position_sgc ))

    # convert from records -> pd.DataFrame
    positions = pd.DataFrame.from_records( positions )
    displacements = pd.DataFrame.from_records( displacements )

    positions = calculate_cartesian_coordinates( positions , 'RA' , 'DEC' , 'Z' )
    displacements = calculate_cartesian_coordinates( displacements , 'RA' , 'Dec' , 'redshift' )

    print('Calculate the cartesian coordinates for the positions, shapes and randoms...')
    position_catalogue = treecorr.Catalog( x=positions['x'], 
                                           y=positions['y'], 
                                           z=positions['z'],
                                           w=positions['WEIGHT'],
                                           npatch =  float( config['treecorr']['npatch'] ) )
        


    print('Calculate the cartesian coordinates for the positions, shapes and randoms...')
    displacement_catalogue = treecorr.Catalog( x=displacements['x'], 
                                                y=displacements['y'], 
                                                z=displacements['z'],
                                                w=displacements['WEIGHT'],
                                                v1 = displacements['dr1'],
                                                v2 = -displacements['dr2'], 
                                                npatch =  float( config['treecorr']['npatch'] ) )


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

        xi_p , xi_x, var_xi = process_nv_rpar_bin( position_catalogue , 
                                                       displacement_catalogue , 
                                                       config ,
                                                       min_rpar,
                                                       max_rpar )

        # Store the results
        xi_gn_p_results.append( xi_p )
        xi_gn_x_results.append( xi_x )
        xi_gn_var_results.append( var_xi )

    xi_gn_p_results = np.array( xi_gn_p_results )
    xi_gn_x_results = np.array( xi_gn_x_results )
    xi_gn_var_results = np.array( xi_gn_var_results )

    # save the correlation functions
    np.save( config['general']['correlation_function_folder'] + 'nv_'+ config['general']['position_tracer'] + '_xi_gn_p.npy', xi_gn_p_results)
    np.save( config['general']['correlation_function_folder'] + 'nv_'+ config['general']['position_tracer'] + '_xi_gn_x.npy', xi_gn_x_results)
    np.save( config['general']['correlation_function_folder'] + 'nv_'+ config['general']['position_tracer'] + '_xi_gn_var.npy', xi_gn_var_results)

    return 


def calculate_shifted_position_correlations( rsd_removed_shape_catalogue, 
                                             recon_shape_catalogue,
                                             rsd_removed_position_catalogue , 
                                             recon_position_catalogue , 
                                             random_catalogue,
                                             tau ):
    return


def calculate_size_correlations( config  ):
    """
    Calculate the correlation between shapes and positions.
    """

    print('Calculating the correlation functions...')

    # load positions, randoms, shapes
    if config['general']['position_type'] == 'observed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/catalogs/'
    elif config['general']['position_type'] == 'reconstructed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results/catalogs_rec/'
    elif config['general']['position_type'] == 'rsd_removed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results_rsd_removal_only/catalogs_rec/'

    position_ngc = fits.open( cat_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data
    position_sgc = fits.open( cat_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    positions = np.concatenate(( position_ngc , position_sgc ))

    randoms_ngc = fits.open( cat_folder + config['general']['position_tracer'] + '_NGC_' + config['general']['random_index'] +'_clustering.ran.fits' )[1].data
    randoms_sgc = fits.open( cat_folder + config['general']['position_tracer'] + '_SGC_' + config['general']['random_index'] +'_clustering.ran.fits' )[1].data

    # Combine the DESI north and south random catalogues
    randoms = np.concatenate(( randoms_ngc , randoms_sgc ))

    shapes = fits.open('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['shape_tracer'] + '_unions_desi_matched_' + config['general']['shape_type'] + '.fits')[1].data

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
    
    size_catalogue = treecorr.Catalog( x=shapes['x'], 
                                        y=shapes['y'], 
                                        z=shapes['z'], 
                                        k = shapes['flux_radius'],
                                        w=shapes['w_iv'] ,
                                        patch_centers = position_catalogue.patch_centers )
        
    # Initialize lists to store results
    xi_results = []
    xi_var_results = []

    rpar_bins = np.linspace(
        float(config['treecorr']['min_rpar']),
        float(config['treecorr']['max_rpar']),
        int(config['treecorr']['n_rpar_bins'])
    )

    # Iterate over rpar bins
    for i in range(len(rpar_bins) - 1):
        min_rpar = rpar_bins[i]
        max_rpar = rpar_bins[i + 1]

        xi, var_xi = process_nr_rpar_bin(
            position_catalogue,
            size_catalogue,
            config,
            min_rpar,
            max_rpar
        )

        # Store the results
        xi_results.append(xi)
        xi_var_results.append(var_xi)

    xi_results = np.array(xi_results)
    xi_var_results = np.array(xi_var_results)

    # save the correlation functions
    np.save(
        config['general']['correlation_function_folder']
        + config['general']['position_tracer']
        + '_xi_nr_'
        + config['general']['shape_type']
        + '_'
        + config['general']['position_type']
        + '.npy',
        xi_results
    )
    np.save(
        config['general']['correlation_function_folder']
        + config['general']['position_tracer']
        + '_xi_nr_var_'
        + config['general']['shape_type']
        + '_'
        + config['general']['position_type']
        + '.npy',
        xi_var_results
    )

    return


def calculate_delta_size_correlations( config  ):
    """
    Calculate the correlation between shapes and positions.
    """

    print('Calculating the correlation functions...')

    # load positions, randoms, shapes
    if config['general']['position_type'] == 'observed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/catalogs/'
    elif config['general']['position_type'] == 'reconstructed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results/catalogs_rec/'
    elif config['general']['position_type'] == 'rsd_removed':
        cat_folder = '/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/results_rsd_removal_only/catalogs_rec/'

    position_ngc = fits.open( cat_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data
    position_sgc = fits.open( cat_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    positions = np.concatenate(( position_ngc , position_sgc ))

    randoms_ngc = fits.open( cat_folder + config['general']['position_tracer'] + '_NGC_' + config['general']['random_index'] +'_clustering.ran.fits' )[1].data
    randoms_sgc = fits.open( cat_folder + config['general']['position_tracer'] + '_SGC_' + config['general']['random_index'] +'_clustering.ran.fits' )[1].data

    # Combine the DESI north and south random catalogues
    randoms = np.concatenate(( randoms_ngc , randoms_sgc ))

    shapes = fits.open('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['shape_tracer'] + '_unions_desi_matched_' + config['general']['shape_type'] + '.fits')[1].data

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
    
    size_catalogue = treecorr.Catalog( x=shapes['x'], 
                                        y=shapes['y'], 
                                        z=shapes['z'], 
                                        k = shapes['delta_size'],
                                        w=shapes['w_iv'] ,
                                        patch_centers = position_catalogue.patch_centers )
        
    # Initialize lists to store results
    xi_results = []
    xi_var_results = []

    rpar_bins = np.linspace(
        float(config['treecorr']['min_rpar']),
        float(config['treecorr']['max_rpar']),
        int(config['treecorr']['n_rpar_bins'])
    )

    # Iterate over rpar bins
    for i in range(len(rpar_bins) - 1):
        min_rpar = rpar_bins[i]
        max_rpar = rpar_bins[i + 1]

        xi, var_xi = process_nr_rpar_bin(
            position_catalogue,
            size_catalogue,
            config,
            min_rpar,
            max_rpar
        )

        # Store the results
        xi_results.append(xi)
        xi_var_results.append(var_xi)

    xi_results = np.array(xi_results)
    xi_var_results = np.array(xi_var_results)

    # save the correlation functions
    np.save(
        config['general']['correlation_function_folder']
        + config['general']['position_tracer']
        + '_delta_size_xi_nr_'
        + config['general']['shape_type']
        + '_'
        + config['general']['position_type']
        + '.npy',
        xi_results
    )
    np.save(
        config['general']['correlation_function_folder']
        + config['general']['position_tracer']
        + '_delta_size_xi_nr_var_'
        + config['general']['shape_type']
        + '_'
        + config['general']['position_type']
        + '.npy',
        xi_var_results
    )

    return


def calculate_displacement_shear_correlations(config):
    """
    Calculate correlations between displacement vectors and shear.
    Since TreeCorr doesn't support this directly, we implement custom pair counting.
    """
    
    print('=== STARTING DISPLACEMENT-SHEAR CORRELATION CALCULATION ===', flush=True)
    print('Step 1: Loading displacement catalogue...', flush=True)
    
    # Load displacement catalogue with displacement vectors
    print('  Opening displacement FITS file...', flush=True)
    fits_path = config['general']['velocity_catalogue_folder'] + config['general']['position_tracer'] + '_displacements.fits'
    print(f'  FITS path: {fits_path}', flush=True)
    print('  Calling fits.open()...', flush=True)
    fits_file = fits.open(fits_path)
    print('  fits.open() completed, accessing [1].data...', flush=True)
    displacements = fits_file[1].data
    print(f'  Loaded {len(displacements):,} displacement vectors', flush=True)
    fits_file.close()
    print('  FITS file closed', flush=True)
    
    print('  About to start Step 2...', flush=True)
    print('Step 2: Loading shape displacements catalogue...', flush=True)
    # Load shape displacements catalogue 
    print('  Opening shape displacement FITS file...', flush=True)
    shape_fits_path = config['general']['velocity_catalogue_folder'] + config['general']['shape_tracer'] + '_shape_displacements.fits'
    print(f'  Shape FITS path: {shape_fits_path}', flush=True)
    print('  Calling fits.open() for shapes...', flush=True)
    shape_fits_file = fits.open(shape_fits_path)
    print('  Shape fits.open() completed, accessing [1].data...', flush=True)
    shape_displacements = shape_fits_file[1].data
    print(f'  Loaded {len(shape_displacements):,} shape displacements', flush=True)
    shape_fits_file.close()
    print('  Shape FITS file closed', flush=True)
    
    print('  Step 2 completed, about to start Step 3...', flush=True)
    print('Step 3: Converting to DataFrames and fixing byte order...', flush=True)
    # Convert to DataFrames and ensure native byte order
    displacements = pd.DataFrame.from_records(displacements)
    shape_displacements = pd.DataFrame.from_records(shape_displacements)
    
    # Fix byte order issues by converting to native byte order
    for col in displacements.columns:
        if hasattr(displacements[col].dtype, 'newbyteorder'):
            displacements[col] = displacements[col].values.astype(displacements[col].dtype.newbyteorder('='))
    
    for col in shape_displacements.columns:
        if hasattr(shape_displacements[col].dtype, 'newbyteorder'):
            shape_displacements[col] = shape_displacements[col].values.astype(shape_displacements[col].dtype.newbyteorder('='))
    
    print('  DataFrame conversion and byte order fix complete')
    
    print('Step 4: Calculating Cartesian coordinates for displacements...')
    # Calculate Cartesian coordinates
    displacements = calculate_cartesian_coordinates(displacements, 'RA', 'Dec', 'redshift')
    print('  Displacement coordinates calculated')
    
    print('Step 5: Calculating Cartesian coordinates for shapes...')
    shape_displacements = calculate_cartesian_coordinates(shape_displacements, 'RA', 'Dec', 'redshift')
    print('  Shape coordinates calculated')
    
    print('Step 5a: Applying spatial pre-filtering...')
    # Find overlapping sky regions to reduce computation
    disp_ra_min, disp_ra_max = displacements['RA'].min(), displacements['RA'].max()
    disp_dec_min, disp_dec_max = displacements['Dec'].min(), displacements['Dec'].max()
    shape_ra_min, shape_ra_max = shape_displacements['RA'].min(), shape_displacements['RA'].max()
    shape_dec_min, shape_dec_max = shape_displacements['Dec'].min(), shape_displacements['Dec'].max()
    
    # Find overlap region
    overlap_ra_min = max(disp_ra_min, shape_ra_min) 
    overlap_ra_max = min(disp_ra_max, shape_ra_max)
    overlap_dec_min = max(disp_dec_min, shape_dec_min)
    overlap_dec_max = min(disp_dec_max, shape_dec_max)
    
    print(f'  Displacement catalogue sky coverage: RA [{disp_ra_min:.2f}, {disp_ra_max:.2f}], Dec [{disp_dec_min:.2f}, {disp_dec_max:.2f}]')
    print(f'  Shape catalogue sky coverage: RA [{shape_ra_min:.2f}, {shape_ra_max:.2f}], Dec [{shape_dec_min:.2f}, {shape_dec_max:.2f}]')
    print(f'  Overlap region: RA [{overlap_ra_min:.2f}, {overlap_ra_max:.2f}], Dec [{overlap_dec_min:.2f}, {overlap_dec_max:.2f}]')
    
    # Filter displacement catalogue to overlap region
    disp_mask = ((displacements['RA'] >= overlap_ra_min) & 
                 (displacements['RA'] <= overlap_ra_max) &
                 (displacements['Dec'] >= overlap_dec_min) & 
                 (displacements['Dec'] <= overlap_dec_max))
    displacements_filtered = displacements[disp_mask].copy()
    
    # Filter shape catalogue to overlap region  
    shape_mask = ((shape_displacements['RA'] >= overlap_ra_min) &
                  (shape_displacements['RA'] <= overlap_ra_max) &
                  (shape_displacements['Dec'] >= overlap_dec_min) &
                  (shape_displacements['Dec'] <= overlap_dec_max))
    shape_displacements_filtered = shape_displacements[shape_mask].copy()
    
    print(f'  Pre-filtering: Displacement catalogue reduced from {len(displacements):,} to {len(displacements_filtered):,} objects')
    print(f'  Pre-filtering: Shape catalogue reduced from {len(shape_displacements):,} to {len(shape_displacements_filtered):,} objects')
    
    # Use filtered catalogues
    displacements = displacements_filtered
    shape_displacements = shape_displacements_filtered
    
    print('Step 6: Setting up binning...')
    # Set up binning
    rpar_bins = np.linspace(
        float(config['treecorr']['min_rpar']),
        float(config['treecorr']['max_rpar']),
        int(config['treecorr']['n_rpar_bins'])
    )
    
    rperp_bins = np.logspace(
        np.log10(float(config['treecorr']['min_rperp'])),
        np.log10(float(config['treecorr']['max_rperp'])),
        int(config['treecorr']['n_rperp_bins']) + 1
    )
    print(f'  rpar bins: {len(rpar_bins)-1} bins from {rpar_bins[0]:.1f} to {rpar_bins[-1]:.1f}')
    print(f'  rperp bins: {len(rperp_bins)-1} bins from {rperp_bins[0]:.1f} to {rperp_bins[-1]:.1f}')
    
    # Initialize result arrays
    xi_plus_results = []
    xi_cross_results = []
    
    print(f'Step 7: Starting rpar bin processing ({len(rpar_bins)-1} bins total)...')
    print('WARNING: This will be very slow with full catalogues!')
    
    # Process each rpar bin
    for i in range(len(rpar_bins) - 1):
        min_rpar = rpar_bins[i]
        max_rpar = rpar_bins[i + 1]
        
        print(f'=== Processing rpar bin {i+1}/{len(rpar_bins)-1}: {min_rpar:.2f} to {max_rpar:.2f} ===')
        
        xi_plus, xi_cross = process_displacement_shear_rpar_bin(
            displacements, shape_displacements, config,
            min_rpar, max_rpar, rperp_bins
        )
        
        xi_plus_results.append(xi_plus)
        xi_cross_results.append(xi_cross)
        print(f'  Completed rpar bin {i+1}/{len(rpar_bins)-1}')
    
    xi_plus_results = np.array(xi_plus_results)
    xi_cross_results = np.array(xi_cross_results)
    
    # Calculate bin centers for output
    rperp_centers = np.sqrt(rperp_bins[:-1] * rperp_bins[1:])
    
    # Save results with appropriate filename
    base_filename = (config['general']['correlation_function_folder'] + 
                    'displacement_shear_' + config['general']['position_tracer'])
    
    # Add region suffix if using region filtering
    if config.has_section('region_filter') and config.getboolean('region_filter', 'use_region_filter'):
        ra_min = config.getfloat('region_filter', 'ra_min')
        ra_max = config.getfloat('region_filter', 'ra_max')
        dec_min = config.getfloat('region_filter', 'dec_min')
        dec_max = config.getfloat('region_filter', 'dec_max')
        base_filename += f'_region_ra{ra_min:.0f}to{ra_max:.0f}_dec{dec_min:.0f}to{dec_max:.0f}'
    
    np.save(base_filename + '_xi_plus.npy', xi_plus_results)
    np.save(base_filename + '_xi_cross.npy', xi_cross_results)
    np.save(base_filename + '_rperp.npy', rperp_centers)
    np.save(base_filename + '_rpar_bins.npy', rpar_bins)
    
    print(f'Results saved with base filename: {base_filename}')
    
    return


def create_region_validation_plots(displacements, shape_displacements, config):
    """
    Create validation plots for the filtered region.
    """
    import matplotlib.pyplot as plt
    
    ra_min = config.getfloat('region_filter', 'ra_min')
    ra_max = config.getfloat('region_filter', 'ra_max')
    dec_min = config.getfloat('region_filter', 'dec_min')
    dec_max = config.getfloat('region_filter', 'dec_max')
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Sky positions
    axes[0,0].scatter(displacements['RA'], displacements['Dec'], 
                      s=1, alpha=0.5, label='Displacements')
    axes[0,0].scatter(shape_displacements['RA'], shape_displacements['Dec'], 
                      s=1, alpha=0.5, label='Shapes')
    axes[0,0].set_xlabel('RA (deg)')
    axes[0,0].set_ylabel('Dec (deg)')
    axes[0,0].set_title('Sky Coverage of Filtered Region')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    # Plot 2: Redshift distributions
    axes[0,1].hist(displacements['redshift'], bins=50, alpha=0.7, 
                   label=f'Displacements (N={len(displacements):,})', density=True)
    axes[0,1].hist(shape_displacements['redshift'], bins=50, alpha=0.7, 
                   label=f'Shapes (N={len(shape_displacements):,})', density=True)
    axes[0,1].set_xlabel('Redshift')
    axes[0,1].set_ylabel('Density')
    axes[0,1].set_title('Redshift Distributions')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # Plot 3: Displacement magnitude distribution
    disp_mag = np.sqrt(displacements['dr1']**2 + displacements['dr2']**2)
    axes[1,0].hist(disp_mag, bins=50, alpha=0.7)
    axes[1,0].set_xlabel('Displacement Magnitude')
    axes[1,0].set_ylabel('Count')
    axes[1,0].set_title('Displacement Vector Magnitudes')
    axes[1,0].grid(True, alpha=0.3)
    
    # Plot 4: Ellipticity distribution
    ellip_mag = np.sqrt(shape_displacements['e1']**2 + shape_displacements['e2']**2)
    axes[1,1].hist(ellip_mag, bins=50, alpha=0.7)
    axes[1,1].set_xlabel('Ellipticity Magnitude')
    axes[1,1].set_ylabel('Count')
    axes[1,1].set_title('Shape Ellipticity Magnitudes')
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot to validation directory
    plot_filename = (f'/home/murray/intrinsic_alignments/validation/correlation_functions/' + 
                    f'displacement_shear_{config["general"]["position_tracer"]}_' +
                    f'region_ra{ra_min:.0f}to{ra_max:.0f}_dec{dec_min:.0f}to{dec_max:.0f}_validation.png')
    plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f'  Validation plot saved: {plot_filename}')


def process_displacement_shear_rpar_bin(displacements, shape_displacements, config, 
                                       min_rpar, max_rpar, rperp_bins):
    """
    Process a single rpar bin for displacement-shear correlations.
    """
    
    print(f'  Extracting coordinates and data...')
    # Extract coordinates and data
    disp_pos = np.column_stack([displacements['x'], displacements['y'], displacements['z']])
    shape_pos = np.column_stack([shape_displacements['x'], shape_displacements['y'], shape_displacements['z']])
    print(f'    Displacement positions: {disp_pos.shape}')
    print(f'    Shape positions: {shape_pos.shape}')
    
    # Displacement vectors (tangential components)
    disp_dr1 = displacements['dr1'].values
    disp_dr2 = displacements['dr2'].values
    
    # Shear components
    shape_e1 = shape_displacements['e1'].values
    shape_e2 = shape_displacements['e2'].values
    
    # Weights
    disp_weights = displacements['WEIGHT'].values
    shape_weights = shape_displacements['w_iv'].values
    print(f'  Data extraction complete')
    
    print(f'  Applying pre-filtering strategies...')
    
    # Strategy 1: Redshift-based filtering for this specific rpar bin
    disp_z = displacements['redshift'].values
    shape_z = shape_displacements['redshift'].values
    
    # Convert rpar range to approximate redshift differences
    z_buffer = abs(max_rpar) / 3000.0  # Very rough conversion
    z_min = disp_z.min() - z_buffer
    z_max = disp_z.max() + z_buffer
    
    z_mask = (shape_z >= z_min) & (shape_z <= z_max)
    
    # Strategy 2: Region-based filtering (configurable)
    region_mask = np.ones(len(shape_z), dtype=bool)  # Default: use all objects
    
    if config.has_section('region_filter') and config.getboolean('region_filter', 'use_region_filter'):
        print('  Applying region-based filtering...')
        ra_min = config.getfloat('region_filter', 'ra_min')
        ra_max = config.getfloat('region_filter', 'ra_max')
        dec_min = config.getfloat('region_filter', 'dec_min')
        dec_max = config.getfloat('region_filter', 'dec_max')
        
        # Apply region filter to both displacement and shape catalogues
        disp_region_mask = ((displacements['RA'] >= ra_min) & 
                           (displacements['RA'] <= ra_max) &
                           (displacements['Dec'] >= dec_min) & 
                           (displacements['Dec'] <= dec_max))
        
        shape_region_mask = ((shape_displacements['RA'] >= ra_min) & 
                            (shape_displacements['RA'] <= ra_max) &
                            (shape_displacements['Dec'] >= dec_min) & 
                            (shape_displacements['Dec'] <= dec_max))
        
        # Filter displacement catalogue to region
        displacements = displacements[disp_region_mask]
        disp_pos = disp_pos[disp_region_mask]
        disp_dr1 = disp_dr1[disp_region_mask]
        disp_dr2 = disp_dr2[disp_region_mask]
        disp_weights = disp_weights[disp_region_mask]
        
        # Update masks for shape catalogue
        region_mask = shape_region_mask
        
        print(f'    Region filter [{ra_min:.1f}, {ra_max:.1f}] × [{dec_min:.1f}, {dec_max:.1f}]')
        print(f'    Displacement objects in region: {np.sum(disp_region_mask):,} / {len(disp_region_mask):,}')
        print(f'    Shape objects in region: {np.sum(region_mask):,} / {len(region_mask):,}')
    
    # Combine all filters
    combined_mask = z_mask & region_mask
    
    # Apply filtering to shape catalogue
    shape_pos_filtered = shape_pos[combined_mask]
    shape_e1_filtered = shape_e1[combined_mask] 
    shape_e2_filtered = shape_e2[combined_mask]
    shape_weights_filtered = shape_weights[combined_mask]
    
    print(f'  Final filtering: Shape catalogue reduced from {len(shape_pos):,} to {len(shape_pos_filtered):,} objects')
    
    # Create validation plots if requested
    if (config.has_section('region_filter') and 
        config.getboolean('region_filter', 'use_region_filter') and
        config.getboolean('region_filter', 'create_validation_plots')):
        print('  Creating validation plots for filtered region...')
        create_region_validation_plots(displacements, shape_displacements[combined_mask], config)
    
    print(f'  Building KD-tree for {len(shape_pos_filtered):,} filtered shape positions...')
    # Build KD-tree for filtered shape positions
    shape_tree = cKDTree(shape_pos_filtered)
    print(f'  KD-tree construction complete')
    
    # Initialize correlation sums
    xi_plus_sum = np.zeros(len(rperp_bins) - 1)
    xi_cross_sum = np.zeros(len(rperp_bins) - 1)
    weight_sum = np.zeros(len(rperp_bins) - 1)
    print(f'  Initialized correlation arrays with {len(rperp_bins)-1} rperp bins')
    
    print(f'  Starting main calculation loop over {len(disp_pos):,} displacement points...')
    print(f'  Progress will be reported every 10,000 points')
    
    # Process each displacement point
    for i in range(len(disp_pos)):
        if i % 10000 == 0:
            print(f'    Processing displacement {i:,}/{len(disp_pos):,} ({100*i/len(disp_pos):.1f}%)')
        
        # Find nearby shape objects
        max_sep = float(config['treecorr']['max_rperp']) * 2  # Search radius
        nearby_indices = shape_tree.query_ball_point(disp_pos[i], r=max_sep)
        
        if len(nearby_indices) == 0:
            continue
            
        # Calculate separations to nearby shapes (using filtered positions)
        nearby_pos = shape_pos_filtered[nearby_indices]
        sep_vectors = nearby_pos - disp_pos[i]
        
        # Calculate parallel and perpendicular separations
        # Use displacement position as reference for line-of-sight
        los_unit = disp_pos[i] / np.linalg.norm(disp_pos[i])
        rpar = np.dot(sep_vectors, los_unit)
        
        # Filter by parallel separation
        rpar_mask = (rpar >= min_rpar) & (rpar < max_rpar)
        if not np.any(rpar_mask):
            continue
            
        # Calculate perpendicular separations
        rpar_filtered = rpar[rpar_mask]
        sep_filtered = sep_vectors[rpar_mask]
        nearby_filtered = np.array(nearby_indices)[rpar_mask]
        
        # Perpendicular component
        rperp_vectors = sep_filtered - np.outer(rpar_filtered, los_unit)
        rperp = np.linalg.norm(rperp_vectors, axis=1)
        
        # Calculate tangential unit vectors
        rperp_unit = rperp_vectors / rperp[:, np.newaxis]
        cos_theta = rperp_unit[:, 0]
        sin_theta = rperp_unit[:, 1]
        
        # Calculate displacement components (vector projection)
        disp_tangent = disp_dr1[i] * cos_theta + disp_dr2[i] * sin_theta
        disp_cross = -disp_dr1[i] * sin_theta + disp_dr2[i] * cos_theta
        
        # Calculate spin-2 rotation angles for ellipticity
        cos_2theta = cos_theta**2 - sin_theta**2
        sin_2theta = 2 * cos_theta * sin_theta
        
        # Calculate ellipticity components (spin-2 projection)
        shape_tangent = -shape_e1_filtered[nearby_filtered] * cos_2theta - shape_e2_filtered[nearby_filtered] * sin_2theta
        shape_cross = -shape_e1_filtered[nearby_filtered] * sin_2theta + shape_e2_filtered[nearby_filtered] * cos_2theta
        
        # Bin by perpendicular separation
        rperp_bin_indices = np.digitize(rperp, rperp_bins) - 1
        valid_bins = (rperp_bin_indices >= 0) & (rperp_bin_indices < len(rperp_bins) - 1)
        
        if not np.any(valid_bins):
            continue
            
        # Add to correlation sums
        for j in range(len(rperp)):
            if not valid_bins[j]:
                continue
                
            bin_idx = rperp_bin_indices[j]
            weight = disp_weights[i] * shape_weights_filtered[nearby_filtered[j]]
            
            xi_plus_sum[bin_idx] += weight * (disp_tangent[j] * shape_tangent[j] + disp_cross[j] * shape_cross[j])
            xi_cross_sum[bin_idx] += weight * (disp_tangent[j] * shape_cross[j] - disp_cross[j] * shape_tangent[j])
            weight_sum[bin_idx] += weight
    
    # Normalize by weights
    xi_plus = np.zeros_like(xi_plus_sum)
    xi_cross = np.zeros_like(xi_cross_sum)
    
    nonzero_mask = weight_sum > 0
    xi_plus[nonzero_mask] = xi_plus_sum[nonzero_mask] / weight_sum[nonzero_mask]
    xi_cross[nonzero_mask] = xi_cross_sum[nonzero_mask] / weight_sum[nonzero_mask]
    
    return xi_plus, xi_cross

##############################

def _ra_calculate_displacement_vectors( config ):
    """
    Calculate the displacement vectors between the RSD-removed and reconstructed catalogues
    for both positions and shapes. Saves the displacement arrays to disk.
    """
    recon_position_ngc = fits.open( desi_recon_catalogue_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data
    recon_position_sgc = fits.open( desi_recon_catalogue_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    recon_positions = np.concatenate(( recon_position_ngc , recon_position_sgc ))

    rsd_removed_position_ngc = fits.open( desi_rsd_removed_catalogue_folder + config['general']['position_tracer'] + '_SGC_clustering.dat.fits' )[1].data
    rsd_removed_position_sgc = fits.open( desi_rsd_removed_catalogue_folder + config['general']['position_tracer'] + '_NGC_clustering.dat.fits' )[1].data

    # Combine the DESI north and south galaxy catalogues
    rsd_removed_positions = np.concatenate(( rsd_removed_position_ngc , rsd_removed_position_sgc ))

    dRA = recon_positions['RA'] - rsd_removed_positions['RA']
    dDEC = recon_positions['DEC'] - rsd_removed_positions['DEC']
    dZ = recon_positions['Z'] - rsd_removed_positions['Z']

    recon_shapes = fits.open('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['shape_tracer'] + '_unions_desi_matched_reconstructed.fits')[1].data
    rsd_removed_shapes = fits.open('/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/shape_catalogs/' + config['general']['shape_tracer'] + '_unions_desi_matched_rsd_removed.fits')[1].data

    dRA_shapes = recon_shapes['RA'] - rsd_removed_shapes['RA']
    dDEC_shapes = recon_shapes['Dec'] - rsd_removed_shapes['Dec']
    dZ_shapes = recon_shapes['redshift'] - rsd_removed_shapes['redshift']

    # Create a new FITS table with RA, Dec, redshift, dRA, dDec, dZ, and weights
    displacement_cols = [
        fits.Column(name='RA', format='E', array=recon_positions['RA']),
        fits.Column(name='Dec', format='E', array=recon_positions['DEC']),
        fits.Column(name='redshift', format='E', array=recon_positions['Z']),
        fits.Column(name='rsd_removed_RA', format='E', array=rsd_removed_positions['RA']),
        fits.Column(name='rsd_removed_Dec', format='E', array=rsd_removed_positions['DEC']),
        fits.Column(name='rsd_removed_redshift', format='E', array=rsd_removed_positions['Z']),
        fits.Column(name='dRA', format='E', array=dRA),
        fits.Column(name='dDec', format='E', array=dDEC),
        fits.Column(name='dZ', format='E', array=dZ),
        fits.Column(name='WEIGHT', format='E', array=recon_positions['WEIGHT']),
        fits.Column(name='rsd_removed_WEIGHT', format='E', array=rsd_removed_positions['WEIGHT'])
    ]
    displacement_hdu = fits.BinTableHDU.from_columns(displacement_cols)
    displacement_hdu.writeto( config['general']['velocity_catalogue_folder'] + config['general']['position_tracer'] + '_displacements.fits', overwrite=True)

    # Create a new FITS table for shape displacements
    shape_displacement_cols = [
        fits.Column(name='RA', format='E', array=recon_shapes['RA']),
        fits.Column(name='Dec', format='E', array=recon_shapes['Dec']),
        fits.Column(name='redshift', format='E', array=recon_shapes['redshift']),
        fits.Column(name='rsd_removed_RA', format='E', array=rsd_removed_shapes['RA']),
        fits.Column(name='rsd_removed_Dec', format='E', array=rsd_removed_shapes['Dec']),
        fits.Column(name='rsd_removed_redshift', format='E', array=rsd_removed_shapes['redshift']),
        fits.Column(name='e1', format='E', array=recon_shapes['e1']),
        fits.Column(name='e2', format='E', array=recon_shapes['e2']),
        fits.Column(name='w_iv', format='E', array=recon_shapes['w_iv']),
        fits.Column(name='dRA', format='E', array=dRA_shapes),
        fits.Column(name='dDec', format='E', array=dDEC_shapes),
        fits.Column(name='dZ', format='E', array=dZ_shapes)
    ]
    shape_displacement_hdu = fits.BinTableHDU.from_columns(shape_displacement_cols)
    shape_displacement_hdu.writeto(
        config['general']['velocity_catalogue_folder'] + config['general']['shape_tracer'] + '_shape_displacements.fits',
        overwrite=True
    )

    return