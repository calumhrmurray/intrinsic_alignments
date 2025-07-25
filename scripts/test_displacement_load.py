#!/usr/bin/env python3
"""
Minimal test script to debug displacement loading issue
"""
import sys
import configparser
from astropy.io import fits

def test_load():
    print("Starting test...", flush=True)
    
    # Load config
    config = configparser.ConfigParser()
    config.read('/home/murray/intrinsic_alignments/scripts/config_files/lrg_observed_displacement_shear.ini')
    
    print("Config loaded", flush=True)
    
    # Test displacement loading
    fits_path = config['general']['velocity_catalogue_folder'] + config['general']['position_tracer'] + '_displacements.fits'
    print(f"FITS path: {fits_path}", flush=True)
    
    print("Opening FITS file...", flush=True)
    fits_file = fits.open(fits_path)
    print("FITS opened", flush=True)
    
    print("Accessing data...", flush=True)
    data = fits_file[1].data
    print(f"Data loaded: {len(data)} records", flush=True)
    
    fits_file.close()
    print("FITS closed", flush=True)
    
    print("Test completed successfully!", flush=True)

if __name__ == "__main__":
    test_load()