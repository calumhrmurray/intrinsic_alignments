#!/usr/bin/env python3
"""
Run displacement-shear correlation function using config file.
"""
import numpy as np
import sys
import os
import argparse
import configparser

from . import measure_correlations
from .measure_correlations import calculate_displacement_shear_correlations

def main():
    parser = argparse.ArgumentParser(description='Calculate displacement-shear correlations')
    parser.add_argument('--config_file', required=True, 
                       help='Path to configuration file')
    
    args = parser.parse_args()
    
    # Load configuration
    config = configparser.ConfigParser()
    config.read(args.config_file)
    
    print(f"Running displacement-shear correlation with config: {args.config_file}")
    
    # Check if required data files exist
    velocity_folder = config['general']['velocity_catalogue_folder']
    position_tracer = config['general']['position_tracer']
    shape_tracer = config['general']['shape_tracer']
    
    disp_file = f"{velocity_folder}{position_tracer}_displacements.fits"
    shape_disp_file = f"{velocity_folder}{shape_tracer}_shape_displacements.fits"
    
    if not os.path.exists(disp_file):
        print(f"Error: Displacement file not found: {disp_file}")
        print("Please run calculate_displacement_vectors first.")
        sys.exit(1)
        
    if not os.path.exists(shape_disp_file):
        print(f"Error: Shape displacement file not found: {shape_disp_file}")
        print("Please run calculate_displacement_vectors first.")
        sys.exit(1)
    
    print("Required data files found. Running displacement-shear correlation...")
    
    try:
        # Run the displacement-shear correlation function
        calculate_displacement_shear_correlations(config)
        print("Displacement-shear correlation calculation completed successfully!")
        
        # Check output files
        output_base = f"{config['general']['correlation_function_folder']}displacement_shear_{position_tracer}"
        output_files = [
            f"{output_base}_xi_plus.npy",
            f"{output_base}_xi_cross.npy", 
            f"{output_base}_rperp.npy",
            f"{output_base}_rpar_bins.npy"
        ]
        
        for outfile in output_files:
            if os.path.exists(outfile):
                data = np.load(outfile)
                print(f"Output file {outfile}: shape {data.shape}")
            else:
                print(f"Warning: Expected output file not found: {outfile}")
                
    except Exception as e:
        print(f"Error running displacement-shear correlation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()