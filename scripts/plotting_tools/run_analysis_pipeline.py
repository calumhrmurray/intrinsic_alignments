#!/usr/bin/env python3
"""
Complete analysis pipeline that runs correlation calculations and plotting in one go.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'measure_correlation_functions'))
import measure_correlations as mcf
import argparse
import subprocess
import os
import configparser

def run_correlation_and_plot(config, correlation_type):
    """Run correlation calculation and plotting for a specific type."""
    
    print(f"\n{'='*60}")
    print(f"Running {correlation_type.upper()} correlation analysis")
    print(f"{'='*60}")
    
    # Run the appropriate correlation function
    if correlation_type == 'displacement_shear':
        mcf.calculate_shear_displacement_correlations(config=config)
        plot_script = 'plot_displacement_shear_results.py'
    elif correlation_type == 'count_displacement':
        mcf.calculate_count_displacement_correlations(config=config)
        plot_script = 'plot_count_displacement_results.py'  # assuming this exists
    elif correlation_type == 'size':
        mcf.calculate_size_correlations(config=config)  
        plot_script = 'plot_size_results.py'  # assuming this exists
    else:
        print(f"Unknown correlation type: {correlation_type}")
        return False
    
    # Run plotting if enabled
    if config.has_option('general', 'auto_plot') and config.getboolean('general', 'auto_plot'):
        print(f"\nCreating {correlation_type} plots...")
        try:
            data_dir = config['general']['correlation_function_folder']
            output_dir = config.get('general', 'figure_output_folder', 
                                  fallback=f'/home/murray/intrinsic_alignments/../results/figures/{correlation_type}')
            
            # Create output directory
            os.makedirs(output_dir, exist_ok=True)
            
            # Run the plotting script if it exists
            plot_script_path = os.path.join(os.path.dirname(__file__), plot_script)
            if os.path.exists(plot_script_path):
                cmd = ['python', plot_script_path, '--data_dir', data_dir, '--output_dir', output_dir]
                print(f"Running: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    print(f"{correlation_type} plotting completed successfully!")
                    if result.stdout:
                        print(result.stdout)
                else:
                    print(f"Plotting failed with return code {result.returncode}")
                    if result.stderr:
                        print("STDERR:", result.stderr)
            else:
                print(f"Plot script {plot_script} not found, skipping plotting")
                
        except Exception as e:
            print(f"Error running plotting: {e}")
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Run complete correlation analysis pipeline.")
    parser.add_argument('--config_file', type=str, required=True, help='Path to the config file')
    parser.add_argument('--correlation_types', type=str, nargs='+', 
                       choices=['displacement_shear', 'count_displacement', 'size', 'all'],
                       default=['displacement_shear'],
                       help='Types of correlations to calculate')
    
    args = parser.parse_args()
    
    config_file = args.config_file
    print(f"Using config file: {config_file}")
    
    # Load config
    config = configparser.ConfigParser()
    config.read(config_file)
    
    # Expand 'all' option
    if 'all' in args.correlation_types:
        correlation_types = ['displacement_shear', 'count_displacement', 'size']
    else:
        correlation_types = args.correlation_types
    
    # Create shape catalogues if needed
    if config.getboolean('general', 'run_shape_matching'):
        print("Creating shape catalogues with DESI position matching...")
        mcf.match_shapes_to_positions(config=config)
    
    # Calculate displacement vectors if needed
    if config.has_option('general', 'run_calculate_displacement_vectors') and config.getboolean('general', 'run_calculate_displacement_vectors'):
        print("Calculating displacement vectors...")
        mcf.calculate_displacement_vectors(config=config)
    
    # Run each correlation type
    success_count = 0
    for correlation_type in correlation_types:
        try:
            if run_correlation_and_plot(config, correlation_type):
                success_count += 1
        except Exception as e:
            print(f"Error processing {correlation_type}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n{'='*60}")
    print(f"Pipeline completed: {success_count}/{len(correlation_types)} correlations successful")
    print(f"{'='*60}")
    
    if config.has_option('general', 'auto_plot') and config.getboolean('general', 'auto_plot'):
        figure_dir = config.get('general', 'figure_output_folder', 
                              fallback='/home/murray/intrinsic_alignments/figures')
        print(f"\nFigures saved to: {figure_dir}")

if __name__ == "__main__":
    main()