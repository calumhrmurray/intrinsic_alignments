from . import measure_correlations as mcf
import argparse
import subprocess
import os

def main():
    parser = argparse.ArgumentParser(description="Run correlation functions with config file.")
    parser.add_argument('--config_file', type=str, required=True, help='Path to the config file')
    args = parser.parse_args()

    config_file = args.config_file
    print(f"Using config file: {config_file}")

    # Now you can load and use the config file as needed
    # For example, if it's an .ini file:
    import configparser
    config = configparser.ConfigParser()
    config.read(config_file)

    print( config )

    # create the shape catalogues, with the matching to DESI positions
    if config.getboolean('general', 'run_shape_matching'):
        mcf.match_shapes_to_positions( config = config )
    
    # calculate the correlation functions
    mcf.calculate_correlations( config = config )
    
    # automatically run plotting if enabled in config
    if config.has_option('general', 'auto_plot') and config.getboolean('general', 'auto_plot'):
        print("\nAuto-plotting enabled, creating figures...")
        try:
            data_dir = config['general']['correlation_function_folder']
            output_dir = config.get('general', 'figure_output_folder', 
                                  fallback='/home/murray/intrinsic_alignments/../results/figures/correlation_functions')
            
            # Run the plotting script
            plot_script = os.path.join(os.path.dirname(__file__), 'plot_displacement_shear_results.py')
            cmd = ['python', plot_script, '--data_dir', data_dir, '--output_dir', output_dir]
            
            print(f"Running: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print("Plotting completed successfully!")
                print(result.stdout)
            else:
                print(f"Plotting failed with return code {result.returncode}")
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)
                
        except Exception as e:
            print(f"Error running plotting script: {e}")
    else:
        print("\nAuto-plotting disabled. To enable, add 'auto_plot = True' to [general] section in config.")
        print("To plot manually, run:")
        data_dir = config['general']['correlation_function_folder']
        print(f"python scripts/plot_displacement_shear_results.py --data_dir {data_dir}")

if __name__ == "__main__":
    main()