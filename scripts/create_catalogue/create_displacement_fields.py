import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'measure_correlation_functions'))
import measure_correlations as mcf
import argparse

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
    if config.getboolean('general', 'run_calculate_displacement_vectors'):
        mcf.calculate_displacement_vectors( config = config )

    # calculate the correlation functions
    #mcf.calculate_count_displacement_correlations( config = config )

    #mcf.calculate_size_correlations( config = config )
    mcf.calculate_delta_size_correlations( config = config )

if __name__ == "__main__":
    main()