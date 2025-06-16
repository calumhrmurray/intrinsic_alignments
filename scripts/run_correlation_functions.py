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

    # create the shape catalogues, with the matching to DESI positions
    if config['general']['run_shape_matching']:
        mcf.match_shapes_to_positions( config = config )
    # calculate the correlation functions
    mcf.calculate_correlations( config = config )

if __name__ == "__main__":
    main()