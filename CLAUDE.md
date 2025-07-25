# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a research project studying intrinsic alignments in galaxy clustering using DESI and UNIONS survey data. The codebase computes correlation functions between galaxy positions/shapes and analyzes intrinsic alignment signals.

## Key Commands

### Running Correlation Functions
- Main entry point: `python scripts/measure_correlation_functions/run_correlation_functions.py --config_file <config_file>`
- Use SLURM script: `sbatch scripts/slurm_scripts/run_measure_correlation_functions.sh`
- Create displacement catalogues: `sbatch scripts/slurm_scripts/run_create_displacement_catalogues.sh`

### Environment Setup
- Activate conda environment: `conda activate calum_conda`
- The project uses Python with dependencies: numpy, treecorr, astropy, matplotlib, pandas, healpy

## Code Architecture

### Core Components

**Scripts Directory Structure**:
- `scripts/measure_correlation_functions/`: Core correlation analysis
  - `run_correlation_functions.py`: Main entry point that loads config and delegates to measure_correlations
  - `measure_correlations.py`: Core correlation function calculations (1000+ lines, main workhorse)
  - `correlation_function_helpers.py`: Helper functions for 2D correlation calculations
  - `run_displacement_shear.py`: Specialized displacement-shear correlations
- `scripts/create_catalogue/`: Catalogue creation and validation
  - `create_displacement_fields.py`: Creates displacement vector catalogues
  - `validate_displacement_catalogues.py`: Validation tools for displacement data
- `scripts/plotting_tools/`: Analysis and visualization
  - `plot_displacement_shear_results.py`: Plotting displacement-shear results
  - `plot_results_simple.py`: Simple correlation function plots
  - `run_analysis_pipeline.py`: Complete analysis pipeline
- `scripts/slurm_scripts/`: HPC job submission scripts
- `scripts/config_files/`: Configuration files for different analyses

**Configuration System**:
- Config files in `scripts/config_files/` use INI format
- Key sections: `[general]`, `[treecorr]`, `[ellipticity_filter]`
- Controls data paths, correlation parameters, filtering options

**Data Flow**:
1. Shape catalogues (UNIONS) matched to position catalogues (DESI)
2. Displacement vectors calculated between reconstructed/RSD-removed catalogues
3. Various correlation functions computed using TreeCorr library
4. Results saved as .npy files in specified output directories

### Key Functions in measure_correlations.py

**Data Processing**:
- `match_shapes_to_positions()`: Cross-matches UNIONS shapes with DESI positions
- `calculate_displacement_vectors()`: Computes position/shape displacements
- `calculate_cartesian_coordinates()`: Converts RA/Dec/z to Cartesian

**Correlation Functions**:
- `calculate_correlations()`: Main galaxy-shear (NG) correlations with ellipticity filtering
- `calculate_shear_displacement_correlations()`: Velocity-shear (VG) correlations  
- `calculate_count_displacement_correlations()`: Count-velocity (NV) correlations
- `calculate_size_correlations()`: Galaxy-size (NR) correlations

**Processing Helpers**:
- `process_ng_rpar_bin()`: Processes single r_parallel bin for NG correlations
- `process_vg_rpar_bin()`: Processes single r_parallel bin for VG correlations
- Functions iterate over r_parallel bins and use TreeCorr for actual calculations

### Notebooks Structure

**Analysis** (`notebooks/analysis/`): Current analysis notebooks
**Preliminary Analysis** (`notebooks/preliminary_analysis/`): Exploratory work on different galaxy samples
**Model** (`notebooks/model/`): Theoretical modeling work

### Data Architecture

**External Data Paths** (hardcoded in code):
- DESI catalogues: `/n17data/murray/desi_data/DESI/catalogs/`
- Reconstructed: `/n17data/murray/desi_data/DESI/results/catalogs_rec/`
- RSD-removed: `/n17data/murray/desi_data/DESI/results_rsd_removal_only/catalogs_rec/`
- Shape catalogues: `/n17data/murray/desi_data/DESI/shape_catalogs/`

**Results Directory Structure** (`results/`):
- `results/catalogues/`: FITS catalogue files (UNIONS-DESI matched catalogues)
- `results/correlation_functions/`: .npy correlation function results organized by type
- `results/figures/`: All plots and visualizations organized by analysis type
- `results/validation/`: Validation plots and diagnostic information

### Configuration Parameters

**Tracer Types**: LRG, BGS, ELG, etc.
**Position Types**: observed, reconstructed, rsd_removed  
**Shape Types**: observed, reconstructed, rsd_removed
**Correlation Types**: NG (galaxy-shear), VG (velocity-shear), NV (count-velocity), NR (galaxy-size)

## Development Notes

- Uses SLURM for job submission on HPC cluster
- TreeCorr library handles correlation function calculations
- Astropy for coordinate transformations and cosmology
- Results are saved as NumPy arrays for analysis in notebooks
- Code assumes specific data directory structure and file naming conventions