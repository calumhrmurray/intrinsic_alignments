#!/usr/bin/env python3
"""
Bin displacement-shear correlation functions to produce monopoles.
Based on bin_correlation_functions.ipynb
"""
import numpy as np
import argparse
import os
import glob

def geometric_weight_function(rper, rpar):
    """Geometric weighting that emphasizes perpendicular separations."""
    return rper**2 / (rpar**2 + rper**2)

def no_weights(rper, rpar):
    """No weighting - uniform weights."""
    return np.ones(rper.shape)

def bin_correlation_function(rper, rpar, xi, var_xi, bins, custom_weight_func, epsilon=1e-10):
    """
    Calculate the binned correlation function xi(r) as a weighted sum, and its error.

    Parameters:
    -----------
    rper : array
        Perpendicular separation bins
    rpar : array  
        Parallel separation bins
    xi : array
        2D correlation function xi(rper, rpar)
    var_xi : array
        Variance of 2D correlation function
    bins : array
        Radial bins for monopole
    custom_weight_func : function
        Weighting function
    epsilon : float
        Small number to avoid division by zero
        
    Returns:
    --------
    r_centers : array
        Bin centers
    xi_r : array
        Weighted xi(r) in each bin
    xi_err : array
        Error on xi(r) in each bin
    weights : array
        Weight grid used
    """

    rper_2d, rpar_2d = np.meshgrid(rper, rpar, indexing='ij')
    r = np.sqrt(rper_2d**2 + rpar_2d**2).T
    w = custom_weight_func(rper_2d, rpar_2d).T / (var_xi + epsilon)

    r_flat = r.flatten()
    xi_flat = xi.flatten()
    w_flat = w.flatten()

    xi_r = []
    xi_var = []
    r_centers = 0.5 * (bins[:-1] + bins[1:])
    
    for i in range(len(bins) - 1):
        mask = (r_flat >= bins[i]) & (r_flat < bins[i+1])
        if np.any(mask):
            xi_r.append(np.average(xi_flat[mask], weights=w_flat[mask]))
            xi_var.append(1.0 / np.sum(w_flat[mask]))
        else:
            xi_r.append(np.nan)
            xi_var.append(np.nan)

    return r_centers, np.array(xi_r), np.array(xi_var), custom_weight_func(rper_2d, rpar_2d)

def calculate_results(xi_p, xi_x, var_xi, bins, custom_weight_func, rper, rpar):
    """Calculate binned results for both plus and cross components."""
    
    r_centers, xi_p_binned, xi_p_var, weights = bin_correlation_function(
        rper, rpar, xi_p, var_xi,
        bins=bins,
        custom_weight_func=custom_weight_func
    )

    r_centers, xi_x_binned, xi_x_var, weights = bin_correlation_function(
        rper, rpar, xi_x, var_xi,
        bins=bins,
        custom_weight_func=custom_weight_func
    )

    return r_centers, xi_p_binned, xi_p_var, xi_x_binned, xi_x_var, weights

def find_displacement_shear_files(data_dir):
    """Find displacement-shear correlation files in data directory."""
    pattern = os.path.join(data_dir, "displacement_shear_*_xi_plus.npy")
    plus_files = glob.glob(pattern)
    
    results = []
    for plus_file in plus_files:
        base = plus_file.replace('_xi_plus.npy', '')
        cross_file = base + '_xi_cross.npy'
        rperp_file = base + '_rperp.npy'
        rpar_file = base + '_rpar_bins.npy'
        
        if all(os.path.exists(f) for f in [cross_file, rperp_file, rpar_file]):
            results.append({
                'base': base,
                'xi_plus': plus_file,
                'xi_cross': cross_file, 
                'rperp': rperp_file,
                'rpar_bins': rpar_file
            })
    
    return results

def process_displacement_shear_correlation(file_info, output_dir):
    """Process a single displacement-shear correlation measurement."""
    
    print(f"Processing: {file_info['base']}")
    
    # Load data
    xi_plus = np.load(file_info['xi_plus'])
    xi_cross = np.load(file_info['xi_cross']) 
    rperp = np.load(file_info['rperp'])
    rpar_bins = np.load(file_info['rpar_bins'])
    
    print(f"  Data shape: xi_plus={xi_plus.shape}, xi_cross={xi_cross.shape}")
    print(f"  rperp bins: {len(rperp)}, rpar bins: {len(rpar_bins)-1}")
    
    # Create rpar centers
    rpar = 0.5 * (rpar_bins[:-1] + rpar_bins[1:])
    
    # Create simple variance estimate (can be improved with proper error propagation)
    var_xi = np.ones_like(xi_plus) * 1e-6  # Placeholder - should use proper variance
    
    # Define radial bins for monopole - logarithmic binning
    #r_bins = np.logspace(np.log10(0.5), np.log10(60.0), 16)  # 10 bins from 0.1 to 60 Mpc/h
    r_bins = np.linspace(0.1, 60, 12)  # 11 bins from 0.5 to 60 Mpc/h

    # Calculate results with different weighting schemes
    print("  Calculating monopole with geometric weighting...")
    r_centers, xi_p_geom, xi_p_var_geom, xi_x_geom, xi_x_var_geom, weights_geom = calculate_results(
        xi_plus, xi_cross, var_xi, r_bins, geometric_weight_function, rperp, rpar
    )
    
    print("  Calculating monopole with no weighting...")
    r_centers, xi_p_flat, xi_p_var_flat, xi_x_flat, xi_x_var_flat, weights_flat = calculate_results(
        xi_plus, xi_cross, var_xi, r_bins, no_weights, rperp, rpar
    )
    
    # Save results
    base_name = os.path.basename(file_info['base'])
    output_base = os.path.join(output_dir, base_name + '_monopole')
    
    # Save geometric weighted results
    np.save(output_base + '_r_centers.npy', r_centers)
    np.save(output_base + '_xi_plus_geometric.npy', xi_p_geom)
    np.save(output_base + '_xi_plus_var_geometric.npy', xi_p_var_geom)
    np.save(output_base + '_xi_cross_geometric.npy', xi_x_geom)
    np.save(output_base + '_xi_cross_var_geometric.npy', xi_x_var_geom)
    
    # Save flat weighted results
    np.save(output_base + '_xi_plus_flat.npy', xi_p_flat)
    np.save(output_base + '_xi_plus_var_flat.npy', xi_p_var_flat)
    np.save(output_base + '_xi_cross_flat.npy', xi_x_flat)
    np.save(output_base + '_xi_cross_var_flat.npy', xi_x_var_flat)
    
    print(f"  Results saved to: {output_base}_*.npy")
    
    return {
        'r_centers': r_centers,
        'xi_plus_geometric': xi_p_geom,
        'xi_plus_var_geometric': xi_p_var_geom,
        'xi_cross_geometric': xi_x_geom,
        'xi_cross_var_geometric': xi_x_var_geom,
        'xi_plus_flat': xi_p_flat,
        'xi_plus_var_flat': xi_p_var_flat,
        'xi_cross_flat': xi_x_flat,
        'xi_cross_var_flat': xi_x_var_flat
    }

def main():
    parser = argparse.ArgumentParser(description='Bin displacement-shear correlations to monopoles')
    parser.add_argument('--data_dir', required=True,
                       help='Directory containing displacement-shear correlation files')
    parser.add_argument('--output_dir', required=True,
                       help='Directory to save monopole results')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Find displacement-shear correlation files
    correlation_files = find_displacement_shear_files(args.data_dir)
    
    if not correlation_files:
        print(f"No displacement-shear correlation files found in {args.data_dir}")
        print("Looking for files matching: displacement_shear_*_xi_plus.npy")
        return
    
    print(f"Found {len(correlation_files)} displacement-shear correlation measurements:")
    for file_info in correlation_files:
        print(f"  {os.path.basename(file_info['base'])}")
    
    # Process each correlation measurement
    all_results = {}
    for file_info in correlation_files:
        try:
            result = process_displacement_shear_correlation(file_info, args.output_dir)
            key = os.path.basename(file_info['base'])
            all_results[key] = result
        except Exception as e:
            print(f"Error processing {file_info['base']}: {e}")
    
    print(f"\nProcessing complete. Results saved to {args.output_dir}")
    print(f"Processed {len(all_results)} correlation measurements successfully.")

if __name__ == "__main__":
    main()