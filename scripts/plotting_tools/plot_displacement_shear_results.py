#!/usr/bin/env python3
"""
Create figures from displacement-shear correlation monopole results.
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import glob

def find_monopole_files(data_dir):
    """Find monopole result files in data directory."""
    pattern = os.path.join(data_dir, "*_monopole_r_centers.npy")
    r_files = glob.glob(pattern)
    
    results = []
    for r_file in r_files:
        base = r_file.replace('_r_centers.npy', '')
        
        required_files = [
            base + '_xi_plus_geometric.npy',
            base + '_xi_plus_var_geometric.npy', 
            base + '_xi_cross_geometric.npy',
            base + '_xi_cross_var_geometric.npy',
            base + '_xi_plus_flat.npy',
            base + '_xi_plus_var_flat.npy',
            base + '_xi_cross_flat.npy',
            base + '_xi_cross_var_flat.npy'
        ]
        
        if all(os.path.exists(f) for f in required_files):
            results.append(base)
    
    return results

def load_monopole_results(base_path):
    """Load monopole results from files."""
    return {
        'r_centers': np.load(base_path + '_r_centers.npy'),
        'xi_plus_geometric': np.load(base_path + '_xi_plus_geometric.npy'),
        'xi_plus_var_geometric': np.load(base_path + '_xi_plus_var_geometric.npy'),
        'xi_cross_geometric': np.load(base_path + '_xi_cross_geometric.npy'),
        'xi_cross_var_geometric': np.load(base_path + '_xi_cross_var_geometric.npy'),
        'xi_plus_flat': np.load(base_path + '_xi_plus_flat.npy'),
        'xi_plus_var_flat': np.load(base_path + '_xi_plus_var_flat.npy'),
        'xi_cross_flat': np.load(base_path + '_xi_cross_flat.npy'),
        'xi_cross_var_flat': np.load(base_path + '_xi_cross_var_flat.npy')
    }

def create_monopole_comparison_plot(results_dict, output_dir):
    """Create comparison plot of monopole results."""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    colors = plt.cm.Set1(np.linspace(0, 1, len(results_dict)))
    
    for i, (name, results) in enumerate(results_dict.items()):
        r = results['r_centers']
        color = colors[i]
        
        # Plot 1: xi_plus with geometric weighting
        axes[0,0].errorbar(r, results['xi_plus_geometric'], 
                          yerr=np.sqrt(results['xi_plus_var_geometric']),
                          label=name, color=color, fmt='o-', alpha=0.8)
        
        # Plot 2: xi_cross with geometric weighting
        axes[0,1].errorbar(r, results['xi_cross_geometric'],
                          yerr=np.sqrt(results['xi_cross_var_geometric']), 
                          label=name, color=color, fmt='s-', alpha=0.8)
        
        # Plot 3: Comparison of weighting schemes for xi_plus
        axes[1,0].errorbar(r, results['xi_plus_geometric'],
                          yerr=np.sqrt(results['xi_plus_var_geometric']),
                          label=f'{name} (geometric)', color=color, fmt='o-', alpha=0.8)
        axes[1,0].errorbar(r*1.05, results['xi_plus_flat'],
                          yerr=np.sqrt(results['xi_plus_var_flat']),
                          label=f'{name} (flat)', color=color, fmt='s--', alpha=0.6)
        
        # Plot 4: r^2 * xi_plus with geometric weighting
        axes[1,1].errorbar(r, r**2 * results['xi_plus_geometric'],
                          yerr=r**2 * np.sqrt(results['xi_plus_var_geometric']),
                          label=name, color=color, fmt='o-', alpha=0.8)
    
    # Format plots
    axes[0,0].set_xlabel('r [Mpc/h]')
    axes[0,0].set_ylabel(r'$\xi_{+}(r)$')
    axes[0,0].set_title('Displacement-Shear Plus Component')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    axes[0,0].set_xscale('log')
    
    axes[0,1].set_xlabel('r [Mpc/h]')
    axes[0,1].set_ylabel(r'$\xi_{\times}(r)$')
    axes[0,1].set_title('Displacement-Shear Cross Component')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    axes[0,1].set_xscale('log')
    
    axes[1,0].set_xlabel('r [Mpc/h]')
    axes[1,0].set_ylabel(r'$\xi_{+}(r)$')
    axes[1,0].set_title('Weighting Scheme Comparison')
    axes[1,0].legend(fontsize=8)
    axes[1,0].grid(True, alpha=0.3)
    axes[1,0].set_xscale('log')
    
    axes[1,1].set_xlabel('r [Mpc/h]')
    axes[1,1].set_ylabel(r'$r^2 \xi_{+}(r)$')
    axes[1,1].set_title('Plus Component × r²')
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)
    axes[1,1].set_xscale('log')
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, 'displacement_shear_monopole_comparison.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Monopole comparison plot saved: {output_file}")

def create_individual_plots(results_dict, output_dir):
    """Create individual plots for each result."""
    
    for name, results in results_dict.items():
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        r = results['r_centers']
        
        # Plot 1: Plus component (linear scale)
        axes[0,0].errorbar(r, results['xi_plus_geometric'],
                          yerr=np.sqrt(results['xi_plus_var_geometric']),
                          fmt='ro-', label='Geometric weight')
        axes[0,0].errorbar(r*1.05, results['xi_plus_flat'],
                          yerr=np.sqrt(results['xi_plus_var_flat']),
                          fmt='bs--', alpha=0.7, label='Flat weight')
        axes[0,0].set_xlabel('r [Mpc/h]')
        axes[0,0].set_ylabel(r'$\xi_{+}(r)$')
        axes[0,0].set_title('Plus Component')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # Plot 2: Cross component (linear scale)
        axes[0,1].errorbar(r, results['xi_cross_geometric'],
                          yerr=np.sqrt(results['xi_cross_var_geometric']),
                          fmt='ro-', label='Geometric weight')
        axes[0,1].errorbar(r*1.05, results['xi_cross_flat'],
                          yerr=np.sqrt(results['xi_cross_var_flat']),
                          fmt='bs--', alpha=0.7, label='Flat weight')
        axes[0,1].set_xlabel('r [Mpc/h]')
        axes[0,1].set_ylabel(r'$\xi_{\times}(r)$')
        axes[0,1].set_title('Cross Component')
        axes[0,1].legend()
        axes[0,1].grid(True, alpha=0.3)
        
        # Plot 3: Log-log plot
        axes[1,0].errorbar(r, np.abs(results['xi_plus_geometric']),
                          yerr=np.sqrt(results['xi_plus_var_geometric']),
                          fmt='ro-', label='Plus (geometric)')
        axes[1,0].errorbar(r*1.05, np.abs(results['xi_cross_geometric']),
                          yerr=np.sqrt(results['xi_cross_var_geometric']),
                          fmt='bs--', alpha=0.7, label='Cross (geometric)')
        axes[1,0].set_xlabel('r [Mpc/h]')
        axes[1,0].set_ylabel(r'$|\xi(r)|$')
        axes[1,0].set_title('Magnitude (Log-Log)')
        axes[1,0].legend()
        axes[1,0].grid(True, alpha=0.3)
        axes[1,0].set_xscale('log')
        axes[1,0].set_yscale('log')
        
        # Plot 4: r^2 scaling
        axes[1,1].errorbar(r, r**2 * results['xi_plus_geometric'],
                          yerr=r**2 * np.sqrt(results['xi_plus_var_geometric']),
                          fmt='ro-', label='Plus × r²')
        axes[1,1].errorbar(r*1.05, r**2 * results['xi_cross_geometric'],
                          yerr=r**2 * np.sqrt(results['xi_cross_var_geometric']),
                          fmt='bs--', alpha=0.7, label='Cross × r²')
        axes[1,1].set_xlabel('r [Mpc/h]')
        axes[1,1].set_ylabel(r'$r^2 \xi(r)$')
        axes[1,1].set_title('Scaled by r²')
        axes[1,1].legend()
        axes[1,1].grid(True, alpha=0.3)
        axes[1,1].set_xscale('log')
        
        plt.suptitle(f'Displacement-Shear Correlation: {name}', fontsize=14)
        plt.tight_layout()
        
        # Save individual plot
        safe_name = name.replace('/', '_').replace(' ', '_')
        output_file = os.path.join(output_dir, f'displacement_shear_{safe_name}.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Individual plot saved: {output_file}")

def create_signal_to_noise_plot(results_dict, output_dir):
    """Create signal-to-noise ratio plots."""
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    colors = plt.cm.Set1(np.linspace(0, 1, len(results_dict)))
    
    for i, (name, results) in enumerate(results_dict.items()):
        r = results['r_centers']
        color = colors[i]
        
        # Calculate S/N ratios
        snr_plus = np.abs(results['xi_plus_geometric']) / np.sqrt(results['xi_plus_var_geometric'])
        snr_cross = np.abs(results['xi_cross_geometric']) / np.sqrt(results['xi_cross_var_geometric'])
        
        axes[0].plot(r, snr_plus, 'o-', color=color, label=name, alpha=0.8)
        axes[1].plot(r, snr_cross, 's-', color=color, label=name, alpha=0.8)
    
    axes[0].set_xlabel('r [Mpc/h]')
    axes[0].set_ylabel('Signal-to-Noise Ratio')
    axes[0].set_title('Plus Component S/N')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_xscale('log')
    
    axes[1].set_xlabel('r [Mpc/h]')
    axes[1].set_ylabel('Signal-to-Noise Ratio')
    axes[1].set_title('Cross Component S/N')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_xscale('log')
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, 'displacement_shear_signal_to_noise.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Signal-to-noise plot saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Plot displacement-shear correlation results')
    parser.add_argument('--data_dir', required=True,
                       help='Directory containing monopole result files')
    parser.add_argument('--output_dir', default='/home/murray/intrinsic_alignments/../results/figures/shear_displacement_correlations',
                       help='Directory to save figures')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Find monopole result files
    monopole_files = find_monopole_files(args.data_dir)
    
    if not monopole_files:
        print(f"No monopole result files found in {args.data_dir}")
        print("Looking for files matching: *_monopole_r_centers.npy")
        return
    
    print(f"Found {len(monopole_files)} monopole result sets:")
    for base_path in monopole_files:
        print(f"  {os.path.basename(base_path)}")
    
    # Load all results
    results_dict = {}
    for base_path in monopole_files:
        try:
            name = os.path.basename(base_path).replace('_monopole', '')
            results = load_monopole_results(base_path)
            results_dict[name] = results
            print(f"Loaded: {name}")
        except Exception as e:
            print(f"Error loading {base_path}: {e}")
    
    if not results_dict:
        print("No results could be loaded successfully.")
        return
    
    print(f"\nCreating figures for {len(results_dict)} result sets...")
    
    # Create plots
    try:
        create_monopole_comparison_plot(results_dict, args.output_dir)
        create_individual_plots(results_dict, args.output_dir)
        create_signal_to_noise_plot(results_dict, args.output_dir)
        
        print(f"\nAll figures saved to: {args.output_dir}")
        
    except Exception as e:
        print(f"Error creating plots: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()