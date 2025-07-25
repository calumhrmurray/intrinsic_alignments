#!/usr/bin/env python3
"""
Simple plotting script for correlation results with current file naming.
"""
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import argparse

def plot_displacement_shear_results(data_dir, output_dir):
    """Plot displacement-shear correlation results."""
    
    # Find displacement-shear result files
    pattern = os.path.join(data_dir, "displacement_shear_*_xi_plus.npy")
    xi_plus_files = glob.glob(pattern)
    
    if not xi_plus_files:
        print(f"No displacement-shear files found in {data_dir}")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    for xi_plus_file in xi_plus_files:
        base_name = xi_plus_file.replace('_xi_plus.npy', '')
        name = os.path.basename(base_name)
        
        # Load data files
        try:
            xi_plus = np.load(xi_plus_file)
            xi_cross = np.load(base_name + '_xi_cross.npy')
            rperp = np.load(base_name + '_rperp.npy')
            
            # Handle 2D results (sum over r_parallel for monopole)
            if xi_plus.ndim == 2:
                xi_plus = np.mean(xi_plus, axis=0)  # Average over rpar bins
                xi_cross = np.mean(xi_cross, axis=0)
            
            # Create plots
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            
            # Plot 1: Plus component
            axes[0,0].plot(rperp, xi_plus, 'ro-', label='Plus component')
            axes[0,0].set_xlabel('r_perp [Mpc/h]')
            axes[0,0].set_ylabel(r'$\xi_{+}(r_\perp)$')
            axes[0,0].set_title('Displacement-Shear Plus Component')
            axes[0,0].grid(True, alpha=0.3)
            axes[0,0].set_xscale('log')
            
            # Plot 2: Cross component  
            axes[0,1].plot(rperp, xi_cross, 'bs-', label='Cross component')
            axes[0,1].set_xlabel('r_perp [Mpc/h]')
            axes[0,1].set_ylabel(r'$\xi_{\times}(r_\perp)$')
            axes[0,1].set_title('Displacement-Shear Cross Component')
            axes[0,1].grid(True, alpha=0.3)
            axes[0,1].set_xscale('log')
            
            # Plot 3: Both components
            axes[1,0].plot(rperp, xi_plus, 'ro-', label='Plus', alpha=0.8)
            axes[1,0].plot(rperp, xi_cross, 'bs-', label='Cross', alpha=0.8)
            axes[1,0].set_xlabel('r_perp [Mpc/h]')
            axes[1,0].set_ylabel(r'$\xi(r_\perp)$')
            axes[1,0].set_title('Both Components')
            axes[1,0].legend()
            axes[1,0].grid(True, alpha=0.3)
            axes[1,0].set_xscale('log')
            
            # Plot 4: r^2 scaling
            axes[1,1].plot(rperp, rperp**2 * xi_plus, 'ro-', label='Plus × r²', alpha=0.8)
            axes[1,1].plot(rperp, rperp**2 * xi_cross, 'bs-', label='Cross × r²', alpha=0.8)
            axes[1,1].set_xlabel('r_perp [Mpc/h]')
            axes[1,1].set_ylabel(r'$r_\perp^2 \xi(r_\perp)$')
            axes[1,1].set_title('Scaled by r²')
            axes[1,1].legend()
            axes[1,1].grid(True, alpha=0.3)
            axes[1,1].set_xscale('log')
            
            plt.suptitle(f'Displacement-Shear Correlation: {name}', fontsize=14)
            plt.tight_layout()
            
            # Save plot
            output_file = os.path.join(output_dir, f'{name}_plot.png')
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"Plot saved: {output_file}")
            
        except Exception as e:
            print(f"Error plotting {name}: {e}")

def plot_other_correlations(data_dir, output_dir):
    """Plot other correlation types."""
    
    # Look for galaxy-shear (NG) correlations
    ng_files = glob.glob(os.path.join(data_dir, "*_xi_gn_p_*.npy"))
    
    if ng_files:
        print(f"Found {len(ng_files)} galaxy-shear correlation files")
        
        for xi_file in ng_files[:3]:  # Plot first 3 as examples
            try:
                base_name = xi_file.replace('_xi_gn_p_', '_').replace('.npy', '')
                name = os.path.basename(base_name)
                
                xi_p = np.load(xi_file)
                xi_x_file = xi_file.replace('_xi_gn_p_', '_xi_gn_x_')
                r_file = xi_file.replace('_xi_gn_p_', '_r_')
                
                if os.path.exists(xi_x_file) and os.path.exists(r_file):
                    xi_x = np.load(xi_x_file)
                    r = np.load(r_file)
                    
                    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
                    
                    axes[0].plot(r, xi_p, 'ro-', label='Plus component')
                    axes[0].set_xlabel('r [Mpc/h]')
                    axes[0].set_ylabel(r'$\xi_{+}(r)$')
                    axes[0].set_title('Galaxy-Shear Plus Component')
                    axes[0].grid(True, alpha=0.3)
                    axes[0].set_xscale('log')
                    
                    axes[1].plot(r, xi_x, 'bs-', label='Cross component')
                    axes[1].set_xlabel('r [Mpc/h]')
                    axes[1].set_ylabel(r'$\xi_{\times}(r)$')
                    axes[1].set_title('Galaxy-Shear Cross Component')
                    axes[1].grid(True, alpha=0.3)
                    axes[1].set_xscale('log')
                    
                    plt.suptitle(f'Galaxy-Shear Correlation: {name}', fontsize=14)
                    plt.tight_layout()
                    
                    output_file = os.path.join(output_dir, f'galaxy_shear_{name}_plot.png')
                    plt.savefig(output_file, dpi=150, bbox_inches='tight')
                    plt.close()
                    
                    print(f"Galaxy-shear plot saved: {output_file}")
                    
            except Exception as e:
                print(f"Error plotting galaxy-shear {xi_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Plot correlation results')
    parser.add_argument('--data_dir', default='/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/correlation_function_measurements_2mpc/',
                       help='Directory containing result files')
    parser.add_argument('--output_dir', default='/home/murray/intrinsic_alignments/../results/figures/correlation_functions',
                       help='Directory to save figures')
    
    args = parser.parse_args()
    
    print(f"Plotting results from: {args.data_dir}")
    print(f"Saving figures to: {args.output_dir}")
    
    # Plot displacement-shear correlations
    plot_displacement_shear_results(args.data_dir, args.output_dir)
    
    # Plot other correlations
    plot_other_correlations(args.data_dir, args.output_dir)
    
    print(f"\nAll plots saved to: {args.output_dir}")

if __name__ == "__main__":
    main()