#!/usr/bin/env python3
"""
Improved validation of LRG displacement catalogues with better visualizations.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import os
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

def load_displacement_catalogues():
    """Load displacement catalogues."""
    print("Loading displacement catalogues...")
    
    velocity_folder = '/n17data/murray/desi_data/DESI/velocity_catalogues/'
    
    # Load position displacements
    disp_file = velocity_folder + 'LRG_displacements.fits'
    if os.path.exists(disp_file):
        displacements = fits.open(disp_file)[1].data
    else:
        print(f"Error: Displacement file not found: {disp_file}")
        return None, None
    
    # Load shape displacements
    shape_disp_file = velocity_folder + 'LRG_shape_displacements.fits'
    if os.path.exists(shape_disp_file):
        shape_displacements = fits.open(shape_disp_file)[1].data
    else:
        print(f"Error: Shape displacement file not found: {shape_disp_file}")
        return displacements, None
    
    return displacements, shape_displacements

def plot_displacement_vectors_on_sky(displacements, shape_displacements):
    """Plot displacement vectors as arrows on a small sky patch."""
    print("\n=== PLOTTING DISPLACEMENT VECTORS ON SKY PATCH ===")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    if displacements is not None:
        # Select a small sky patch (~20 deg^2)
        # Find a region with good coverage
        ra_center = 180.0  # Example center
        dec_center = 30.0
        patch_size = 4.5  # degrees (roughly 20 deg^2)
        
        # Filter to small patch
        ra_mask = (displacements['RA'] > ra_center - patch_size/2) & (displacements['RA'] < ra_center + patch_size/2)
        dec_mask = (displacements['Dec'] > dec_center - patch_size/2) & (displacements['Dec'] < dec_center + patch_size/2)
        patch_mask = ra_mask & dec_mask
        
        if np.sum(patch_mask) < 50:
            # Try different patch centers if first one has too few galaxies
            for ra_test, dec_test in [(150, 20), (200, 40), (220, 10), (160, 50)]:
                ra_mask = (displacements['RA'] > ra_test - patch_size/2) & (displacements['RA'] < ra_test + patch_size/2)
                dec_mask = (displacements['Dec'] > dec_test - patch_size/2) & (displacements['Dec'] < dec_test + patch_size/2)
                test_mask = ra_mask & dec_mask
                if np.sum(test_mask) > 50:
                    patch_mask = test_mask
                    ra_center, dec_center = ra_test, dec_test
                    break
        
        # Extract patch data
        ra_patch = displacements['RA'][patch_mask]
        dec_patch = displacements['Dec'][patch_mask]
        dr1_patch = displacements['dr1'][patch_mask]
        dr2_patch = displacements['dr2'][patch_mask]
        
        # Further downsample if too many points (keep every Nth)
        if len(ra_patch) > 500:
            step = len(ra_patch) // 300  # Keep ~300 vectors
            subset = slice(None, None, step)
            ra_patch = ra_patch[subset]
            dec_patch = dec_patch[subset]
            dr1_patch = dr1_patch[subset]
            dr2_patch = dr2_patch[subset]
        
        # Scale displacement for visibility (larger scale for small patch)
        scale_factor = 0.05  # Increased scale for better visibility
        dra = dr1_patch * scale_factor
        ddec = dr2_patch * scale_factor
        
        # Plot displacement vectors
        axes[0].quiver(ra_patch, dec_patch, dra, ddec, 
                      angles='xy', scale_units='xy', scale=1,
                      alpha=0.8, width=0.003, headwidth=4, color='red', 
                      label=f'Displacement vectors (×{scale_factor:.2f})')
        axes[0].scatter(ra_patch, dec_patch, s=8, alpha=0.6, c='blue', label='Galaxy positions')
        axes[0].set_xlabel('RA (deg)')
        axes[0].set_ylabel('Dec (deg)')
        axes[0].set_title(f'Position Displacements - Sky Patch ({len(ra_patch)} galaxies)\nCenter: RA={ra_center}°, Dec={dec_center}°')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        axes[0].set_aspect('equal')
        
        print(f"Position displacements: {len(ra_patch)} vectors in {patch_size:.1f}° × {patch_size:.1f}° patch")
    
    if shape_displacements is not None:
        # Same patch selection for shapes
        ra_mask_s = (shape_displacements['RA'] > ra_center - patch_size/2) & (shape_displacements['RA'] < ra_center + patch_size/2)
        dec_mask_s = (shape_displacements['Dec'] > dec_center - patch_size/2) & (shape_displacements['Dec'] < dec_center + patch_size/2)
        patch_mask_s = ra_mask_s & dec_mask_s
        
        ra_patch_s = shape_displacements['RA'][patch_mask_s]
        dec_patch_s = shape_displacements['Dec'][patch_mask_s]
        dr1_patch_s = shape_displacements['dr1'][patch_mask_s]
        dr2_patch_s = shape_displacements['dr2'][patch_mask_s]
        
        # Downsample if needed
        if len(ra_patch_s) > 200:
            step_s = len(ra_patch_s) // 150
            subset_s = slice(None, None, step_s)
            ra_patch_s = ra_patch_s[subset_s]
            dec_patch_s = dec_patch_s[subset_s]
            dr1_patch_s = dr1_patch_s[subset_s]
            dr2_patch_s = dr2_patch_s[subset_s]
        
        # Scale for visibility
        scale_factor = 0.05
        dra_s = dr1_patch_s * scale_factor
        ddec_s = dr2_patch_s * scale_factor
        
        # Plot shape displacement vectors
        axes[1].quiver(ra_patch_s, dec_patch_s, dra_s, ddec_s,
                      angles='xy', scale_units='xy', scale=1,
                      alpha=0.8, width=0.003, headwidth=4, color='green',
                      label=f'Displacement vectors (×{scale_factor:.2f})')
        axes[1].scatter(ra_patch_s, dec_patch_s, s=10, alpha=0.6, c='orange', label='Shape positions')
        axes[1].set_xlabel('RA (deg)')
        axes[1].set_ylabel('Dec (deg)')
        axes[1].set_title(f'Shape Displacements - Sky Patch ({len(ra_patch_s)} shapes)\nCenter: RA={ra_center}°, Dec={dec_center}°')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        axes[1].set_aspect('equal')
        
        print(f"Shape displacements: {len(ra_patch_s)} vectors in {patch_size:.1f}° × {patch_size:.1f}° patch")
    
    plt.tight_layout()
    plt.savefig('/home/murray/intrinsic_alignments/validation/create_catalogues/displacement_vectors_sky_patch.png', dpi=150, bbox_inches='tight')
    print("Saved: /home/murray/intrinsic_alignments/validation/create_catalogues/displacement_vectors_sky_patch.png")

def plot_physically_accurate_displacement_vectors(displacements, shape_displacements):
    """Plot displacement vectors with correct physical angular lengths in redshift bins."""
    print("\n=== PLOTTING PHYSICALLY ACCURATE DISPLACEMENT VECTORS BY REDSHIFT ===")
    
    # Set up cosmology (same as in measure_correlations.py)
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    if displacements is not None:
        # Much larger patch for better statistics
        ra_center = 180.0
        dec_center = 30.0
        patch_size = 9.5  # Much larger patch (~90 deg²)
        
        # Filter to larger patch
        ra_mask = (displacements['RA'] > ra_center - patch_size/2) & (displacements['RA'] < ra_center + patch_size/2)
        dec_mask = (displacements['Dec'] > dec_center - patch_size/2) & (displacements['Dec'] < dec_center + patch_size/2)
        patch_mask = ra_mask & dec_mask
        
        # Try different patch centers if needed
        if np.sum(patch_mask) < 200:
            for ra_test, dec_test in [(150, 20), (200, 40), (220, 10), (160, 50)]:
                ra_mask = (displacements['RA'] > ra_test - patch_size/2) & (displacements['RA'] < ra_test + patch_size/2)
                dec_mask = (displacements['Dec'] > dec_test - patch_size/2) & (displacements['Dec'] < dec_test + patch_size/2)
                test_mask = ra_mask & dec_mask
                if np.sum(test_mask) > 200:
                    patch_mask = test_mask
                    ra_center, dec_center = ra_test, dec_test
                    break
        
        # Extract patch data
        ra_patch = displacements['RA'][patch_mask]
        dec_patch = displacements['Dec'][patch_mask]
        z_patch = displacements['redshift'][patch_mask]
        dr1_patch = displacements['dr1'][patch_mask]  # Mpc/h
        dr2_patch = displacements['dr2'][patch_mask]  # Mpc/h
        
        print(f"Found {len(ra_patch)} galaxies in {patch_size:.1f}° × {patch_size:.1f}° patch")
        
        # Create 4 redshift bins
        z_bins = np.percentile(z_patch, [0, 25, 50, 75, 100])
        print(f"Redshift bins: {z_bins}")
        
        for bin_idx in range(4):
            z_min, z_max = z_bins[bin_idx], z_bins[bin_idx + 1]
            z_mask = (z_patch >= z_min) & (z_patch < z_max if bin_idx < 3 else z_patch <= z_max)
            
            ra_bin = ra_patch[z_mask]
            dec_bin = dec_patch[z_mask]
            z_bin = z_patch[z_mask]
            dr1_bin = dr1_patch[z_mask]
            dr2_bin = dr2_patch[z_mask]
            
            if len(ra_bin) == 0:
                continue
                
            # Downsample if too many points
            if len(ra_bin) > 200:
                step = len(ra_bin) // 150
                subset = slice(None, None, step)
                ra_bin = ra_bin[subset]
                dec_bin = dec_bin[subset]
                z_bin = z_bin[subset]
                dr1_bin = dr1_bin[subset]
                dr2_bin = dr2_bin[subset]
        
            # Convert physical displacements to angular displacements for this redshift bin
            da_bin = cosmo.angular_diameter_distance(z_bin).to(u.Mpc).value  # Mpc (not h^-1)
            
            # Convert Mpc/h to Mpc (assuming h=0.7 as in the cosmology)
            h = 0.7
            dr1_mpc_bin = dr1_bin / h  # Convert Mpc/h to Mpc for each galaxy
            dr2_mpc_bin = dr2_bin / h
            
            # Convert to angular coordinates: θ = physical_size / angular_diameter_distance
            dra_radians_bin = dr1_mpc_bin / da_bin  # radians per galaxy
            ddec_radians_bin = dr2_mpc_bin / da_bin  # radians per galaxy
            
            # Convert to degrees for matplotlib plotting
            dra_physical_bin = np.degrees(dra_radians_bin)  # degrees per galaxy  
            ddec_physical_bin = np.degrees(ddec_radians_bin)  # degrees per galaxy
            
            # Calculate typical sizes in arcminutes
            disp_mag_arcmin_bin = np.sqrt(dra_radians_bin**2 + ddec_radians_bin**2) * 180/np.pi * 60
            
            # Plot displacement as simple lines from galaxy position to displaced position
            for i in range(len(ra_bin)):
                axes[bin_idx].plot([ra_bin[i], ra_bin[i] + dra_physical_bin[i]], 
                            [dec_bin[i], dec_bin[i] + ddec_physical_bin[i]], 
                            'r-', alpha=0.7, linewidth=1)
            
            # Plot galaxy positions as markers
            axes[bin_idx].scatter(ra_bin, dec_bin, s=12, alpha=0.7, c='blue', label='Galaxy positions')
            
            # Add a few example displacement endpoints for clarity
            n_examples = min(15, len(ra_bin))
            axes[bin_idx].scatter(ra_bin[:n_examples] + dra_physical_bin[:n_examples], 
                           dec_bin[:n_examples] + ddec_physical_bin[:n_examples], 
                           s=8, alpha=0.5, c='red', marker='x', label=f'Displacement endpoints (first {n_examples})')
            
            axes[bin_idx].set_xlabel('RA (deg)')
            axes[bin_idx].set_ylabel('Dec (deg)')
            axes[bin_idx].set_title(f'Redshift Bin {bin_idx+1}: {z_min:.2f} < z < {z_max:.2f}\n'
                              f'{len(ra_bin)} galaxies, ~{np.mean(disp_mag_arcmin_bin):.1f} arcmin, D_A={np.mean(da_bin):.0f} Mpc')
            axes[bin_idx].legend(fontsize=8)
            axes[bin_idx].grid(True, alpha=0.3)
            axes[bin_idx].set_aspect('equal')
            
            print(f"Bin {bin_idx+1} (z={z_min:.2f}-{z_max:.2f}): {len(ra_bin)} galaxies, "
                  f"mean displacement {np.mean(disp_mag_arcmin_bin):.1f} arcmin")
        
    
    plt.tight_layout()
    plt.savefig('/home/murray/intrinsic_alignments/validation/create_catalogues/displacement_vectors_physical_scale.png', dpi=150, bbox_inches='tight')
    print("Saved: /home/murray/intrinsic_alignments/validation/create_catalogues/displacement_vectors_physical_scale.png")

def plot_displacement_distributions(displacements, shape_displacements):
    """Plot displacement distributions as line plots."""
    print("\n=== DISPLACEMENT DISTRIBUTION ANALYSIS ===")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    if displacements is not None:
        # Position displacement components
        dr1 = displacements['dr1']
        dr2 = displacements['dr2'] 
        drp = displacements['drp']
        
        # Create line plots instead of histograms
        bins = np.linspace(-50, 50, 100)
        
        hist_dr1, _ = np.histogram(dr1, bins=bins, density=True)
        hist_dr2, _ = np.histogram(dr2, bins=bins, density=True)
        hist_drp, _ = np.histogram(drp, bins=bins, density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        
        axes[0,0].plot(bin_centers, hist_dr1, linewidth=2, color='blue')
        axes[0,0].set_xlabel('dr1 (Mpc/h)')
        axes[0,0].set_ylabel('Density')
        axes[0,0].set_title('Position Displacement dr1')
        axes[0,0].axvline(0, color='red', linestyle='--', alpha=0.7)
        axes[0,0].grid(True, alpha=0.3)
        
        axes[0,1].plot(bin_centers, hist_dr2, linewidth=2, color='blue')
        axes[0,1].set_xlabel('dr2 (Mpc/h)')
        axes[0,1].set_ylabel('Density')
        axes[0,1].set_title('Position Displacement dr2')
        axes[0,1].axvline(0, color='red', linestyle='--', alpha=0.7)
        axes[0,1].grid(True, alpha=0.3)
        
        axes[0,2].plot(bin_centers, hist_drp, linewidth=2, color='blue')
        axes[0,2].set_xlabel('drp (Mpc/h)')
        axes[0,2].set_ylabel('Density')
        axes[0,2].set_title('Position Displacement drp (parallel)')
        axes[0,2].axvline(0, color='red', linestyle='--', alpha=0.7)
        axes[0,2].grid(True, alpha=0.3)
        
        print(f"Position displacements:")
        print(f"  dr1 - mean: {np.mean(dr1):.6f}, std: {np.std(dr1):.6f}")
        print(f"  dr2 - mean: {np.mean(dr2):.6f}, std: {np.std(dr2):.6f}")
        print(f"  drp - mean: {np.mean(drp):.6f}, std: {np.std(drp):.6f}")
    
    if shape_displacements is not None:
        # Shape displacement components
        shape_dr1 = shape_displacements['dr1']
        shape_dr2 = shape_displacements['dr2']
        shape_drp = shape_displacements['drp']
        
        hist_sdr1, _ = np.histogram(shape_dr1, bins=bins, density=True)
        hist_sdr2, _ = np.histogram(shape_dr2, bins=bins, density=True)
        hist_sdrp, _ = np.histogram(shape_drp, bins=bins, density=True)
        
        axes[1,0].plot(bin_centers, hist_sdr1, linewidth=2, color='green')
        axes[1,0].set_xlabel('dr1 (Mpc/h)')
        axes[1,0].set_ylabel('Density')
        axes[1,0].set_title('Shape Displacement dr1')
        axes[1,0].axvline(0, color='red', linestyle='--', alpha=0.7)
        axes[1,0].grid(True, alpha=0.3)
        
        axes[1,1].plot(bin_centers, hist_sdr2, linewidth=2, color='green')
        axes[1,1].set_xlabel('dr2 (Mpc/h)')
        axes[1,1].set_ylabel('Density')
        axes[1,1].set_title('Shape Displacement dr2')
        axes[1,1].axvline(0, color='red', linestyle='--', alpha=0.7)
        axes[1,1].grid(True, alpha=0.3)
        
        axes[1,2].plot(bin_centers, hist_sdrp, linewidth=2, color='green')
        axes[1,2].set_xlabel('drp (Mpc/h)')
        axes[1,2].set_ylabel('Density')
        axes[1,2].set_title('Shape Displacement drp (parallel)')
        axes[1,2].axvline(0, color='red', linestyle='--', alpha=0.7)
        axes[1,2].grid(True, alpha=0.3)
        
        print(f"Shape displacements:")
        print(f"  dr1 - mean: {np.mean(shape_dr1):.6f}, std: {np.std(shape_dr1):.6f}")
        print(f"  dr2 - mean: {np.mean(shape_dr2):.6f}, std: {np.std(shape_dr2):.6f}")
        print(f"  drp - mean: {np.mean(shape_drp):.6f}, std: {np.std(shape_drp):.6f}")
    
    plt.tight_layout()
    plt.savefig('/home/murray/intrinsic_alignments/validation/create_catalogues/displacement_distributions_lines.png', dpi=150, bbox_inches='tight')
    print("Saved: /home/murray/intrinsic_alignments/validation/create_catalogues/displacement_distributions_lines.png")

def plot_displacement_magnitude_analysis(displacements, shape_displacements):
    """Plot displacement magnitude analysis."""
    print("\n=== DISPLACEMENT MAGNITUDE ANALYSIS ===")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    if displacements is not None:
        # Calculate 2D and 3D displacement magnitudes
        dr_2d = np.sqrt(displacements['dr1']**2 + displacements['dr2']**2)
        dr_3d = np.sqrt(displacements['dr1']**2 + displacements['dr2']**2 + displacements['drp']**2)
        
        # Line plots for magnitudes
        bins_mag = np.linspace(0, 40, 80)
        hist_2d, _ = np.histogram(dr_2d, bins=bins_mag, density=True)
        hist_3d, _ = np.histogram(dr_3d, bins=bins_mag, density=True)
        bin_centers_mag = (bins_mag[:-1] + bins_mag[1:]) / 2
        
        axes[0,0].plot(bin_centers_mag, hist_2d, linewidth=2, label='2D magnitude', color='blue')
        axes[0,0].plot(bin_centers_mag, hist_3d, linewidth=2, label='3D magnitude', color='red')
        axes[0,0].set_xlabel('Displacement magnitude (Mpc/h)')
        axes[0,0].set_ylabel('Density')
        axes[0,0].set_title('Position Displacement Magnitudes')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # 2D displacement field visualization
        # Show displacement direction preferences
        angles = np.arctan2(displacements['dr2'], displacements['dr1'])
        hist_angles, angle_bins = np.histogram(angles, bins=36, density=True)
        angle_centers = (angle_bins[:-1] + angle_bins[1:]) / 2
        
        axes[0,1].plot(angle_centers * 180/np.pi, hist_angles, linewidth=2, color='purple')
        axes[0,1].set_xlabel('Displacement angle (degrees)')
        axes[0,1].set_ylabel('Density')
        axes[0,1].set_title('Position Displacement Direction Distribution')
        axes[0,1].grid(True, alpha=0.3)
        
        print(f"Position displacement magnitudes:")
        print(f"  2D magnitude - mean: {np.mean(dr_2d):.3f}, std: {np.std(dr_2d):.3f}")
        print(f"  3D magnitude - mean: {np.mean(dr_3d):.3f}, std: {np.std(dr_3d):.3f}")
    
    if shape_displacements is not None:
        # Shape displacement magnitudes
        sdr_2d = np.sqrt(shape_displacements['dr1']**2 + shape_displacements['dr2']**2)
        sdr_3d = np.sqrt(shape_displacements['dr1']**2 + shape_displacements['dr2']**2 + shape_displacements['drp']**2)
        
        hist_s2d, _ = np.histogram(sdr_2d, bins=bins_mag, density=True)
        hist_s3d, _ = np.histogram(sdr_3d, bins=bins_mag, density=True)
        
        axes[1,0].plot(bin_centers_mag, hist_s2d, linewidth=2, label='2D magnitude', color='green')
        axes[1,0].plot(bin_centers_mag, hist_s3d, linewidth=2, label='3D magnitude', color='orange')
        axes[1,0].set_xlabel('Displacement magnitude (Mpc/h)')
        axes[1,0].set_ylabel('Density')
        axes[1,0].set_title('Shape Displacement Magnitudes')
        axes[1,0].legend()
        axes[1,0].grid(True, alpha=0.3)
        
        # Shape displacement directions
        s_angles = np.arctan2(shape_displacements['dr2'], shape_displacements['dr1'])
        hist_s_angles, _ = np.histogram(s_angles, bins=36, density=True)
        
        axes[1,1].plot(angle_centers * 180/np.pi, hist_s_angles, linewidth=2, color='brown')
        axes[1,1].set_xlabel('Displacement angle (degrees)')
        axes[1,1].set_ylabel('Density')
        axes[1,1].set_title('Shape Displacement Direction Distribution')
        axes[1,1].grid(True, alpha=0.3)
        
        print(f"Shape displacement magnitudes:")
        print(f"  2D magnitude - mean: {np.mean(sdr_2d):.3f}, std: {np.std(sdr_2d):.3f}")
        print(f"  3D magnitude - mean: {np.mean(sdr_3d):.3f}, std: {np.std(sdr_3d):.3f}")
    
    plt.tight_layout()
    plt.savefig('/home/murray/intrinsic_alignments/validation/create_catalogues/displacement_magnitude_analysis.png', dpi=150, bbox_inches='tight')
    print("Saved: /home/murray/intrinsic_alignments/validation/create_catalogues/displacement_magnitude_analysis.png")

def main():
    print("=== IMPROVED LRG DISPLACEMENT CATALOGUE VALIDATION ===")
    
    # Create output directory
    os.makedirs('../figures/validation', exist_ok=True)
    
    # Load displacement catalogues
    displacements, shape_displacements = load_displacement_catalogues()
    
    if displacements is None:
        print("Error: Could not load displacement catalogues!")
        return
    
    # Generate improved visualizations
    plot_displacement_vectors_on_sky(displacements, shape_displacements)
    plot_physically_accurate_displacement_vectors(displacements, shape_displacements)
    plot_displacement_distributions(displacements, shape_displacements)
    plot_displacement_magnitude_analysis(displacements, shape_displacements)
    
    print("\n=== IMPROVED VALIDATION COMPLETE ===")
    print("Generated plots:")
    print("- ../figures/validation/displacement_vectors_sky_patch.png")
    print("- ../figures/validation/displacement_vectors_physical_scale.png")
    print("- ../figures/validation/displacement_distributions_lines.png") 
    print("- ../figures/validation/displacement_magnitude_analysis.png")

if __name__ == "__main__":
    main()