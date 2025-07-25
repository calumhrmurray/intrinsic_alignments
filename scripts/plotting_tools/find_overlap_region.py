#!/usr/bin/env python3
"""
Find regions where displacement and shape catalogues overlap well.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import configparser

def find_overlap_region():
    # Load config
    config = configparser.ConfigParser()
    config.read('/home/murray/intrinsic_alignments/scripts/config_files/lrg_observed_displacement_shear.ini')
    
    print("Loading displacement catalogue...")
    displacements = fits.open(config['general']['velocity_catalogue_folder'] + 
                             config['general']['position_tracer'] + '_displacements.fits')[1].data
    
    print("Loading shape displacements catalogue...")
    shape_displacements = fits.open(config['general']['velocity_catalogue_folder'] + 
                                   config['general']['shape_tracer'] + '_shape_displacements.fits')[1].data
    
    print(f"Displacement catalogue: {len(displacements):,} objects")
    print(f"Shape catalogue: {len(shape_displacements):,} objects")
    
    # Plot sky coverage separately
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot 1: Displacement catalogue
    axes[0].scatter(displacements['RA'], displacements['Dec'], s=0.1, alpha=0.5, c='blue')
    axes[0].set_xlabel('RA (deg)')
    axes[0].set_ylabel('Dec (deg)')
    axes[0].set_title(f'Displacement Catalogue (N={len(displacements):,})')
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Shape catalogue
    axes[1].scatter(shape_displacements['RA'], shape_displacements['Dec'], s=0.1, alpha=0.5, c='red')
    axes[1].set_xlabel('RA (deg)')
    axes[1].set_ylabel('Dec (deg)')
    axes[1].set_title(f'Shape Catalogue (N={len(shape_displacements):,})')
    axes[1].grid(True, alpha=0.3)
    
    # Plot 3: Overlay to see overlap
    axes[2].scatter(displacements['RA'], displacements['Dec'], s=0.1, alpha=0.3, c='blue', label='Displacements')
    axes[2].scatter(shape_displacements['RA'], shape_displacements['Dec'], s=0.1, alpha=0.3, c='red', label='Shapes')
    axes[2].set_xlabel('RA (deg)')
    axes[2].set_ylabel('Dec (deg)')
    axes[2].set_title('Overlap of Both Catalogues')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/home/murray/intrinsic_alignments/validation/correlation_functions/sky_coverage_analysis.png', dpi=150, bbox_inches='tight')
    print("Sky coverage plot saved: /home/murray/intrinsic_alignments/validation/correlation_functions/sky_coverage_analysis.png")
    plt.close()
    
    # Find regions with good overlap
    print("\n=== Sky Coverage Analysis ===")
    print(f"Displacement RA range: [{displacements['RA'].min():.1f}, {displacements['RA'].max():.1f}]")
    print(f"Displacement Dec range: [{displacements['Dec'].min():.1f}, {displacements['Dec'].max():.1f}]")
    print(f"Shape RA range: [{shape_displacements['RA'].min():.1f}, {shape_displacements['RA'].max():.1f}]")
    print(f"Shape Dec range: [{shape_displacements['Dec'].min():.1f}, {shape_displacements['Dec'].max():.1f}]")
    
    # Find overlap region
    overlap_ra_min = max(displacements['RA'].min(), shape_displacements['RA'].min())
    overlap_ra_max = min(displacements['RA'].max(), shape_displacements['RA'].max())
    overlap_dec_min = max(displacements['Dec'].min(), shape_displacements['Dec'].min())
    overlap_dec_max = min(displacements['Dec'].max(), shape_displacements['Dec'].max())
    
    print(f"\nOverlap region: RA [{overlap_ra_min:.1f}, {overlap_ra_max:.1f}], Dec [{overlap_dec_min:.1f}, {overlap_dec_max:.1f}]")
    
    # Test different subregions and count objects
    test_regions = [
        # Small dense regions for testing
        (140, 160, 20, 40),   # 20×20 deg region
        (160, 180, 30, 50),   # 20×20 deg region  
        (200, 220, 0, 20),    # 20×20 deg region
        (120, 150, 10, 40),   # 30×30 deg region
        (180, 210, -10, 20),  # 30×30 deg region
    ]
    
    print("\n=== Testing subregions for good overlap ===")
    for i, (ra_min, ra_max, dec_min, dec_max) in enumerate(test_regions):
        disp_in_region = np.sum((displacements['RA'] >= ra_min) & 
                               (displacements['RA'] <= ra_max) &
                               (displacements['Dec'] >= dec_min) & 
                               (displacements['Dec'] <= dec_max))
        
        shapes_in_region = np.sum((shape_displacements['RA'] >= ra_min) & 
                                 (shape_displacements['RA'] <= ra_max) &
                                 (shape_displacements['Dec'] >= dec_min) & 
                                 (shape_displacements['Dec'] <= dec_max))
        
        area = (ra_max - ra_min) * (dec_max - dec_min)
        disp_density = disp_in_region / area if area > 0 else 0
        shape_density = shapes_in_region / area if area > 0 else 0
        
        print(f"Region {i+1}: RA [{ra_min}, {ra_max}] × Dec [{dec_min}, {dec_max}]")
        print(f"  Displacements: {disp_in_region:,} objects (density: {disp_density:.1f}/deg²)")
        print(f"  Shapes: {shapes_in_region:,} objects (density: {shape_density:.1f}/deg²)")
        print(f"  Area: {area:.0f} deg²")
        
        if disp_in_region > 1000 and shapes_in_region > 100:
            print(f"  *** GOOD CANDIDATE: Both catalogues have good coverage ***")
        print()

if __name__ == "__main__":
    find_overlap_region()