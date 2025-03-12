import numpy as np
from collections import defaultdict
import csv


import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MaxNLocator


import vtk


def read_vtk(file_path):
    # Read the VTK file
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(file_path)
    reader.Update()

    # Get the output data
    polydata = reader.GetOutput()

    # Extract points
    points_vtk = polydata.GetPoints()
    points = []
    for i in range(points_vtk.GetNumberOfPoints()):
        points.append(points_vtk.GetPoint(i))

    # Extract radii (assuming it's stored as point data, e.g., scalar field)
    radii = []
    radius_data = polydata.GetPointData().GetScalars("radius")
    if radius_data:
        for i in range(radius_data.GetNumberOfTuples()):
            radii.append(radius_data.GetValue(i))
    else:
        raise ValueError("No 'radius' scalar field found in the VTK file.")

    return points, radii


# Function to calculate the volume of a sphere
def sphere_volume(radius):
    return (4 / 3) * np.pi * radius**3


# Function to calculate packing density with height exclusion
def calculate_packing_density(points, radii, domain_volume, max_height):
    total_particle_volume = 0
    for point, radius in zip(points, radii):
        # Calculate the height based on the y-coordinate and the radius
        height = point[2] + radius
        # Exclude the particles above the max height
        if height > max_height:
            continue

        total_particle_volume += sphere_volume(radius)

    packing_density = total_particle_volume / domain_volume
    return packing_density


def sphere_segment_volume(radius, z1, z2, sphere_center_z):
    if z1 > z2:
        z1, z2 = z2, z1  # Ensure z1 <= z2

    # Limits of the spherical segment
    z_lower = max(sphere_center_z - radius, z1)
    z_upper = min(sphere_center_z + radius, z2)

    # Check for no overlap
    if z_upper <= z_lower:
        return 0

    # Calculate the volume of the spherical segment
    R = radius
    zc = sphere_center_z

    term1 = R**2 * (z_upper - z_lower)
    term2 = (1 / 3) * ((z_upper - zc) ** 3 - (z_lower - zc) ** 3)
    volume = np.pi * (term1 - term2)
    return volume


# Function to calculate packing densities for N regions, from min_z to max_z
def calculate_packing_density_regions(points, radii, domain_volume, num_regions, max_z):
    region_volumes = np.zeros(num_regions)

    # Initialize particle counts grouped by radius for each region
    particle_counts = {i: defaultdict(float) for i in range(num_regions)}

    min_z = 0

    total_height = max_z - min_z
    region_height = total_height / num_regions
    print(f"Region height: {region_height}")

    for point, radius in zip(points, radii):
        # Get the top and bottom z-coordinates of the particle
        z_bottom = point[2] - radius
        z_top = point[2] + radius

        # Determine which regions the particle overlaps
        start_region = int(np.floor((z_bottom - min_z) / region_height))
        end_region = int(np.floor((z_top - min_z) / region_height))

        # Ensure start_region and end_region are within bounds
        start_region = max(0, start_region)
        end_region = min(num_regions - 1, end_region)

        # Loop over all regions that the particle overlaps with
        for region in range(start_region, end_region + 1):
            region_bottom = min_z + region * region_height
            region_top = min_z + (region + 1) * region_height

            # Calculate the particle's height contribution within this region
            region_contribution_bottom = max(region_bottom, z_bottom)
            region_contribution_top = min(region_top, z_top)

            if region_contribution_top > region_contribution_bottom:
                # Calculate the volume of the sphere segment in this region
                volume_contribution = sphere_segment_volume(
                    radius,
                    region_contribution_bottom,
                    region_contribution_top,
                    point[2],
                )

                # Add to the corresponding region's volume
                region_volumes[region] += volume_contribution

                # Calculate the proportional height contribution
                height_contribution = (
                    region_contribution_top - region_contribution_bottom
                ) / (z_top - z_bottom)

                # Add proportional contribution to particle count for this radius
                particle_counts[region][radius] += height_contribution

    # Calculate packing density for each region
    region_packing_densities = region_volumes / (domain_volume / num_regions)

    return region_packing_densities, dict(particle_counts)


# Function to calculate the domain volume
def calculate_domain_volume(height):
    # Box dimensions as per the provided vertices
    width = 0.075  # x-direction
    depth = 0.075  # z-direction

    # Volume of the bounding box
    return width * height * depth


"""
all_files = [
    "final_results-ILs/1-1/particles_330000.vtk",
    "final_results-ILs/1-2/particles_330000.vtk",
    "final_results-ILs/1-3/particles_330000.vtk",
    "final_results-ILs/1-4/particles_330000.vtk",
    "final_results-ILs/2-1/particles_320000.vtk",
    "final_results-ILs/2-2/particles_320000.vtk",
    "final_results-ILs/2-3/particles_320000.vtk",
    "final_results-ILs/2-4/particles_320000.vtk",
    "final_results-ILs/3-1/particles_420000.vtk",
    "final_results-ILs/3-2/particles_340000.vtk",
    "final_results-ILs/3-3/particles_320000.vtk",
]


all_files = [
    "final_results-ILMs/1-1/particles_330000.vtk",
    "final_results-ILMs/1-2/particles_380000.vtk",
    "final_results-ILMs/1-3/particles_480000.vtk",
    "final_results-ILMs/2-1/particles_330000.vtk",
    "final_results-ILMs/2-2/particles_380000.vtk",
    "final_results-ILMs/2-3/particles_480000.vtk",
    "final_results-ILMs/3-1/particles_330000.vtk",
    "final_results-ILMs/3-2/particles_380000.vtk",
    "final_results-ILMs/3-3/particles_480000.vtk",
]
"""

all_files = [
    "ILs/Simulation_4-1/post/particles_450000.vtk",
    "ILs/Simulation_4-2/post/particles_450000.vtk",
    "ILs/Simulation_4-3/post/particles_330000.vtk",
    "ILs/Simulation_4-4/post/particles_330000.vtk",
]


# Main Execution

# Ensure the results directory exists
results_dir = "results_ILs"
os.makedirs(results_dir, exist_ok=True)

for vtk_file in all_files:
    points, radii = read_vtk(vtk_file)

    # Create a dictionary to count the occurrences of each radius
    radius_count = defaultdict(int)

    # Count how many particles have the same radius
    for radius in radii:
        radius_count[radius] += 1

    # Convert defaultdict to a regular dict (optional, for display purposes)
    radius_count = dict(radius_count)

    # Print the dictionary of radii and their respective counts
    print("Radius counts:")
    for radius, count in radius_count.items():
        print(f"Radius: {radius}, Count: {count}")

    # Calculate the max height for height exclusion which is the highest z-coordinate
    max_height = max(point[2] - radius for point, radius in zip(points, radii))

    print(f"Max Height: {max_height * 542.8} nm")

    # Overrule the max height if needed
    # max_height = 0.075

    # Calculate domain volume
    domain_volume = calculate_domain_volume(max_height)

    # Calculate packing density with height exclusion
    packing_density = calculate_packing_density(
        points, radii, domain_volume, max_height
    )
    print(f"Packing Density: {packing_density}")

    # Number of regions to divide the domain
    num_regions = 40

    # Calculate packing densities for each region
    region_packing_densities, particle_counts = calculate_packing_density_regions(
        points, radii, domain_volume, num_regions, max_height
    )

    # Generate a base name for plots using the file name
    simulation_name = os.path.basename(os.path.dirname(os.path.dirname(vtk_file)))

    # Define distances (in nm) along the domain for each region
    domain_length_nm = max_height * 542.8
    distances = np.linspace(0, domain_length_nm, num_regions)

    # Plotting the packing densities
    plt.figure(figsize=(8, 6))
    plt.plot(distances, region_packing_densities, marker="o")

    # Add labels and title
    plt.xlabel("Distance (nm)")
    plt.ylabel("Packing Density")
    plt.title(f"Packing Density by Distance for {simulation_name}")

    # Limit x-axis to show a maximum of 6 ticks
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=False, nbins=10))

    # Save the plot
    packing_density_plot_path = os.path.join(
        results_dir, f"{simulation_name}_packing_density.png"
    )
    plt.tight_layout()
    plt.savefig(packing_density_plot_path)
    plt.close()

    # Plot the particle counts for each region
    plt.figure(figsize=(8, 6))

    # Loop through each radius and plot the counts in each region
    for radius in radius_count.keys():
        counts_per_region = [
            particle_counts[region].get(radius, 0) for region in range(num_regions)
        ]
        plt.plot(
            distances, counts_per_region, marker="o", label=f"Radius: {radius:.8f}"
        )

    # Add labels and title
    plt.xlabel("Distance (nm)")
    plt.ylabel("Particle Count Contribution")
    plt.title(f"Particle Count Contribution by Radius for {simulation_name}")

    # Limit x-axis to show a maximum of 6 ticks
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=False, nbins=10))
    plt.legend(title="Radius", loc="upper right")

    # Save the plot
    particle_count_plot_path = os.path.join(
        results_dir, f"{simulation_name}_particle_counts.png"
    )
    plt.tight_layout()
    plt.savefig(particle_count_plot_path)
    plt.close()

    # ===== Save CSV Files =====

    # Extract a clean simulation name from the VTK file path
    simulation_name = os.path.basename(os.path.dirname(os.path.dirname(vtk_file)))
    particle_name = os.path.basename(os.path.dirname(vtk_file))
    formatted_simulation_name = f"{simulation_name}_{particle_name}"

    # Save CSV file for packing density vs. distance
    packing_density_csv_path = os.path.join(
        results_dir, f"{formatted_simulation_name}_packing_density.csv"
    )
    with open(packing_density_csv_path, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Distance (nm)", "Packing Density"])
        for d, pd_val in zip(distances, region_packing_densities):
            csv_writer.writerow([d, pd_val])
    print(f"Saved CSV: {packing_density_csv_path}")

    # Save CSV file for particle counts per region (with a column for each radius)
    particle_counts_csv_path = os.path.join(
        results_dir, f"{formatted_simulation_name}_particle_counts.csv"
    )
    # Sort the radii keys to ensure consistent column ordering
    sorted_radii = sorted(radius_count.keys())
    with open(particle_counts_csv_path, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        # Header: Distance column followed by one column per unique radius
        header = ["Distance (nm)"] + [f"Radius {r:.8f}" for r in sorted_radii]
        csv_writer.writerow(header)
        # Write one row per region
        for region in range(num_regions):
            row = [distances[region]]
            for r in sorted_radii:
                row.append(particle_counts[region].get(r, 0))
            csv_writer.writerow(row)
    print(f"Saved CSV: {particle_counts_csv_path}")
