import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import freud
from scipy.spatial import ConvexHull, QhullError


def compute_tetrahedron_centroid(v0, v1, v2, v3):
    return 0.25 * (v0 + v1 + v2 + v3)


def compute_tetrahedron_signed_volume(v0, v1, v2, v3):
    return np.dot(v1 - v0, np.cross(v2 - v0, v3 - v0)) / 6.0


def compute_polyhedron_centroid_from_vertices(cell_vertices, interior_point):
    
    hull = ConvexHull(cell_vertices)

    total_volume = 0.0
    weighted_centroid = np.zeros(3)

    for simplex in hull.simplices:
        v1, v2, v3 = cell_vertices[simplex]

        signed_volume = compute_tetrahedron_signed_volume(interior_point, v1, v2, v3)

        volume = abs(signed_volume)

        if (volume <= 1.0e-15):
            continue

        centroid = compute_tetrahedron_centroid(interior_point, v1, v2, v3)

        weighted_centroid += volume * centroid
        total_volume += volume

    if (total_volume <= 1.0e-15):
        raise ValueError("Degenerate Voronoi cell with zero volume")

    return weighted_centroid / total_volume, total_volume


def lloyd_relaxation_3d(initial_points, box, w=1.0, iterations=10):
    points = box.wrap(np.asarray(initial_points, dtype=np.float64).copy())

    for _ in range(iterations):
        voro = freud.locality.Voronoi()
        voro.compute((box, points))

        cell_vertices_all = voro.polytopes
        new_points = points.copy()

        for i in range(len(points)):
            site = points[i, :]
            rel_vertices = box.wrap(cell_vertices_all[i] - site)
            cell_vertices = site + rel_vertices

            try:
                centroid, cell_volume = compute_polyhedron_centroid_from_vertices(cell_vertices, site)

            except (QhullError, ValueError):
                print(f"Warning: could not compute centroid for cell {i}")
                continue

            # periodic displacement from site to centroid.
            disp = box.wrap(centroid - site)

            new_points[i, :] = site + w * disp

        points = box.wrap(new_points)

    return points


if (__name__ == '__main__'): 
    print('running 3D...')

    # setup 
    phi = 0.05

    D = 0.1
    L = 10 * D

    Lx = L
    Ly = 2 * L
    Lz = L

    output_dir = '../runs/moving_particle_array'
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    N_sphere = int( 6*phi*Lx*Ly*Lz / (np.pi*D**3) )
    print(f'volume fraction phi: {phi}, number of spheres: {N_sphere}')
    print(f'actual phi value: {N_sphere*4/3*np.pi*(D/2)**3/(Lx*Ly*Lz)}')

    x_i = Lx/2 * np.random.uniform(-1, 1, N_sphere)
    y_i = Ly/2 * np.random.uniform(-1, 1, N_sphere)
    z_i = Lz/2 * np.random.uniform(-1, 1, N_sphere)

    initial_points = np.stack((x_i, y_i, z_i), axis=1)
    box = freud.box.Box(Lx, Ly, Lz)
    
    relaxed_points = lloyd_relaxation_3d(initial_points, box, iterations=60)
    print(np.shape(relaxed_points))

    np.savetxt(output_dir+'/sphere_array_locations.txt', relaxed_points)

    # check no spheres are overlapping, including periodic images
    min_dist = L
    overlap_tol = 1.05 * D

    for i in range(N_sphere):
        for j in range(i + 1, N_sphere):

            # displacement
            dr = relaxed_points[i, :] - relaxed_points[j, :]

            # minimum-image periodic displacement
            dr_periodic = box.wrap(dr)

            dist = np.linalg.norm(dr_periodic)
            min_dist = min(min_dist, dist)

            if dist <= overlap_tol:
                print(f'spheres overlapping, periodic dist={dist}, spheres #: {i}, {j}')
                print(f'locations: ({relaxed_points[i, :]}), ({relaxed_points[j, :]})')
                print(f'raw displacement:      {dr}')
                print(f'periodic displacement: {dr_periodic}')
                print(f'periodic image shift:  {dr_periodic - dr}')

    print(f'closest neighbors physical distance: {min_dist}, distance in particle diameters: {min_dist/D}, separation in diameters: {min_dist/D - 1}')

    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.scatter(initial_points[:, 0], initial_points[:, 1], initial_points[:, 2], color='blue', s=10)
    ax1.set_title('initial points')
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.scatter(relaxed_points[:, 0], relaxed_points[:, 1], relaxed_points[:, 2], color='red', s=10)
    ax2.set_title('relaxed points')
    plt.show()
    plt.close()
