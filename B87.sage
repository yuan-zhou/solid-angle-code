# FILEPATH: /Users/allisonfitisone/solid-angle-code/B87.sage
import logging
import csv
deg = 50000000000000000
eps = 1e-3
decompose_to_tridiag = False
base_ring = RR
verbose = False
output_file_path='b.csv'
# Define your column titles
column_titles = ["Vertex", "Total_Num_Cones", "Epsilon", "Solid_Angle", "CPU_Time"]

# Open the CSV file in write mode and write the column titles
with open(output_file_path, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(column_titles)
    
    vect = (3/5, 2/3, 1/5, 4/5, 1/3, 2/5)
    rays = [[-1, 0, -2, 0, 0, 0], [-1, 0, 0, 0, 0, -1], [0, -1, 0, 0, -1, 0], [0, 0, -1, -1, 0, 0], [0, 0, -1, 0, 0, -2], [0, 0, -1, 0, 0, 0], [0, 0, 0, 0, -1, 0]]
    triangulation = [[0, 1, 2, 3, 5, 6], [1, 2, 3, 4, 5, 6]]
    T = [matrix([rays[i] for i in simplex]) for simplex in triangulation]
    Start_Time = time.process_time()
    sa_vect = solid_angle_measure(T, deg=deg, eps=eps, decompose_to_tridiag=decompose_to_tridiag, base_ring=base_ring, verbose=verbose)
    Execution_Time = time.process_time() - Start_Time
    row_to_append = [vect, total_num_cones(T), eps, sa_vect, Execution_Time]
    writer.writerow(row_to_append)
    print('{} has solid angle {}'.format(vect, sa_vect))
    print(' ')

    vect = (7/9, 2/3, 5/9, 4/9, 1/3, 2/9)
    rays = [[-1, 0, 0, 0, 0, -1], [0, -1, 0, 0, -1, 0], [0, 0, -1, -1, 0, 0], [0, 0, -1, 0, 0, -2], [0, 0, 0, -1, -1, -1], [0, 0, 0, 0, -1, -3], [0, 0, 0, 0, -1, 0]]
    triangulation = [[0, 1, 2, 3, 5, 6], [0, 1, 2, 4, 5, 6]]
    T = [matrix([rays[i] for i in simplex]) for simplex in triangulation]
    Start_Time = time.process_time()
    sa_vect = solid_angle_measure(T, deg=deg, eps=eps, decompose_to_tridiag=decompose_to_tridiag, base_ring=base_ring, verbose=verbose)
    Execution_Time = time.process_time() - Start_Time
    row_to_append = [vect, total_num_cones(T), eps, sa_vect, Execution_Time]
    writer.writerow(row_to_append)
    print('{} has solid angle {}'.format(vect, sa_vect))
    print(' ')

    vect = (1/7, 2/7, 3/7, 4/7, 5/7, 6/7)
    rays = [[-3, 0, 0, -1, 0, 0], [-2, 0, 0, 0, -1, 0], [-1, -3, 0, 0, 0, 0], [-1, -1, 0, -1, 0, 0], [-1, 0, -2, 0, 0, 0], [-1, 0, 0, 0, 0, -1], [-1, 0, 0, 0, 0, 0], [0, -2, -1, 0, 0, 0], [0, -1, 0, 0, -1, 0], [0, 0, -1, -1, 0, 0]]
    triangulation = [[0, 1, 5, 6, 8, 9], [0, 3, 5, 6, 8, 9], [1, 4, 5, 6, 8, 9], [2, 3, 5, 6, 8, 9], [2, 5, 6, 7, 8, 9], [4, 5, 6, 7, 8, 9]]
    T = [matrix([rays[i] for i in simplex]) for simplex in triangulation]
    Start_Time = time.process_time()
    sa_vect = solid_angle_measure(T, deg=deg, eps=eps, decompose_to_tridiag=decompose_to_tridiag, base_ring=base_ring, verbose=verbose)
    Execution_Time = time.process_time() - Start_Time
    row_to_append = [vect, total_num_cones(T), eps, sa_vect, Execution_Time]
    writer.writerow(row_to_append)
    print('{} has solid angle {}'.format(vect, sa_vect))
    print(' ')

    vect = (1/3, 2/3, 1, 0, 1/3, 2/3)
    rays = [[-3, 0, 0, -1, 0, 0], [-2, 0, 0, 0, -1, 0], [-1, -1, 0, -1, 0, 0], [-1, 0, 0, 0, 0, -1], [0, -1, 0, 0, -1, 0], [0, 0, -1, -1, 0, 0], [0, 0, 0, -1, -1, -1], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0]]
    triangulation = [[0, 1, 3, 4, 5, 8], [0, 2, 3, 4, 5, 7], [0, 3, 4, 5, 7, 8], [3, 4, 5, 6, 7, 8]]
    T = [matrix([rays[i] for i in simplex]) for simplex in triangulation]
    Start_Time = time.process_time()
    sa_vect = solid_angle_measure(T, deg=deg, eps=eps, decompose_to_tridiag=decompose_to_tridiag, base_ring=base_ring, verbose=verbose)
    Execution_Time = time.process_time() - Start_Time
    row_to_append = [vect, total_num_cones(T), eps, sa_vect, Execution_Time]
    writer.writerow(row_to_append)
    print('{} has solid angle {}'.format(vect, sa_vect))
    print(' ')

    vect = (1/3, 2/3, 1/3, 2/3, 1/3, 2/3)
    rays = [[-2, 0, 0, 0, -1, 0], [-1, 0, -2, 0, 0, 0], [-1, 0, 0, 0, 0, -1], [0, -1, 0, 0, -1, 0], [0, 0, -1, -1, 0, 0], [0, 0, 0, 0, -1, 0]]
    triangulation = [[0, 1, 2, 3, 4, 5]]
    T = [matrix([rays[i] for i in simplex]) for simplex in triangulation]
    Start_Time = time.process_time()
    sa_vect = solid_angle_measure(T, deg=deg, eps=eps, decompose_to_tridiag=decompose_to_tridiag, base_ring=base_ring, verbose=verbose)
    Execution_Time = time.process_time() - Start_Time
    row_to_append = [vect, total_num_cones(T), eps, sa_vect, Execution_Time]
    writer.writerow(row_to_append)
    print('{} has solid angle {}'.format(vect, sa_vect))
    print(' ')

    vect = (1, 0, 1, 0, 1, 0)
    rays = [[-1, -3, 0, 0, 0, 0], [-1, -1, 0, -1, 0, 0], [-1, 0, 0, 0, 0, -1], [0, -2, -1, 0, 0, 0], [0, -1, 0, 0, -1, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, -1, 0, 0], [0, 0, -1, 0, 0, -2], [0, 0, 0, -1, -1, -1], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, -3], [0, 0, 0, 0, 0, -1]]
    triangulation = [[0, 1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 6], [1, 2, 4, 5, 6, 7], [1, 2, 4, 5, 7, 9], [1, 2, 4, 6, 7, 9], [1, 4, 5, 6, 7, 9], [2, 3, 4, 5, 6, 7], [2, 4, 5, 7, 9, 10], [2, 4, 6, 7, 9, 10], [2, 4, 6, 8, 9, 10], [2, 5, 7, 9, 10, 11]]
    T = [matrix([rays[i] for i in simplex]) for simplex in triangulation]
    Start_Time = time.process_time()
    sa_vect = solid_angle_measure(T, deg=deg, eps=eps, decompose_to_tridiag=decompose_to_tridiag, base_ring=base_ring, verbose=verbose)
    Execution_Time = time.process_time() - Start_Time
    row_to_append = [vect, total_num_cones(T), eps, sa_vect, Execution_Time]
    writer.writerow(row_to_append)
    print('{} has solid angle {}'.format(vect, sa_vect))
    print(' ')
    
    vect = (3/5, 2/5, 1/5, 4/5, 3/5, 2/5)
    rays = [[-1, 0, -2, 0, 0, 0], [-1, 0, 0, 0, 0, -1], [0, -2, -1, 0, 0, 0], [0, -1, 0, 0, -1, 0], [0, 0, -1, -1, 0, 0], [0, 0, -1, 0, 0, -2], [0, 0, -1, 0, 0, 0]]
    triangulation = [[0, 1, 2, 3, 4, 6], [1, 2, 3, 4, 5, 6]]
    T = [matrix([rays[i] for i in simplex]) for simplex in triangulation]
    Start_Time = time.process_time()
    sa_vect = solid_angle_measure(T, deg=deg, eps=eps, decompose_to_tridiag=decompose_to_tridiag, base_ring=base_ring, verbose=verbose)
    Execution_Time = time.process_time() - Start_Time
    row_to_append = [vect, total_num_cones(T), eps, sa_vect, Execution_Time]
    writer.writerow(row_to_append)
    print('{} has solid angle {}'.format(vect, sa_vect))
