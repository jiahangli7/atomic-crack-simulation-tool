import os

def is_atom_to_fix(x, y, z):

    if 15 < x < 16 or 15 < x < 16 or 36 < x < 37 or 36 < y < 37:
        return True
    return False

for i in range(0, 21):
    directory_path = f"./simulation/{i}"

    if not os.path.exists(directory_path):
        continue

    stru_file_path = os.path.join(directory_path, "POSCAR")  


    with open(stru_file_path, 'r') as file:
        lines = file.readlines()


    cartesian_line_index = lines.index("Cartesian\n")
    lines.insert(cartesian_line_index, "Selective Dynamics\n")


    coordinates_start_line = cartesian_line_index + 1
    coordinates_end_line = coordinates_start_line + 98  
    coordinates_lines = lines[coordinates_start_line:coordinates_end_line]

 
    modified_coordinates = []
    atom_indices_to_fix = [] 

    for index, line in enumerate(coordinates_lines):
        x, y, z, q, w, e = line.split()

    
        if is_atom_to_fix(float(x), float(y), float(z)):
            modified_coordinates.append(f"{x} {y} {z} F F F\n")
            atom_indices_to_fix.append(index)
        else:
            modified_coordinates.append(f"{x} {y} {z} T T T\n")


    lines[coordinates_start_line:coordinates_end_line] = modified_coordinates

 
    with open(stru_file_path, 'w') as file:
        file.writelines(lines)


    print(f"文件夹 {i} 中要固定的原子序号列表：{atom_indices_to_fix}")