import numpy as np
import sys

try:
    vaspfile_ref = sys.argv[1]
    outfile = sys.argv[2]
except IndexError:
    print("Please enther the input and output file name!")
    print("****************************")
    print("Example usage: python convert.py inputfile outputfile")
    print("****************************")
    exit()

# Read data from vasp file
def readxyz(poscar):
    data = open(poscar,'r')
    lines = data.readlines()

    xsize = float(lines[2].split()[0])
    ysize = float(lines[3].split()[1])
    zsize = float(lines[4].split()[2])

    natoms = int(lines[6])
    boundaries = ([xsize, ysize, zsize])
    print("boundary:",boundaries)
    position_all = []
    id_all = []
    for i, line in enumerate(lines[8:8+natoms]):
        line = line.split()
        position = [float(line[0]), float(line[1]), float(line[2])]
        position_all.append(position)
        id_all.append(i + 1)
            
    return natoms, id_all, position_all, boundaries

# Convert the region to a PBC cell.
def convert_pbc(atominfo, xlatspacing, ylatspacing):
    origin_coords = []
    id_all = atominfo[1]
    print("id_all",id_all)
    position_all = atominfo[2]
    boundaries = atominfo[3]
    for i in range(len(id_all)):
        new_coords = (id_all[i], position_all[i][0], position_all[i][1], position_all[i][2])
        origin_coords.append(new_coords)
    # measure the distance along X, Y and Z
    origin_coords = np.array(origin_coords)
    rbt_coords = []
    # move the lower left atom to origin (0,0)
    for everyele in origin_coords:
        rbt_x = everyele[1] - min(origin_coords[:, 1])
        #print("rbt_x:",rbt_x)
        rbt_y = everyele[2] - min(origin_coords[:, 2])
        #print("rbt_y:",rbt_y)
        rbt_z = everyele[3]
        #print("rbt_z:",rbt_z)
        rbt = ([everyele[0], rbt_x, rbt_y, rbt_z])
        #print("Before conversion:", origin_coords)
        rbt_coords.append(rbt)
        #print("After rbt conversion:", rbt_coords)

    rbt_coords = np.array(rbt_coords)
    #print("After rbt conversion:", rbt_coords)
    # copy image and rotate
    # create an symmetry image, rotation symmetry, counterclockwise by pi.
    # can be represented by R=[-1,0;0,-1]
    # [x,y] --> R.*[x,y]=[-x,-y]
    # displace the image rigidly: move up by the distance of (Ymax-Ymin-one_lattice)
    
    pbc_coords_img = []
    image_coords = []
    for coords in rbt_coords:
        atomid = int(id_all[0])
        rot_x = -coords[1] 
        rot_y = -coords[2] 
        rot_z = coords[3]
        pbc_coords_img.append([len(rbt_coords) + atomid,rot_x, rot_y, rot_z])

    pbc_coords_img = np.array(pbc_coords_img)
    print("After pbc1 conversion:", pbc_coords_img)
        #print("pbc_image_coords:",pbc_coords_img)
    # the left bottom corner atom of rotated image is
    # the right upper of origin before rotation
    # find the coordinate of that atom
    for everyele in pbc_coords_img:
        rbt_x = everyele[1] + (max(pbc_coords_img[:,1]) - min(pbc_coords_img[:,1]))
        rbt_y = everyele[2] + 2 * (max(pbc_coords_img[:,2]) - min(pbc_coords_img[:,2])) + ylatspacing
        rbt_z = everyele[3]
        rbt = ([everyele[0], rbt_x, rbt_y, rbt_z])
        image_coords.append(rbt)
    image_coords = np.array(image_coords)
    print("After image conversion:", image_coords)
    # find the rigid displacement of the rotated image
    # in order to align with the original one
    # displacement along X
    all_coords = np.vstack((rbt_coords, image_coords))
    print(all_coords.shape)
    all_coords = np.array(all_coords)
    print("all coords:", all_coords)
    pbc_coords = all_coords[:,1:4]
    # add Gaussin noisy (0,0.05)
    noise = np.random.normal(0,0.05,[len(pbc_coords),3])
    pbc_coords = (pbc_coords + noise).tolist()
    pbc_coords = np.array(pbc_coords)

    # Tune the size of boundaries
    xmax = np.max(pbc_coords[:,1])
    xlo = 0
    xhi = xmax + xlatspacing
    print("xhi:", xhi)

    ymax = np.max(pbc_coords[:,2])
    ylo = 0
    yhi = ymax + ylatspacing
    print("yhi:", yhi)

    zlo = 0
    zhi = 3.2705400
    boundaries = ([xlo,xhi,ylo,yhi,zlo,zhi])
    print("boundary:",boundaries)
    
    return boundaries, pbc_coords


# Write structure information to vasp input file
def write_vasp(final_data, out_filename):
    
    # determine the 3d matrix
    xlength = final_data[0][1] - final_data[0][0]
    ylength = final_data[0][3] - final_data[0][2]
    zlength = final_data[0][5] - final_data[0][4]
    # write vasp input file
    poscar = open(out_filename,'w')
    poscar.write('#periodic structure of crack\n')
    poscar.write("1.000000\n")
    poscar.write(f'     {xlength}        0.00000000        0.00000000\n')
    poscar.write(f'     0.00000000       {ylength}        0.00000000\n')
    poscar.write(f'     0.00000000       0.00000000        {zlength}\n')
    poscar.write(" Nb\n")
    poscar.write("    82\n")
    poscar.write("Cartesian\n")
    for coords in final_data[1]:
        poscar.write('%12.8f %12.8f %12.8f\n' % (coords[0], coords[1], coords[2]))
    print("****************************")
    print("Write vasp input file successful!")
    print("vasp input file name: %s" % outfile)
    print("****************************")
    
    return

###################################################################
# Execute the code 
###################################################################
# define some variables related to the crystallography
latconst = 3.27
xlatspacing = latconst * 0.5 * np.sqrt(2)
ylatspacing = latconst * 0.5 * np.sqrt(1)
#----------------------
# Read current configuration
cur_config = readxyz(vaspfile_ref)
# Convert crack tip to a PBC cell
pbc_final = convert_pbc(cur_config, xlatspacing, ylatspacing)
# Write input file for vasp
write_vasp(pbc_final, outfile)