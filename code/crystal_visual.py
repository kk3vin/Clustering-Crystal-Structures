from CifFile import ReadCif 
import numpy as np
import re
import math
import csv 
import os 
import plotly.graph_objects as go
import numpy as np
import plotly.express as px
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.neighbors import NearestNeighbors
import clustering

class motif_point:
    def __init__(self, coordinate, labels, types):
        self.coordinate = coordinate
        self.labels = labels
        self.types = types

    def __eq__(self,other):
        if self.coordinate == other.coordinate and self.labels == other.labels and self.types == other.types:
            return True

class UnitCellAxis:
    def __init__(self, vector, motif_point):
        self.vector = vector
        self.motif_point = motif_point

class list_cell:
    def __init__(self, cells, center_index):
        self.cells = cells
        self.center_index = center_index




def getUnitCellAxis(ciffile) :

    namespace = (ciffile.keys())[1]
    
    x = ciffile[namespace]['_cell_length_a']
    y = ciffile[namespace]['_cell_length_b']
    z = ciffile[namespace]['_cell_length_c']
    alpha = ciffile[namespace]['_cell_angle_alpha']
    beta = ciffile[namespace]['_cell_angle_beta']
    gamma = ciffile[namespace]['_cell_angle_gamma']
    m_xs = ciffile[namespace]['_atom_site_fract_x']
    m_ys = ciffile[namespace]['_atom_site_fract_y']
    m_zs = ciffile[namespace]['_atom_site_fract_z']
    m_label = ciffile[namespace]['_atom_site_label']
    m_type = ciffile[namespace]['_atom_site_type_symbol']
    motif_points = [None]*len(m_xs)
    m_x = [None]*len(m_xs)
    m_y = [None]*len(m_ys)
    m_z = [None]*len(m_zs)
    m_labels = [None]*len(m_label)
    m_types = [None]*len(m_type)
    for i in range(len(m_xs)):
        m_xs[i] = re.sub("[\(\[].*?[\)\]]", "",m_xs[i])
        m_ys[i] = re.sub("[\(\[].*?[\)\]]", "",m_ys[i])
        m_zs[i] = re.sub("[\(\[].*?[\)\]]", "",m_zs[i])
        m_label[i] = re.sub("[\(\[].*?[\)\]]", "",m_label[i])
        m_type[i] = re.sub("[\(\[].*?[\)\]]", "",m_type[i])
        m_x[i] = float(m_xs[i])
        m_y[i] = float(m_ys[i])
        m_z[i] = float(m_zs[i])
        m_labels[i] = str(m_label[i])
        m_types[i] = str(m_type[i])
        motif_points[i] = motif_point(np.array([m_x[i],m_y[i],m_z[i]]),m_labels[i],m_types[i])

    alpha_d = float(alpha)
    beta_d = float(beta)
    gamma_d = float(gamma)

    alpha = alpha_d*math.pi/180
    beta = beta_d*math.pi/180
    gamma = gamma_d*math.pi/180

    x = float(re.sub("[\(\[].*?[\)\]]", "", x))
    a = np.linalg.norm(x)
    a_v = [a,0,0]
    y = float(re.sub("[\(\[].*?[\)\]]", "", y))
    b = np.linalg.norm(y)
    b_v = [b*math.cos(gamma), b*math.sin(gamma), 0]
    z = float(re.sub("[\(\[].*?[\)\]]", "", z))
    c = np.linalg.norm(z)
    k = math.pow((math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma), 2)

    c_z = c*(math.sqrt(abs(1-math.pow(math.cos(beta), 2)- k )))
    c_v = [c*math.cos(beta), ((c*(math.cos(alpha)-math.cos(beta)*math.cos(gamma)))/math.sin(gamma)), c_z]

    abc_m = np.array([a_v, b_v, c_v])
    transpose_abc_m = np.transpose(abc_m)
    for i in range(len(motif_points)):
        motif_points[i].coordinate = transpose_abc_m.dot(motif_points[i].coordinate).tolist()
    r_UnitCellAxis = UnitCellAxis(abc_m, motif_points)
    return r_UnitCellAxis


def build_cell(axis_matrix,lattice_p, motif_points):
    
    # get vectors a_v,b_v,c_v
    a_v = axis_matrix[0]
    b_v = axis_matrix[1]
    c_v = axis_matrix[2]
    new_motif_points = [None]*len(motif_points)

    
    for i in range(len(new_motif_points)):
        new_motif_points[i] = motif_point(list(map(sum,zip(motif_points[i].coordinate,lattice_p))), motif_points[i].labels, motif_points[i].types)
    a_b_p = list(map(sum, zip(a_v, b_v)))
    a_c_p = list(map(sum, zip(a_v, c_v)))
    a_b_c_p = list(map(sum, zip(a_v, b_v, c_v)))
    b_c_p = list(map(sum, zip(c_v, b_v)))
    Z = np.array([lattice_p, list(map(sum,zip(a_v,lattice_p))), 
    list(map(sum,zip(a_b_p,lattice_p))), list(map(sum,zip(b_v,lattice_p))),
    list(map(sum,zip(c_v,lattice_p))),list(map(sum,zip(a_c_p,lattice_p))), 
    list(map(sum,zip(a_b_c_p,lattice_p))),list(map(sum,zip(b_c_p,lattice_p)))])
  
    cell = [Z,new_motif_points]
    return cell


def draw_cell(cell_list, motif_list):
    labels= []
    xs = []
    ys = []
    zs = []
    Z_list = []
    color_scale = px.colors.qualitative.Light24
    fig = go.Figure()
    cell_l = cell_list.cells
    for i in range(len(cell_l)):
        Z = cell_l[i][0]
        Z_list.append(Z.tolist())

        for j in range(len(cell_l[i][1])):
            labels.append(cell_l[i][1][j].labels)
            xs.append((cell_l[i][1][j].coordinate)[0])
            ys.append((cell_l[i][1][j].coordinate)[1])
            zs.append((cell_l[i][1][j].coordinate)[2])
    fig.add_trace(go.Scatter3d(x = xs, y = ys, z = zs, mode = 'markers', hovertext=labels, hoverinfo='text', marker_color = 'rgba(102,102,102,0.7)'))
    for z_element in Z_list:
        verts = [[z_element[0],z_element[1],z_element[2],z_element[3],z_element[0]],
        [z_element[0],z_element[4],z_element[5],z_element[1],z_element[0]], 
        [z_element[0],z_element[4],z_element[7],z_element[3]], 
        [z_element[3],z_element[7],z_element[6],z_element[2]], 
        [z_element[2],z_element[6],z_element[5],z_element[5],z_element[1]]]
        xs = []
        ys=[]
        zs=[]
        for p in verts:
            for x in p:
                xs.append(x[0])
                ys.append(x[1])
                zs.append(x[2])

        fig.add_trace(go.Scatter3d(x = xs, y = ys, z = zs, mode = 'lines', marker_color = 'rgba(0, 0, 0, 1)'))
    
    for motif_p in motif_list:
        xs=[]
        ys=[]
        zs=[]
        labels = []
        xs.append(motif_p[1].coordinate[0])
        xs.append(motif_p[2].coordinate[0])
        ys.append(motif_p[1].coordinate[1])
        ys.append(motif_p[2].coordinate[1])
        zs.append(motif_p[1].coordinate[2])
        zs.append(motif_p[2].coordinate[2])
        labels.append(motif_p[1].labels)
        labels.append(motif_p[2].labels)
        distance = "Nomuber." + str(motif_list.index(motif_p)+1) + " distance: " + str(motif_p[0]) + " between " + str(motif_p[1].labels) + " and " + str(motif_p[2].labels)
        color_num = motif_list.index(motif_p)
        color_num = color_num%len(color_scale)
        fig.add_trace(go.Scatter3d(x = xs, y = ys, z = zs, mode='lines', hovertext=[distance,distance], hoverinfo = 'text', marker_color = color_scale[color_num]))
        fig.add_trace(go.Scatter3d(x = xs, y = ys, z = zs, mode='markers', hovertext=labels, marker_color = color_scale[color_num]))
    fig.update_layout(showlegend = False)
    fig.show()



def compute_distance_with_labels(centerIndex, cell_list):
    dist_list = []
    for j in range(len(cell_list[centerIndex][1])):
        for i in range(len(cell_list)):
                for k in range(len(cell_list[i][1])):
                    center_modif_p = cell_list[centerIndex][1][j]
                    center_modif_p_coor = np.array(center_modif_p.coordinate)
                    p2 = cell_list[i][1][k]
                    p2_coor = np.array(p2.coordinate)
                    squard_dist = np.sum((center_modif_p_coor - p2_coor)**2, axis = 0)
                    dist = np.sqrt(squard_dist)
                    if dist != 0:
                        dist_list.append((dist,center_modif_p,p2))
    sorted_dist_list = sorted(dist_list,key=lambda x: x[0])
    return sorted_dist_list



def compute_distance_kdtree(centerIndex, cell_list, num_dist):
    #get coordinates of all points
    kd_tree_coordinates = []
    for cell in cell_list:
        cell_coordinates = get_coordinate_matrix(cell[1])
        for cell_coordinate in cell_coordinates:
            kd_tree_coordinates.append(cell_coordinate)
    X = np.array(kd_tree_coordinates)
    #buile kd tree
    num_dist = num_dist+1
    nbrs = NearestNeighbors(n_neighbors=num_dist, algorithm='kd_tree').fit(X)

    center_modif_p_coor = get_coordinate_matrix(cell_list[centerIndex][1])
    X = np.array(center_modif_p_coor)
    distance, indices = nbrs.kneighbors(X)
    return distance

def get_coordinate_matrix(cell):
    coordinate_matrix = []
    for point in cell:
        coordinate_matrix.append(point.coordinate)
    return coordinate_matrix

def compute_avrage_dist(dist_list):
    data = np.array(dist_list)
    avrage_list = np.average(data, axis = 0).tolist()
    return avrage_list

def write_csv_file(cif_filename, result, result_file_path):
    fields = ['Name']
    dist = 'DIST'
    for i in range(len(result)):
        fields.append(dist+str(i+1))
    str_filename = cif_filename.split('.')
    csv_name = str_filename[0] + ".csv"
    file_exist = os.path.isfile(result_file_path)
    with open(result_file_path,"a+", newline='') as result_file:
        #check if csv file has header
        wr = csv.DictWriter(result_file,fieldnames = fields)
        if file_exist == False:
            wr.writeheader()

        wr = csv.writer(result_file)
        wr.writerow([csv_name]+result)

def write_bond_file(results,files):
    head_name = ['Bond', 'Number']
    filepath = os.path.split(files)[0]
    type_bonds = {}
    result_path = os.path.join(filepath, 'bond_num.csv')
    for result in results:
        type_1 = result[1].types
        type_2 = result[2].types
        type_bond = type_1 + '-' + type_2
        if type_bond in type_bonds:
            type_bonds[type_bond] = type_bonds[type_bond] + 1
        else:
            type_bonds[type_bond] = 1

    with open(result_path, "w") as bond_file:
        wr = csv.DictWriter(bond_file, fieldnames=head_name)
        wr.writeheader()
        for key in type_bonds.keys():
            bond_file.write("%s, %s\n" % (key, type_bonds[key]))

def get_cellList(cif_filename, num_expand, atom_type):
    cf = ReadCif(cif_filename)
    list_c = []
    e=int(num_expand/2)
    # axis, motif_points = getUnitCellAxis(cf)
    UnitCellAxis = getUnitCellAxis(cf)
    selected_motif_points = []
    a_v = UnitCellAxis.vector[0]
    b_v = UnitCellAxis.vector[1]
    c_v = UnitCellAxis.vector[2]
    if atom_type:
        for motif_point in UnitCellAxis.motif_point:
            if motif_point.types in atom_type and motif_point.types != 'H':
                selected_motif_points.append(motif_point)
    else:
        for motif_point in UnitCellAxis.motif_point:
            if motif_point.types != 'H':
                selected_motif_points.append(motif_point)
    for i in range(-e,e+1):
        for j in range(-e,e+1):
            for l in range(-e,e+1):
                lattice_p = i * a_v + j * b_v + l * c_v
                cell = build_cell(UnitCellAxis.vector, lattice_p, selected_motif_points)
                list_c.append(cell)
                # save (i,j,l) in a list or dictionary
    center_index = 0
    for i in range(len(list_c)):
        if list_c[i][1] == selected_motif_points:
            center_index = i
    cellList = list_cell(list_c, center_index)
    return cellList

def get_DistList(cellList,num_dist):
    center_index = cellList.center_index
    cellList = cellList.cells 
    dist_list = compute_distance_with_labels(center_index, cellList)
    new_dist_list = remove_repeat_item(dist_list, num_dist)
    return  new_dist_list


def remove_repeat_item(dist_list,num_dist):
    new_dist_list = [dist_list[0]]
    for i in range(len(dist_list)):
        repeat = False
        for j in range(len(new_dist_list)):
            if new_dist_list[j][0] == dist_list[i][0]:
                repeat = True
        if repeat == False:
            new_dist_list.append(dist_list[i])
        if len(new_dist_list) == num_dist:
            break
    
    return new_dist_list



def get_AvgDist(cellList, num_dist):    
    center_index = cellList.center_index
    cellList = cellList.cells 
    dist_list_2d = compute_distance_kdtree(center_index, cellList, num_dist)
    avrage_list = compute_avrage_dist(dist_list_2d)[1:]
    return avrage_list

def final_step(files, atom_type, num_file,num_expand=5, num_dist=200):
    if os.path.isfile(files) and files.endswith(".cif"):
        cellList = get_cellList(files,num_expand,atom_type)
        dist_list= get_DistList(cellList, num_dist)
        write_bond_file(dist_list, files)
        draw_cell(cellList,dist_list)
    else:
        result_files_path = os.path.join(files, 'result_dir.csv')
        if os.path.isfile(result_files_path):
            os.remove(result_files_path)
        for cif_file in os.listdir(files):
            if cif_file.endswith(".cif"):
                cif_fullname = os.path.join(files, cif_file)
                cellList = get_cellList(cif_fullname, num_expand, atom_type)
                result = get_AvgDist(cellList, num_dist)
                write_csv_file(cif_file, result,result_files_path)
        clustering.do_clustering(files, num_file)
        
 




