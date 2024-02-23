import mesh
import math
import numpy as np
from collections import defaultdict
import trimesh
import collections
dictt = collections.defaultdict(list) 
class Edge:
    def __init__(self, vertex1, vertex2):
        self.vertices = (vertex1, vertex2)
        self.vertices = tuple(sorted(self.vertices))

    def __eq__(self, other):
        return self.vertices == other.vertices

    def __hash__(self):
        return hash(self.vertices)

    def __repr__(self):
        return f"Edge({self.vertices[0]}, {self.vertices[1]})"

def get_coordinate(vertex_index):
    return vertices[vertex_index].get_vertex()

def sharedFace(my_edge, face_points):
    remote_points = [0] * 2
    for fp in face_points:
        if my_edge.vertices[0] != fp and my_edge.vertices[1] != fp:
            remote_points = fp

    return remote_points

def isShareEdge(my_edge, facets):
    count = 0
    remote_point_list = []
    for face in facets:
        a = get_coordinate(face.a)
        b = get_coordinate(face.b)
        c = get_coordinate(face.c)

        face_points = [a,b,c]
        
        if my_edge.vertices[0] in face_points and my_edge.vertices[1] in face_points:
            remote_point_list.append(sharedFace(my_edge, face_points))
            count += 1
        if count == 2:
            return remote_point_list
    if count == 1:
        return remote_point_list
    else:
        return

def get_interior_odd_vertex(edge, remote1, remote2):
    v0 = edge.vertices[0]
    v1 = remote1
    v2 = edge.vertices[1]
    v3 = remote2

    new_x = 3/8 * (v0[0] + v2[0]) + 1/8 * (v3[0] + v1[0])
    new_y = 3/8 * (v0[1] + v2[1]) + 1/8 * (v3[1] + v1[1])
    new_z = 3/8 * (v0[2] + v2[2]) + 1/8 * (v3[2] + v1[2])
    return [new_x, new_y, new_z]

def get_boundary_odd_vertex(edge):
    v0 = edge.vertices[0]
    v2 = edge.vertices[1]

    new_x = 1/2 * (v0[0] + v2[0])
    new_y = 1/2 * (v0[1] + v2[1])
    new_z = 1/2 * (v0[2] + v2[2])
    
    return [new_x, new_y, new_z]


mp_dict = defaultdict(list)
boundary = []
def get_odd_vertice_and_new_facets(vertices, facets):
    new_facets = []
    new_vertices = []
    
    for vertex in vertices:
        new_vertices.append(vertex.get_vertex())

    for face in facets:

        a = get_coordinate(face.a)
        b = get_coordinate(face.b)
        c = get_coordinate(face.c)

        edge1 = Edge(a,b)
        edge2 = Edge(b,c)
        edge3 = Edge(c,a)

        remote_points_edge_1 = isShareEdge(edge1, facets)
        if len(remote_points_edge_1) == 2:
            odd_vertex_ab = get_interior_odd_vertex(edge1, remote_points_edge_1[0], remote_points_edge_1[1])
        elif len(remote_points_edge_1) == 1: 
            # boundary mapping
            odd_vertex_ab = get_boundary_odd_vertex(edge1)
            boundary.append(face.a)
            mp_dict[face.a].append(face.b)
            mp_dict[face.b].append(face.a)
    
        remote_points_edge_2 = isShareEdge(edge2, facets)
        if len(remote_points_edge_2) == 2:
            odd_vertex_bc = get_interior_odd_vertex(edge2, remote_points_edge_2[0], remote_points_edge_2[1])
        elif len(remote_points_edge_2) == 1: 
            odd_vertex_bc = get_boundary_odd_vertex(edge2)
            boundary.append(face.b)
            mp_dict[face.c].append(face.b)
            mp_dict[face.b].append(face.c)

        remote_points_edge_3 = isShareEdge(edge3, facets)
        if len(remote_points_edge_3) == 2:
            odd_vertex_ca = get_interior_odd_vertex(edge3, remote_points_edge_3[0], remote_points_edge_3[1])
        elif len(remote_points_edge_3) == 1:  
            odd_vertex_ca = get_boundary_odd_vertex(edge3)
            boundary.append(face.c)
            mp_dict[face.c].append(face.a)
            mp_dict[face.a].append(face.c)
 
        if odd_vertex_ab in new_vertices:
            index_ab = new_vertices.index(odd_vertex_ab)
        else:
            new_vertices.append(odd_vertex_ab)
            index_ab = new_vertices.index(odd_vertex_ab)

        if odd_vertex_bc in new_vertices:
            index_bc = new_vertices.index(odd_vertex_bc)
        else:
            new_vertices.append(odd_vertex_bc)
            index_bc = new_vertices.index(odd_vertex_bc)

        if odd_vertex_ca in new_vertices:
            index_ca = new_vertices.index(odd_vertex_ca)
        else:
            new_vertices.append(odd_vertex_ca)
            index_ca = new_vertices.index(odd_vertex_ca)     
        
        new_face_0 = [face.a, index_ab, index_ca]
        new_face_1 = [index_ab, face.b, index_bc]
        new_face_2 = [index_bc, face.c, index_ca]
        new_face_3 = [index_ab, index_bc, index_ca]

        new_facets.append(new_face_0)
        new_facets.append(new_face_1)
        new_facets.append(new_face_2)
        new_facets.append(new_face_3)

    return new_vertices, new_facets

def calculate_neighbor(vertices, facets):
    neighbor_index = defaultdict(list) 
    neighbor_num = {} 

    for vertex in vertices:

        for f in facets:
            fff = [f.a, f.b, f.c]
            if vertex.index in fff:
                neighbor_num[vertex.index] = neighbor_num.get(vertex.index, 0) + 1
                for i in fff:
                    if vertex.index != i and i not in neighbor_index[vertex.index]:
                        neighbor_index[vertex.index].append(i)
    return neighbor_index, neighbor_num

def calculate_beta(neighbors_count):
    if neighbors_count == 3:
        return 3/16
    else:
        return (1.0 / neighbors_count) * (5.0 / 8.0 - ((3.0 / 8.0) + (1.0 / 4.0) * math.cos(2 * math.pi / neighbors_count))** 2)
    

def update_old_vertice(vertices, facets):
    neighbor_index, neighbor_num = calculate_neighbor(vertices, facets)

    bound = [] # store all the boundary index into list
    count_interval = 0
    count_allsmall = 0
    count_allbig = 0

    bbo = set(boundary)
   
    update_old = []
    for i in range(len(neighbor_num)):

        total_vi_x = 0
        total_vi_y = 0
        total_vi_z = 0

        if i in bbo:
            closest_num1 = mp_dict[i][0]
            closest_num2 = mp_dict[i][1]
            dictt[i].append(closest_num1)
            dictt[i].append(closest_num2)
            total_vi_x = vertices[closest_num1].get_vertex()[0] + vertices[closest_num2].get_vertex()[0] 
            total_vi_y = vertices[closest_num1].get_vertex()[1] + vertices[closest_num2].get_vertex()[1]
            total_vi_z = vertices[closest_num1].get_vertex()[2] + vertices[closest_num2].get_vertex()[2]

            v_update_x = (3.0 / 4.0) * vertices[i].get_vertex()[0] + (1.0 / 8.0) * total_vi_x
            v_update_y = (3.0 / 4.0) * vertices[i].get_vertex()[1] + (1.0 / 8.0) * total_vi_y
            v_update_z = (3.0 / 4.0) * vertices[i].get_vertex()[2] + (1.0 / 8.0) * total_vi_z
            update_old.append([round(v_update_x, 7), round(v_update_y, 7), round(v_update_z, 7)])
        else:
            for j in neighbor_index[i]:
                total_vi_x += calculate_beta(neighbor_num[i]) * vertices[j].get_vertex()[0] 
            for j in neighbor_index[i]:
                total_vi_y += calculate_beta(neighbor_num[i]) * vertices[j].get_vertex()[1] 
            for j in neighbor_index[i]:
                total_vi_z += calculate_beta(neighbor_num[i]) * vertices[j].get_vertex()[2] 
            v_update_x = (1-neighbor_num[i]*calculate_beta(neighbor_num[i])) * vertices[i].get_vertex()[0] + total_vi_x
            v_update_y = (1-neighbor_num[i]*calculate_beta(neighbor_num[i])) * vertices[i].get_vertex()[1] + total_vi_y
            v_update_z = (1-neighbor_num[i]*calculate_beta(neighbor_num[i])) * vertices[i].get_vertex()[2] + total_vi_z

            update_old.append([round(v_update_x, 6), round(v_update_y, 6), round(v_update_z, 6)])
            #update_old.append([v_update_x, v_update_y, v_update_z])
    # print(neighbor_num)
    # print(neighbor_index)
    # print(count_interval)
    # print(count_allsmall)
    # print(count_allbig)

    return update_old


def save_halfmesh_as_obj(all_vertex, new_facet, file_name):
    with open(file_name, 'w') as open_file:
        for vertex in all_vertex:
            x = vertex[0]
            y = vertex[1]
            z = vertex[2]
            open_file.write("v {} {} {}\n".format(x, y, z))

        for face in new_facet:
            f0 = face[0]
            f1 = face[1]
            f2 = face[2]
            open_file.write("f {} {} {}\n".format(f0+1, f1+1, f2+1))




data_path = 'assets/input/cube.obj'
iterations = 3
tmesh = trimesh.load_mesh(data_path)


for i in range(iterations):
    m = mesh.HalfedgeMesh(tmesh)
    vertices, facets = m.vertices, m.facets
    new_vertice, new_facets = get_odd_vertice_and_new_facets(vertices, facets)

    old_updated_vertice = update_old_vertice(vertices, facets)

    new_vertice[:(len(old_updated_vertice))] = old_updated_vertice

    save_halfmesh_as_obj(new_vertice, new_facets, 'assets/output/loop_subdivision/cube_subdivided' + str(i) + '.obj')
    tmesh = trimesh.load_mesh('assets/output/loop_subdivision/cube_subdivided' + str(i) + '.obj')

