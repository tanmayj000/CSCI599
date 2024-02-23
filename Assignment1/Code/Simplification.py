import numpy as np
import scipy as sp
import heapq
import copy
import os

class Mesh:
    def __init__(self, path, build_code=False, build_mat=False):
        self.path = path
        self.vertices, self.faces = self.read_obj(path)
        self.compute_face_normals()
        self.compute_face_center()
        self.simp = False      
        self.build_adj_matrices() 
        self.build_v2v()
        self.build_vf()

    def read_obj(self, path):
        vertices, faces = [], []
        f = open(path)
        for line in f:
            line = line.strip()
            splitted_line = line.split()
            if not splitted_line:
                continue
            elif splitted_line[0] == 'v':
                vertices.append([float(v) for v in splitted_line[1:4]])
            elif splitted_line[0] == 'f':
                face_vertex_ids = [int(c.split('/')[0]) for c in splitted_line[1:]]
                assert len(face_vertex_ids) == 3
                face_vertex_ids = [(ind - 1) if (ind >= 0) else (len(vertices) + ind) for ind in face_vertex_ids]
                faces.append(face_vertex_ids)
        f.close()
        vertices = np.asarray(vertices)
        faces = np.asarray(faces, dtype=int)

        assert np.logical_and(faces >= 0, faces < len(vertices)).all()
        return vertices, faces

    def build_adj_matrices(self):
        self.vertex_edge = [[] for _ in self.vertices]
        self.vertex_edge_idx = [[] for _ in self.vertices]
        edge_nb = []
        sides = []
        edge2key = dict()
        edges = []
        edges_count = 0
        nb_count = []
        for face_id, face in enumerate(self.faces):
            faces_edges = []
            for i in range(3):
                cur_edge = (face[i], face[(i + 1) % 3])
                faces_edges.append(cur_edge)
            for idx, edge in enumerate(faces_edges):
                edge = tuple(sorted(list(edge)))
                faces_edges[idx] = edge
                if edge not in edge2key:
                    edge2key[edge] = edges_count
                    edges.append(list(edge))
                    edge_nb.append([-1, -1, -1, -1])
                    sides.append([-1, -1, -1, -1])
                    self.vertex_edge[edge[0]].append(edges_count)
                    self.vertex_edge[edge[1]].append(edges_count)
                    self.vertex_edge_idx[edge[0]].append(0)
                    self.vertex_edge_idx[edge[1]].append(1)
                    nb_count.append(0)
                    edges_count += 1
            for idx, edge in enumerate(faces_edges):
                edge_key = edge2key[edge]
                edge_nb[edge_key][nb_count[edge_key]] = edge2key[faces_edges[(idx + 1) % 3]]
                edge_nb[edge_key][nb_count[edge_key] + 1] = edge2key[faces_edges[(idx + 2) % 3]]
                nb_count[edge_key] += 2
            for idx, edge in enumerate(faces_edges):
                edge_key = edge2key[edge]
                sides[edge_key][nb_count[edge_key] - 2] = nb_count[edge2key[faces_edges[(idx + 1) % 3]]] - 1
                sides[edge_key][nb_count[edge_key] - 1] = nb_count[edge2key[faces_edges[(idx + 2) % 3]]] - 2
        self.edges = np.array(edges, dtype=np.int32)
        self.gemm_edges = np.array(edge_nb, dtype=np.int64)
        self.sides = np.array(sides, dtype=np.int64)
        self.edges_count = edges_count

    def compute_face_normals(self):
        face_normals = np.cross(self.vertices[self.faces[:, 1]] - self.vertices[self.faces[:, 0]], self.vertices[self.faces[:, 2]] - self.vertices[self.faces[:, 0]])
        norm = np.linalg.norm(face_normals, axis=1, keepdims=True) + 1e-24
        face_areas = 0.5 * np.sqrt((face_normals**2).sum(axis=1))
        face_normals /= norm
        self.face_normals, self.fa = face_normals, face_areas

    
    def compute_face_center(self):
        faces = self.faces
        vertices = self.vertices
        self.face_centers = np.sum(vertices[faces], 1) / 3.0

    def build_vf(self):
        vf = [set() for _ in range(len(self.vertices))]
        for i, f in enumerate(self.faces):
            vf[f[0]].add(i)
            vf[f[1]].add(i)
            vf[f[2]].add(i)
        self.vf = vf

    def build_v2v(self):
        v2v = [[] for _ in range(len(self.vertices))]
        for i, e in enumerate(self.edges):
            v2v[e[0]].append(e[1])
            v2v[e[1]].append(e[0])
        self.v2v = v2v

        edges = self.edges
        v2v_inds = edges.T
        v2v_inds = np.concatenate([v2v_inds, v2v_inds[[1, 0]]], axis=1).astype(np.int64)
        v2v_vals = np.ones(v2v_inds.shape[1], dtype=np.float32)
        self.v2v_mat = sp.sparse.csr_matrix((v2v_vals, v2v_inds), shape=(len(self.vertices), len(self.vertices)))
        self.v_dims = np.sum(self.v2v_mat.toarray(), axis=1)

    def simplification(self, target_v, valence_aware=True, midpoint=False):
        vertices, vf, fn, fc, edges = self.vertices, self.vf, self.face_normals, self.face_centers, self.edges

        """ 1. compute Q for each vertex """
        Q_s = [[] for _ in range(len(vertices))]
        E_s = [[] for _ in range(len(vertices))]
        for i, v in enumerate(vertices):
            f_s = np.array(list(vf[i]))
            fc_s = fc[f_s]
            fn_s = fn[f_s]
            d_s = - 1.0 * np.sum(fn_s * fc_s, axis=1, keepdims=True)
            abcd_s = np.concatenate([fn_s, d_s], axis=1)
            Q_s[i] = np.matmul(abcd_s.T, abcd_s)
            v4 = np.concatenate([v, np.array([1])])
            E_s[i] = np.matmul(v4, np.matmul(Q_s[i], v4.T))

        """ 2. compute E for every possible pairs and create heapq """
        E_heap = []
        for i, e in enumerate(edges):
            v_0, v_1 = vertices[e[0]], vertices[e[1]]
            Q_0, Q_1 = Q_s[e[0]], Q_s[e[1]]
            Q_new = Q_0 + Q_1

            Q_lp = np.eye(4)
            Q_lp[:3] = Q_new[:3]
            try:
                Q_lp_inv = np.linalg.inv(Q_lp)
                v4_new = np.matmul(Q_lp_inv, np.array([[0,0,0,1]]).reshape(-1,1)).reshape(-1)
            except:
                v_new = 0.5 * (v_0 + v_1)
                v4_new = np.concatenate([v_new, np.array([1])])
            
            
            E_new = np.matmul(v4_new, np.matmul(Q_new, v4_new.T))
            heapq.heappush(E_heap, (E_new, (e[0], e[1])))
        
        """ 3. collapse minimum-error vertex """
        simp_mesh = copy.deepcopy(self)

        vi_mask = np.ones([len(simp_mesh.vertices)]).astype(np.bool_)
        fi_mask = np.ones([len(simp_mesh.faces)]).astype(np.bool_)

        vert_map = [{i} for i in range(len(simp_mesh.vertices))]
        # pbar = tqdm(total=np.sum(vi_mask)-target_v, desc="Processing")
        while np.sum(vi_mask) > target_v:
            if len(E_heap) == 0:
                break

            E_0, (vi_0, vi_1) = heapq.heappop(E_heap)

            if (vi_mask[vi_0] == False) or (vi_mask[vi_1] == False):
                continue

            """ edge collapse """
            shared_vv = list(set(simp_mesh.v2v[vi_0]).intersection(set(simp_mesh.v2v[vi_1])))
            merged_faces = simp_mesh.vf[vi_0].intersection(simp_mesh.vf[vi_1])

            if len(shared_vv) != 2:
                continue

            elif len(merged_faces) != 2:
                continue

            else:
                self.edge_collapse(simp_mesh, vi_0, vi_1, merged_faces, vi_mask, fi_mask, vert_map, Q_s, E_heap)
                # print(np.sum(vi_mask), np.sum(fi_mask))
        
        self.rebuild_mesh(simp_mesh, vi_mask, fi_mask, vert_map)
        simp_mesh.simp = True        
        return simp_mesh
    
    
    def edge_collapse(self, simp_mesh, vi_0, vi_1, merged_faces, vi_mask, fi_mask, vert_map, Q_s, E_heap):
        shared_vv = list(set(simp_mesh.v2v[vi_0]).intersection(set(simp_mesh.v2v[vi_1])))
        new_vi_0 = set(simp_mesh.v2v[vi_0]).union(set(simp_mesh.v2v[vi_1])).difference({vi_0, vi_1})
        simp_mesh.vf[vi_0] = simp_mesh.vf[vi_0].union(simp_mesh.vf[vi_1]).difference(merged_faces)
        simp_mesh.vf[vi_1] = set()
        simp_mesh.vf[shared_vv[0]] = simp_mesh.vf[shared_vv[0]].difference(merged_faces)
        simp_mesh.vf[shared_vv[1]] = simp_mesh.vf[shared_vv[1]].difference(merged_faces)

        simp_mesh.v2v[vi_0] = list(new_vi_0)
        for v in simp_mesh.v2v[vi_1]:
            if v != vi_0:
                simp_mesh.v2v[v] = list(set(simp_mesh.v2v[v]).difference({vi_1}).union({vi_0}))
        simp_mesh.v2v[vi_1] = []
        vi_mask[vi_1] = False

        vert_map[vi_0] = vert_map[vi_0].union(vert_map[vi_1])
        vert_map[vi_0] = vert_map[vi_0].union({vi_1})
        vert_map[vi_1] = set()
        
        fi_mask[np.array(list(merged_faces)).astype(np.int32)] = False

        simp_mesh.vertices[vi_0] = 0.5 * (simp_mesh.vertices[vi_0] + simp_mesh.vertices[vi_1])

        """ recompute E """
        Q_0 = Q_s[vi_0]
        for vv_i in simp_mesh.v2v[vi_0]:
            v_mid = 0.5 * (simp_mesh.vertices[vi_0] + simp_mesh.vertices[vv_i])
            Q_1 = Q_s[vv_i]
            Q_new = Q_0 + Q_1
            v4_mid = np.concatenate([v_mid, np.array([1])])

            E_new = np.matmul(v4_mid, np.matmul(Q_new, v4_mid.T)) 
            heapq.heappush(E_heap, (E_new, (vi_0, vv_i)))

    @staticmethod
    def rebuild_mesh(simp_mesh, vi_mask, fi_mask, vert_map):
        face_map = dict(zip(np.arange(len(vi_mask)), np.cumsum(vi_mask)-1))
        simp_mesh.vertices = simp_mesh.vertices[vi_mask]
        
        vert_dict = {}
        for i, vm in enumerate(vert_map):
            for j in vm:
                vert_dict[j] = i

        for i, f in enumerate(simp_mesh.faces):
            for j in range(3):
                if f[j] in vert_dict:
                    simp_mesh.faces[i][j] = vert_dict[f[j]]

        simp_mesh.faces = simp_mesh.faces[fi_mask]
        for i, f in enumerate(simp_mesh.faces):
            for j in range(3):
                simp_mesh.faces[i][j] = face_map[f[j]]
        
        simp_mesh.compute_face_normals()
        simp_mesh.compute_face_center()
        simp_mesh.build_adj_matrices()
        simp_mesh.build_v2v()
        simp_mesh.build_vf()
            
    def save(self, filename):
        assert len(self.vertices) > 0
        vertices = np.array(self.vertices, dtype=np.float32).flatten()
        indices = np.array(self.faces, dtype=np.uint32).flatten()

        with open(filename, 'w') as fp:
            # Write positions
            for i in range(0, vertices.size, 3):
                x = vertices[i + 0]
                y = vertices[i + 1]
                z = vertices[i + 2]
                fp.write('v {0:.8f} {1:.8f} {2:.8f}\n'.format(x, y, z))

            # Write indices
            for i in range(0, len(indices), 3):
                i0 = indices[i + 0] + 1
                i1 = indices[i + 1] + 1
                i2 = indices[i + 2] + 1
                fp.write('f {0} {1} {2}\n'.format(i0, i1, i2))


    
def main():
    inp_base = 'assets/input/'
    op_base = 'assets/output/quadratic_error_simplification/'
    probability_list = [0.2, 0.4, 0.6, 0.9]

    files = os.listdir(inp_base)

    obj_files = [file for file in files if file.endswith('.obj')]

    print(obj_files)

    for file in obj_files:
        directory = inp_base + file
        mesh = Mesh(directory)
        for probability in probability_list:
            target_v = int(len(mesh.vertices) * probability)
            simp_mesh = mesh.simplification(target_v=target_v)
            simp_mesh.save(op_base + file.split('.')[0] + '_' + str(int(probability*100)) + 'perc.obj')
    # print("[FIN] Simplification Completed!")

if __name__ == "__main__":
    main()