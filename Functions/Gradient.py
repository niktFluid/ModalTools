import math
import numpy as np
from scipy import linalg

from Functions.BoundaryCondition import BoundaryCondition as BDcond
from Functions.Variables import Variables


class Gradient(Variables):
    def __init__(self, mesh, grad_type='GLSQ', is2d=False):
        super(Gradient, self).__init__(mesh, n_return=3)
        self.bd_cond = BDcond(mesh, is2d)

        self._grad_type = grad_type
        self._beta = None

        # Least Square matrix
        n_cell = mesh.n_cell
        self.matLU1 = np.zeros((n_cell, 3, 3), dtype=np.float64)
        self.matLU2 = np.zeros((n_cell, 3), dtype=np.float64)
        self._set_left_mat()

    def _return_ref_cells(self, id_cell):
        cell_list = [id_cell] + self.mesh.cell_neighbours(id_cell)
        ref_cells = [i_cell for i_cell in cell_list if not self.mesh.is_boundary_cell(i_cell)]

        return list(set(ref_cells))

    def formula(self, data, id_cell, id_val=0):
        vec_rhs = self._set_rhs(data, id_cell, id_val)
        return linalg.lu_solve((self.matLU1[id_cell], self.matLU2[id_cell]), vec_rhs)

    def _set_left_mat(self):
        self._beta = np.empty(self.mesh.n_cell, dtype=np.float64)
        if self._grad_type == 'LS':
            self._beta = 1.0
        elif self._grad_type == 'GG':
            self._beta = 0.0
        elif self._grad_type == 'GLSQ':
            for i_cell in range(self.mesh.n_cell):
                face_dis_max = 1.0e-30
                face_area_max = 1.0e-30

                for i_face in self.mesh.cell_faces[i_cell]:
                    r_vec_0 = self.mesh.face_centers[i_face] - self.mesh.centers[i_cell]
                    dis0 = math.sqrt(r_vec_0[0] * r_vec_0[0] + r_vec_0[1] * r_vec_0[1] + r_vec_0[2] * r_vec_0[2])
                    if dis0 > face_dis_max:
                        face_dis_max = dis0

                    area = self.mesh.face_area[i_face]
                    if area > face_area_max:
                        face_area_max = area

                self._beta[i_cell] = min(1.0, self.mesh.volumes[i_cell] / (face_dis_max * face_area_max))

        for id_cell in range(self.mesh.n_cell):
            matA = self._set_left_mat_by_cell(id_cell)
            self.matLU1[id_cell], self.matLU2[id_cell] = linalg.lu_factor(matA)

    def _set_left_mat_by_cell(self, id_cell):
        # Left side matrix for gradient calculation.
        mat_cell = np.zeros((3, 3), dtype=np.float64)

        if self._grad_type == 'LS' or self._grad_type == 'GLSQ':
            mat_cell = self._ls_left_mat(mat_cell, id_cell)
        elif self._grad_type == 'GG' or self._grad_type == 'GLSQ':
            mat_cell = self._gg_left_mat(mat_cell, id_cell)

        return mat_cell

    def _ls_left_mat(self, mat_cell, id_cell):
        nb_cells = self.mesh.cell_neighbours(id_cell)
        faces = self.mesh.cell_faces[id_cell]

        for id_nb, id_face in zip(nb_cells, faces):
            flip = self.mesh.get_face_direction(id_cell, id_nb, id_face)
            vec_lr = self.mesh.vec_lr[id_face] * flip

            weight_ls = self._get_weight_ls(id_cell, id_nb, id_face)
            mat_cell += vec_lr.reshape(3, 1) @ vec_lr.reshape(1, 3) * weight_ls * self._beta[id_cell]
        return mat_cell

    def _get_weight_ls(self, id_cell, id_nb, id_face):
        vec_lr = self.mesh.vec_lr[id_face]
        lr_inv = 1.0 / math.sqrt(vec_lr[0] * vec_lr[0] + vec_lr[1] * vec_lr[1] + vec_lr[2] * vec_lr[2])

        r_vec_0 = self.mesh.face_centers[id_face] - self.mesh.centers[id_cell]
        dis0 = math.sqrt(r_vec_0[0] * r_vec_0[0] + r_vec_0[1] * r_vec_0[1] + r_vec_0[2] * r_vec_0[2])
        r_vec_1 = self.mesh.face_centers[id_face] - self.mesh.centers[id_nb]
        dis1 = math.sqrt(r_vec_1[0] * r_vec_1[0] + r_vec_1[1] * r_vec_1[1] + r_vec_1[2] * r_vec_1[2])
        dis_inv = 1.0 / (dis0 + dis1)

        area = self.mesh.face_area[id_face]

        return (2.0 * dis0 * dis_inv) ** 2 * area * lr_inv

    def _gg_left_mat(self, mat_cell, id_cell):
        volume = self.mesh.volumes[id_cell]
        mat_eye = np.eye(3, dtype=np.float64)

        mat_cell += 2.0 * volume * mat_eye * (1.0 - self._beta[id_cell])
        return mat_cell

    def _set_rhs(self, data, id_cell, id_val):
        nb_cells = self.mesh.cell_neighbours(id_cell)
        faces = self.mesh.cell_faces[id_cell]

        val_vec = np.array([data[id_cell, i_val] for i_val in range(data.shape[1])])

        rhs_vec = np.zeros(3, dtype=np.float64)
        for id_nb, id_face in zip(nb_cells, faces):
            val_diff = self._get_val_diff(data, val_vec, id_nb, id_face, id_val)

            if self._grad_type == 'LS' or self._grad_type == 'GLSQ':
                rhs_vec = self._ls_rhs(rhs_vec, val_diff, id_cell, id_nb, id_face)
            elif self._grad_type == 'GG' or self._grad_type == 'GLSQ':
                rhs_vec = self._gg_rhs(rhs_vec, val_diff, id_cell, id_nb, id_face)

        return rhs_vec

    def _ls_rhs(self, rhs_vec, val_diff, id_cell, id_nb, id_face):
        flip = self.mesh.get_face_direction(id_cell, id_nb, id_face)
        vec_lr = self.mesh.vec_lr[id_face] * flip

        weight_ls = self._get_weight_ls(id_cell, id_nb, id_face)
        rhs_vec += val_diff * vec_lr * weight_ls * self._beta[id_cell]

        return rhs_vec

    def _gg_rhs(self, rhs_vec, val_diff, id_cell, id_nb, id_face):
        flip = self.mesh.get_face_direction(id_cell, id_nb, id_face)
        vec_n = self.mesh.face_mat[id_face, 0] * flip
        area = self.mesh.face_area[id_face]

        rhs_vec += val_diff * area * vec_n * (1.0 - self._beta[id_cell])
        return rhs_vec

    def _get_val_diff(self, data, vec_0, id_k, id_k_face, id_val):
        if not self.mesh.is_boundary_face(id_k_face):  # For inner cells
            return data[id_k, id_val] - vec_0[id_val]
        else:  # For boundary cells
            val_bd = self.bd_cond.get_bd_val(vec_0, id_k_face)
            return val_bd[id_val] - vec_0[id_val]
