# import math
import numpy as np
from scipy import linalg
# from scipy.sparse import csr_matrix

from Functions.BoundaryCondition import BoundaryCondition as BDcond
from Functions.Variables import Variables


class Gradient(Variables):
    def __init__(self, mesh, is2d=False):
        super(Gradient, self).__init__(mesh, n_return=3)
        self.bd_cond = BDcond(mesh, is2d)
        # self.axis = axis

        # self._val_vec = np.empty(5, dtype=np.float64)

        # Least Square matrix
        n_cell = mesh.n_cell
        self.matLU1 = np.zeros((n_cell, 3, 3), dtype=np.float64)
        self.matLU2 = np.zeros((n_cell, 3), dtype=np.float64)
        self._set_left_mat()

    def _return_ref_cells(self, id_cell):
        cell_list = [id_cell] + self.mesh.cell_neighbours(id_cell)
        ref_cells = [i_cell for i_cell in cell_list if i_cell >= 0]

        return list(set(ref_cells))

    def formula(self, data, id_cell, id_val=0):
        vec_rhs = self._set_rhs(data, id_cell, id_val)
        return linalg.lu_solve((self.matLU1[id_cell], self.matLU2[id_cell]), vec_rhs)

    def _set_left_mat(self):
        for id_cell in range(self.mesh.n_cell):
            matA = self._set_left_mat_by_cell(id_cell)
            self.matLU1[id_cell], self.matLU2[id_cell] = linalg.lu_factor(matA)

    def _set_left_mat_by_cell(self, id_cell):
        # Left side matrix for the least square method.
        nb_cells = self.mesh.cell_neighbours(id_cell)
        faces = self.mesh.cell_faces[id_cell]

        mat_cell = np.zeros((3, 3), dtype=np.float64)
        for id_nb, id_face in zip(nb_cells, faces):
            vec_lr = self._get_pos_diff(id_cell, id_nb, id_face)
            mat_cell += vec_lr.reshape(3, 1) @ vec_lr.reshape(1, 3)

        return mat_cell

    def _get_pos_diff(self, id_cell, nb_cell, nb_face):
        flip = -1.0 + 2.0 * float(nb_cell > id_cell or nb_cell < 0)
        return self.mesh.vec_lr[nb_face] * flip

    def _set_rhs(self, data, id_cell, id_val):
        nb_cells = self.mesh.cell_neighbours(id_cell)
        faces = self.mesh.cell_faces[id_cell]

        val_vec = np.array([data[id_cell, i_val] for i_val in range(data.shape[1])])
        # val_vec = np.zeros(data.shape[1], dtype=np.float64)
        # for i_val in range(data.shape[1]):
        #     val_vec[i_val] = data[id_cell, i_val]

        rhs_vec = np.zeros(3, dtype=np.float64)
        for id_nb, id_face in zip(nb_cells, faces):
            vec_lr = self._get_pos_diff(id_cell, id_nb, id_face)
            val_diff = self._get_val_diff(data, val_vec, id_nb, id_face, id_val)
            rhs_vec += val_diff * vec_lr
        return rhs_vec

    def _get_val_diff(self, data, vec_0, id_k, id_k_face, id_val):
        if not self.mesh.is_boundary_face(id_k_face):  # For inner cells
            return data[id_k, id_val] - vec_0[id_val]
        else:  # For boundary cells
            val_bd = self.bd_cond.get_bd_val(vec_0, id_k_face)
            return val_bd[id_val] - vec_0[id_val]
