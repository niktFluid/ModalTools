import math
import numpy as np
from scipy import linalg

from Functions.Variables import Variables


class Gradient(Variables):
    def __init__(self, mesh, bd_cond, grad_type='GLSQ'):
        """
        Calculate gradient of the variables on an unstructured grid.
        The implementation is based on 'Shima, E., Kitamura, K., and Haga, T. , AIAA Journal,
         Vol. 51, No. 11, 2013, pp. 2740–2747. doi:10.2514/1.J052095.'

        :param mesh: Class contains computational grid and related methods.
        :param bd_cond: Class for boundary condition.
        :param grad_type:   GG - Green gauss method.
                            LS - Weighted least square method.
                            GLSQ - Hybrid method of GG and LS.
        """

        super(Gradient, self).__init__(mesh, n_return=3)
        self.bd_cond = bd_cond

        self._grad_type = grad_type
        self._beta = None

        # Left matrix for gradient calculation.
        n_cell = mesh.n_cell
        self.matLU1 = np.zeros((n_cell, 3, 3), dtype=np.float64)
        self.matLU2 = np.zeros((n_cell, 3), dtype=np.float64)
        self._set_left_mat()

    def _return_ref_cells(self, id_cell):
        """
        Return list of cells which are 'directory' referenced from this class for calculation.
        For this class, the list contains the target cell and the adjoined cells.

        :param id_cell: Target cell.
        :return: List of cell.
        """

        cell_list = [id_cell] + self.mesh.cell_neighbours(id_cell)
        ref_cells = [i_cell for i_cell in cell_list if not self.mesh.is_boundary_cell(i_cell)]

        return list(set(ref_cells))

    def formula(self, data, id_cell, id_val=0):
        vec_rhs = self._set_rhs(data, id_cell, id_val)
        return linalg.lu_solve((self.matLU1[id_cell], self.matLU2[id_cell]), vec_rhs)

    def _set_left_mat(self):
        """
        Create the left matrix for gradient calculation.
        This method also calculate LU factorization for solving GLSQ equations.
        """

        # Calculate switching factors beta for GLSQ method.
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
                    # Vector from cell center to face center.
                    r_vec_0 = self.mesh.face_centers[i_face] - self.mesh.centers[i_cell]
                    dis0 = math.sqrt(r_vec_0[0] * r_vec_0[0] + r_vec_0[1] * r_vec_0[1] + r_vec_0[2] * r_vec_0[2])
                    if dis0 > face_dis_max:
                        face_dis_max = dis0

                    area = self.mesh.face_area[i_face]
                    if area > face_area_max:
                        face_area_max = area

                self._beta[i_cell] = min(1.0, self.mesh.volumes[i_cell] / (face_dis_max * face_area_max))

        # Calculate the left matrix and its LU factorization.
        for id_cell in range(self.mesh.n_cell):
            matA = self._set_left_mat_by_cell(id_cell)
            self.matLU1[id_cell], self.matLU2[id_cell] = linalg.lu_factor(matA)

    def _set_left_mat_by_cell(self, id_cell):
        """
        Calculate the left matrix for each cells.

        :param id_cell: Target cell.
        :return: The left matrix
        """
        mat_cell = np.zeros((3, 3), dtype=np.float64)

        if self._grad_type == 'LS' or self._grad_type == 'GLSQ':
            mat_cell = self._ls_left_mat(mat_cell, id_cell)
        elif self._grad_type == 'GG' or self._grad_type == 'GLSQ':
            mat_cell = self._gg_left_mat(mat_cell, id_cell)

        return mat_cell

    def _ls_left_mat(self, mat_cell, id_cell):
        """
        Calculate the left matrix for the weighted least square method.

        :param mat_cell: Ndarray for the matrix.
        :param id_cell: Target cell.
        :return: The matrix
        """

        nb_cells = self.mesh.cell_neighbours(id_cell)  # Get a list of the adjoined cells.
        faces = self.mesh.cell_faces[id_cell]  # List of faces on the target cell.

        for id_nb, id_face in zip(nb_cells, faces):
            # Coefficient which flips the vectors so that they always point outside of the target cell.
            flip = self.mesh.get_face_direction(id_cell, id_nb, id_face)
            vec_lr = self.mesh.vec_lr[id_face] * flip  # Vector from center of the target cell to the adjoined cell.

            weight_ls = self._get_weight_ls(id_cell, id_nb, id_face)  # Weight for LS method.
            mat_cell += vec_lr.reshape(3, 1) @ vec_lr.reshape(1, 3) * weight_ls * self._beta[id_cell]
        return mat_cell

    def _get_weight_ls(self, id_cell, id_nb, id_face):
        """
        Calculate weight for the weighted least square method.

        :param id_cell: Target cells.
        :param id_face: Face between the target and adjoined cell.
        :return: Weight.
        """

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
        """
        Calculate left matrix for Green-Gauss method.

        :param mat_cell: Ndarray for the matrix.
        :param id_cell: Target cell.
        :return: Matrix.
        """

        volume = self.mesh.volumes[id_cell]
        mat_eye = np.eye(3, dtype=np.float64)

        mat_cell += 2.0 * volume * mat_eye * (1.0 - self._beta[id_cell])
        return mat_cell

    def _set_rhs(self, data, id_cell, id_val):
        """
        Calculate right hand side vector for GLSQ method.

        :param data: Data for calculating gradient.
        :param id_cell: Target cell.
        :param id_val: Targed variables.
        :return: RHS vector.
        """

        nb_cells = self.mesh.cell_neighbours(id_cell)  # List of the adjoined cells.
        faces = self.mesh.cell_faces[id_cell]  # List of the face on the cell.

        val_vec = np.array([data[id_cell, i_val] for i_val in range(data.shape[1])])

        rhs_vec = np.zeros(3, dtype=np.float64)
        for id_nb, id_face in zip(nb_cells, faces):
            # Calculate difference of variables between cells.
            val_diff = self._get_val_diff(data, val_vec, id_nb, id_face, id_val)

            if self._grad_type == 'LS' or self._grad_type == 'GLSQ':
                rhs_vec = self._ls_rhs(rhs_vec, val_diff, id_cell, id_nb, id_face)
            elif self._grad_type == 'GG' or self._grad_type == 'GLSQ':
                rhs_vec = self._gg_rhs(rhs_vec, val_diff, id_cell, id_nb, id_face)

        return rhs_vec

    def _ls_rhs(self, rhs_vec, val_diff, id_cell, id_nb, id_face):
        """
        Calculate right hand side vector for the weighted least square part.

        :param rhs_vec: RHS vector.
        :param val_diff: Difference of variables between cells.
        :param id_cell: Target cell ID.
        :param id_nb: Adjoined cell ID.
        :param id_face: Face ID between the target and adjoined cell.
        :return: RHS vector.
        """

        flip = self.mesh.get_face_direction(id_cell, id_nb, id_face)
        vec_lr = self.mesh.vec_lr[id_face] * flip

        weight_ls = self._get_weight_ls(id_cell, id_nb, id_face)
        rhs_vec += val_diff * vec_lr * weight_ls * self._beta[id_cell]

        return rhs_vec

    def _gg_rhs(self, rhs_vec, val_diff, id_cell, id_nb, id_face):
        """
        Calculate right hand side vector for the green gauss part.

        :param rhs_vec: RHS vector.
        :param val_diff: Difference of variables between cells.
        :param id_cell: Target cell ID.
        :param id_nb: Adjoined cell ID.
        :param id_face: Face ID between the target and adjoined cell.
        :return: RHS vector.
        """

        flip = self.mesh.get_face_direction(id_cell, id_nb, id_face)
        vec_n = self.mesh.face_mat[id_face, 0] * flip
        area = self.mesh.face_area[id_face]

        rhs_vec += val_diff * area * vec_n * (1.0 - self._beta[id_cell])
        return rhs_vec

    def _get_val_diff(self, data, vec_0, id_k, id_k_face, id_val):
        """
        Calculate difference of variable value between the target and adjoined cell.

        :param data: Data vector.
        :param vec_0: Vector of value for the target cell.
        :param id_k: Cell ID for the adjoined cell.
        :param id_k_face: Face ID between the target and adjoined cell.
        :param id_val: ID of the variables.
        :return: Difference of the value.
        """

        if not self.mesh.is_boundary_face(id_k_face):  # For inner cells
            return data[id_k, id_val] - vec_0[id_val]
        else:  # For boundary cells
            val_bd = self.bd_cond.get_bd_val(vec_0, id_k_face)
            return val_bd[id_val] - vec_0[id_val]
