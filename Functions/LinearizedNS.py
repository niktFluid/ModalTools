from itertools import product
import numpy as np
import math

from Functions.Variables import Variables
from Functions.Gradient import Gradient


class LNS(Variables):
    def __init__(self, mesh, bd_cond, ave_field, mu, pr, grad_type='GLSQ'):
        """
        This class provides formulation of Linearized Navier-Stokes equations.
        The equations are discretized by finite volume method (FVM).
        A class method named 'formula' calculates right hand side of the equation
        with specific input array and time-averaged data.

        :param mesh: Class for grid information and related methods.
        :param bd_cond: Boundary condition corresponding solvers.
        :param ave_field: Time averaged field.
        :param mu: Nondimensionalized viscosity that is equal to [Mach number]/[Reynolds number].
        :param pr: Prandtl number.
        :param grad_type: Type of gradient. Default: 'GLSQ'
        """

        self.gamma = 1.4
        self.gamma_1 = 1.0 / (1.4 - 1.0)
        self.gamma_inv = 1.0 / 1.4
        self.mu = mu
        self.pr = pr
        self.coef_heat_flux = mu / ((self.gamma - 1) * pr)

        self.n_cell = mesh.n_cell
        self.n_val = ave_field.n_val

        self.mesh = mesh
        self.bd_cond = bd_cond

        self._grad = Gradient(mesh, bd_cond, grad_type=grad_type)
        sub_list = [self._grad]

        super(LNS, self).__init__(mesh, n_return=5, sub_list=sub_list)

        self._data = None
        self.ave_data = ave_field.data
        self._grad_ave = self._grad_ave_field()

        self._grad_refs = self._GradData(self._grad)

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

    def formula(self, data, id_cell, **kwargs):
        """
        Calculate right hand side of the equations.

        :param data: Input vector of the equation.
        :param id_cell: Target cell.
        :return: Right hand side value of the equation.
        """

        self._data = data
        # Calculate gradient of the variables for the target and adjoined cells.
        self._grad_refs.set_grad(data, self._return_ref_cells(id_cell))

        nb_cells = self.mesh.cell_neighbours(id_cell)  # Get list of the adjoined cells.
        faces = self.mesh.cell_faces[id_cell]  # Get list of the face of the target cell.

        rhs_vec = np.zeros(5, dtype=np.float64)  # Right hand side vector.

        # Calculate surface integration of the flux.
        for nb_cell, nb_face in zip(nb_cells, faces):
            area = self.mesh.face_area[nb_face]  # Area of the surface.

            # Coefficient which flips the surface normal vectors so that they always point outside of the target cell.
            flip = self.mesh.get_face_direction(id_cell, nb_cell, nb_face)

            rhs_vec -= self._calc_inviscid_flux(id_cell, nb_cell, nb_face) * area * flip  # Inviscid flux.
            rhs_vec += self._calc_viscous_flux(id_cell, nb_cell, nb_face) * area * flip  # Viscous flux.
        return self._conv2prime(rhs_vec, id_cell) / self.mesh.volumes[id_cell]

    class _GradData:
        def __init__(self, grad):
            """
            Calculate the gradient of variables at the target cells.
            :param grad: Gradient module.
            """
            self.n_data = 7
            self.data = None
            self.ref_cells = None
            self.grad_refs = None

            self.grad_list = [1, 2, 3, 4]  # Target vaiables. (u-vel, v-vel, w-vel, temperature)
            self.grad = grad

        def set_grad(self, data, ref_cells):
            self.data = data
            self.ref_cells = ref_cells

            def grad(i_cell):
                grad_array = np.zeros((self.n_data, 3), dtype=np.float64)
                for i_val in self.grad_list:
                    grad_array[i_val] = self.grad.formula(self.data, i_cell, i_val)
                return grad_array

            self.grad_refs = [grad(id_cell) for id_cell in ref_cells]

        def __getitem__(self, i_cell):
            return self.grad_refs[self.ref_cells.index(i_cell)]

    def _grad_ave_field(self):
        grad_ave = np.zeros((self.n_cell, self.n_val, 3), dtype=np.float64)
        for i_cell, i_val in product(range(self.n_cell), range(self.n_val)):
            grad_ave[i_cell, i_val] = self._grad.formula(self.ave_data, i_cell, i_val)
        return grad_ave

    def _calc_inviscid_flux(self, id_cell, nb_cell, nb_face):
        """
        Calculate the inviscid flux of the equations.
        This method uses KEP scheme like discretization.

        "Jameson, A., Journal of Scientific Computing,
        Vol. 34, No. 2, 2008, pp. 188–208. doi:10.1007/s10915-007-9172-6."

        :param id_cell: Target cells.
        :param nb_cell: Adjoined cell.
        :param nb_face: Face between the target and adjoined cell.
        :return: Inviscid flux on the face.
        """

        def flux(vec_face, ave_face):
            """
            Inviscid flux.

            :param vec_face: Basic variables vector at the cell face.
            :param ave_face: Time averaged variables vector.
            :return: Inviscid flux.
            """

            rho = vec_face[0]
            u = vec_face[1]
            v = vec_face[2]
            w = vec_face[3]
            p = vec_face[5]
            e = vec_face[6]

            rho_ave = ave_face[0]
            u_ave = ave_face[1]
            v_ave = ave_face[2]
            w_ave = ave_face[3]
            p_ave = ave_face[5]
            e_ave = ave_face[6]

            f = np.empty(5, dtype=np.float64)
            f[0] = rho * u_ave + rho_ave * u
            f[1] = 2.0 * rho_ave * u_ave * u + rho * u_ave * u_ave + p
            f[2] = rho_ave * u_ave * v + rho_ave * u * v_ave + rho * u_ave * v_ave
            f[3] = rho_ave * u_ave * w + rho_ave * u * w_ave + rho * u_ave * w_ave
            f[4] = (e_ave + p_ave) * u + (e + p) * u_ave
            return f

        # The value of the variables at the face is calculated as the average value of the target and adjoined cell.
        # The velocity vector are converted from global coordinate to face oriented coordinate.
        vec_0, vec_nb = self._get_cell_vals(self._data, id_cell, nb_cell, nb_face)
        vec_f = self.mesh.g2l_vel(0.5 * (vec_0 + vec_nb), nb_face)

        ave_0, ave_nb = self._get_cell_vals(self.ave_data, id_cell, nb_cell, nb_face)
        ave_f = self.mesh.g2l_vel(0.5 * (ave_0 + ave_nb), nb_face)

        # Calculate flux and reconvert the flux vector.
        return self.mesh.l2g_vel(flux(vec_f, ave_f), nb_face)

    def _calc_viscous_flux(self, id_cell, nb_cell, nb_face):
        """
        Calculate the viscous flux of the equations.

        :param id_cell: Target cells.
        :param nb_cell: Adjoined cell.
        :param nb_face: Face between the target and adjoined cell.
        :return: Viscous flux on the face.
        """

        flux = np.zeros(5, dtype=np.float64)
        face_normal_vec = self.mesh.face_mat[nb_face, 0]  # Normal vector of the face.

        vec_a, vec_b = self._get_cell_vals(self._data, id_cell, nb_cell, nb_face)
        vec_f = 0.5 * (vec_a + vec_b)  # Estimated basic value at the face.
        u_vel = vec_f[1:4]

        g_face = self._get_face_grad(self._data, self._grad_refs, id_cell, nb_cell, nb_face)  # Gradient at the face.
        tau = self._get_stress_tensor(g_face)  # Get stress tensor at the face.

        ave_a, ave_b = self._get_cell_vals(self.ave_data, id_cell, nb_cell, nb_face)
        ave_f = 0.5 * (ave_a + ave_b)
        u_ave = ave_f[1:4]

        g_face_ave = self._get_face_grad(self.ave_data, self._grad_ave, id_cell, nb_cell, nb_face)
        tau_ave = self._get_stress_tensor(g_face_ave)

        flux[1:4] = tau @ face_normal_vec  # Viscous flux for the kinetic part.

        # Viscous flux for the energy part.
        energy_flux = tau @ u_ave + tau_ave @ u_vel + self.coef_heat_flux * g_face[4, :]
        flux[4] = energy_flux @ face_normal_vec

        return flux

    def _conv2prime(self, vec_conv, id_cell):
        """
        Convert the variation of conservative variables to the variation of prime variables.
        [rho', rho-u', rho-v', rho-w', e'] -> [rho', u', v', w', T']
        """
        rho = vec_conv[0]
        ru = vec_conv[1]
        rv = vec_conv[2]
        rw = vec_conv[3]
        e = vec_conv[4]

        rho_ave = self.ave_data[id_cell, 0]
        u_ave = self.ave_data[id_cell, 1]
        v_ave = self.ave_data[id_cell, 2]
        w_ave = self.ave_data[id_cell, 3]
        e_ave = self.ave_data[id_cell, 6]
        ra_inv = 1.0 / rho_ave

        vec_pr = np.empty(5, dtype=np.float64)
        vec_pr[0] = rho
        vec_pr[1] = ra_inv * ru - u_ave * ra_inv * rho
        vec_pr[2] = ra_inv * rv - v_ave * ra_inv * rho
        vec_pr[3] = ra_inv * rw - w_ave * ra_inv * rho

        # Concert to temperature
        u = vec_pr[1]
        v = vec_pr[2]
        w = vec_pr[3]
        vec_pr[4] = 1.4 * 0.4 * (e * ra_inv - e_ave * rho * ra_inv * ra_inv)
        vec_pr[4] += - 1.4 * 0.4 * (u * u_ave + v * v_ave + w * w_ave)

        return vec_pr

    def _get_cell_vals(self, data, id_cell, nb_cell, nb_face):
        """
        Get basic values for the target and adjoined cell.
        This method also calculate boundary values if it is needed.

        :param data: Input vector.
        :param id_cell: Target cell.
        :param nb_cell: Adjoined cell.
        :param nb_face: The face between the target and adjoined cell.
        :return: Cell values.
        """

        def get_vals(i_cell):
            val_vec = np.empty(self.n_val, dtype=np.float64)
            for i_val in range(self.n_val):
                val_vec[i_val] = data[i_cell, i_val]
            return val_vec

        val_vec_0 = get_vals(id_cell)  # Values for the target cell.

        # Calculate values for the adjoined cell.
        if not self.mesh.is_boundary_face(nb_face):  # For inner cells
            val_vec_nb = get_vals(nb_cell)
        else:  # For boundary cells
            val_vec_nb = self.bd_cond.get_bd_val(val_vec_0, nb_face)

        return val_vec_0, val_vec_nb

    def _get_face_vals(self, data, grad_data, id_cell, nb_cell, nb_face):
        """
        Estimate value of the basic variables at the cell face.
        MUSCL interpolation.

        :param data: Input vector.
        :param grad_data: Gradient matrix at the target cell￿.
        :param id_cell: Target cell.
        :param nb_cell: Adjoined cell.
        :param nb_face: Face between the target and adjoined cell.
        :return: Interpolated values for the both side of the face.
        """
        val_0, val_nb = self._get_cell_vals(data, id_cell, nb_cell, nb_face)

        def reconstruct(val_vec, ind):
            """
            Reconstruct distribution of the variables inside the cell

            :param val_vec: Variables vector.
            :param ind: Target cell.
            :return: Interpolated values.
            """

            grad = grad_data[ind]
            r_vec = self.mesh.face_centers[nb_face] - self.mesh.centers[ind]  # Vector from cell center to face center.
            return val_vec + grad @ r_vec

        if not self.mesh.is_boundary_face(nb_face):  # For inner cells.
            val_vec_0 = reconstruct(val_0, id_cell)
            val_vec_nb = reconstruct(val_nb, nb_cell)
        else:  # For boundary cells.
            val_vec_0, val_vec_nb = val_0, val_nb

        return val_vec_0, val_vec_nb

    def _get_face_grad(self, data, grad_data, id_cell, nb_cell, nb_face):
        """
        Calculate the gradient at the cell faces.

        :param data: Input vector.
        :param grad_data: Gradient data.
        :param id_cell: Target cell.
        :param nb_cell: Adjoined cell.
        :param nb_face: Face between the target and adjoined cell.
        :return: Gradient at the face.
        """

        grad_id = grad_data[id_cell]  # Gradient of the variables at the face.
        vol_id = self.mesh.volumes[id_cell]  # Volume of the target cell

        if not self.mesh.is_boundary_face(nb_face):  # For inner faces
            grad_nb = grad_data[nb_cell]  # Gradient of the adjoined cell.
            vol_nb = self.mesh.volumes[nb_cell]
            grad_face = (grad_id * vol_nb + grad_nb * vol_id) / (vol_id + vol_nb)  # Estimated gradient at the face.
        else:  # For boundary faces.
            grad_face = grad_id

        # For prevent even-odd instability.
        # Nishikawa, H., 40th Fluid Dynamics Conference and Exhibit, 2010, AIAA 2010-5093. See 'Appendix C'.

        flip = self.mesh.get_face_direction(id_cell, nb_cell, nb_face)

        vec_lr = self.mesh.vec_lr[nb_face]  # Vector from the center of target cell to adjoined cell.
        inv_lr = 1.0 / math.sqrt(vec_lr[0]*vec_lr[0] + vec_lr[1]*vec_lr[1] + vec_lr[2]*vec_lr[2])

        vec_a, vec_b = self._get_cell_vals(data, id_cell, nb_cell, nb_face)
        coef = (grad_face @ vec_lr - (vec_b - vec_a) * flip) * inv_lr  # Coefficient for correction.

        return grad_face - coef.reshape(self.n_val, 1) @ vec_lr.reshape(1, 3) * inv_lr

    def _get_stress_tensor(self, grad):
        """
        Calculate stress tensor.￿

        :param grad: Gradient of the data.
        :return: Stress tensor.
        """
        tensor = np.empty((3, 3), dtype=np.float64)
        mu = self.mu

        dudx = grad[1, 0]
        dudy = grad[1, 1]
        dudz = grad[1, 2]
        dvdx = grad[2, 0]
        dvdy = grad[2, 1]
        dvdz = grad[2, 2]
        dwdx = grad[3, 0]
        dwdy = grad[3, 1]
        dwdz = grad[3, 2]

        div_u = (dudx + dvdy + dwdz) / 3.0
        tensor[0, 0] = 2.0 * mu * (dudx - div_u)
        tensor[1, 1] = 2.0 * mu * (dvdy - div_u)
        tensor[2, 2] = 2.0 * mu * (dwdz - div_u)

        tensor[0, 1] = mu * (dudy + dvdx)
        tensor[0, 2] = mu * (dudz + dwdx)
        tensor[1, 2] = mu * (dwdy + dvdz)

        tensor[1, 0] = tensor[0, 1]
        tensor[2, 0] = tensor[0, 2]
        tensor[2, 1] = tensor[1, 2]

        return tensor
