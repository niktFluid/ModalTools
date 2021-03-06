import os
import time
from itertools import product
import pickle

import numpy as np
import scipy as sp
from scipy.sparse import linalg
from scipy import sparse

from Functions.FieldData import FieldData
from Functions.LinearizedNS import LNS
from Functions.Gradient import Gradient


class ModalData(FieldData):
    def __init__(self, mesh, operator_name=None, n_val=5, k=10, mpi_comm=None, **kwargs):
        """
        Abstract class for modal analysis.

        :param mesh: Mesh class.
        :param operator_name: File name of the liner operator.
        :param n_val: Number of variables. Default: 5.
        :param k: Number of modes for each calculations.
        :param mpi_comm: MPI communicator.
        :param kwargs: Options.
        """

        # kwargs for set up the operator.
        self._k = k
        self._n_q = n_val

        super(ModalData, self).__init__(mesh, n_val=self._data_num(), data_list=self._data_name_list())

        if operator_name is not None:
            self.operator = self._set_operator(sparse.load_npz(operator_name), **kwargs)
        else:
            self.operator = None

        self.vec_data = None

        self._arpack_options = {
            'k': self._k,
            'sigma': None,
            'which': 'LM',
            'tol': 1.0e-8,
        }

        self._comm = mpi_comm
        if mpi_comm is not None:
            self._size = mpi_comm.Get_size()
            self._rank = mpi_comm.Get_rank()
            self._mpi = True
        else:
            self._size = 1
            self._rank = 0
            self._mpi = False
        self._is_root = self._rank == 0

    def _init_field(self, *args, **kwargs):
        self.data = np.empty((self.n_cell, self._data_num()), dtype=np.float64)

    def _data_num(self):
        """
        Calculate number of data stored self.data.

        :return: Number of data.
        """

        raise NotImplementedError

    def _data_name_list(self):
        """
        Generate list of data name for visualization.

        :return: Name list.
        """

        raise NotImplementedError

    def _set_operator(self, operator, **kwargs):
        raise NotImplementedError

    def solve(self, **kwargs):
        """
        Solve problems.

        :param kwargs:
        :return:
        """

        self.vec_data = self._calculate(**kwargs)
        self._set_data(self.vec_data)  # Set self.data and self._vec_data for the visualization.

    def save_data(self, filename='modalData.pickle'):
        with open(filename, 'wb') as file_obj:
            pickle.dump(self.vec_data, file_obj)

    def load_data(self, filename='modalData.pickle'):
        with open(filename, 'rb') as file_obj:
            self.vec_data = pickle.load(file_obj)
        self._set_data(self.vec_data)

    def _calculate(self, **kwargs):
        raise NotImplementedError

    def _set_data(self, data):
        """
        Set data for self.data.

        :param data: Moe data.
        """

        raise NotImplementedError


class LinearStabilityMode(ModalData):
    def __init__(self, mesh, operator, n_val=5, k=5, n_grid=10, add_pres=False, **kwargs):
        """
        Class for global stability analysis.
        This class does not support MPI computations.

        :param mesh: Mesh class.
        :param operator_name: File name of the liner operator.
        :param n_val: Number of variables. Default: 5.
        :param k: Number of modes for each calculations.
        :param n_grid: Number of grid points used for shift-inverted mode.
        :param add_pres: Flag for visualized data of pressure mode.
        """

        self._n_grid = n_grid
        self._add_pres = add_pres

        if add_pres:
            n_val += 1

        super(LinearStabilityMode, self).__init__(mesh, operator, n_val, k, **kwargs)
        self._arpack_options.update(**kwargs)

        if self._mpi:
            raise Exception('This class does not support MPI computations.')

    def _data_num(self):
        return self._n_q * self._k * self._n_grid

    def _data_name_list(self):
        data_list_base = ['rho', 'u', 'v', 'w', 'T']

        if self._add_pres:
            data_list_base.append('p')

        data_list = []
        for i_mode in range(self._k * self._n_grid):
            data_list += ['mode{:0=4}_'.format(i_mode) + x for x in data_list_base]
        return data_list

    def solve(self, grid_list):
        if len(grid_list) != self._n_grid:
            raise Exception

        eig_list = []
        vec_list = []
        for i_grid, (sigma_real, sigma_imag) in enumerate(grid_list):
            sigma = sigma_real + 1.0j * sigma_imag
            self._arpack_options.update({'sigma': sigma})
            print('Sigma =', sigma)

            eigs, vecs = self._calculate()
            vec_list.append(vecs)
            eig_list.append(eigs)

        self.vec_data = (np.hstack(eig_list), np.hstack(vec_list))
        self._set_data(self.vec_data)

    def _set_operator(self, operator, **kwargs):
        return operator

    def _calculate(self):
        return linalg.eigs(self.operator, **self._arpack_options)

    def _set_data(self, data):
        _, vecs = data

        for i_mode, vec in enumerate(vecs.T):
            i_start = self._n_q * i_mode
            i_end = self._n_q * (i_mode + 1)

            if self._add_pres:
                w_vec = self._calc_pres_mode(vec.reshape((self.n_cell, int(vec.size / self.n_cell)), order='F'))
            else:
                w_vec = vec.reshape((self.n_cell, self._n_q), order='F')

            self.data[:, i_start:i_end] = np.real(w_vec)

    def _calc_pres_mode(self, data):
        rho = data[:, 0]
        t = data[:, 4]
        p = rho * t / 1.4
        return np.hstack((data, p.reshape((self.n_cell, 1))))

    def save_data(self, filename='modalData.pickle'):
        save_name, _ = os.path.splitext(filename)

        eigs, _ = self.vec_data
        # noinspection PyTypeChecker
        np.savetxt(save_name + '_eigs.txt', np.vstack((np.real(eigs), np.imag(eigs))).T)

        with open(filename, 'wb') as file_obj:
            pickle.dump(self.vec_data, file_obj)


class ResolventMode(ModalData):
    def __init__(self, mesh, ave_field, operator, n_val=5, k=6, mode=None, mpi_comm=None, **kwargs):
        """
        Class for resolvent analysis.

        :param mesh: Mesh class.
        :param ave_field: Time averaged flow field. Used to calculate Tu norm.
        :param operator: File name of the operator.
        :param n_val: Number of variables.
        :param k: Number of mode for each calculations.
        :param mode: Mode to be calculated. 'Forcing' or 'Response' or 'Both'.
        :param mpi_comm: MPI communicator.
        :param kwargs: Options.
        """

        self._ave_field = ave_field

        self._mode = mode  # 'F' for the forcing mode or 'R' for the response mode. 'None' will get both.
        self._mode_f = mode == 'Both' or mode == 'Forcing'
        self._mode_r = mode == 'Both' or mode == 'Response'

        # Matrix for the numerical quadrature.
        self.qi = None
        self.qo = None
        self._resolvent = None

        super(ResolventMode, self).__init__(mesh, operator, n_val, k, mpi_comm, **kwargs)
        self._arpack_options.update(sigma=0.0, which='LM', **kwargs)

    def solve(self, grid_list, save_dir):
        os.makedirs(save_dir, exist_ok=True)
        gain_file = save_dir + '/gains.dat'

        f = open(gain_file, 'w', encoding='utf-8')
        f.close()

        if self._is_root:
            print('Start resolvent operations.')

        t_start = time.time()
        for i_grid, grid in self._grid_queue(grid_list):  # Calculate resolvent modes for each omega and alpha values.
            if grid is not None:
                omega, alpha = grid

                resolvent = self._make_resolvent(omega, alpha)
                gain, mode_r, mode_f = self._calculate(resolvent)

                self.vec_data = (omega, alpha, gain, mode_r, mode_f)
                self._set_data(self.vec_data)

                save_name = save_dir + '/modes_{:0=5}'.format(i_grid)
                self.save_data(save_name + '.pickle')
                self.vis_tecplot(save_name + '.dat')

                w_list = [i_grid, omega, alpha] + list(gain)
            else:
                w_list = None

            if self._mpi:
                gains_list = self._comm.gather(w_list, root=0)  # Gather list of gain from threads to root.
            else:
                gains_list = [w_list]

            if self._is_root:
                # Write gains list to text file.
                if not all(item is None for item in gains_list):
                    data_list = [item for item in gains_list if item is not None]
                    data_list.sort()

                    with open(gain_file, mode='a') as f_obj:
                        for w_list in data_list:
                            f_obj.writelines(list(map(lambda x: str(x) + ' ', w_list)) + ['\n'])

                            omega = w_list[1]
                            alpha = w_list[2]
                            gain = w_list[3:]
                            print('Omega = {:.6f}'.format(omega) + ', Alpha = {:.6f}'.format(alpha))
                            print('Gains: ', gain)

                    t_current = time.time() - t_start
                    print('Elapsed time: {:.1f} [sec.].'.format(t_current))

                    stop_operation = False
                else:
                    print('Done resolvent operations.')
                    stop_operation = True
            else:
                stop_operation = None

            if self._mpi:
                stop_operation = self._comm.bcast(stop_operation, root=0)

            # Continue operation or not.
            if stop_operation:
                break
            else:
                continue

    def _grid_queue(self, grid_list):
        """
        Generate grid queue for parallel computations.

        :param grid_list: List of grid contains omega and alpha.
        :return: Grid data.
        """

        i_step = self._size
        i_ind = self._rank % self._size

        while i_ind < len(grid_list):
            yield i_ind, grid_list[i_ind]
            i_ind += i_step

        while True:
            yield None, None

    def _data_num(self):
        if self._mode == 'Both':
            return self._n_q * self._k * 2
        else:
            return self._n_q * self._k

    def _data_name_list(self):
        data_list_base = ['rho', 'u', 'v', 'w', 'T']
        data_list = []
        for i_mode in range(self._k):
            if self._mode_f:
                data_list += ['forcing{:0=4}_'.format(i_mode) + x for x in data_list_base]
            if self._mode_r:
                data_list += ['response{:0=4}_'.format(i_mode) + x for x in data_list_base]
        return data_list

    def _set_operator(self, operator, **kwargs):
        self.qi, self.qo = self._get_norm_quadrature()
        return operator

    def _make_resolvent(self, omega, alpha):
        eye = sparse.eye(self.operator.shape[0], dtype=np.complex128, format='csc')
        omegaI = 1.0j * (omega + 1.0j * alpha) * eye
        return self.qo * (-omegaI - self.operator) * self.qi

    def _calculate(self, resolvent):
        svs = None
        if self._mode_f:
            matF = resolvent * resolvent.H
            svs, mode_f = linalg.eigsh(matF, **self._arpack_options)
            print('Eigenvalues for forcing: ', svs)
        else:
            mode_f = None

        if self._mode_r:
            matR = resolvent.H * resolvent
            svs, mode_r = linalg.eigsh(matR, **self._arpack_options)
            print('Eigenvalues for response: ', svs)
        else:
            mode_r = None

        print('Singular values: ', np.sqrt(svs))
        # print('Gains: ', 1.0 / np.sqrt(svs))
        return 1.0 / np.sqrt(svs), self.qi @ mode_r, self.qi @ mode_f

    def _set_data(self, data):
        """
        Set data arrays for resolvent analysis.

        :param data: Mode vectors.
        """

        _, _, _, r_vecs, f_vecs = data

        coef_ind_1 = 1 + int(self._mode == 'Both')
        coef_ind_2 = self._n_q * int(self._mode == 'Both')
        for i_mode in range(self._k):
            if self._mode_f:
                f_vec = f_vecs[:, i_mode]
                fw_vec = f_vec.reshape((self.n_cell, self._n_q), order='F')
                i_start = coef_ind_1 * self._n_q * i_mode
                i_end = i_start + self._n_q
                self.data[:, i_start:i_end] = np.real(fw_vec)

            if self._mode_r:
                r_vec = r_vecs[:, i_mode]
                rw_vec = r_vec.reshape((self.n_cell, self._n_q), order='F')
                i_start = coef_ind_1 * self._n_q * i_mode + coef_ind_2
                i_end = i_start + self._n_q
                self.data[:, i_start:i_end] = np.real(rw_vec)

    def _get_norm_quadrature(self):
        """
        Calculate numerical quadrature by Chu's energy norm.

        :return:
        """

        ave_data = self._ave_field.data
        rho_data = ave_data[:, 0]
        t_data = ave_data[:, 4]

        gamma = 1.4  # heat ratio
        r_gas = 1.0 / 1.4  # Non-dimensionalized gas constant.

        vols = self.mesh.volumes / np.linalg.norm(self.mesh.volumes)  # Weights by cell volumes.

        diag_rho = vols * r_gas * t_data / rho_data
        diag_u = vols * rho_data
        diag_t = vols * r_gas * rho_data / ((gamma - 1) * t_data)
        diags = np.hstack((diag_rho, diag_u, diag_u, diag_u, diag_t))

        qi = sparse.diags(1.0 / np.square(diags), format='csc')
        qo = sparse.diags(np.square(diags), format='csc')

        return qi, qo


class RandomizedResolventMode(ResolventMode):
    def __init__(self, mesh, bd_cond, ave_field, operator, n_val=5, k=6, mode='Both', **kwargs):
        """
        Class for calculating Randomized resolvent analysis.
        Jean Hélder Marques Ribeiro, Chi-An Yeh, Kunihiko Taira, 2020, https://arxiv.org/abs/1902.01458.

        :param mesh: Mesh class.
        :param ave_field: Time averaged flow field. Used to calculate Tu norm.
        :param operator: File name of the operator.
        :param n_val: Number of variables.
        :param k: Number of mode for each calculations.
        :param mode: Mode to be calculated. 'Forcing' or 'Response' or 'Both'.
        :param mpi_comm: MPI communicator.
        :param kwargs: Options.
        """

        super(RandomizedResolventMode, self).__init__(mesh, ave_field, operator, n_val, k, mode, **kwargs)

        self._scaling = self._get_scaling_factor(ave_field.data, bd_cond)

    def _get_scaling_factor(self, ave_data, bd_cond):
        grad = Gradient(self.mesh, bd_cond)

        grad_vel = np.zeros((self.n_cell, 3, 3), dtype=np.float64)
        for i_cell, i_val in product(range(self.n_cell), [1, 2, 3]):
            grad_vel[i_cell, i_val-1] = grad.formula(ave_data, i_cell, i_val)

        # noinspection PyTypeChecker
        phi = np.tile(np.linalg.norm(grad_vel, axis=(1, 2)), 5)

        return sparse.diags(phi, format='csc')

    def _calculate(self, resolvent):
        resolvent_lu = linalg.splu(resolvent)

        m = self.n_cell * 5  # = resolvent.shape[0]
        k = self._k  # Number of mode
        matO = self._scaling @ np.random.normal(0.0, 0.1, (m, k))
        matY = resolvent_lu.solve(matO)

        matQ, _ = sp.linalg.qr(matY, mode='economic')
        matB = resolvent_lu.solve(matQ, trans='H')

        _, _, V = sp.linalg.svd(matB.T.conj(), full_matrices=False)
        matUS = resolvent_lu.solve(V.T.conj())

        U, Sigma, Vapp = sp.linalg.svd(matUS, full_matrices=False)
        V = V.T.conj() @ Vapp.T.conj()

        return Sigma, self.qi @ U, self.qi @ V


class RHS(FieldData):
    def __init__(self, mesh, field, mu, pr, is2d=False):
        """
        Calculate right hand side of LNS operator.

        :param mesh: Mesh class.
        :param field: Field data.
        :param mu: viscosity.
        :param pr: Prandtl number.
        :param is2d: Flag for 2D computation.
        """

        super(RHS, self).__init__(mesh, n_val=5, data_list=['Rho', 'U', 'V', 'W', 'T'])

        self.rhs_ns = LNS(mesh, field, mu, pr, is2d)

    def _init_field(self, *args, **kwargs):
        self.data = np.empty((self.n_cell, self.n_val), dtype=np.float64)

    def calculate(self, data):
        for i_cell in range(self.n_cell):
            self.data[i_cell] = self.rhs_ns.formula(data, i_cell)
