import time

import itertools
import numpy as np
from scipy import sparse
from scipy.sparse import lil_matrix, csc_matrix

from Functions.Variables import Variables


# from Functions.Identity import Identity


class MatMaker:
    def __init__(self, target, n_cell, n_val=5, ave_field=None, mpi_comm=None):
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

        self.n_cell = n_cell
        self.n_val = n_val
        self.n_size_in = n_cell * n_val

        self._target = target
        self.n_return = self._target.n_return
        self.n_size_out = n_cell * self.n_return

        self._check_target()

        self.operator = None

        if ave_field is None:
            n_val_ph = 5  # rho, u, v, w, temperature
        else:
            n_val_ph = 7  # add energy and pressure
        self._ph = PlaceHolder(n_cell, n_val_ph, ave_field)

    def _check_target(self):
        if not isinstance(self._target, Variables):
            raise TypeError('"target" of "MatMaker" must be a sub class of "Variables".')

    def make_mat(self):
        if self._is_root:
            print('Calculation start.')
        t_start = time.time()

        i_start = self._rank % self._size
        i_step = self._size

        val_list = []
        # val_array = np.empty((0, 5), dtype=np.float64)
        for id_cell in range(i_start, self.n_cell, i_step):
            for val in self._calc_values(id_cell):
                val_list.append(val)
            if self._is_root:
                self._print_progress(id_cell, t_start)
        val_array = np.array(val_list)

        if self._is_root:
            t_end = time.time() - t_start
            print('Calculation done. Elapsed time: {:.0f} [sec.].'.format(t_end))

        if self._mpi:
            self._comm.barrier()
            array_list = self._comm.gather(val_array, root=0)
        else:
            array_list = val_array

        if self._is_root:
            self._set_mat(np.vstack(array_list))

    def _calc_values(self, id_cell):
        ref_cells = self._target.get_leaves(id_cell)
        for ref_cell, ref_val in itertools.product(ref_cells, range(self.n_val)):
            self._ph.set_ph(ref_cell, ref_val)
            func_val = self._target.formula(self._ph, id_cell)
            for id_val, val in enumerate(func_val):
                if val != 0.0:
                    yield np.array([id_cell, id_val, ref_cell, ref_val, val], dtype=np.float64)

    def _print_progress(self, id_cell, t_start):
        interval = 10
        prog_1 = int(100.0 * id_cell / self.n_cell)
        prog_0 = int(100.0 * (id_cell - self._size) / self.n_cell)

        if int(prog_1 / interval) > int(prog_0 / interval):
            t_elapse = time.time() - t_start
            print(str(prog_1) + ' %, Elapsed time: {:.0f}'.format(t_elapse) + ' [sec.]')

    def _set_mat(self, val_array):
        # print(val_array.shape)
        operator = lil_matrix((self.n_size_out, self.n_size_in), dtype=np.float64)
        for id_cell, id_val, ref_cell, ref_val, val in val_array:
            i_row = int(id_cell) + int(id_val) * self.n_cell
            i_col = int(ref_cell) + int(ref_val) * self.n_cell
            operator[i_row, i_col] = val
        self.operator = csc_matrix(operator)

    def save_mat(self, filename='matL.npz'):
        if self._is_root:
            sparse.save_npz(filename, self.operator)


class PlaceHolder:
    def __init__(self, n_cell, n_val, ave_field=None):
        self._gamma = 1.4
        self._gamma_2 = 1.0 / (1.4 - 1.0)
        self._gamma_3 = 1.0 / 1.4

        self.n_cell = n_cell
        self.n_val = n_val
        self.shape = (self.n_cell, self.n_val)

        self._val_vec = None

        self.i_cell = -1
        self.i_val = -1

        self._ave_field = ave_field

        if n_val > 5 and ave_field is None:
            # For calculating energy and pressure
            raise Exception('Time averaged field is necessary for computing energy and pressure fluctuations.')

    def set_ph(self, i_cell, i_val):
        self.i_cell = i_cell
        self.i_val = i_val

        self._val_vec = np.zeros(self.n_val)
        self._val_vec[i_val] = 1.0
        self._val_vec[5] = self._calc_pressure(i_cell)
        self._val_vec[6] = self._calc_energy(i_cell)

    def __getitem__(self, x):
        i_cell = x[0]
        i_val = x[1]
        if i_cell == self.i_cell:
            return self._val_vec[i_val]
        else:
            return 0.0

    def _calc_energy(self, i_cell):
        # for energy variation
        g2 = self._gamma_2
        g3 = self._gamma_3

        rho = self[i_cell, 0]
        u = self[i_cell, 1]
        v = self[i_cell, 2]
        w = self[i_cell, 3]
        # p = self[i_cell, 4]
        t = self[i_cell, 4]

        ave = self._ave_field.data
        rho_ave = ave[i_cell, 0]
        u_ave = ave[i_cell, 1]
        v_ave = ave[i_cell, 2]
        w_ave = ave[i_cell, 3]
        t_ave = ave[i_cell, 4]

        # e = g2 * p
        e = g2 * g3 * (rho * t_ave + rho_ave * t)
        e += 0.5 * rho * (u_ave * u_ave + v_ave * v_ave + w_ave * w_ave)
        e += rho_ave * (u * u_ave + v * v_ave + w * w_ave)
        return e

    def _calc_pressure(self, i_cell):
        # for calculate pressure
        g3 = self._gamma_3

        rho = self[i_cell, 0]
        t = self[i_cell, 4]

        ave = self._ave_field.data
        rho_ave = ave[i_cell, 0]
        t_ave = ave[i_cell, 4]

        return g3 * (rho * t_ave + rho_ave * t)

    # def _calc_temperature(self, i_cell):
    #     for temperature variation
        # g1 = self._gamma
        #
        # rho = self[i_cell, 0]
        # p = self[i_cell, 4]
        #
        # ave = self._ave_field.data
        # rho_ave = ave[i_cell, 0]
        # p_ave = ave[i_cell, 4]
        #
        # t = g1 * p / rho_ave
        # t += g1 * p_ave * rho / (rho_ave * rho_ave)
        # return t
