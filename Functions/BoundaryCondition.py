import numpy as np


class BoundaryCondition:
    def __init__(self, mesh, is2d=False):
        self.mesh = mesh
        self.is2d = is2d

        self.wall = 'wall'
        self.patch = 'patch'
        self.empty = 'empty'

    def get_bd_val(self, val_vec, id_face):
        bd_data = self.mesh.get_bd_tuple(id_face)
        bd_type = bd_data.type

        if self.is2d and bd_type == self.empty:
            val_vec[3] *= -1.0
            return val_vec
        else:
            vec_loc = self.mesh.g2l_vel(val_vec, id_face)
            bd_func = self._select_bd_func(bd_type)
            bd_vec_loc = bd_func(vec_loc, id_face)

            return self.mesh.l2g_vel(bd_vec_loc, id_face)

    def _select_bd_func(self, bd_type):
        if bd_type == self.wall:
            return self._bd_wall
        elif bd_type == self.empty:
            return self._bd_symmetry
        elif bd_type == self.patch:
            return self._bd_0th_ex
        else:
            raise Exception

    @staticmethod
    def _bd_0th_ex(val_vec, _):
        return val_vec

    @staticmethod
    def _bd_symmetry(val_vec, _):
        val_vec[1] *= -1.0
        return val_vec

    # @staticmethod
    def _bd_wall(self, val_vec, id_face):
        bd_data = self.mesh.get_bd_tuple(id_face)
        vel_wall_g = bd_data.u_val

        if vel_wall_g is None or not np.any(vel_wall_g):
            vel_wall = 0.0
        else:
            vel = np.zeros(5, dtype=np.float64)
            vel[1:4] = vel_wall_g
            vel_loc = self.mesh.g2l_vel(vel, id_face)
            vel_wall = vel_loc[1:4]

        val_vec[1:4] = 2.0 * vel_wall - val_vec[1:4]
        # val_vec[1] *= -1.0
        # val_vec[2] *= -1.0
        # val_vec[3] *= -1.0
        return val_vec
