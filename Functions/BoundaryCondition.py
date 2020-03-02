import numpy as np


class BoundaryCondition:
    def __init__(self, mesh, bd_dict, is2d=False, key_2d=None):
        """
        Abstract class for boundary conditions.

        :param mesh: Mesh class.
        :param bd_dict: Dictionary contains boundary name and functions.
        :param is2d: Flag to apply 2D mode.
        :param key_2d: Specify name of boundary name which is used for 2D boundaries.
        """

        self.mesh = mesh
        self.bd_dict = bd_dict

        self.is2d = is2d
        self.key_2d = key_2d

        if is2d and key_2d is None:
            raise Exception('BD key for 2D boundary must be set if is2d is True.')

    def get_bd_val(self, val_vec, id_face):
        bd_data = self.mesh.get_bd_tuple(id_face)
        bd_type = bd_data.type

        if self.is2d and bd_type == self.key_2d:
            return self._bd_2d(val_vec, id_face)
        else:
            # Convert velocity vector from global coordinate to face oriented coordinate.
            vec_loc = self.mesh.g2l_vel(val_vec, id_face)
            bd_func = self.bd_dict[bd_type]
            bd_vec_loc = bd_func(vec_loc, id_face)

            return self.mesh.l2g_vel(bd_vec_loc, id_face)

    def _bd_2d(self, val_vec, id_face):
        raise NotImplementedError


class OFBC(BoundaryCondition):
    def __init__(self, mesh, is2d=False):
        """
        This class provides boundary conditions for OpenFOAM data.

        :param mesh: Mesh class.
        :param is2d: Flag to apply 2D mode.
        """

        bd_dict = {
            'wall': self._bd_wall,
            'patch': self._bd_0th_ex,
            'empty': self._bd_symmetry
                   }

        super(OFBC, self).__init__(mesh, bd_dict, is2d, key_2d='empty')

    def _bd_2d(self, val_vec, _):
        val_vec[3] *= -1.0
        return val_vec

    @staticmethod
    def _bd_0th_ex(val_vec, _):
        return val_vec

    @staticmethod
    def _bd_symmetry(val_vec, _):
        val_vec[1] *= -1.0
        return val_vec

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
        return val_vec
