import numpy as np
import Ofpp
from collections import namedtuple


class Mesh:
    def __init__(self):
        """
        This class hold and calculate geometry data of calculation grid.

        """

        dummy = np.empty(0)
        # Grid sizes
        self.n_node = 0  # Number of node.
        self.n_face = 0  # Number of face.
        self.n_bdface = 0  # Number of boundary face.
        self.n_cell = 0  # Number of cell.

        # Connectivity data
        self.nodes = dummy  # Node positions.
        self.face_nodes = dummy  # List of nodes constructs each faces.
        self.cell_faces = dummy  # List of faces constructs each cells.
        self.cell_neighbour = dummy  # List of adjoined cells for each cells.

        # Cell geometry
        self.centers = dummy  # Position of cell centres.
        self.volumes = dummy  # Volume of cells.

        # Face geometry
        self.face_area = dummy  # Area of faces.
        self.face_centers = dummy  # Center position of faces.

        # Matrix for global -> local coordinate transformation
        self.face_mat = dummy  # Normal and tangential vectors of each faces.

        # Vectors: Owner cell -> neighbour cell
        self.vec_lr = dummy  # Vectors from cell center to center of adjoined cell.

        # BC information
        self._bd_tuple = namedtuple('Boundary', ['name', 'type', 'id', 'i_start', 'num', 'u_val', 'p_val'])
        self.boundary = [self._bd_tuple(None, None, None, None, None, None, None)]
        self.bd_nums = [0]
        self.bd_faces = [dummy]
        self.bd_cells = [dummy]

        # Initialization
        self._init_mesh()

    def _init_mesh(self):
        """
        Initialize grid geometry.
        """
        raise NotImplementedError

    def cell_neighbours(self, id_cell):
        """
        Return a list of neighbour cell IDs including boundary (ghost) cells.

        :param id_cell: ID for target cell.
        :return: Adjoined cell list.
        """
        raise NotImplementedError

    def is_boundary_face(self, id_face):
        """
        # Detect boundary face.

        :param id_face: ID for the target face.
        :return: True or False
        """
        raise NotImplementedError

    def is_boundary_cell(self, id_cell):
        """
        Detect boundary cell.

        :param id_cell: ID for the target cell.
        :return: True or False
        """
        raise NotImplementedError

    def get_face_direction(self, id_cell, nb_cell, id_face):
        """
        Return a coefficient which flips the vectors so that they always point outside of the target cell.
        This method must return '1.0' or '-1.0' according as the order of input cell ID.

        :param id_cell: ID for the target cell
        :param nb_cell: ID for the adjoined cell
        :param id_face: ID for the face between the target and adjoined cell
        :return:
        """
        raise NotImplementedError

    def get_bd_tuple(self, id_face):
        """
        Return bd_tuple that the target face is be affiliated with.

        :param id_face: Target face.
        :return: bd_tuple
        """

        for bd_data in self.boundary:
            if bd_data.i_start <= id_face < bd_data.i_start + bd_data.num:
                return bd_data
        return None

    def g2l_vel(self, val_vec, id_face):
        """
        Convert velocity vector from global coordinate to face oriented coordinate.

        :param val_vec: Variable vector.
        :param id_face: ID for the target face.
        :return: Converted vector.
        """

        face_mat = self.face_mat[id_face]

        u_vel = val_vec[1:4]
        vec_loc = np.copy(val_vec)
        vec_loc[1:4] = face_mat @ u_vel
        return vec_loc

    def l2g_vel(self, val_vec, id_face):
        """
        Convert velocity vector from face oriented to global coordinate.

        :param val_vec: Variable vector.
        :param id_face: ID for the target face.
        :return: Converted vector.
        """

        face_mat = self.face_mat[id_face].T

        u_vel = val_vec[1:4]
        vec_loc = np.copy(val_vec)
        vec_loc[1:4] = face_mat @ u_vel
        return vec_loc


class OfMesh(Mesh):
    def __init__(self, path_dir, path_centres, path_vols, path_bd_u, path_bd_p):
        """
        Mesh class for OpenFOAM grids. This class uses Ofpp to parse data.

        :param path_dir: Case directory of the OpenFOAM
        :param path_centres: File path including data of cell centers.
        :param path_vols: File path including cell volumes.
        :param path_bd_u: File path for the velocity data.
        :param path_bd_p: File path for the pressure data.
        """

        self.path_dir = path_dir
        self.path_centres = path_centres
        self.path_vols = path_vols
        self.path_bd_u = path_bd_u
        self.path_bd_p = path_bd_p

        dummy = np.empty(0)
        self.owner = dummy  # List of the owner cells.
        self.neighbour = dummy  # List of the neighbour cells.

        super(OfMesh, self).__init__()

        # Find minimum face ID for the boundary faces.
        i_bd_face_min = 1e64
        for bd_data in self.boundary:
            if bd_data.i_start <= i_bd_face_min:
                i_bd_face_min = bd_data.i_start
        self._i_start_bd = i_bd_face_min

    def cell_neighbours(self, id_cell):
        return [self.owner[x] + self.neighbour[x] - id_cell for x in self.cell_faces[id_cell]]

    def is_boundary_face(self, id_face):
        if id_face >= self._i_start_bd:
            return True
        else:
            return False

    def is_boundary_cell(self, id_cell):
        if id_cell < 0:
            return True
        else:
            return False

    def get_face_direction(self, id_cell, nb_cell, id_face):
        return -1.0 + 2.0 * float(nb_cell > id_cell or self.is_boundary_face(id_face))

    def _init_mesh(self):
        mesh = Ofpp.FoamMesh(self.path_dir)
        mesh.read_cell_centres(self.path_dir + self.path_centres)
        mesh.read_cell_volumes(self.path_dir + self.path_vols)

        self.nodes = mesh.points
        self.face_nodes = mesh.faces
        self.cell_faces = mesh.cell_faces

        self.n_node = len(self.nodes)
        self.n_face = mesh.num_face
        self.n_bdface = mesh.num_face - mesh.num_inner_face
        self.n_cell = len(self.cell_faces)

        self.owner = mesh.owner
        self.neighbour = mesh.neighbour

        self.centers = mesh.cell_centres
        self.volumes = mesh.cell_volumes

        self._set_boundary(mesh)
        self._calc_face_centers()
        self._calc_face_vec()
        self._calc_vec_lr()

    def _calc_face_vec(self):
        """
        Calculate face geometries.
        """

        self.face_area = np.empty(self.n_face, dtype=np.float64)
        self.face_mat = np.empty((self.n_face, 3, 3), dtype=np.float64)

        for i_face in range(self.n_face):
            face_nodes = np.array(self.face_nodes[i_face])

            vec_a = self.nodes[face_nodes[2]] - self.nodes[face_nodes[0]]
            vec_b = self.nodes[face_nodes[3]] - self.nodes[face_nodes[1]]
            vec_c = np.cross(vec_a, vec_b)  # Calculate face normal vector.

            abs_norm = np.linalg.norm(vec_c)
            self.face_area[i_face] = 0.5 * abs_norm  # Area of face.

            vec_n = vec_c / abs_norm  # Face normal vector
            flip = self._check_normal_vec(i_face, vec_n)
            vec_n *= flip

            vec_t1 = vec_a / np.linalg.norm(vec_a)  # Face tangential vector 1
            vec_t2 = np.cross(vec_n, vec_t1)  # Face tangential vector 2.

            self.face_mat[i_face, 0] = vec_n
            self.face_mat[i_face, 1] = vec_t1
            self.face_mat[i_face, 2] = vec_t2

    def _check_normal_vec(self, i_face, vec_n):
        """
        Make sure the face normal vector points to neighbour cells.

        :param i_face: ID for the target face.
        :param vec_n: Normal vector.
        :return: Flip factor.
        """

        id_o = self.owner[i_face]
        id_n = self.neighbour[i_face]

        if id_n >= 0:  # For inner faces
            vec_lr = self.centers[id_n] - self.centers[id_o]
        else:  # For boundary faces
            vec_lr = self.face_centers[i_face] - self.centers[id_o]

        if np.dot(vec_lr, vec_n) < 0.0:
            print('Flip normal vector!')
            return -1.0
        else:
            return 1.0

    def _calc_face_centers(self):
        """
        Calculate center position for each faces.￿
        """

        points = self.nodes[self.face_nodes]
        self.face_centers = np.mean(points, axis=1)

    def _calc_vec_lr(self):
        """
        Calculate vector from cell center to cell center of adjoined cell.
        """

        self.vec_lr = np.zeros((self.n_face, 3), dtype=np.float64)
        centers = self.centers

        for i_face in range(self.n_face):
            id_o = self.owner[i_face]
            id_n = self.neighbour[i_face]

            if id_n >= 0:  # For inner faces
                self.vec_lr[i_face] = centers[id_n] - centers[id_o]
            else:  # For boundary faces
                face_vec_n = self.face_mat[i_face, 0]
                face_centers = self.face_centers

                dist_fc = np.linalg.norm(face_centers[i_face] - centers[id_o])
                self.vec_lr[i_face] = 2.0 * dist_fc * face_vec_n

    def _set_boundary(self, mesh):
        """
        Set boundary tuples.

        :param mesh: Ofpp object.
        """
        bd_u = Ofpp.parse_boundary_field(self.path_dir + self.path_bd_u)
        bd_p = Ofpp.parse_boundary_field(self.path_dir + self.path_bd_p)

        def make_bd_tuple(bd_key):
            source_dic = mesh.boundary[bd_key]

            if b'value' in bd_u[bd_key]:
                u_val = bd_u[bd_key][b'value']
            else:
                u_val = None

            if b'value' in bd_p[bd_key]:
                p_val = bd_p[bd_key][b'value']
            else:
                p_val = None

            bd_tuple = self._bd_tuple(bd_key.decode(),
                                      source_dic.type.decode(),
                                      source_dic.id,
                                      source_dic.start,
                                      source_dic.num,
                                      u_val, p_val)
            return bd_tuple

        bd = [make_bd_tuple(key) for key in mesh.boundary.keys()]
        owner = np.array(mesh.owner)

        self.boundary = bd
        self.bd_nums = [bd[i].num for i in range(len(bd))]
        self.bd_faces = [np.arange(bd[i].i_start, bd[i].i_start + bd[i].num) for i in range(len(bd))]
        self.bd_cells = [owner[self.bd_faces[i]] for i in range(len(bd))]
