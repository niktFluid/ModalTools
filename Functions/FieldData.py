import numpy as np
import Ofpp


class FieldData:
    def __init__(self, mesh, n_val=5, data_list=None):
        if data_list is None:
            self.data_list = ['rho', 'u', 'v', 'w', 'T']
        else:
            self.data_list = data_list

        self.mesh = mesh
        self.n_val = n_val
        self.n_cell = mesh.n_cell
        self.data = np.empty(0)

        self._init_field()

    def _init_field(self, *args, **kwargs):
        raise NotImplementedError

    def vis_tecplot(self, file_name='Default.dat'):
        header = self._make_tec_header()
        with open(file_name, mode='w') as file_obj:
            def write_data(w_data):
                for w_list in self._make_write_list(w_data):
                    file_obj.writelines(list(map(lambda x: str(x) + ' ', w_list)) + ['\n'])

            file_obj.write(header)
            for ind in range(3):
                write_data(self.mesh.nodes[:, ind])
            for ind in range(self.n_val):
                write_data(self.data[:, ind])
            for connectivity in self._make_connectivity():
                c_list = list(map(lambda x: str(x) + ' ', connectivity))
                file_obj.writelines(c_list + ['\n'])

    def _make_tec_header(self):
        header = 'Title = "Data" \n'
        header += 'Variables = "x"\n"y"\n"z"\n' + '\n'.join(['"' + x + '"' for x in self.data_list]) + '\n'
        header += 'Zone T = "Fluid area" \n'
        # header += 'StrandID=1, SolutionTime=' + str(self.sol_time) + ' \n'
        header += 'Nodes=' + str(self.mesh.n_node) + ' \n'
        header += 'Elements=' + str(self.mesh.n_cell) + ' \n'
        header += 'DATAPACKING=BLOCK \nZONETYPE=FEBrick \n'
        header += 'VARLOCATION=([4-' + str(self.n_val + 3) + ']=CELLCENTERED) \n'
        return header

    def _make_connectivity(self):
        # This subroutine should be modified for clarify.
        mesh = self.mesh
        for faces in mesh.cell_faces:
            nodes = set(sum([mesh.face_nodes[x] for x in faces], []))

            id_0 = faces[0]  # Face ID of a base face.
            fn_0 = set(mesh.face_nodes[id_0])  # Node IDs of the base face.
            fn_o = nodes - fn_0  # Node IDs of a face which is on the opposite side of the cell.
            id_o = self._find_opposite(faces, fn_o)  # Find the face ID of the opposite face.

            face_nodes_0 = mesh.face_nodes[id_0]  # Node list of the base face.
            # Arrange the node list of the opposite face to meet connectivity.
            face_nodes_op = self._find_connection(id_0, id_o, faces)

            # Node number starts from '1' in the TecPlot format.
            yield list(map(lambda x: x + 1, face_nodes_0)) + list(map(lambda x: x + 1, face_nodes_op))

    def _find_opposite(self, face_list, fn_set):
        for i_fn in face_list:
            i_fn_set = set(self.mesh.face_nodes[i_fn])
            if fn_set == i_fn_set:
                return i_fn

    def _find_connection(self, id_base, id_opposite, faces):
        nodes_0 = self.mesh.face_nodes[id_base]
        side_0 = set(nodes_0[0:2])

        i_nb = -1
        for i_face in set(faces) - {id_base, id_opposite}:
            node_set = set(self.mesh.face_nodes[i_face])
            if len(side_0 & node_set) == 2:
                i_nb = i_face
        nodes_nb = self.mesh.face_nodes[i_nb]

        nb_0 = nodes_nb.index(nodes_0[0])
        i_corner = nodes_nb[nb_0 - 2]
        last_node = set(nodes_nb) - side_0 - {i_corner}
        i_next_0 = last_node.pop()
        side_opposite = [i_next_0, i_corner]

        nodes_op = self.mesh.face_nodes[id_opposite]
        for i_rotation in range(3):
            for i_shift in range(len(nodes_op)):
                shifted_list = nodes_op[-i_shift:] + nodes_op[:-i_shift]
                if shifted_list[0:2] == side_opposite:
                    return shifted_list
            nodes_op.reverse()

    @staticmethod
    def _make_write_list(data, n_line=8):
        n_data = data.size
        w_data = data.reshape(n_data)

        i_end = 0
        n_row = int(n_data / n_line)
        for i_row in range(n_row):
            i_start = n_line * i_row
            i_end = n_line * (i_row + 1)
            row_data = w_data[i_start:i_end]
            yield row_data

        last_data = w_data[i_end:]
        if not len(last_data) == 0:
            yield last_data


class FlowData(FieldData):
    def __init__(self, mesh, add_e=False, add_pres=False, data_list=None):
        super(FlowData, self).__init__(mesh, data_list=data_list)

        if add_pres:
            self._add_pressure()
        if add_e:
            self._add_energy()
        # if add_temp:
        #     self._add_temperature()

    def _init_field(self, *args, **kwargs):
        raise NotImplementedError

    def _add_pressure(self):
        self.n_val += 1
        self.data_list.append('p')

        rho = self.data[:, 0]
        t = self.data[:, 4]
        p = rho * t / 1.4

        self.data = np.hstack((self.data, p.reshape((self.n_cell, 1))))

    def _add_energy(self):
        self.n_val += 1
        self.data_list.append('e')

        rho = self.data[:, 0]
        u = self.data[:, 1]
        v = self.data[:, 2]
        w = self.data[:, 3]
        p = self.data[:, 4]

        g2 = 1.0 / (1.4 - 1.0)
        e = g2 * p + 0.5 * rho * (u * u + v * v + w * w)

        self.data = np.hstack((self.data, e.reshape((self.n_cell, 1))))

    def _add_temperature(self):
        self.n_val += 1
        self.data_list.append('T')

        rho = self.data[:, 0]
        p = self.data[:, 4]
        temp = 1.4 * p / rho

        self.data = np.hstack((self.data, temp.reshape((self.n_cell, 1))))


class OfData(FlowData):
    def __init__(self, mesh, path_dir, path_u, path_p, path_rho, add_e=False, add_pres=False):
        self.path_dir = path_dir
        self.path_u = path_u
        self.path_p = path_p
        self.path_rho = path_rho

        super(OfData, self).__init__(mesh, add_e=add_e, add_pres=add_pres)

    def _init_field(self):
        self.update_from_file()

    def update_from_file(self, path_u=None, path_p=None, path_rho=None):
        if path_u is None:
            path_u = self.path_u

        if path_p is None:
            path_p = self.path_p

        if path_rho is None:
            path_rho = self.path_rho

        u_data = Ofpp.parse_internal_field(self.path_dir + path_u)
        p_data = Ofpp.parse_internal_field(self.path_dir + path_p)
        rho_data = Ofpp.parse_internal_field(self.path_dir + path_rho)

        # Calculate temperature
        t_data = 1.4 * p_data / rho_data

        self.n_cell = u_data.shape[0]
        # self.data = np.hstack((rho_data[:, np.newaxis], u_data, p_data[:, np.newaxis]))
        self.data = np.hstack((rho_data[:, np.newaxis], u_data, t_data[:, np.newaxis]))


class OfConData(FlowData):
    def __init__(self, mesh, path_dir, path_u, path_p, path_rho):
        self.path_dir = path_dir
        self.path_u = path_u
        self.path_p = path_p
        self.path_rho = path_rho
        d_list = ['Rho', 'RhoU', 'RhoV', 'RhoW', 'E']

        super(OfConData, self).__init__(mesh, add_e=False, add_pres=False, data_list=d_list)

    def _init_field(self):
        self.update_from_file()

    def update_from_file(self, path_u=None, path_p=None, path_rho=None):
        if path_u is None:
            path_u = self.path_u

        if path_p is None:
            path_p = self.path_p

        if path_rho is None:
            path_rho = self.path_rho

        u_data = Ofpp.parse_internal_field(self.path_dir + path_u)
        p_data = Ofpp.parse_internal_field(self.path_dir + path_p)
        rho_data = Ofpp.parse_internal_field(self.path_dir + path_rho)

        # Calculate conservative variables
        u = u_data[:, 0]
        v = u_data[:, 1]
        w = u_data[:, 2]

        ru = rho_data * u
        rv = rho_data * v
        rw = rho_data * w
        e = p_data / (1.4 - 1.0) + 0.5 * rho_data * (u * u + v * v + w * w)

        self.n_cell = u_data.shape[0]
        # self.data = np.hstack((rho_data[:, np.newaxis], u_data, p_data[:, np.newaxis]))
        self.data = np.vstack((rho_data, ru, rv, rw, e)).T
