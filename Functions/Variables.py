import itertools
# import numpy as np


class Variables:
    def __init__(self, mesh, n_return=1, sub_list=None):
        self.mesh = mesh
        # self.flow_data = flow_data

        self.leave_list = None
        self._sub_list = sub_list
        # if sub_list is not None:
        #     self._check_sub_list()

        self.n_return = n_return

    def _check_sub_list(self):
        for sub in self._sub_list:
            if not issubclass(type(sub), type(self)):
                raise TypeError('"sub_list" must contains only sub classes of "Variables".')

    def get_leaves(self, id_cell):
        my_ref_cells = self._return_ref_cells(id_cell)

        sub_ref_cells = []
        if self._sub_list is not None:
            for sub_eqs, ref_cell in itertools.product(self._sub_list, my_ref_cells):
                sub_ref_cells += sub_eqs.get_leaves(ref_cell)

        self.leave_list = list(set(my_ref_cells + sub_ref_cells))

        return self.leave_list

    def formula(self, data, id_cell, **kwargs):
        raise NotImplementedError

    def _return_ref_cells(self, id_cell):
        raise NotImplementedError
