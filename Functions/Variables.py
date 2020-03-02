import itertools
# import numpy as np


class Variables:
    def __init__(self, mesh, n_return=1, sub_list=None):
        """
        Abstract class for equations.

        :param mesh: Mesh class.
        :param n_return: Number of variables that is returned for calculations on each cells.
        :param sub_list: List of subclass of Variables that is used for calculating.
        """

        self.mesh = mesh

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
        """
        Return list of cells which are referenced to calculate function values.

        :param id_cell: Target cell ID.
        :return: List of cells.
        """

        my_ref_cells = self._return_ref_cells(id_cell)

        sub_ref_cells = []
        if self._sub_list is not None:
            for sub_eqs, ref_cell in itertools.product(self._sub_list, my_ref_cells):
                sub_ref_cells += sub_eqs.get_leaves(ref_cell)

        self.leave_list = list(set(my_ref_cells + sub_ref_cells))

        return self.leave_list

    def formula(self, data, id_cell, **kwargs):
        """
        Specific formulation for equations as a function of data and cell ID.

        :param data: Input vector of the functions.
        :param id_cell: Target cell ID.
        :param kwargs: Options
        :return: Calculated values.
        """
        raise NotImplementedError

    def _return_ref_cells(self, id_cell):
        """
        Return list of cells which are 'directory' referenced from this class for calculation.
        For this class, the list contains the target cell and the adjoined cells.

        :param id_cell: Target cell.
        :return: List of cell.
        """

        raise NotImplementedError
