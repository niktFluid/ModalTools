from Functions.Variables import Variables


class Identity(Variables):
    # Identity.py function
    def __init__(self, mesh):
        super(Identity, self).__init__(mesh)

    def _return_ref_cells(self, id_cell):
        return [id_cell]

    def formula(self, ph, id_cell, id_val=0, **kwargs):
        return ph[id_cell, id_val]
