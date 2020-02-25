from Functions.Mesh import OfMesh
from Functions.FieldData import OfConData
from Functions.ModalAnalysis import RHS


def main(case_dir, time, mu, pr):
    mesh = OfMesh(case_dir, time + 'C', time + 'V', time + 'U', time + 'p')
    field = OfConData(mesh, case_dir + time, 'U', 'p', 'rho')

    rhs = RHS(mesh, field, mu, pr, is2d=True)
    rhs.vis_tecplot('rhs_data.dat')


if __name__ == '__main__':
    main('CylinderFlow/', '1/', 1.33333e-3, 0.7)
