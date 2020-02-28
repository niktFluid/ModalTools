import argparse
import numpy as np

from Functions.Mesh import OfMesh
from Functions.FieldData import OfConData
from Functions.ModalAnalysis import RHS


def main(case_dir, time, mu, pr):
    mesh = OfMesh(case_dir, time + 'C', time + 'V', time + 'U', time + 'p')
    field = OfConData(mesh, case_dir + time, 'U', 'p', 'rho')

    rhs = RHS(mesh, field, mu, pr, is2d=True)

    data = np.ones((mesh.n_cell, 7), dtype=np.float64)
    rhs.calculate(data)

    rhs.vis_tecplot('rhs_data.dat')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate Right-Hand-Side of the equations.')

    parser.add_argument('-d', '--dir', default='CylinderFlow/', help='Target directory for flow data.')
    parser.add_argument('-t', '--time', default='5000/', help='Time directory.')
    parser.add_argument('-mu', '--viscosity', type=float, default=1.33333e-3, help='Viscosity.')
    parser.add_argument('-pr', '--Pr', type=float, default=0.7, help='Prandtl number.')
    args = parser.parse_args()

    main(args.dir, args.time, args.viscosity, args.Pr)
