import os
import errno

import argparse
import configparser

from Functions.Mesh import OfMesh
from Functions.FieldData import OfData

from Functions.MatMaker import MatMaker
from Functions.LinearizedNS import LNS


def main(param_file='Parameter.dat', profile='Default'):
    if not os.path.exists(param_file):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), param_file)

    if profile == 'Default':
        profile = 'GenerateOperatorDefault'

    profiles = configparser.ConfigParser()
    profiles.read(param_file, encoding='utf-8')
    params = profiles[profile]

    case_dir = params['CaseDir']
    time_dir = params['TimeDir']
    operator_name = params['Operator']
    mu = float(params['Viscosity'])
    pr = float(params['PrandtlNumber'])

    MakeOperator(case_dir, time_dir, operator_name, mu, pr)


def MakeOperator(case_dir, time, filename, mu, pr):
    mesh = OfMesh(case_dir, time + 'C', time + 'V', time + 'U', time + 'p')
    ave_field = OfData(mesh, case_dir + time, 'UMean', 'pMean', 'rhoMean', add_e=True, add_pres=True)

    linear_ns = LNS(mesh, ave_field, mu=mu, pr=pr, is2d=True)  # viscosity and Prandtl number

    mat_maker = MatMaker(linear_ns, mesh.n_cell, ave_field=ave_field)
    mat_maker.make_mat()
    mat_maker.save_mat(filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make a liner operator from OpenFOAM results. MPI version.')

    parser.add_argument('-f', '--filename', default='Parameter.dat', help='Parameter file for the calculation.')
    parser.add_argument('-p', '--profile', default='Default', help='Profile for the parameters.')
    args = parser.parse_args()

    main(args.filename, args.profile)
