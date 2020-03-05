import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import os
import errno

import argparse
import configparser

import numpy as np

from Functions.Mesh import OfMesh
from Functions.FieldData import OfData

from Functions.ModalAnalysis import ResolventMode as Resolvent


def main(param_file='Parameter.dat', profile='Default', length=None, velocity=None):
    if (length is None and velocity is not None) or (length is not None and velocity is None):
        print('CAUTION: If you want to scale the angular frequency to Strouhal number, '
              'you have to specify both characteristic length and velocity. '
              'The graph will be plotted with using the angular frequencies.')

    if not os.path.exists(param_file):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), param_file)

    profiles = configparser.ConfigParser()
    profiles.read(param_file, encoding='utf-8')
    params = profiles[profile]

    print('Calculation start.\nParameters:')
    for key, value in params.items():
        print(key + ' = ' + value)

    # Set up
    mode = params['Mode']
    if not mode == 'RandomizedResolvent' or mode == 'Resolvent':
        raise Exception('You made a deviation!')

    case_dir = params['CaseDir']
    time_dir = params['TimeDir']
    operator = params['Operator']
    save_name = params['SaveName']
    k = int(params['ModeNum'])

    PlotResolventGain(case_dir, time_dir, operator, save_name, k=k, length=length, velocity=velocity)


def PlotResolventGain(case_dir, time, operator_name, save_name, k=3, length=None, velocity=None):
    flag_st_conv = length is not None and velocity is not None

    mesh = OfMesh(case_dir, time + 'C', time + 'V', time + 'U', time + 'p')
    ave_field = OfData(mesh, case_dir + time, 'UMean', 'pMean', 'rhoMean')
    resolvent_mode = Resolvent(mesh, ave_field, operator_name, k=k, mode='Both')

    gains_filename = save_name + '/gains.dat'
    file_num, omega, alpha, gains = get_gains(gains_filename)

    print('Calculating gains...')
    gain_grads = get_gain_grad(resolvent_mode, save_name, file_num, gains, k)
    omega_plot, gain_plot = get_plot_data(omega, gains[:, 0:k], gain_grads)
    print('Done.')

    figure, axis = make_figure()

    if flag_st_conv:
        axis.set_xlabel(r'$St = \omega L / 2 \pi U_\infty$')
    else:
        axis.set_xlabel(r'$\omega$')

    axis.set_ylabel(r'Gain: $\sigma$')

    for gain_array in gain_plot:
        if flag_st_conv:
            axis.plot(omega_plot * length / (2.0 * np.pi * velocity), gain_array)
        else:
            axis.plot(omega_plot, gain_array)

    # for gain in gains[:, 0:k].T:
        # axis.scatter(omega * L / (2.0 * np.pi * U), gain)

    if not flag_st_conv:
        file_name = '/omega-gains.pdf'
    else:
        file_name = '/st-gains.pdf'

    figure.savefig(save_name + file_name, bbox_inches='tight')


def get_gains(filename):
    data_array = np.loadtxt(filename)

    file_num = data_array[:, 0]
    omega = data_array[:, 1]
    alpha = data_array[:, 2]
    gains = data_array[:, 3:]

    return file_num, omega, alpha, gains


def get_gain_grad(mode_data, save_name, file_num_array, gain_array, k):

    grad_list = []
    for file_num, gain in zip(file_num_array, gain_array):
        filename = save_name + '/modes_{:0=5}.pickle'.format(int(file_num))

        mode_data.load_data(filename)
        _, _, _, response, forcing = mode_data.vec_data

        w_r = mode_data.qo @ response[:, 0:k]
        w_f = mode_data.qo @ forcing[:, 0:k]

        grad_gain = -gain[0:k]**2 * np.imag(np.diag(w_f.conj().T @ w_r))

        grad_list.append(grad_gain)

    return np.vstack(grad_list)


def get_plot_data(omega_array, gain_array, grad_array, n_section=25):
    omega_list = []
    gain_list = []
    for i_section in range(len(omega_array) - 1):
        omega0 = omega_array[i_section]
        omega1 = omega_array[i_section + 1]

        gain0 = gain_array[i_section]
        gain1 = gain_array[i_section + 1]

        grad0 = grad_array[i_section]
        grad1 = grad_array[i_section + 1]

        matL = np.array([
            [omega0**3, omega0**2, omega0, 1.0],
            [omega1**3, omega1**2, omega1, 1.0],
            [3.0*omega0**2, 2.0*omega0, 1.0, 0.0],
            [3.0*omega1**2, 2.0*omega1, 1.0, 0.0]
        ])

        rhs = np.array([gain0, gain1, grad0, grad1])

        coef_array = np.linalg.solve(matL, rhs)

        omega_plot = np.linspace(omega0, omega1, n_section)
        gain_plot = np.array([np.polyval(coef, omega_plot) for coef in coef_array.T])

        omega_list.append(omega_plot)
        gain_list.append(gain_plot)

    return np.hstack(omega_list), np.hstack(gain_list)


def make_figure():
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams['mathtext.fontset'] = 'cm'  # Computer Modern
    plt.rcParams['font.family'] = 'cmr10'  # Computer Modern
    plt.rcParams['font.size'] = 14
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    figure_1 = plt.figure(figsize=(5, 3.5))
    ax1 = figure_1.add_subplot(111)

    ax1.set_yscale('log')

    return figure_1, ax1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Modal analysis.')

    parser.add_argument('-f', '--filename', default='Parameter.dat', help='Parameter file for the calculation.')
    parser.add_argument('-p', '--profile', default='Default', help='Profile for the parameters.')
    parser.add_argument('-L', '--length', default=None,
                        help='Characteristic length. If this parameter is specified, '
                             'frequency will be normalized as Strouhal number.', type=float)
    parser.add_argument('-U', '--velocity', default=None,
                        help='Characteristic velocity. If this parameter is specified, '
                             'frequency will be normalized as Strouhal number.', type=float)
    args = parser.parse_args()

    main(param_file=args.filename, profile=args.profile, length=args.length, velocity=args.velocity)
