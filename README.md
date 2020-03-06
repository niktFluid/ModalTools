# ModalTools
This codes constructs a linear operator from flow field data of compressible fluid flows.
(日本語版の説明は[こちら](README_JP.md)．)

## Description
ModalTools is a tool for constructing and analyzing a linearized Navier-Stokes (LNS) operator in a modal manner. The project contains necessary codes for estimating the Right-Hand-Side of the LNS equation, as well as a matrix generator that uses the so-called zero-one trick for making a sparse matrix representation of the LNS operator. The operator can be used for modal analysis including linear stability analysis and resolvent analysis (see [here](https://doi.org/10.2514/1.J056060)). We also provide a low Reynolds number compressible cylinder flow data, obtained via [OpenFOAM](https://www.openfoam.com/), for demonstration purposes.


## Demo
This section shows a demonstration of ModalTools for a two-dimensional compressible cylinder flow at <img src="https://latex.codecogs.com/gif.latex?Re=150"> and <img src="https://latex.codecogs.com/gif.latex?M_\infty=0.2">. The simulation is conducted using OpenFOAM. The flow field shows vortex shedding behind the cylinder at a frequency of <img src="https://latex.codecogs.com/gif.latex?St=0.178">. The LNS operator is constructed about the time-averaged flow.

![Figure1-2](https://user-images.githubusercontent.com/47338366/76111253-cc8fdc00-5f94-11ea-880f-b9a0953f5330.png)

### Stability Analysis
Solving an eigenvalue problem of the operator provides a spectrum that represents the characteristics of the dynamic response of a linear system. Eigenvectors corresponding to each eigenvalue, depict structures of the fluid flow. In our formulation, the real part of each eigenvalue indicates the growth rate of the structure, while the imaginary part correspond to its frequency. The following figure shows the results of stability analysis for the cylinder flow. We can find some eigenvalues on the unstable plane (<img src="https://latex.codecogs.com/gif.latex?{\rm&space;Re}(\lambda)&space;>&space;0">). The imaginary part of the higher frequency one (<img src="https://latex.codecogs.com/gif.latex?\lambda&space;=&space;7.25&space;\times&space;10^{-3}&space;&plus;&space;0.208i">) gives a frequency of <img src="https://latex.codecogs.com/gif.latex?St=0.165"> which is close to the shedding frequency obtained from the simulation. Indeed, the eigenvector that corresponds to this eigenvalue shows structures like vortex shedding behind the cylinder and can be seen on the right side of the following figure. Some unstable eigenvalues that are just on the real axis may come from insufficient time-averaging for the data (try to visualize these modes).

![Figure2-2](https://user-images.githubusercontent.com/47338366/76111278-d580ad80-5f94-11ea-966a-903eed397f2a.png)

### Resolvent Analysis
The other feature implemented in ModalTools involves the pseudospectrum of a linear operator. Resolvent analysis can be understood in the context of an analysis of a harmonic input forcing its output structure, through the linear system. Singular Value Decomposition (SVD) of the operator delivers pairs of orthogonal basis vectors correspond to the forcing and response modes of the operator, with the singular value giving the energy amplitude ratio (gain) between them. The following picture shows the distribution of the three largest gains over frequency, and the resolvent modes corresponding to the shedding frequency. The leading gain shows its highest value at the shedding frequency since the phenomenon comes from fluid dynamic unsteadiness. The response modes show a vortex shedding like structure while the forcing is looks like some sort of velocity fluctuations in the vicinity of the cylinder walls. We also provide an implementation of [randomized resolvent analysis](https://arxiv.org/abs/1902.01458) that makes your assessments faster.

![Figure3-2](https://user-images.githubusercontent.com/47338366/76111281-d7e30780-5f94-11ea-9f57-79398097fc97.png)


## Requirements
* Python 3.*
* NumPy 
* SciPy
* [Ofpp](https://github.com/dayigu/ofpp) - Used to parse OpenFOAM data.
* [mpi4py](https://mpi4py.readthedocs.io/en/stable/) - Optional but recommended.

## Usage
The author will provide instructions about how to use this codes when using OpenFOAM to solve the compressible Navier-Stokes equations. Please note that the current version only supports a hexahedra cells mesh. ModalTools assumes the variables used in the simulation and codes are normalized as following.

<img src="https://latex.codecogs.com/gif.latex?x&space;=&space;\frac{\widetilde{x}}{L},&space;\:&space;y&space;=&space;\frac{\widetilde{y}}{L},&space;\:&space;z&space;=&space;\frac{\widetilde{z}}{L}">, 

<img src="https://latex.codecogs.com/gif.latex?\rho&space;=&space;\frac{\widetilde{\rho}}{\rho_\infty},&space;\:&space;u&space;=&space;\frac{\widetilde{u}}{a_\infty},&space;\:&space;T&space;=&space;\frac{\widetilde{T}}{T_\infty}">.

Here, x, y, z are position, <img src="https://latex.codecogs.com/gif.latex?\rho"> is density, <img src="https://latex.codecogs.com/gif.latex?u"> is velocity, <img src="https://latex.codecogs.com/gif.latex?T"> is temperature and <img src="https://latex.codecogs.com/gif.latex?a"> is the speed of sound. The subscript <img src="https://latex.codecogs.com/gif.latex?\infty"> indicates far-field variables. The easiest way to get this normalized field is to set constant pressure specific heat <img src="https://latex.codecogs.com/gif.latex?c_p"> and molar mass <img src="https://latex.codecogs.com/gif.latex?m"> as <img src="https://latex.codecogs.com/gif.latex?c_p=2.5">, <img src="https://latex.codecogs.com/gif.latex?m=11640.3"> as you can find in [thermophysical settings of the CylinderFlow case](CylinderFlow/constant/thermophysicalProperties). 

Please make sure that the simulation case has time averaged data of each variables (`CylinderFlow/1000/*Mean`) as in [this directory](CylinderFlow/1000). ModalTools also needs data of center position (`CylinderFlow/1000/C`) and the volume (`CylinderFlow/1000/V`) of each computational cell. You can get this data by using the following commands in the case directory.
```
postProcess -funcs writeCellCentres 
postProcess -funcs writeCellVolumes
```

The demonstration case ([CylinderFlow](CylinderFlow)) has a script named [Allrun](CylinderFlow/Allrun) for running simulation and post processing. So all you have to do for obtaining demonstration data is just to execute [Allrun](CylinderFlow/Allrun).

We should mention that in our trial, the Ofpp module which is used to parse OpenFOAM data worked well only if the `writeFormat` option in [controlDict](CylinderFlow/system/controlDict) is set as `ascii`. If you have a problem regarding the parsing process in Ofpp, please try to convert your data into binary format by using `foamFormatConvert` with set `writeFormat` as `binary` and reconvert them into ascii format.


### Generating Operator
You can use [GenerateOperator.py](GenerateOperator.py) (for serial computation) or [GenerateOperatorMPI.py](GenerateOperatorMPI.py) (for parallel computation) for generating the operator. 
```
python3 GenerateOperator.py -f [Parameter file name] -p [Parameter name]
mpiexec -np [Number of thread] python3 GenerateOperatorMPI.py -f [Parameter file name] -p [Parameter name]
```

You can find an example of parameter in [Parameters.dat](Parameters.dat). `CaseDir` and` TimeDir` set the case directory and the directory containing the flow field data at the time represented by the directory name, respectively.
```
[GenerateOperator]
CaseDir = CylinderFlow/
TimeDir = 1000/
Operator = CylinderFlow/matL.npz
Viscosity = 1.333333e-3
PrandtlNumber = 0.7
```


### Stability or Resolvent Analysis
[CalcMode.py](CalcMode.py) or [CalcModeMPI.py](CalcModeMPI.py) provide basic functions for stability or resolvent analysis. The code outputs TecPlot ASCII format file (.dat) for visualization and binary files (.pickle) for data saving. You will also get text files that have data of eigenvalues or gains, and a figure of resolvent gain plots which are cubic interpolated by using gains value and its gradients ([reference1](https://web.stanford.edu/group/ctr/Summer/SP14/08_Transition_and_turbulence/11_fosas.pdf), [reference2](https://spiral.imperial.ac.uk/handle/10044/1/72876)).
```
python3 CalcMode.py -f [Parameter file name] -p [Parameter name]
mpiexec -np [Number of thread] python3 CalcModeMPI.py -f [Parameter file name] -p [Parameter name]
```

The ModalTools uses shift-invert mode for eigenvalue problem and sweeps complex plane by parameters you set. Parameters regarding the sigma value (`Sigma*` in the following example) define grid points that are used for shift-invert mode. `ModeNum` indicates how many eigenvalues are calculated for each grid point. In the following example, you will get 104 (= 4 * 26) modes after the computation. `Which` controls priority order to find eigenvectors and eigenvalues (see [here](https://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html) for details).
```
[Stability]
Mode = Stability
CaseDir = CylinderFlow/
TimeDir = 1000/
Operator = CylinderFlow/matL.npz
SaveName = CylinderFlow/stability
ModeNum = 4
SigmaRealStart = 0.01
SigmaRealEnd = 0.01
SigmaRealNum = 1
SigmaImagStart = 0.0
SigmaImagEnd = 1.0
SigmaImagNum = 26
Which = LM
```

The author recommends using [the randomized method](https://arxiv.org/abs/1902.01458) for resolvent analysis since this method is much faster than normal SVD. Parameters regarding `Omega` define omega values which is used in the sweeping operation on frequency. `Alpha` parameters are involving with 'time discounting' (see [this paper](https://doi.org/10.1017/jfm.2019.163) for the detail). `ResolventMode` controls which mode is saved. You can use `Both`, `Response`, `Forcing` for the parameter.
```
[Resolvent]
Mode = RandomizedResolvent
CaseDir = CylinderFlow/
TimeDir = 1000/
Operator = CylinderFlow/matL.npz
SaveName = CylinderFlow/resolvent
ModeNum = 3
OmegaStart = 0.0
OmegaEnd = 1.0
OmegaNum = 101
AlphaStart = 0.0
AlphaEnd = 0.0
AlphaNum = 1
ResolventMode = Both
```


## Author
* **Yoimi Kojima** - [niktFluid](https://github.com/niktFluid)


## Licence

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE) file for details


## Acknowledgments

The author gratefully acknowledges members of [Taira laboratory](http://www.seas.ucla.edu/fluidflow/index.html) for insightful comments and discussions on the implementations.
