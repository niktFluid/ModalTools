# ModalTools
This codes construct a linear operator from a field data of compressible fluid flows.


## Description
ModalTools is a tool for constructing and analyzing a linearized Navier-Stokes (LNS) operator in a modal manner. The project contains necessary codes for estimating the Right-Hand-Side of LNS equation and a matrix generator that uses so-called zero-one trick for making a sparse matrix representing LNS operator. The operator can be used for modal analysis including linear stability analysis and resolvent analysis (See [here](https://doi.org/10.2514/1.J056060)). We also provide a low Reynolds number compressible cylinder flow data via [OpenFOAM](https://www.openfoam.com/) for the demonstration.


## Demo
This section shows a demonstration of the tool for compressible cylinder flow field at <img src="https://latex.codecogs.com/gif.latex?Re=150">, <img src="https://latex.codecogs.com/gif.latex?M_\infty=0.2">. The simulation is conducted using OpenFOAM. The flow field shows vortex shedding behind the cylinder at frequency of <img src="https://latex.codecogs.com/gif.latex?St=0.178">. 

![Figure1-1](https://user-images.githubusercontent.com/47338366/75612946-26982980-5add-11ea-97c5-cc24cd3bc953.png)

### Stability Analysis
Solving an eigenvalue problem of the operator provides a spectrum that represents the characteristic of the dynamic response of a linear system. Eigenvectors corresponding to each eigenvalue, depicts structures of fluid flow. In our formulation, the real parts of eigenvalues indicate growth rates of the structure, while the imaginary parts correspond to its frequency. The following figure shows the results of stability analysis for the cylinder flow. We can find some eigenvalue on the unstable plane (<img src="https://latex.codecogs.com/gif.latex?{\rm&space;Re}(\lambda)&space;>&space;0">). The imaginary part of the higher frequency one (<img src="https://latex.codecogs.com/gif.latex?\lambda&space;=&space;7.25&space;\times&space;10^{-3}&space;&plus;&space;0.208i">) gives a frequency of <img src="https://latex.codecogs.com/gif.latex?St=0.165"> which is close to the shedding frequency obtained from the simulation. Indeed, the eigenvector corresponds to this eigenvalue apparently shows structures like vortex shedding behind the cylinder as it can be seen on the right side of the following figure. Some unstable eigenvalues that are just on the real axis may come from insufficient time-averaging for the data (try to visualize these modes).

![Figure2-1](https://user-images.githubusercontent.com/47338366/75613193-b2ab5080-5adf-11ea-94dc-2ba803502a39.png)

### Resolvent Analysis
The other feature implemented on the codes involves with a pseudospectrum of a liner operator. The resolvent analysis can be understood in the context that an analysis of harmonic input forcing and output structure of the linear systems. Singular Value Decomposition (SVD) of the operator delivers pairs of orthogonal basis correspond to the forcing, response modes of the operator and energy amplitude ratio (gain) between them. The following picture shows distribution of three largest gains over frequency and resolvent modes correspond to the shedding frequency. The leading gain shows its highest value at the shedding frequency since the phenomenon comes from fluid dynamic unsteadiness. The response modes shows vortex shedding like structure while the forcing is looked like some sort of velocity fluctuations in the vicinity of the cylinder walls. We also provides an implementation of [randomized resolvent analysis](https://arxiv.org/abs/1902.01458) that makes your assessments faster.

![Figure3-1](https://user-images.githubusercontent.com/47338366/75613466-e340b980-5ae2-11ea-960e-13f0339524d6.png)


## Requirements
* Python 3.*
* NumPy 
* SciPy
* [Ofpp](https://github.com/dayigu/ofpp) - Used to parse OpenFOAM data.
* [mpi4py](https://mpi4py.readthedocs.io/en/stable/) - Optional but recommended.

## Usage
We will provide instructions about how to use this codes when we use OpenFOAM and solving compressible Navier-Stokes equations. Please note that the current version only support hexahedra cells mesh. We assumes the variables used in the simulation and codes are normalized as following equations.

<img src="https://latex.codecogs.com/gif.latex?x&space;=&space;\frac{\widetilde{x}}{L},&space;\:&space;y&space;=&space;\frac{\widetilde{y}}{L},&space;\:&space;z&space;=&space;\frac{\widetilde{z}}{L}">
<img src="https://latex.codecogs.com/gif.latex?\rho&space;=&space;\frac{\widetilde{\rho}}{\rho_\infty},&space;\:&space;u&space;=&space;\frac{\widetilde{u}}{a_\infty},&space;\:&space;T&space;=&space;\frac{\widetilde{T}}{T_\infty}">

Here, x, y, z are position, <img src="https://latex.codecogs.com/gif.latex?\rho"> is density, <img src="https://latex.codecogs.com/gif.latex?u"> is velocity, <img src="https://latex.codecogs.com/gif.latex?T"> is temperature and <img src="https://latex.codecogs.com/gif.latex?a"> is the speed of sound. Subscript <img src="https://latex.codecogs.com/gif.latex?\infty"> indicates far-field variables. The easiest way to get this normalized field is to set constant pressure specific heat <img src="https://latex.codecogs.com/gif.latex?c_p"> and molar mass <img src="https://latex.codecogs.com/gif.latex?M"> as <img src="https://latex.codecogs.com/gif.latex?c_p=2.5">, <img src="https://latex.codecogs.com/gif.latex?M=11640.3"> as you can find in [thermophysical settings of CylinderFlow case](CylinderFlow/constant/thermophysicalProperties). 

Please make sure that the simulation case has time averaged data of each variables (`CylinderFlow/1000/*Mean`) as like [this directory](CylinderFlow/1000). The codes also need data of center  position (`CylinderFlow/1000/C`) and volume (`CylinderFlow/1000/V`) of each computational cells. You can get these data by using following commands in the case directory.
```
postProcess -funcs writeCellCentres 
postProcess -funcs writeCellVolumes
```

We guess we should mention that Ofpp module which is used to parse OpenFOAM data works well only `writeFormat` option in [controlDict](CylinderFlow/system/controlDict) is set as `ascii` in our trial. If you have a problem regarding the parsing process in Ofpp, please try to convert your data into binary format by using `foamFormatConvert` with set `writeFormat` as `binary` and reconvert them into ascii format.


### Generating Operator
You can use [GenerateOperator.py](GenerateOperator.py) (for serial computation) or [GenerateOperatorMPI.py](GenerateOperatorMPI.py) (for parallel computation) for generating the operator. 
```
python3 GenerateOperator.py -f [Parameter file name] -p [Parameter name]
mpiexec -np [Number of thread] python3 GenerateOperatorMPI.py -f [Parameter file name] -p [Parameter name]
```

You can find an example of parameter in [Parameters.dat](Parameters.dat).
```
[GenerateOperator]
CaseDir = CylinderFlow/
TimeDir = 1000/
Operator = CylinderFlow/matL.npz
Viscosity = 1.333333e-3
PrandtlNumber = 0.7
```


### Stability or Resolvent Analysis
[CalcMode.py](CalcMode.py) or [CalcModeMPI.py](CalcModeMPI.py) provide basic functions for stability or resolvent analysis. The code output TecPlot ASCII format file (*.dat) for visualization and binary files (*.pickle) for data saving. You will also get text files which have data of eigenvalues or gains.
```
python3 CalcMode.py -f [Parameter file name] -p [Parameter name]
mpiexec -np [Number of thread] python3 CalcModeMPI.py -f [Parameter file name] -p [Parameter name]
```

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

```
[ResolventDefault]
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

This project is licensed under the [MIT License](https://github.com/tcnksm/tool/blob/master/LICENCE) - see the [LICENSE.md](LICENSE) file for details


## Acknowledgments

The author gratefully acknowledge members of [Taira laboratory](http://www.seas.ucla.edu/fluidflow/group.html) for insightful comments and discussions on the implementations.
