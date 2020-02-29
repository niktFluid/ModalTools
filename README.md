# ModalTools
This codes construct a linear operator from a field data of compressible fluid flows.


## Description
ModalTools is a tool for constructing and analyzing a linearized Navier-Stokes (LNS) operator in a modal manner. The project contains necessary codes for estimating the Right-Hand-Side of LNS equation and a matrix generator that uses so-called zero-one trick for making a sparse matrix representing LNS operator. The operator can be used for modal analysis including linear stability analysis and resolvent analysis (See [here](https://doi.org/10.2514/1.J056060)). We also provide a low Reynolds number compressible cylinder flow data via [OpenFOAM](https://www.openfoam.com/) for the demonstration.


## Demo
This section shows a demonstration of the tool for compressible cylinder flow field at Re=150, M=0.2. The simulation is conducted using OpenFOAM. The flow field shows vortex shedding behind the cylinder at frequency of St=0.178. 

![Figure1](https://user-images.githubusercontent.com/47338366/75595128-06665d00-5a40-11ea-8ac9-87b0e0d344b9.png)

* **Stability Analysis**

Solving an eigenvalue problem of the operator provides a spectrum that represents the characteristic of the dynamic response of a linear system. Eigenvectors corresponding to each eigenvalue, depicts structures of fluid flow. In our formulation, the real parts of eigenvalues indicate growth rates of the structure, while the imaginary parts correspond to its frequency. The following figure shows the results of stability analysis for the cylinder flow. We can find some eigenvalue on the unstable plane Real(lambda) > 0. The imaginary part of the higher value one (lambda = 2.89e-3 + 0.218i) gives a frequency of St = 0.173 which is fairly close to the shedding frequency obtained from the simulation. Indeed, the eigenvector corresponds to this eigenvalue apparently shows structures like vortex shedding behind the cylinder as it can be seen on the right side of the following figure. Some unstable eigenvalues that are just on the imaginary axis may come from insufficient time-averaging for the data.

![Figure2](https://user-images.githubusercontent.com/47338366/75592649-ad92c680-5a37-11ea-84fd-17f2ed069b0f.png)

* **Resolvent Analysis**

The other feature implemented on the codes involves with a pseudospectrum of a liner operator. The resolvent analysis can be understood in the context that an analysis of harmonic input forcing and output structure of the linear systems. Singular Value Decomposition (SVD) of the operator delivers pairs of orthogonal basis correspond to the forcing,￿ response modes of the operator and energy amplitude ratio (gain) between them. The following picture shows distribution of three largest gains over frequency and resolvent modes correspond to the shedding frequency. The leading gain shows its highest value at the shedding frequency since the phenomenon comes from fluid dynamic unsteadiness. The response modes shows vortex shedding like structure while the forcing is looked like some sort of velocity fluctuations in the vicinity of the cylinder walls. We also provides an implementation of [randomized resolvent analysis](https://arxiv.org/abs/1902.01458) that makes your assessments faster.

![Figure3](https://user-images.githubusercontent.com/47338366/75595182-2bf36680-5a40-11ea-8814-381d05360ef8.png)


## Requirements
* Python 3.*
* NumPy 
* SciPy
* [Ofpp](https://github.com/dayigu/ofpp) - Used to parse OpenFOAM data.
* [mpi4py](https://mpi4py.readthedocs.io/en/stable/) - Optional but recommended.

## Usage
* **Generating Operator**


* **Stability or Resolvent Analysis**



## Author
* **Yoimi Kojima** - [niktFluid](https://github.com/niktFluid)


## Licence

This project is licensed under the [MIT License](https://github.com/tcnksm/tool/blob/master/LICENCE) - see the [LICENSE.md](LICENSE) file for details


## Acknowledgments

The author gratefully acknowledge members of [Taira laboratory](http://www.seas.ucla.edu/fluidflow/group.html) for insightful comments and discussions on the implementations.
