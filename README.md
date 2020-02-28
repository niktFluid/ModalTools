# ModalTools

This codes construct a linear operator from a field data of compressible fluid flows.

## Description

ModalTools is a tool for constructing and analyzing a linearized Navier-Stokes (LNS) operator in a modal manner. The project contains necessary codes for estimating the Right-Hand-Side of LNS equation and a matrix generator that uses so-called zero-one trick for making a sparse matrix representing LNS operator. The operator can be used for modal analysis including linear stability analysis and resolvent analysis (See [here](https://doi.org/10.2514/1.J056060)). We also provide a low Reynolds number compressible cylinder flow data via OpenFOAM for the demonstration.

## Demo

## Requirements
* Python 3.*
* NumPy 
* SciPy
* [Ofpp](https://github.com/dayigu/ofpp) - Used to parse OpenFOAM data.
* [mpi4py](https://mpi4py.readthedocs.io/en/stable/) - Optional but recommended.

## Usage

## Author

* **Yoimi Kojima** - [niktFluid](https://github.com/niktFluid)

## Licence

This project is licensed under the [MIT License]((https://github.com/tcnksm/tool/blob/master/LICENCE)) - see the [LICENSE.md](LICENSE) file for details


## Acknowledgments

The author gratefully acknowledge members of [Taira laboratory](http://www.seas.ucla.edu/fluidflow/group.html) for insightful comments and discussions on the implementations.
