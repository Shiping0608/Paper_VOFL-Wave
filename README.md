# Fourier pseudospectral methods for the  variable-order space fractional wave equations

**Paper DOI:** to be updated  
**Authors:** Yanzhi Zhang, Xiaofei Zhao, and Shiping Zhou  
**Publication:** to be updated

[![Paper](https://img.shields.io/badge/Paper-PDF-red)](https://arxiv.org/abs/2311.13049)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Abstract

In this paper, we propose Fourier pseudospectral methods to solve the variable-order space fractional wave equation and develop an accelerated matrix-free approach for its effective implementation.
In constant-order cases, fast algorithms can be designed via the fast Fourier transforms (FFTs), and the computational cost at each time step is ${\mathcal O}(N\log N)$ with $N$  the total number of spatial points.
In variable-order cases, however, the spatial dependence in the power $s(\mathbf{x})$ leads to the failure of inverse FFTs.
While the direct matrix-vector multiplication approach becomes impractical due to excessive memory requirements.
Hence, we propose an accelerated matrix-free approach for effective implementation in variable-order cases.
The computational and storage costs are ${\mathcal O}(MN\log N)$ and ${\mathcal O}(MN)$, respectively, with $M \ll N$.
Moreover, our method can be easily parallelized to further enhance efficiency.
Numerical studies show that our methods are effective in solving the variable-order space fractional wave equations, especially in high-dimensional cases.
Wave propagation in heterogeneous media is studied in comparison to homogeneous counterparts. We find that wave dynamics in fractional cases become more intricate due to nonlocal interactions. Particularly, dynamics in heterogeneous media are more complex than those in homogeneous media.


## Repository Structure

```
├── Section 3.1/         #
├── Section 3.2/         #
├── Section 3.3/         #
├── Section 4.1/         #
├── Section 4.2/         #
├── Section 4.3/         #
├── Example 1/           #
├── Example 2/           #
├── Example 3/           #
├── requirements.md      # Dependencies
└── README.md            # This file
```
## How to Use
- All code is executed through **text.m**, which generate the simulation data.
- Separate plotting scripts are provided to produce the corresponding figures.

## Citation

```bibtex
@article{Zhang2025,
  title = {Fourier pseudospectral methods for the  variable-order space fractional wave equations},
  author = {Zhang, Y., Zhao, X., and Zhou, S.},
  year = {2025},
  journal = {arXiv},
  doi = {10.48550/arXiv.2311.13049}
}
```

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.
