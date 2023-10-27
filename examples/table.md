Generated using pyGCS (Grid Convergence Study)
- https://github.com/tomrobin-teschner/pyGCS
- https://pypi.org/project/pygcs/

Table 1: Grid convergence study over 3 grids. phi represents the {INSERT MEANING OF PHI HERE} and phi_extrapolated its extrapolated value. N_cells is the number of grid elements, r the refinement ration between two successive grids. GCI is the grid convergence index in percent and its asymptotic value is provided by GCI_asymptotic, where a value close to unity indicates a grid independent solution. The order achieved in the simulation is given by p.

|        |  phi      |   N_cells   |  r  |  GCI  | GCI_asymptotic |  p   | phi_extrapolated |
|--------|:---------:|:-----------:|:---:|:-----:|:--------------:|:----:|:----------------:|
|        |           |             |     |       |                |      |                  |
| Grid 1 | 1.030e+00 |        4500 | 1.3 | 7.48% |                |      |                  |
| Grid 2 | 1.020e+00 |        2500 | 1.3 | 8.78% |      1.010     | 0.51 |     1.09e+00     |
| Grid 3 | 1.010e+00 |        1500 | -   | -     |                |      |                  |
|        |           |             |     |       |                |      |                  |