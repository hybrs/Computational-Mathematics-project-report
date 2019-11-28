## Primal-Dual Interior-Point method for Constrained CQP problem
#### Computational Matemathics project report | a.y. 2018/19 | Prof. Poloni, Prof. Frangioni
##### Authors : [Giuliano Cornacchia](https://github.com/GiulianoCornacchia), [Mario Leonardo Salinas](https://github.com/hybrs)

### Intro
This repository contains the implemented soluition for the noML project #3 of *[Computational Mathematics for Learning and Data Analysis @ UniPi MD Computer Science](https://elearning.di.unipi.it/enrol/index.php?id=131)*: that is to design deploy and evaluate a Primal-Dual Interior-Point algorithm to solve a specific convex quadratic optimization problem with linear constraints. 

In ![this report](main.pdf) we describe the problem, our implementation choices and experimental results. All the code is developed in ``MATLAB`` and can be found in the ``src`` folder of this repository.

#### Usage

To run an experiment with our implemented primal-dual algorithm you have to:
- initialize a constrained CQP problem with **proper dimensions** (see [report](main.pdf) for details)
- run the algorithm with the desired linear solver (``'gmres' | 'ldl'``)

```
>> n = 1000; m = 200; delta = 0.4;
>> p = genProblem(n, m, delta);
>> [x, fval, lambda, exit_code, tm, it, gaps] = PDIP(p, 100, 1e-14, "ldl");

iter     fval     gap         dualf     primalf     sTx
-----------------------------------------------------------------------------
1     2.402e+01 1.431e+01   4.694e+01 5.661e-16 3.438e+02
2    -4.094e+00 3.757e+01   7.478e+00 7.764e-15 1.538e+02
                            .
                            .
                            .
19    6.074e+01 7.405e-14   4.846e-13 1.404e-15 4.500e-12
20   -6.074e+01 7.019e-15   3.892e-14 1.369e-15 4.255e-13

Execution terminated because duality gap reduced under the threshold
Primal-Dual Interior Point method terminated in 20 iterations
elapsed time is 3.3502 seconds
fval = -6.074e+01 and complementary gap = 7.019e-15
```

