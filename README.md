# 说明


-  openmp + eigen
- openmp omp_num_threads=2 最快
- eigen  修改了nuts中，尤其是buildtree中，p_sharp_right在递归中无法输入，因而

```c++

      p_sharp_dummy = p_sharp_plus;
       int n1 = BuildTree(Data, phi, theta, phi_propose, theta_propose, p_sharp_minus, p_sharp_plus, rho_left, util, depth - 1, epsilon, lambda, z);
       p_sharp_plus = p_sharp_dummy;


       if (!util.criterion) {
                return 0; // early stopping
        }

        p_sharp_dummy = p_sharp_minus;
                int n2 = BuildTree(Data, phi, theta, phi_propose_new, theta_propose_new,
                                   p_sharp_minus, p_sharp_plus, rho_right, util, depth - 1, epsilon, lambda, z);
                p_sharp_minus = p_sharp_dummy;

```

在计算rhat时，调用bayesplot包计算rhat

- eigen比arma包要快得多，80000次原来需要23D,现在160000只需要24h
