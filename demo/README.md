# Large K

- 种子
113584957
487306701
812614905
- 3类，每一类100条曲线，
<img src =  "simudata.png" width=100%>

参数如下
 w 1 0.5 10
 v 0.2 1 0.2
 sigma2 0.0025 0.001 0.0005

## 迭代40000步，结果如下
- 对数似然
<img src =  "log lik.png" width=100%>

成分数的变化
<img src =  "K.png" width=100%>

- 参数估计结果

取后一半迭代结果，前一半作为warm-up, 并每隔100次保留一次样本。
对参数进行排序。
计算参数的均值如下
w 1.05 0.49 10.14
v 0.18 1.06 0.19
sigma2 0.0025 0.001 0.0005


相对误差%

w 4.77775 2.22120 1.38735
v 8.7320 5.5924 3.9915
sigma2 0 0 0


- 参数随着迭代的变化
<img src =  "para trace.png" width=100%>
<img src =  "para panel.png" width=100%>

- 参数收敛的指标

 迭代20000步

    Rhat    ESS
pi  1      122       96
1      125       79
 1.03      200       78
w 1.22       15       23
 1.77        4       32
 1.14       11       24
v 1.32       13       45
 1.21        9       19
  1.37        5       15
sigma2  1      100      100
1      100      100
1      100      100

迭代40000步
Rhat    ESS
pi  1      203      184
1      179      182
0.99      223      226
w 1.43        5       25
1.03       22       22
1.03       29       81
v 1.06       32       68
1.1       34       75
 1.02       27       61
sigma2 1      196      196
1      196      196
1      196      196

For each parameter, Bulk_ESS and Tail_ESS are crude measures of
effective sample size for bulk and tail quantities respectively (an ESS > 100
per chain is considered good), and Rhat is the potential scale reduction
factor on rank normalized split chains (at convergence, Rhat <= 1.05).

这说明40000步迭代还不够收敛。

- 参数的描述性分析
<img src =  "K box.png" width=100%>

summary(as.factor(datha$K))
    2     3     4     5     6
    5 39676   297    19     3

- 聚类的结果
对验证集的误差0.67%

- 分类结果
 新生成600条曲线，用参数的平均值进行分类，错误率0%
