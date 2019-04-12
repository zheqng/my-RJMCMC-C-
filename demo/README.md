# Fix k moves

- 预测10000步，结果如下
<img src =  "fix k move 10000.png" width=100%>

- z:

1  1  1  1  1  1  0  0  0
类别数目不正确。考虑用原来的数据做。

- rmse:
0.07275178

- 每一类的rmse

0.0127 0.0064 0.1252

- 每一类的相关性

0.9996  0.99997    0.9698

-对比原来的方法：

EM算法

<img src =  "mix-EM total rmse.png" width=50%>

fix MC和变化的MC

<img src =  "fix-MC rmse.png" width>

只有第三类的rmse大
