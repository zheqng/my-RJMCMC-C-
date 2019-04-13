<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>
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

<img src =  "mix-EM total rmse.png" width>

fix MC和变化的MC

<img src =  "fix-MC rmse.png" width>

只有第三类的rmse大

- log likelihood

<img src =  "fix k move loglik.png" width=50%>

整体下降趋势，而且偶尔上升。不正确、最大似然应该朝着似然上升的方向发展，检查程序，发现nuts中的hamilton能量用的$ log L(\theta) - K(\phi) $ 而Neal用的是$ -log L(\theta) + K(\phi)$
改正

- fix k move中不用death空成分。
