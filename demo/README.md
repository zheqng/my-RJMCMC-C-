# varied k moves
没有对参数进行排序

每一类选择2条曲线，共8条曲线。
每条曲线减去自身的均值。

- 先验采用平坦先验
<img src = "prior.png" width = 100%>
- 20000步
分类的结果如下：

<img src =  "register.png" width=100%>

聚类的结果
<img src =  "classification.png" width=100%>
但聚类的结果和专家统计结果不同
有两类两类划分为一类。
- 考虑在RJMCMC中加入简单的均值。
