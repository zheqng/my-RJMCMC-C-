# Fix k moves
- 修改了xixj函数，因为预测时，xixj(traindata$X,testdata$X)维度不同，所以矩阵不对称
- 预测9000多步，结果如下
<img src =  "fix k move 9229 iter.png" width=100%>
- 发现z没有变化，在主程序中加入gibbs sample z
<img src =  "fix k move 2.png" width=100%>

z:

2  2  2  2  1  2  0  0  0,
2  2  2  2  1  2  0  0  0

预测结果比9000多步好．
