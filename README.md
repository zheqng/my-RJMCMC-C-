# fix k moves

- 修改logsumexp
- 验证了solve,虽然计算结果和R不同，但都满足A*x = b
  这样，即使logP不同，但都可以，基本差不多
- 修改nuts中迭代iter=2,t0=2
