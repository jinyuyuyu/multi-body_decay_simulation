# Ne16

* 模拟代码在[sims目录](https://github.com/jinyuyuyu/multi-body_decay_simulation/tree/master/Ne16/sims)下。

* 模拟的主要物理过程

  * [主程序 sim.cpp](https://jinyuyuyu.github.io/multi-body_decay_simulation/Ne16/man/sim.html)
  
* 模拟中调用的类的说明

  * [CPlf](https://jinyuyuyu.github.io/multi-body_decay_simulation/Ne16/man/CPlf.html)：projectile-like fragment，描述2p衰变母核。
  
  * CFrag：描述末态粒子的信息、在物质中的反应和在探测器中的探测。
  
  * [CFrame](https://jinyuyuyu.github.io/multi-body_decay_simulation/Ne16/man/CFrame.html)：描述原子核的基本性质和运动学特征，进行速度、能量、动量的相互转换和不同坐标系的速度变换 (可选相对论/非相对论)。
  
  * CDecay：描述衰变过程、不变质量谱和关联信息的重建。
  
  * CLoss：描述能损过程。
  
  * CRange：计算质子在CsI中的射程和横向歧离。
  
  * CMulScat：描述多重散射过程。
  
  * CMomDist：描述敲出反应的动量分布。
  
  * histo：存储histogram。
  
  * constants.h：存储物理学常数和原子核质量信息。
