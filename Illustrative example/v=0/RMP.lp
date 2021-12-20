\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  z_0 + z_1
Subject To
 budget: y_0_3 + y_1_2 <= 1
 CUT_0_0: y_0_3 + z_0 >= 3
 CUT_0_1: 2 y_1_2 + z_1 >= 3
 CUT_1_0: y_0_3 + z_0 >= 3
 CUT_1_1: z_1 >= 1
Bounds
Binaries
 y_0_3 y_1_2
End
