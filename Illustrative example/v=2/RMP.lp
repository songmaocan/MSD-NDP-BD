\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  z_0 + z_1
Subject To
 budget: y_0_3 + y_1_2 <= 1
 CUT_0_0: 4 y_0_3 + z_0 >= 7
 CUT_0_1: 3 y_1_2 + z_1 >= 7
 CUT_1_0: z_0 >= 4
 CUT_1_1: 3 y_1_2 + z_1 >= 7
 CUT_2_0: z_0 >= 4
 CUT_2_1: 3 y_1_2 + z_1 >= 7
Bounds
Binaries
 y_0_3 y_1_2
End
