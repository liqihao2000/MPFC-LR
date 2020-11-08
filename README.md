# MPFC-LR
Matlab code of paper "Efficient second-order unconditionally stable numerical schemes for the modified phase field crystal model with long-range interaction" by Qi Li, Liquan Mei, and Yibao Li.

# Steps:
You can follow the steps below to run this program:

1. Temporal accuracy test
```matlab
cd example06_2D_CN_in_Time_exact_o_2
run main06  % need long time
cd ../example06_2D_CN_in_Time_exact_ok
run main06  % need long time
cd ../example09_2D_BDF_in_Time_exact_o_2
run main09  % need long time
cd ../example09_2D_BDF_in_Time_exact_ok
run main09  % need long time
cd ../
run error_plot
```
Get three figures in folder **figure_MPFC_SAV**. This step will take a long time.

2. Comparison with IEQ approach
```matlab
cd example06_2D_CN_in_time_exact_SAVvsIEQ
run main06_SAV
run main06_IEQ
cd ../example09_2D_CN_in_time_exact_SAVvsIEQ
run main09_SAV
run main09_IEQ
```
Get the numerical results displayed in **Table 1**.

3. Energy stability test
```matlab
cd ../example14_2D_CN_energy_and_mass_ok
run main14_1_1
run main14_1_2
run main14_2_1
run main14_2_2
run main14_2_3
run main14_fig3
cd ../example15_2D_BDF_energy_and_mass_ok
run main15_fig3
```

4. The evolution of the phase transition behavior in 2D
```matlab
cd ../example16_2D_random_init_v2
run main16_2_1
run main16_2_2
run main16_2_3
run main16_plot
```
Note: need to change the variable **dirname**.

5. The evolution of the phase transition behavior in 3D
```matlab
cd ../example19_3D_random_init_ok_L50
run main19_2_1
run main19_2_2
run main19_2_3
% Plot
run main19_2_plot
run main19_2_plot_slice
run main19_2_plot_slice3
```
Note: need to change the variable **dirname**.

6. Crystal growth in 2D
```matlab
cd ../example17_2D_3Difference_Crystral_growth_ok
run main17_2
run main17_plot
cd ../example21_2D_single_init
run main21_1
run main21_plot
```

7. Crystal growth in 3D
```matlab
cd ../example18_3D_random_init_ok
run main18_2
run main18_2_plot_slice
```

After done the above steps, some figures will be generated in folder **figure_MPFC_SAV**, and then using the application **JPG2EPS.exe** to convert these figures into EPS format.



If you have any questions, please contact liqihao2000@126.com.

