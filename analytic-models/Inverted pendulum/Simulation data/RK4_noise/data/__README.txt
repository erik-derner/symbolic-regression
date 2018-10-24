Parameters in the filename:

1dof_N1L#1#_RK4_fri00_5V_s#2#_#3#.txt

#1# ... number of data entries
#2# ... noise: 00 = no noise, 01 = noise_std is 0.01 * range, 05 = noise_std is 0.05 * range, 10 = noise_std is 0.1 * range
        In particular:
		#2# value    | 00 | 01        | 05        | 10
        noise_std_x1 |  0 | 0.01 * pi | 0.05 * pi | 0.1 * pi
        noise_std_x2 |  0 | 0.01 * 40 | 0.05 * 40 | 0.1 * 40
#3# ... output: 4 = x1 (angle), 5 = x2 (angular velocity)
