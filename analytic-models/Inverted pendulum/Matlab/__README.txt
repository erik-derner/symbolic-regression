To run a real experiment:
1. Open Matlab
2. switch_sym_models('../Real data/models/random_x1/model_580006.m', '../Real data/models/random_x2/model_580530.m'); % initial model
OR switch_sym_models('../Real data/models/refined_x1/model_583030.m', '../Real data/models/refined_x2/model_583509.m'); % refined model
3. pend_disc_VI();
4. rtswingup();

To run a simulation:
1. Calculate analytic models (use the bat files in the subdirectories of the 'Simulation data' directory)
2. Open Matlab
3. Use switch_sym_models() to switch to analytic models generated in step 1 (see above for example usage)
4. pend_disc_VI();
5. simpi_movie_V_pi();
