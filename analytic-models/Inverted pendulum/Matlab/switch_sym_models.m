function switch_sym_models(pathname_model_x1, pathname_model_x2)
%SWITCH_SYM_MODELS 1DOF - Set symbolic models to be used to calculate next state
% Input:
%  pathname_model_x1 - [string] pathname to the symbolic model for angle 
%                      in the next state
%  pathname_model_x2 - [string] pathname to the symbolic model for angular 
%                      velocity in the next state

eval(sprintf('copyfile ''%s'' fsm_model_x1.m', pathname_model_x1));
eval(sprintf('copyfile ''%s'' fsm_model_x2.m', pathname_model_x2));

end

