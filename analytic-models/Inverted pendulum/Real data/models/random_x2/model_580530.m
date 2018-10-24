function V = model_580530(x)

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

V = -57.88307806 + ...
-0.17708351 *  ( ( ( sin(((x2*1.0 * 1.0)))  .*  sign(((x3*1.0 * 1.0))) )  +  ( (18.0 -  (((x1*1.0 * 1.0)) - -302.89091258) )  -  (((x1*1.0 * 1.0)) - -16.0) ) ) )  + ...
4.57603252 *  ( ( sin( sin(((x1*1.0 * 1.0))) )  +  (-0.42995636 .*  (((x1*1.0 * 1.0)) - -16.0) ) ) )  + ...
-0.40684514 *  ( ( sign( sin(((x1*1.0 * 1.0))) )  +  sin( ( sign(((x2*1.0 * 1.0)))  + ((x1*1.0 * 1.0))) ) ) )  + ...
35.30200124 *  ( ( sin( ( sin(((x1*1.0 * 1.0)))  + -0.125) )  .* -0.96302607) )  + ...
26.83578032 *  ( sin( sin( ( sin(((x1*1.0 * 1.0)))  + -0.125) ) ) )  + ...
0.00414278 *  ( ( ( (((x1*1.0 * 1.0)) - -302.89091258)  +  sin( ( (-22.59718444 -  (((x3*1.0 * 1.0)) + ((x1*1.0 * 1.0))) )  - ((x1*1.0 * 1.0))) ) )  .*  (24.04540769 +  (((x3*1.0 * 1.0)) + ((x1*1.0 * 1.0))) ) ) )  + ...
-0.07077078 *  ( ( sign( cos(((x3*1.0 * 1.0))) )  +  cos( ( sin(((x2*1.0 * 1.0)))  .* ((x2*1.0 * 1.0))) ) ) )  + ...
-0.09238959 *  ( (((x2*1.0 * 1.0)) .*  cos(((x1*1.0 * 1.0))) ) )  + ...
-0.56505364 *  ( sign( (((x3*1.0 * 1.0)) - ((x2*1.0 * 1.0))) ) )  + ...
0.2999161 *  ( ( (((x2*1.0 * 1.0)) -  sign(((x2*1.0 * 1.0))) )  +  ( ( sin(((x1*1.0 * 1.0)))  +  (((x1*1.0 * 1.0)) - -16.0) )  +  (((x2*1.0 * 1.0)) + ((x2*1.0 * 1.0))) ) ) )  + ...
0;

% MSE = 0.016346788832507196

% complexity = 127


% Configuration:
%         seed: 580501
%         nbOfRuns: 30
%         dataset: data/rtrandom2018_5s_8_x2_train.txt
%         maxGenerations: 30000
% Default nbOfThreads: 2
%         epochLength: 1000
%         maxEpochs: 30
%         populationSize: 500
%         nbOfTransformedVar: -1
%         lsIterations: 100
% Default maxNodeEvaluations: 9223372036854775807
%         depthLimit: 7
% Default probHeadNode: 0.1
% Default probTransformedVarNode: 0.2
%         useIdentityNodes: true
% Default saveModelsForPython: false
%         optMisplacementPenalty: 0.0
%         desiredOptimum: 0 0 0
%         regressionClass: LeastSquaresFit
%         populationClass: PartitionedPopulation
%         resultsDir: rtrandom2018_5s_8_x2/
%         tailFunctionSet: Multiply, Plus, Minus, Sine, Cosine, Sgn
% Default solverName: SolverMultiThreaded
%         nbRegressors: 10
%         nbPredictors: 10
% Default improvementThreshold: 0.0
% Default maxNonImprovingEpochs: 2147483647
% Default identityNodeType: identity
