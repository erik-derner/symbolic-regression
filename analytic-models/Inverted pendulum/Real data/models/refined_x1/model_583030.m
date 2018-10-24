function V = model_583030(x)

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

V = -0.03416977 + ...
-0.00270509 *  ( ( cos(((x1*1.0 * 1.0)))  .*  (((x2*1.0 * 1.0)) -  sin( ( cos(((x3*1.0 * 1.0)))  + -0.85154002) ) ) ) )  + ...
-0.00342441 *  ( sign( ( ( (-0.58578644 - ((x2*1.0 * 1.0)))  .*  (((x2*1.0 * 1.0)) - 1.22507908) )  .*  cos(((x3*1.0 * 1.0))) ) ) )  + ...
0.05169766 *  ( (((x3*1.0 * 1.0)) .* 0.64963694) )  + ...
-0.13978685 *  ( sin(((x1*1.0 * 1.0))) )  + ...
-0.00469305 *  ( sign( sign( (((x2*1.0 * 1.0)) +  (((x3*1.0 * 1.0)) .* 2.24016786) ) ) ) )  + ...
0.00290958 *  ( ( cos( ( (((x1*1.0 * 1.0)) .* 4.0)  +  (((x1*1.0 * 1.0)) .* 4.0) ) )  .*  ( sign(((x2*1.0 * 1.0)))  +  (1.31607401 -  cos(((x1*1.0 * 1.0))) ) ) ) )  + ...
0.01208414 *  ( ( sign( (0.84089642 - ((x2*1.0 * 1.0))) )  .*  cos( cos(((x1*1.0 * 1.0))) ) ) )  + ...
-0.00149642 *  ( ( ( sign(((x2*1.0 * 1.0)))  +  (((x3*1.0 * 1.0)) .* 2.24016786) )  .*  sign( ( (0.87758256 +  cos(((x2*1.0 * 1.0))) )  - ((x3*1.0 * 1.0))) ) ) )  + ...
0.04942012 *  ( ( (0.82335966 .* ((x1*1.0 * 1.0)))  -  (-0.58578644 - ((x2*1.0 * 1.0))) ) )  + ...
0.95917665 *  (((x1*1.0 * 1.0)))  + ...
0;

% MSE = 7.282985884809281E-5

% complexity = 110


% Configuration:
%         seed: 583016
%         nbOfRuns: 15
%         dataset: data/refined2018_5s_8_x1_train.txt
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
%         resultsDir: refined2018_5s_8_x1/
%         tailFunctionSet: Multiply, Plus, Minus, Sine, Cosine, Sgn
% Default solverName: SolverMultiThreaded
%         nbRegressors: 10
%         nbPredictors: 10
% Default improvementThreshold: 0.0
% Default maxNonImprovingEpochs: 2147483647
% Default identityNodeType: identity
