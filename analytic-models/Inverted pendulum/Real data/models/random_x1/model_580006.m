function V = model_580006(x)

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

V = 0.18242701 + ...
-0.0167668 *  ( ( ( sign( (((x2*1.0 * 1.0)) + 0.10276921) )  +  cos( cos(((x3*1.0 * 1.0))) ) )  + ((x3*1.0 * 1.0))) )  + ...
0.00156189 *  ( ( sign(((x1*1.0 * 1.0)))  +  ( ( cos(((x3*1.0 * 1.0)))  - ((x2*1.0 * 1.0)))  .*  cos(((x1*1.0 * 1.0))) ) ) )  + ...
-0.00401201 *  ( ( cos(((x2*1.0 * 1.0)))  .*  sign(((x3*1.0 * 1.0))) ) )  + ...
0.00417499 *  ( cos( ( cos(((x3*1.0 * 1.0)))  .* 7.81426766) ) )  + ...
2.3016E-4 *  ( ( ( ( ( cos(((x2*1.0 * 1.0)))  + ((x3*1.0 * 1.0)))  -  cos( (((x2*1.0 * 1.0)) .* 1.73205081) ) )  .*  ( ( ( cos(((x1*1.0 * 1.0)))  +  cos(((x1*1.0 * 1.0))) )  -  (((x1*1.0 * 1.0)) .* -0.42227491) )  - 2.87886793) )  .*  sign( (((x2*1.0 * 1.0)) + 0.10276921) ) ) )  + ...
-0.01164439 *  ( cos( (1.73205081 +  cos(((x1*1.0 * 1.0))) ) ) )  + ...
-0.1381032 *  ( sin(((x1*1.0 * 1.0))) )  + ...
0.20312861 *  ( (-0.24044433 .*  ( (3.57754615 - ((x3*1.0 * 1.0)))  - ((x2*1.0 * 1.0))) ) )  + ...
0.00549741 *  ( sin( ( ( (((x2*1.0 * 1.0)) - 4.96196185)  .* 15.49991819)  .* 2.09635705) ) )  + ...
-3.99880096 *  ( (((x1*1.0 * 1.0)) .* -0.25) )  + ...
0;

% MSE = 4.97015317374297E-6

% complexity = 112


% Configuration:
%         seed: 580001
%         nbOfRuns: 30
%         dataset: data/rtrandom2018_5s_8_x1_train.txt
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
%         resultsDir: rtrandom2018_5s_8_x1/
%         tailFunctionSet: Multiply, Plus, Minus, Sine, Cosine, Sgn
% Default solverName: SolverMultiThreaded
%         nbRegressors: 10
%         nbPredictors: 10
% Default improvementThreshold: 0.0
% Default maxNonImprovingEpochs: 2147483647
% Default identityNodeType: identity
