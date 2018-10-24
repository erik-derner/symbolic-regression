function V = model_583509(x)

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

V = 0.06669847 + ...
-0.16882929 *  ( ( cos( sin( ( cos(((x1*1.0 * 1.0)))  -  (((x2*1.0 * 1.0)) + 1.0) ) ) )  .*  cos( ( sin(((x3*1.0 * 1.0)))  .*  ( sin( sin(((x1*1.0 * 1.0))) )  .*  ( (((x2*1.0 * 1.0)) + 1.0)  - 25.23151858) ) ) ) ) )  + ...
-5.64918586 *  ( sin(((x1*1.0 * 1.0))) )  + ...
-0.07441436 *  ( ( ( (-0.5 +  cos(((x1*1.0 * 1.0))) )  +  cos(((x1*1.0 * 1.0))) )  .*  cos( ( (4.13320934 -  cos(((x1*1.0 * 1.0))) )  .*  ( sign(((x3*1.0 * 1.0)))  -  (-18.8495559 + ((x1*1.0 * 1.0))) ) ) ) ) )  + ...
0.92217377 *  ( ( ( sin( sin(((x1*1.0 * 1.0))) )  +  (((x3*1.0 * 1.0)) .* -8.17575769) )  + ((x2*1.0 * 1.0))) )  + ...
1.25204962 *  ( ( ( sin(((x1*1.0 * 1.0)))  -  sign( sin(((x1*1.0 * 1.0))) ) )  .*  sin(((x1*1.0 * 1.0))) ) )  + ...
-0.0469503 *  ( sign( cos( (((x2*1.0 * 1.0)) + 1.0) ) ) )  + ...
0.01531207 *  ( ( (12.66514798 - ((x2*1.0 * 1.0)))  +  ( (0.18628876 .* ((x3*1.0 * 1.0)))  .*  (2998.41000769 + ((x3*1.0 * 1.0))) ) ) )  + ...
0.01108634 *  ( ( ( ( (((x3*1.0 * 1.0)) -  cos(((x3*1.0 * 1.0))) )  -  ( cos(((x1*1.0 * 1.0)))  -  (((x2*1.0 * 1.0)) + 1.0) ) )  .*  ( sign( (((x2*1.0 * 1.0)) +  sign(((x3*1.0 * 1.0))) ) )  .* ((x3*1.0 * 1.0))) )  .*  (1.0 +  ( cos( sin(((x1*1.0 * 1.0))) )  -  cos(((x1*1.0 * 1.0))) ) ) ) )  + ...
0.17777743 *  ( ( ( sign( (0.47189234 +  (((x2*1.0 * 1.0)) -  sign(((x3*1.0 * 1.0))) ) ) )  -  cos( sin( (((x3*1.0 * 1.0)) .* -8.17575769) ) ) )  +  ( ( sin( sin(((x1*1.0 * 1.0))) )  .*  ( (((x2*1.0 * 1.0)) + 1.0)  - 25.23151858) )  .* 0.01249568) ) )  + ...
-0.12184723 *  ( ( ( cos( (-18.8495559 + ((x1*1.0 * 1.0))) )  .*  (0.91392937 .* ((x2*1.0 * 1.0))) )  +  ( sin( cos(((x3*1.0 * 1.0))) )  -  sin( sin( (0.18628876 .* ((x3*1.0 * 1.0))) ) ) ) ) )  + ...
0;

% MSE = 0.0662481102491341

% complexity = 180


% Configuration:
%         seed: 583501
%         nbOfRuns: 15
%         dataset: data/refined2018_5s_8_x2_train.txt
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
%         resultsDir: refined2018_5s_8_x2/
%         tailFunctionSet: Multiply, Plus, Minus, Sine, Cosine, Sgn
% Default solverName: SolverMultiThreaded
%         nbRegressors: 10
%         nbPredictors: 10
% Default improvementThreshold: 0.0
% Default maxNonImprovingEpochs: 2147483647
% Default identityNodeType: identity
