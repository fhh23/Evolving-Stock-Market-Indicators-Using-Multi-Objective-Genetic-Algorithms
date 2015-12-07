% CSE  848: Evolutionary Computation Semester Project
% Authors: Farhan Hormasji and Bonnie Reiff

%% plotParetoEvolution_exponentialSamplingScale: TODO
% Input: inputDir, the name of the folder of results to use (one folder
%           level about unified_nsga3)
%        indicatorName, the name of the indicator to be used for the plot

function [] = plotParetoEvolution_exponentialSamplingScale( inputDir, indicatorName )

% Assumption that the first run (run000) will be used for the plot
relativePath = strcat(inputDir, '/unified_nsga3/generation_wise_run000/');

%% Read in the results of the five desired generations
% Values are multiplied by -1 to show the maximum of the original objective
% functions
g001 = dlmread(strcat(relativePath, 'gen_0001_obj.dat'));
g001 = g001 .* (-1);
% g004 = dlmread(strcat(relativePath, 'gen_0004_obj.dat'));
% g004 = g004 .* (-1);
g009 = dlmread(strcat(relativePath, 'gen_0009_obj.dat'));
g009 = g009 .* (-1);
% g016 = dlmread(strcat(relativePath, 'gen_0016_obj.dat'));
% g016 = g016 .* (-1);
% g025 = dlmread(strcat(relativePath, 'gen_0025_obj.dat'));
% g025 = g025 .* (-1);
g036 = dlmread(strcat(relativePath, 'gen_0036_obj.dat'));
g036 = g036 .* (-1);
% g049 = dlmread(strcat(relativePath, 'gen_0049_obj.dat'));
% g049 = g049 .* (-1);
g081 = dlmread(strcat(relativePath, 'gen_0081_obj.dat'));
g081 = g081 .* (-1);
g144 = dlmread(strcat(relativePath, 'gen_0144_obj.dat'));
g144 = g144 .* (-1);
g199 = dlmread(strcat(relativePath, 'gen_0199_obj.dat'));
g199 = g199 .* (-1);

%% Plot all read in generations on the same graph
% Use different colors and marker types for each generation

figure;
hold on

plot(g001(:,1), g001(:,2), 'ro');
plot(g009(:,1), g009(:,2), 'md');
plot(g036(:,1), g036(:,2), 'b+');
plot(g081(:,1), g081(:,2), 'c^');
plot(g144(:,1), g144(:,2), 'g*');
plot(g199(:,1), g199(:,2), 'kx');

legend('G0001', 'G0009', 'G0036', 'G0081', 'G0144', 'G0199');

hold off

%% Set various graph properties
xlabel('Objective 1 Value', 'FontSize', 12);
ylabel('Objective 2 Value', 'FontSize', 12);
title(strcat('Evolution of the Pareto Front for the ', indicatorName, ' Indicator for 800 Individuals'), 'FontSize', 12);
% xlim([])
% ylim([])

%% Summary statistics for the individual technical indicator
% Average, maximum, and minimum value of the objective in the final
% population (using g199 data)
averageObj1 = mean(g199(:,1));
averageObj2 = mean(g199(:,2));
minValues = min(g199);
maxValues = max(g199);
% Display the summary statistics:
fprintf('Average Objective 1 Value: %.4f, Range: [ %.4f, %.4f ]\n', 100 * averageObj1, 100 * minValues(1), 100 * maxValues(1) );
fprintf('Average Objective 2 Value: %.4f, Range: [ %.4f, %.4f ]\n', averageObj2, minValues(2), maxValues(2) );

end