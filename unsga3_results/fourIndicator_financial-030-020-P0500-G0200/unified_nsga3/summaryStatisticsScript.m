% CSE  848: Evolutionary Computation Semester Project
% Authors: Farhan Hormasji and Bonnie Reiff
% Filename: summaryStatisticsScript.m
% Description: Given a number of runs of U-NSGA-III computes summary
% statistics and outputs graphs. The statistical analyses performed can be
% controlled by binary values at the beginnging of the script.
% Assumption: The script is placed at the same level as the
% "generation_wise_runXXX" folders and therefore requires no input.

clear all; % Clear the variable workspace

%% Summary Statistics to be Performed %%
% Modify here to turn statistics/graphs on or off
avgFitnessPerGenerationGraphs             = 1;
paretoFrontFinalPopulationGraphAllRuns    = 1;
paretoFrontFinalPopulationGraphSingleRun  = 0;
objValueSummaryStatisticsAllRuns          = 1;
objValueSummaryStatisticsSingleRun        = 0;
parameterSummaryStatisticsAllRuns         = 1;
parameterSummaryStatisticsSingleRun       = 0;

if ( avgFitnessPerGenerationGraphs || ...
     paretoFrontFinalPopulationGraphAllRuns || ...
     objValueSummaryStatisticsAllRuns || ...
     parameterSummaryStatisticsAllRuns )
 
    %% Put all run data into a single file
    objectiveValuesPerGenerationData = [];
    finalPopulationObjectiveValues = [];
    finalPopulationParameterValues = [];
    numRuns = 20;
    for run = 1:numRuns
        % Construct the folder name
        zeroBasedRunNum = run - 1;
        if (zeroBasedRunNum < 10)
            folderPrefix = 'generation_wise_run00';
        elseif ((zeroBasedRunNum >= 10) && (zeroBasedRunNum < 100))
            folderPrefix = 'generation_wise_run0';
        elseif (zeroBasedRunNum >= 100)
            folderPrefix = 'generation_wise_run';
        end
        folderName = strcat(folderPrefix, int2str(zeroBasedRunNum));
        fprintf('Entering folder %s...\n', folderName); % DEBUG. REMOVE!
        
        callingFolder = cd(folderName);
        
        % Get a list of all objective files in the folder
        objFiles = dir('*_obj.dat');
        objFileFilenames = {objFiles.name};
     
        if ( avgFitnessPerGenerationGraphs )            
            for gen = 1:length(objFileFilenames)
                fileData = importdata(objFileFilenames{gen});
                if (run == 1)
                    objectiveValuesPerGenerationData(gen, 1) = size(fileData, 2); % number of columns = number of individuals
                    objectiveValuesPerGenerationData(gen, 2) = sum(fileData(:,1));
                    objectiveValuesPerGenerationData(gen, 3) = sum(fileData(:,2));
                    % objectiveValuesPerGenerationData(gen, 4) = sum(fileData(:,3));
                else
                    objectiveValuesPerGenerationData(gen, 1) = objectiveValuesPerGenerationData(gen, 1) + size(fileData, 2); % number of columns = number of individuals
                    objectiveValuesPerGenerationData(gen, 2) = objectiveValuesPerGenerationData(gen, 2) + sum(fileData(:,1));
                    objectiveValuesPerGenerationData(gen, 3) = objectiveValuesPerGenerationData(gen, 3) + sum(fileData(:,2));
                    % objectiveValuesPerGenerationData(gen, 4) = objectiveValuesPerGenerationData(gen, 4) + sum(fileData(:,3));
                end
            end
        end
        if ( paretoFrontFinalPopulationGraphAllRuns || objValueSummaryStatisticsAllRuns )
            % Note: final generation number hardcoded. change this for
            % different runs!
            fileData = importData('gen_0199_obj.dat');
            finalPopulationObjectiveValues = [finalPopulationObjectiveValues; fileData];
        end
        if ( parameterSummaryStatisticsAllRuns )
            % Note: final generation number hardcoded. change this for
            % different runs!
            fileData = importData('gen_0199_var.dat');
            finalPopulationParameterValues = [finalPopulationParameterValues; fileData];
        end
        
        cd(callingFolder);
    end
    
    if (avgFitnessPerGenerationGraphs)
        generations = [1:size(objectiveValuesPerGenerationData, 1)];
        figure;
        yvalues = ((-1) .* objectiveValuesPerGenerationData(:, 2)) ./ objectiveValuesPerGenerationData(:, 1);
        plot(generations, yvalues, '-o');
        title('Average Objective 1 Value Per Generation', 'FontSize', 12);
        xlabel('Generation Number', 'FontSize', 12);
        ylabel('Objective 1 Value', 'FontSize', 12);
        figure;
        yvalues = ((-1) .* objectiveValuesPerGenerationData(:, 3)) ./ objectiveValuesPerGenerationData(:, 1);
        plot(generations, yvalues, '-o');
        title('Average Objective 2 Value Per Generation', 'FontSize', 12);
        xlabel('Generation Number', 'FontSize', 12);
        ylabel('Objective 2 Value', 'FontSize', 12);
    end
    
    if ( paretoFrontFinalPopulationGraphAllRuns )
        finalPopulationObjectiveValues = finalPopulationObjectiveValues .* (-1);
        figure;
        scatter(finalPopulationObjectiveValues(:,1), finalPopulationObjectiveValues(:, 2));
        title('Final Population Objective 2 versus Objective 1 Values for All Runs', 'FontSize', 12);
        xlabel('Objective 1 Value', 'FontSize', 12);
        ylabel('Objective 2 Value', 'FontSize', 12);
    end
    
    if ( objValueSummaryStatisticsAllRuns )
        finalPopulationObjectiveValues = finalPopulationObjectiveValues .* (-1);
        fprintf('[ Objective Value Summary Statistics Over All Runs ]\n');
        meanValues = mean(finalPopulationObjectiveValues);
        minValues = min(finalPopulationObjectiveValues);
        maxValues = max(finalPopulationObjectiveValues);
        fprintf('\tObjective 1: Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', 100 * meanValues(1), 100 * minValues(1), 100 * maxValues(1)); 
        fprintf('\tObjective 2: Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(2), minValues(2), maxValues(2));
    end
    
    if ( parameterSummaryStatisticsAllRuns )
        fprintf('[ Parameter Value Summary Statistics Over All Runs ]\n');
        meanValues = mean(finalPopulationParameterValues);
        minValues = min(finalPopulationParameterValues);
        maxValues = max(finalPopulationParameterValues);
        fprintf('\tDEMAC Short Lookback: Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(1), minValues(1), maxValues(1)); 
        fprintf('\tDEMAC Long Lookback : Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(2), minValues(2), maxValues(2));
        fprintf('\tMACD Short Lookback : Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(3), minValues(3), maxValues(3));
        fprintf('\tMACD Long Lookback  : Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(4), minValues(4), maxValues(4));
        fprintf('\tMACD Signal Lookback: Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(5), minValues(5), maxValues(5));
        fprintf('\tRSI Lookback        : Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(6), minValues(6), maxValues(6));
        fprintf('\tRSI Lower Boundary  : Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(7), minValues(7), maxValues(7));
        fprintf('\tRSI Upper Boundary  : Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(8), minValues(8), maxValues(8));
        fprintf('\tMARSI RSI Lookback  : Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(9), minValues(9), maxValues(9));
        fprintf('\tMARSI Lower Boundary: Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(10), minValues(10), maxValues(10));
        fprintf('\tMARSI Upper Boundary: Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(11), minValues(11), maxValues(11));
        fprintf('\tMARSI SMA Lookback  : Average = %.4f, Minimum = %.4f, Maximum = %.4f\n', meanValues(12), minValues(12), maxValues(12));
    end
end

if ( paretoFrontFinalPopulationGraphSingleRun || ...
     objValueSummaryStatisticsSingleRun || ...
     parameterSummaryStatisticsSingleRun )
end

