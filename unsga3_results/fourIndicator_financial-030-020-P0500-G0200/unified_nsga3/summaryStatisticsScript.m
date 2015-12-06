% CSE  848: Evolutionary Computation Semester Project
% Authors: Farhan Hormasji and Bonnie Reiff
% Filename: summaryStatisticsScript.m
% Description: Given a number of runs of U-NSGA-III computes summary
% statistics and outputs graphs. The statistical analyses performed can be
% controlled by binary values at the beginnging of the script.
% Assumption: The script is placed at the same level as the
% "generation_wise_runXXX" folders and therefore requires no input.

clar all; % Clear the variable workspace

%% Summary Statistics to be Performed %%
% Modify here to turn statistics/graphs on or off
paretoFrontFinalPopulationGraphAllRuns    = 0;
paretoFrontFinalPopulationGraphSingleRun  = 0;
objValueSummaryStatisticsAllRuns          = 0;
objValueSummaryStatisticsSingleRun        = 0;
avgFitnessPerGenetationGraphs             = 1;
finalPopulationObjectiveStatsAllRuns      = 0;
finalPopulationObjectiveStatsSingleRun    = 0;
parameterSummaryStatisticsAllRuns         = 0;
parameterSummaryStatisticsSingleRun       = 0;

if ( avgFitnessPerGenetationGraphs || ...
     paretoFrontFinalPopulationGraphAllRuns || ...
     objValueSummaryStatisticsAllRuns || ...
     finalPopulationObjectiveStatsAllRuns || ...
     parameterSummaryStatisticsAllRuns )
 
    %% Put all run data into a single file
    objectiveValuesPerGenerationData = [];
    finalPopulationObjectiveValues = [];
    numRuns = 20;
    for run = 1:numRuns
        % Construct the folder name
        if (run < 10)
            folderPrefix = 'generation_wise_run00';
        elseif ((run >= 10) && (run < 100))
            folderPrefix = 'generation_wise_run0';
        elseif (run >= 100)
            folderPrefix = 'generation_wise_run';
        end
        folderName = strcat(folderPrefix, int2str(run));
        fprtinf('Using folder %s...\n', folderName); % DEBUG. REMOVE!
        
        callingFolder = cd(folderName);
     
        % Get a list of all objective files in the folder
        objFiles = dir('*_obj.dat');
        objFileFilenames = {objFiles.name};
     
        if ( avgFitnessPerGenetationGraphs )
            for gen = 1:length(objFileFilenames)
                fileData = importdata(objFileFilenames{idx});
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
        if ( paretoFrontFinalPopulationGraphAllRuns )
            if (run == 1)
                fileData = importData('gen_0199_odj.dat');
                finalPopulationObjectiveValues = [finalPopulationObjectiveValues; fileData];
            end
        end
        
        cd(callingFolder);
    end
    
    if (avgFitnessPerGenerationGraphs)
        generations = [1:size(objectiveValuesPerGenerationData, 2)];
        figure;
        yvalues = ((-1) * objectiveValuesPerGenerationData(:,2)) / objectiveValuesPerGenerationData(:,1);
        plot(generations, yvalues);
        title('Average Objective 1 Value Per Generation', 'FontSize', 12);
        xlabel('Generation Number', 'FontSize', 12);
        ylabel('Objective 1 Value', 'FontSize', 12);
        figure;
        yvalues = ((-1) * objectiveValuesPerGenerationData(:,3)) / objectiveValuesPerGenerationData(:,1);
        plot(generations, yvalues);
        title('Average Objective 2 Value Per Generation', 'FontSize', 12);
        xlabel('Generation Number', 'FontSize', 12);
        ylabel('Objective 2 Value', 'FontSize', 12);
    end
    
    if ( paretoFrontFinalPopulationGraphAllRuns )
        finalPopulationObjectiveValues = finalPopulationObjectiveValues .* (-1);
        figure;
        plot(finalPopulationObjectiveValues(:,1), finalPopulationObjectiveValues(:, 2));
        title('Final Population Objective 2 versus Objective 1 Values for All Runs', 'FontSize', 12);
        xlabel('Objective 1 Value', 'FontSize', 12);
        ylabel('Objective 2 Value', 'FontSize', 12);
    end
end

if ( paretoFrontFinalPopulationGraphSingleRun || ...
     objValueSummaryStatisticsSingleRun || ...
     finalPopulationObjectiveStatsSingleRun || ...
     parameterSummaryStatisticsSingleRun )
end

