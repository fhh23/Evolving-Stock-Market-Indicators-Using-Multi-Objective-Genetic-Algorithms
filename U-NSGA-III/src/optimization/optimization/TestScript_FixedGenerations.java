/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization;

import Utilities.GeneralUtilities;
import Utilities.InputOutput;
import Utilities.Mathematics;
import Utilities.PerformanceMetrics;
import Utilities.RandomNumberGenerator;
import emo.DoubleAssignmentException;
import emo.Individual;
import emo.OptimizationProblem;
import emo.VirtualIndividual;
import engines.AbstractGeneticEngine;
import engines.UnifiedNSGA3Engine;
import evaluators.*;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import javax.xml.stream.XMLStreamException;
import net.sourceforge.jeval.EvaluationException;
import parsing.IndividualEvaluator;
import parsing.InvalidOptimizationProblemException;
import parsing.StaXParser;

/**
 *
 * @author toshiba
 */
public class TestScript_FixedGenerations extends TestScript {

    // According to the values of the following flags the script calculates
    // Hypervolume(for two objectives only), IGD and GD.
    static boolean calculateHV = false;
    static boolean calculateGD = false;
    static boolean calculateIGD = false;
    // Number of auto generated Pareto points (used to calculate GD and/or IGD)
    static int paretoOptimalFrontPointsCount = 10000;
    // Number of runs performed to take averages
    public static int runsCount = 1; // TODO: change the value to get multiple runs of the algorithm
    public static int currentRunIndex = 0;

    /**
     * @param args the command line arguments
     * @throws javax.xml.stream.XMLStreamException
     * @throws parsing.InvalidOptimizationProblemException
     * @throws net.sourceforge.jeval.EvaluationException
     * @throws java.io.IOException
     * @throws emo.DoubleAssignmentException
     */
    public static void main(String[] args)
            throws
            XMLStreamException,
            InvalidOptimizationProblemException,
            EvaluationException,
            IOException,
            DoubleAssignmentException {
        InputStream in = null;
        try {
            // ********** MODIFY ********** MODIFY ********** MODIFY **********
            // *                         MODIFY START                         *
            // ********** MODIFY ********** MODIFY ********** MODIFY **********
            // Read Problem
            URL url = UnifiedNSGA3Engine.class.getResource("../samples/rsi_financial.xml");
            in = url.openStream();
            OptimizationProblem optimizationProblem = StaXParser.readProblem(in);
            // Create Evaluator
            IndividualEvaluator individualEvaluator = new SingleIndicator_FinancialEvaluator();
            // Uncomment the following line only if you need the scaled version of the problem (Scaling is supported only for DTLZ1 and DTLZ2 until now(30Oct.2014))
            //((GeneralDTLZ1Evaluator)individualEvaluator).setScaled(true);
            // Create random seed
            double seed = 0.2;
            RandomNumberGenerator.setSeed(seed);
            // ********** MODIFY ********** MODIFY ********** MODIFY **********
            // *                          MODIFY END                          *
            // ********** MODIFY ********** MODIFY ********** MODIFY **********
            VirtualIndividual[] paretoFront = null;
            if (calculateGD || calculateIGD) {
                System.out.format("* Generating reference Pareto front...%n");
                paretoFront = individualEvaluator.getParetoFront(
                        optimizationProblem.objectives.length,
                        paretoOptimalFrontPointsCount);
                System.out.format("%d Pareto points generated.%n", paretoFront.length);
            }
            double[] REFERENCE_POINT = null, IDEAL_POINT = null;
            if (calculateHV) {
                IDEAL_POINT = individualEvaluator.getIdealPoint();
                REFERENCE_POINT = individualEvaluator.getReferencePoint();
                // Introduce marginal displacements in ref. point position (used to
                // calculate HV) away from Nadir point position in both directions.
                for (int i = 0; i < REFERENCE_POINT.length; i++) {
                    REFERENCE_POINT[i] += REFERENCE_POINT[i] * 0.01;
                }
            }

            // Overriding parameters settings (for testing purposes)
            // etaC Parameter
            int etaCStart = 30;
            int etaCEnd = 30;
            // etaM Parameter
            int etaMStart = 20;
            int etaMEnd = 20;
            // Population Parameter
            int popSizeStart = 200;
            int popSizeEnd = 200;
            // Generations Parameter
            int genCountStart = 200;
            int genCountEnd = 200;
            // Parameters Steps
            int etaCStep = 10;
            int etaMStep = 10;
            int popSizeStep = 20;
            int genCountStep = 10;
            // Set your epsilon for epsilon-domination
            double epsilon = 0.0;

            //List<String> hypervolumesList = new ArrayList<String>();
            for (int etaC = etaCStart; etaC <= etaCEnd; etaC += etaCStep) {
                optimizationProblem.setRealCrossoverDistIndex(etaC);
                for (int etaM = etaMStart; etaM <= etaMEnd; etaM += etaMStep) {
                    optimizationProblem.setRealMutationDistIndex(etaM);
                    for (int popSize = popSizeStart; popSize <= popSizeEnd; popSize += popSizeStep) {
                        optimizationProblem.setPopulationSize(popSize);
                        // The following line is very important. It overrides the number of steps in
                        // order to guarantee a number of directions equal to the population size
                        // ONLY IN TWO OBJECTIVES.. ONLY IN TWO OBJECTIVES.. ONLY IN TWO OBJECTIVES..
                        if (optimizationProblem.objectives.length == 2) {
                            optimizationProblem.setSteps(popSize - 1);
                        }
                        for (int genCount = genCountStart; genCount <= genCountEnd; genCount += genCountStep) {
                            optimizationProblem.setGenerationsCount(genCount);
                            // ********** MODIFY ********** MODIFY ********** MODIFY **********
                            // *                          MODIFY START                          *
                            // ********** MODIFY ********** MODIFY ********** MODIFY **********
                            // Create the engine
                            AbstractGeneticEngine geneticEngine = new UnifiedNSGA3Engine(optimizationProblem, individualEvaluator);
                            //int[] divisions = {3,2};
                            //AbstractGeneticEngine geneticEngine = new UnifiedNSGA3Engine(optimizationProblem, individualEvaluator, divisions);
                            //AbstractGeneticEngine geneticEngine = new UnifiedNSGA3Engine(optimizationProblem, individualEvaluator, "D:/Extra/R2.ref");
                            // ********** MODIFY ********** MODIFY ********** MODIFY **********
                            // *                          MODIFY END                          *
                            // ********** MODIFY ********** MODIFY ********** MODIFY **********

                            String outputDir = outDir + String.format("%s-%03d-%03d-P%04d-G%04d/%s/",
                                    optimizationProblem.getProblemID(),
                                    optimizationProblem.getRealCrossoverDistIndex(),
                                    optimizationProblem.getRealMutationDistIndex(),
                                    optimizationProblem.getPopulationSize(),
                                    optimizationProblem.getGenerationsCount(),
                                    geneticEngine.getAlgorithmName()
                            );
                            // Make directories
                            File dirs = new File(outputDir);
                            if (!dirs.exists()) {
                                dirs.mkdirs();
                            }
                            // perform several runs (runsCount) and collect metrics
                            double[] hyperVolume = null, gd = null, igd = null;
                            if (calculateHV) {
                                hyperVolume = new double[runsCount];
                            }
                            if (calculateGD) {
                                gd = new double[runsCount];
                            }
                            if (calculateIGD) {
                                igd = new double[runsCount];
                            }
                            int totalEvaluationsCount = 0;
                            for (int runIndex = 0; runIndex < runsCount; runIndex++) {
                                System.out.println("*****************");
                                System.out.format("    Run(%03d)  %n", runIndex);
                                System.out.println("*****************");
                                // Create a genetic engine for the problem
                                Individual[] finalPopulation = geneticEngine.start(outputDir, runIndex, epsilon);

                                // Generate output data file
                                String dataFileName = String.format(outputDir + "%s-G%03d-P%03d-run%03d-data.dat",
                                        optimizationProblem.getProblemID(),
                                        optimizationProblem.getGenerationsCount(),
                                        optimizationProblem.getPopulationSize(),
                                        runIndex);
                                InputOutput.dumpPopulation("Final Population",
                                        optimizationProblem,
                                        finalPopulation,
                                        dataFileName);
                                if (optimizationProblem.objectives.length == 2 || optimizationProblem.objectives.length == 3) {
                                    // Prepare Matlab script file name
                                    String matlabScriptFilePath
                                            = String.format("Matlab_%s_G%03d_P%03d_run%03d_data.m",
                                                    optimizationProblem.getProblemID(),
                                                    optimizationProblem.getGenerationsCount(),
                                                    optimizationProblem.getPopulationSize(),
                                                    runIndex);
                                    matlabScriptFilePath = GeneralUtilities.replaceDashesWithUnderscores(matlabScriptFilePath);
                                    matlabScriptFilePath = outputDir + matlabScriptFilePath;
                                    // Dump the Matlab plotting script
                                    if (optimizationProblem.objectives.length == 2) {
                                        //InputOutput.dumpMatlabPlottinScriptFor2dPoints(matlabScriptFilePath, individualEvaluator, finalPopulation);
                                    } else {
                                        //InputOutput.dumpMatlabPlottinScriptFor3dPoints(matlabScriptFilePath, individualEvaluator, finalPopulation);
                                    }
                                }
                                // Reset the number of function evaluations to start
                                // counting from Zero again in the next iteration
                                individualEvaluator.resetEvaluationsCount();
                                // Hypervolume
                                if (optimizationProblem.objectives.length == 2) {
                                    if (calculateHV) {
                                        hyperVolume[runIndex]
                                                = PerformanceMetrics.calculateHyperVolumeForTwoObjectivesOnly(
                                                        geneticEngine,
                                                        finalPopulation,
                                                        REFERENCE_POINT,
                                                        IDEAL_POINT,
                                                        epsilon);
                                    }
                                }
                                // GD
                                if (calculateGD) {
                                    gd[runIndex] = PerformanceMetrics.calculateGenerationalDistance(
                                            optimizationProblem.objectives.length,
                                            finalPopulation,
                                            paretoFront,
                                            2);
                                }
                                // IGD
                                if (calculateIGD) {
                                    igd[runIndex] = PerformanceMetrics.calculateInvertedGenerationalDistance(
                                            optimizationProblem.objectives.length,
                                            finalPopulation,
                                            paretoFront,
                                            2);
                                }
                                // Increment run index
                                currentRunIndex++;
                            }
                            System.out.format("Total Evaluations Count = %d%n", totalEvaluationsCount);
                            System.out.format("Avg Evaluations Count(per run) = %d%n", (totalEvaluationsCount / runsCount));
                            if (calculateHV || calculateGD || calculateIGD) {
                                String metricsFileName = String.format(outputDir + "%s-G%03d-P%03d-metrics.dat",
                                        optimizationProblem.getProblemID(),
                                        optimizationProblem.getGenerationsCount(),
                                        optimizationProblem.getPopulationSize());
                                InputOutput.dumpPerformanceMetrics(
                                        optimizationProblem,
                                        hyperVolume,
                                        gd,
                                        igd,
                                        metricsFileName);
                            }
                        }
                    }
                }
            }
        } finally {
            if (in != null) {
                in.close();
            }
        }
    }
}
