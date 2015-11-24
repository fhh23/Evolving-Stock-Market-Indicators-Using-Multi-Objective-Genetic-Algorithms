/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package engines;

import Utilities.GeneralUtilities;
import Utilities.InputOutput;
import Utilities.Mathematics;
import Utilities.RandomNumberGenerator;
import emo.DoubleAssignmentException;
import emo.Individual;
import emo.IndividualsSet;
import emo.OptimizationProblem;
import emo.OptimizationUtilities;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import net.sourceforge.jeval.EvaluationException;
import optimization.TestScript;
import optimization.TestScript_FixedGenerations;
import optimization.TestScript_SingleObjective;
import parsing.IndividualEvaluator;
import reference_directions.ReferenceDirection;

/**
 *
 * @author toshiba
 */
public abstract class AbstractGeneticEngine {

    public static final double MIN_DOUBLE_VALUE = Math.pow(10, -6);
    public static final double MAX_DOUBLE_VALUE = Math.pow(10, 12);
    public final OptimizationProblem optimizationProblem;
    public final IndividualEvaluator individualEvaluator;
    //Individual[] population;
    public final static boolean DEBUG_ALL = false;
    public final static boolean DEBUG_REFERENCE_DIRECTIONS = false;
    public final static boolean DEBUG_POPULATIONS = true;
    public final static boolean DEBUG_RANKING = false;
    public final static boolean DEBUG_IDEAL_POINT = false;
    public final static boolean DEBUG_TRANSLATION = false;
    public final static boolean DEBUG_INTERCEPTS = false;
    public final static boolean DEBUG_ASSOCIATION = false;
    public final static boolean EXTREME_POINTS_DEEP_DEBUG = false;

    public final static boolean DUMP_ALL_GENERATIONS_DECISION_SPACE = true;
    public final static boolean DUMP_ALL_GENERATIONS_OBJECTIVE_SPACE = true;
    public final static boolean DUMP_ALL_GENERATIONS_MATLAB_SCRIPTS = true;
    // If DUMP_ALL_GENERATIONS_NORMALIZED_MATLAB_SCRIPTS is true if and only if
    // DUMP_ALL_GENERATIONS_OBJECTIVE_SPACE is also true,
    // otherwise DUMP_ALL_GENERATIONS_NORMALIZED_MATLAB_SCRIPTS value will be
    // ignored.
    public final static boolean DUMP_ALL_GENERATIONS_NORMALIZED_MATLAB_SCRIPTS
            = false;
    // If DUMP_ANIMATED_MATLAB_SCRIPT is true if and only if
    // DUMP_ALL_GENERATIONS_OBJECTIVE_SPACE is also true,
    // otherwise DUMP_ANIMATED_MATLAB_SCRIPT value will be ignored.
    public final static boolean DUMP_ANIMATED_MATLAB_SCRIPT
            = false;

    public AbstractGeneticEngine(
            OptimizationProblem optimizationProblem,
            IndividualEvaluator individualEvaluator) throws EvaluationException {
        this.optimizationProblem = optimizationProblem;
        this.individualEvaluator = individualEvaluator;
    }

    public Individual[] generateInitialPopulation() throws EvaluationException {
        Individual[] population = new Individual[optimizationProblem.getPopulationSize()];
        for (int i = 0; i < population.length; i++) {
            population[i] = new Individual(optimizationProblem, individualEvaluator);
        }
        return population;
    }

    protected Individual tournamentSelect(IndividualsSet subset) {
        Individual individual1 = subset.getIndividual1();
        Individual individual2 = subset.getIndividual2();
        // If the problem is constrained and at least one of the
        // individuals under investigation is infeasible return the feasible
        // one (which is the dominating individual).
        if (/*individual1.dominates(individual2)*/individual1.getRank() < individual2.getRank()) {
            // Individual-1 dominates Individual-2
            return individual1;
        } else if (/*individual2.dominates(individual1)*/individual1.getRank() > individual2.getRank()) {
            // Individual-2 dominates Individual-1
            return individual2;
        } else if (Mathematics.compare(individual1.getNsga2crowdingDistance(), individual2.getNsga2crowdingDistance()) == 1) {
            // Both Individual-1 & Individual-2 are non-dominated with respect 
            // to each other but Individual-1 has a larger crowding distance.
            return individual1;
        } else if (Mathematics.compare(individual1.getNsga2crowdingDistance(), individual2.getNsga2crowdingDistance()) == -1) {
            // Both Individual-1 & Individual-2 are non-dominated with respect 
            // to each other but Individual-2 has a larger crowding distance.
            return individual2;
        } else {
            // Both Individual-1 & Individual-2 are non-dominated with respect 
            // to each other and have the same crodingDistance (a rare case).
            // In this case select one of the two individuals randomly.
            if (RandomNumberGenerator.randomperc() <= 0.5) {
                return individual1;
            } else {
                return individual2;
            }
        }
    }

    protected Individual[] getOffspringPopulation(Individual[] oldPopulation) throws EvaluationException {
        Individual[] newPopulation = new Individual[optimizationProblem.getPopulationSize()];
        int[] a1 = new int[optimizationProblem.getPopulationSize()];
        int[] a2 = new int[optimizationProblem.getPopulationSize()];
        int temp;
        int i;
        int rand;
        Individual parent1, parent2;
        IndividualsSet childrenSet;
        for (i = 0; i < optimizationProblem.getPopulationSize(); i++) {
            a1[i] = a2[i] = i;
        }
        for (i = 0; i < optimizationProblem.getPopulationSize(); i++) {
            rand = RandomNumberGenerator.rnd(i, optimizationProblem.getPopulationSize() - 1);
            temp = a1[rand];
            a1[rand] = a1[i];
            a1[i] = temp;
            rand = RandomNumberGenerator.rnd(i, optimizationProblem.getPopulationSize() - 1);
            temp = a2[rand];
            a2[rand] = a2[i];
            a2[i] = temp;
        }
        for (i = 0; i < optimizationProblem.getPopulationSize(); i += 4) {
            parent1 = tournamentSelect(new IndividualsSet(oldPopulation[a1[i]], oldPopulation[a1[i + 1]]));
            parent2 = tournamentSelect(new IndividualsSet(oldPopulation[a1[i + 2]], oldPopulation[a1[i + 3]]));
            childrenSet = crossover(new IndividualsSet(parent1, parent2));
            newPopulation[i] = childrenSet.getIndividual1();
            newPopulation[i + 1] = childrenSet.getIndividual2();

            parent1 = tournamentSelect(new IndividualsSet(oldPopulation[a2[i]], oldPopulation[a2[i + 1]]));
            parent2 = tournamentSelect(new IndividualsSet(oldPopulation[a2[i + 2]], oldPopulation[a2[i + 3]]));
            childrenSet = crossover(new IndividualsSet(parent1, parent2));
            newPopulation[i + 2] = childrenSet.getIndividual1();
            newPopulation[i + 3] = childrenSet.getIndividual2();
        }
        // Nullify (set as null) all the reference directions of the offspring.
        // This is a protective step to make the code run into an exception
        // if an individual from the offspring was forgotten in the association
        // step. Remember that offsprings are results of crossovers and mutation
        // so most probably the direction they have (copied from their parents)
        // are no longer valid according to their new objective values.
        for (Individual individual : newPopulation) {
            individual.setReferenceDirection(null);
        }
        return newPopulation;
    }

    protected IndividualsSet crossover(IndividualsSet parents) throws EvaluationException {
        IndividualsSet children = new IndividualsSet(
                new Individual(optimizationProblem, parents.getIndividual1(), individualEvaluator),
                new Individual(optimizationProblem, parents.getIndividual2(), individualEvaluator));
        //realCrossover_bad(parents, children);
        realCrossover_good(parents, children);
        binaryCrossover(parents, children);
        return children;
    }

    //private void crossover(IndividualsSet parents, IndividualsSet children) {
    //    binaryCrossover(parents, children);
    //    realCrossover(parents, children);
    //}
    public void binaryCrossover(IndividualsSet parents, IndividualsSet children) {
        Individual parent1 = parents.getIndividual1();
        Individual parent2 = parents.getIndividual2();
        Individual child1 = children.getIndividual1();
        Individual child2 = children.getIndividual2();
        // Loop over all available binary variables
        int binaryVarCount = parent1.binary.length;
        for (int i = 0; i < binaryVarCount; i++) {
            int bitsCount = parent1.binary[i].getSpecs().getNumberOfBits();
            // Perform binary crossover according to the user defined probability
            if (/*RandomNumberGenerator.nextDouble()*/RandomNumberGenerator.randomperc() < optimizationProblem.getBinaryCrossoverProbability()) {
                // Randomly pick two cut positions
                //int startCutIndex = RandomNumberGenerator.nextIntegerWithin(0, bitsCount);
                //int endCutIndex = RandomNumberGenerator.nextIntegerWithin(0, bitsCount);
                int startCutIndex = RandomNumberGenerator.rnd(0, bitsCount);
                int endCutIndex = RandomNumberGenerator.rnd(0, bitsCount);
                // Reorder cut points (ascendingly) if required
                if (startCutIndex > endCutIndex) {
                    int temp = startCutIndex;
                    startCutIndex = endCutIndex;
                    endCutIndex = temp;
                }
                // Binary crossover
                for (int bitIndex = 0; bitIndex < bitsCount; bitIndex++) {
                    if (bitIndex < startCutIndex || bitIndex > endCutIndex) {
                        // Copy bits from P1 to C1
                        if (parent1.binary[i].getValueOfBit(bitIndex) == 0) {
                            child1.binary[i].setBitToZero(bitIndex);
                        } else {
                            child1.binary[i].setBitToOne(bitIndex);
                        }
                        // Copy bits from P2 to C2
                        if (parent2.binary[i].getValueOfBit(bitIndex) == 0) {
                            child2.binary[i].setBitToZero(bitIndex);
                        } else {
                            child2.binary[i].setBitToOne(bitIndex);
                        }
                    } else {
                        // Copy bits from P1 to C2
                        if (parent1.binary[i].getValueOfBit(bitIndex) == 0) {
                            child2.binary[i].setBitToZero(bitIndex);
                        } else {
                            child2.binary[i].setBitToOne(bitIndex);
                        }
                        // Copy bits from P2 to C1
                        if (parent2.binary[i].getValueOfBit(bitIndex) == 0) {
                            child1.binary[i].setBitToZero(bitIndex);
                        } else {
                            child1.binary[i].setBitToOne(bitIndex);
                        }
                    }
                }
            } else {
                // Copy all bits from P1 to C1 & from P2 to C2
                for (int bitIndex = 0; bitIndex < bitsCount; bitIndex++) {
                    if (parent1.binary[i].getValueOfBit(bitIndex) == 0) {
                        child1.binary[i].setBitToZero(bitIndex);
                    } else {
                        child1.binary[i].setBitToOne(bitIndex);
                    }
                    if (parent2.binary[i].getValueOfBit(bitIndex) == 0) {
                        child2.binary[i].setBitToZero(bitIndex);
                    } else {
                        child2.binary[i].setBitToOne(bitIndex);
                    }
                }
            }
        }
    }

    public void realCrossover_bad(IndividualsSet parents, IndividualsSet children) {
        Individual parent1 = parents.getIndividual1();
        Individual parent2 = parents.getIndividual2();
        Individual child1 = children.getIndividual1();
        Individual child2 = children.getIndividual2();
        // Get the number of real variables
        int realVarCount = parent1.real.length;
        // Perform real crossover according to the user defined probability
        double rrnndd = RandomNumberGenerator.randomperc();
        if (rrnndd < optimizationProblem.getRealCrossoverProbability()) {
            for (int i = 0; i < realVarCount; i++) {

                // Get the specs object of that real variable
                OptimizationProblem.RealVariableSpecs specs = null;
                int realVarIndex = -1;
                for (int j = 0; j < optimizationProblem.variablesSpecs.length; j++) {
                    if (optimizationProblem.variablesSpecs[j] instanceof OptimizationProblem.RealVariableSpecs) {
                        realVarIndex++;
                        if (realVarIndex == i) {
                            specs = (OptimizationProblem.RealVariableSpecs) optimizationProblem.variablesSpecs[j];
                        }
                    }
                }

                // ******************************************************
                // IMPORTANT FROM HERE
                // ******************************************************
                // Perform real crossover per variable with 50% probability
                if (RandomNumberGenerator.randomperc() < 0.5) {
                    if (Math.abs(parent1.real[i] - parent2.real[i]) < Mathematics.EPSILON) {
                        // Copy the current real variable values from P1 to C1 & from P2 to C2
                        child1.real[i] = parent1.real[i];
                        child2.real[i] = parent2.real[i];
                    } else {
                        // Perform real crossover on the current real variable
                        double y1 = Math.min(parent1.real[i], parent2.real[i]);
                        double y2 = Math.max(parent1.real[i], parent2.real[i]);
                        double yl = specs.getMinValue();
                        double yu = specs.getMaxValue();
                        double alphaL = 2 - Math.pow((1 + 2 * (y1 - yl) / (y2 - y1)), -(optimizationProblem.getRealCrossoverDistIndex() + 1));
                        double alphaR = 2 - Math.pow((1 + 2 * (yu - y2) / (y2 - y1)), -(optimizationProblem.getRealCrossoverDistIndex() + 1));
                        //double uL = RandomNumberGenerator.nextDouble();
                        //double uR = RandomNumberGenerator.nextDouble();
                        double uL = RandomNumberGenerator.randomperc();
                        double uR = RandomNumberGenerator.randomperc();
                        double betaL;
                        if (uL <= 1 / alphaL) {
                            betaL = Math.pow((uL * alphaL), 1.0 / (optimizationProblem.getRealCrossoverDistIndex() + 1));
                        } else {
                            betaL = Math.pow((1.0 / (2 - uL * alphaL)), 1.0 / (optimizationProblem.getRealCrossoverDistIndex() + 1));
                        }
                        child1.real[i] = (y1 + y2 - betaL * (y2 - y1)) / 2;
                        double betaR;
                        if (uR <= 1 / alphaR) {
                            betaR = Math.pow((uR * alphaR), 1.0 / (optimizationProblem.getRealCrossoverDistIndex() + 1));
                        } else {
                            betaR = Math.pow((1.0 / (2 - uR * alphaR)), 1.0 / (optimizationProblem.getRealCrossoverDistIndex() + 1));
                        }
                        child2.real[i] = (y1 + y2 + betaR * (y2 - y1)) / 2;
                        // If children's values went beyond bounds pull them back in
                        // Child1
                        if (child1.real[i] < yl) {
                            child1.real[i] = yl;
                        } else if (child1.real[i] > yu) {
                            child1.real[i] = yu;
                        }
                        // Child2
                        if (child2.real[i] < yl) {
                            child2.real[i] = yl;
                        } else if (child2.real[i] > yu) {
                            child2.real[i] = yu;
                        }
                        // With 50% probability, swap values between the two children
                        if (RandomNumberGenerator.randomperc() <= 0.5) {
                            double temp = child1.real[i];
                            child1.real[i] = child2.real[i];
                            child2.real[i] = temp;
                        }
                    }
                    continue;
                }
                // Copy the current real variable values from P1 to C1 & from P2 to C2
                child1.real[i] = parent1.real[i];
                child2.real[i] = parent2.real[i];
                // ******************************************************
                // IMPORTANT UNTIL HERE
                // ******************************************************

            }
        } else {
            // Copy all real variables values from P1 to C1 & from P2 to C2
            for (int i = 0; i < realVarCount; i++) {
                child1.real[i] = parent1.real[i];
                child2.real[i] = parent2.real[i];
            }
        }
    }

    ///*
    public void realCrossover_good(IndividualsSet parents, IndividualsSet children) throws EvaluationException {
        Individual parent1 = parents.getIndividual1();
        Individual parent2 = parents.getIndividual2();
        Individual child1 = children.getIndividual1();
        Individual child2 = children.getIndividual2();
        // Get the number of real variables
        int realVarCount = parent1.real.length;

        double par1Value, par2Value, child1Value, child2Value, betaq, beta, alpha;
        double y1, y2, yu, yl, expp;

        //Check Whether the cross-over to be performed
        double rnd = RandomNumberGenerator.randomperc();
        if (rnd <= optimizationProblem.getRealCrossoverProbability()) {

            //Loop over no of variables
            for (int j = 0; j < realVarCount; j++) {
                par1Value = parent1.real[j];
                par2Value = parent2.real[j];

                // Get the specs object of that real variable
                OptimizationProblem.RealVariableSpecs specs = null;
                int realVarIndex = -1;
                for (int k = 0; k < optimizationProblem.variablesSpecs.length; k++) {
                    if (optimizationProblem.variablesSpecs[k] instanceof OptimizationProblem.RealVariableSpecs) {
                        realVarIndex++;
                        if (realVarIndex == j) {
                            specs = (OptimizationProblem.RealVariableSpecs) optimizationProblem.variablesSpecs[k];
                        }
                    }
                }

                yl = specs.getMinValue();
                yu = specs.getMaxValue();

                rnd = RandomNumberGenerator.randomperc();

                // Check whether variable is selected or not
                if (rnd <= 0.5) {
                    //Variable selected
                    if (Math.abs(par1Value - par2Value) > 0.000001) // changed by Deb (31/10/01)
                    {
                        if (par2Value > par1Value) {
                            y2 = par2Value;
                            y1 = par1Value;
                        } else {
                            y2 = par1Value;
                            y1 = par2Value;
                        }

                        //Find beta value
                        if ((y1 - yl) > (yu - y2)) {
                            beta = 1 + (2 * (yu - y2) / (y2 - y1));
                            //printf("beta = %f\n",beta);
                        } else {
                            beta = 1 + (2 * (y1 - yl) / (y2 - y1));
                            //printf("beta = %f\n",beta);
                        }

                        //Find alpha
                        expp = optimizationProblem.getRealCrossoverDistIndex() + 1.0;

                        beta = 1.0 / beta;

                        alpha = 2.0 - Math.pow(beta, expp);

                        if (alpha < 0.0) {
                            System.out.format("ERRRROR %f %f %f", alpha, par1Value, par2Value);
                            System.exit(-1);
                        }

                        rnd = RandomNumberGenerator.randomperc();

                        if (rnd <= 1.0 / alpha) {
                            alpha = alpha * rnd;
                            expp = 1.0 / (optimizationProblem.getRealCrossoverDistIndex() + 1.0);
                            betaq = Math.pow(alpha, expp);
                        } else {
                            alpha = alpha * rnd;
                            alpha = 1.0 / (2.0 - alpha);
                            expp = 1.0 / (optimizationProblem.getRealCrossoverDistIndex() + 1.0);
                            if (alpha < 0.0) {
                                System.out.println("ERRRORRR");
                                System.exit(-1);
                            }
                            betaq = Math.pow(alpha, expp);
                        }

                        //Generating two children
                        child1Value = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                        child2Value = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

                    } else {

                        betaq = 1.0;
                        y1 = par1Value;
                        y2 = par2Value;

                        //Generation two children
                        child1Value = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                        child2Value = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

                    }
                    // added by deb (31/10/01)
                    if (child1Value < yl) {
                        child1Value = yl;
                    }
                    if (child1Value > yu) {
                        child1Value = yu;
                    }
                    if (child2Value < yl) {
                        child2Value = yl;
                    }
                    if (child2Value > yu) {
                        child2Value = yu;
                    }
                } else {
                    //Copying the children to parents
                    child1Value = par1Value;
                    child2Value = par2Value;
                }
                child1.real[j] = child1Value;
                child2.real[j] = child2Value;
            }
        } else {
            for (int j = 0; j < realVarCount; j++) {
                par1Value = parent1.real[j];
                par2Value = parent2.real[j];
                child1Value = par1Value;
                child2Value = par2Value;
                child1.real[j] = child1Value;
                child2.real[j] = child2Value;
            }
        }
    }
    //*/

    public Individual[] merge(Individual[] oldPopulation, Individual[] newPopulation) {
        Individual[] combinedPopulation
                = new Individual[oldPopulation.length + newPopulation.length];
        System.arraycopy(oldPopulation, 0, combinedPopulation, 0, oldPopulation.length);
        System.arraycopy(newPopulation, 0, combinedPopulation, oldPopulation.length, newPopulation.length);
        return combinedPopulation;
    }

    int checkDominance(Individual individual1, Individual individual2, double epsilon) {
        if (individual1.dominates(individual2, epsilon)) {
            return 1;
        } else if (individual2.dominates(individual1, epsilon)) {
            return -1;
        } else {
            return 0;
        }
    }

    int count = 0;

    public List<List<Individual>> assign_rank(Individual[] individuals, double epsilon) {
        int flag;
        int i;
        int end;
        int front_size;
        int rank = 1;
        ArrayList<Individual> orig = new ArrayList<Individual>();
        ArrayList<Individual> cur = new ArrayList<Individual>();
        int temp1;
        int temp2;
        for (i = 0; i < individuals.length; i++) {
            orig.add(individuals[i]);
        }
        do {
            if (orig.size() == 1) {
                orig.get(0).setRank(rank);
                orig.get(0).validRankValue = true;
                break;
            }
            cur.add(0, orig.remove(0));
            temp1 = 0;
            temp2 = 0;
            front_size = 1;
            do {
                temp2 = 0;
                do {
                    end = 0;
                    flag = checkDominance(orig.get(temp1), cur.get(temp2), epsilon);
                    if (flag == 1) {
                        orig.add(0, cur.remove(temp2));
                        temp1++;
                        front_size--;
                    } else if (flag == -1) {
                        end = 1;
                    } else {
                        temp2++;
                    }
                } while (end != 1 && temp2 < cur.size());
                if (flag == 0 || flag == 1) {
                    cur.add(0, orig.get(temp1));
                    temp2++;
                    front_size++;
                    orig.remove(temp1);
                }
                if (flag == -1) {
                    temp1++;
                }
            } while (temp1 < orig.size());
            for (Individual individual : cur) {
                individual.setRank(rank);
                individual.validRankValue = true;
            }
            cur.clear();
            rank++;
        } while (!orig.isEmpty());
        // Prepare return List of Lists
        List<List<Individual>> fronts = new ArrayList<List<Individual>>();
        // Get Max rank
        int maxRank = 1;
        for (Individual individual : individuals) {
            if (individual.getRank() > maxRank) {
                maxRank = individual.getRank();
            }
        }
        // Initialize the fronts lists
        for (int f = 0; f < maxRank; f++) {
            fronts.add(new ArrayList<Individual>());
        }
        // Add each individual to its corresponding rank
        for (Individual individual : individuals) {
            fronts.get(individual.getRank() - 1).add(individual);
        }
        // Return the final list of fronts
        return fronts;
    }

    public double[] getInitialIdealPoint(Individual[] individuals) {
        double[] idealPoint = new double[optimizationProblem.objectives.length];
        // Let the ideal point be the first point in the population (just as a start)
        for (int i = 0; i < idealPoint.length; i++) {
            idealPoint[i] = /*individuals[0].getObjective(i)*/ Math.pow(10, 12);
        }
        // Update the value of each objective in the population if a smaller value
        // is found in any subsequent population member
        for (int i = 0; i < individuals.length; i++) {
            if (Mathematics.compare(individuals[i].getTotalConstraintViolation(), 0) == 0) {
                for (int j = 0; j < idealPoint.length; j++) {
                    if (individuals[i].getObjective(j) < idealPoint[j]) {
                        idealPoint[j] = individuals[i].getObjective(j);
                    }
                }
            }
        }
        // Return the ideal point
        return idealPoint;
    }

    public void mutate(Individual[] individuals) {
        for (int i = 0; i < optimizationProblem.getPopulationSize(); i++) {
            mutation_ind(individuals[i]);
        }
    }

    /* Function to perform mutation of an individual */
    void mutation_ind(Individual individual) {
        if (individual.real.length != 0) {
            realMutateIndividual(individual);
        }
        if (individual.binary.length != 0) {
            binaryMutateIndividual(individual);
        }
    }

    /* Routine for binary mutation of an individual */
    public void binaryMutateIndividual(Individual individual) {
        for (int j = 0; j < individual.binary.length; j++) {
            for (int k = 0; k < individual.binary[j].getSpecs().getNumberOfBits(); k++) {
                double prob = RandomNumberGenerator.randomperc();
                if (prob <= optimizationProblem.getBinaryMutationProbabilty()) {
                    if (individual.binary[j].getValueOfBit(k) == 0) {
                        individual.binary[j].setBitToOne(k);
                    } else {
                        individual.binary[j].setBitToZero(k);
                    }
                }
            }
        }
    }

    //public static List<Double> originalY = new ArrayList<Double>();
    //public static List<Double> newY = new ArrayList<Double>();
    /* Routine for real polynomial mutation of an individual */
    public void realMutateIndividual(Individual individual) {
        int j;
        double rnd, delta1, delta2, mut_pow, deltaq;
        double y, yl, yu, val, xy;
        for (j = 0; j < individual.real.length; j++) {
            if (RandomNumberGenerator.randomperc() <= optimizationProblem.getRealMutationProbability()) {

                // Get the specs object of that real variable
                OptimizationProblem.RealVariableSpecs specs = null;
                int realVarIndex = -1;
                for (int k = 0; k < optimizationProblem.variablesSpecs.length; k++) {
                    if (optimizationProblem.variablesSpecs[k] instanceof OptimizationProblem.RealVariableSpecs) {
                        realVarIndex++;
                        if (realVarIndex == k) {
                            specs = (OptimizationProblem.RealVariableSpecs) optimizationProblem.variablesSpecs[j];
                        }
                    }
                }

                // ***********************************************************
                // IMPORTANT FROM HERE
                // ***********************************************************
                y = individual.real[j];
                //originalY.add(y);
                yl = specs.getMinValue();
                yu = specs.getMaxValue();
                delta1 = (y - yl) / (yu - yl);
                delta2 = (yu - y) / (yu - yl);
                rnd = RandomNumberGenerator.randomperc();
                mut_pow = 1.0 / (optimizationProblem.getRealMutationDistIndex() + 1.0);
                if (rnd <= 0.5) {
                    xy = 1.0 - delta1;
                    val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (Math.pow(xy, (optimizationProblem.getRealMutationDistIndex() + 1.0)));
                    deltaq = Math.pow(val, mut_pow) - 1.0;
                } else {
                    xy = 1.0 - delta2;
                    val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (Math.pow(xy, (optimizationProblem.getRealMutationDistIndex() + 1.0)));
                    deltaq = 1.0 - (Math.pow(val, mut_pow));
                }
                y = y + deltaq * (yu - yl);
                if (y < yl) {
                    y = yl;
                }
                if (y > yu) {
                    y = yu;
                }
                individual.real[j] = y;
                //newY.add(y);
                // ***********************************************************
                // IMPORTANT UNTIL HERE
                // ***********************************************************

            }
        }
    }

    public double[] getUpdatedIdealPoint(Individual[] mergedPopulation, double[] idealPoint) {
        double[] updatedIdealPoint = new double[optimizationProblem.objectives.length];
        for (int i = 0; i < updatedIdealPoint.length; i++) {
            updatedIdealPoint[i] = MAX_DOUBLE_VALUE;
        }
        for (int i = 0; i < optimizationProblem.objectives.length; i++) {
            for (Individual individual : mergedPopulation) {
                if (Mathematics.compare(individual.getTotalConstraintViolation(), 0) == 0) {
                    double minValue;
                    if (individual.getObjective(i) < idealPoint[i]) {
                        minValue = individual.getObjective(i);
                    } else {
                        minValue = idealPoint[i];
                    }
                    if (minValue < updatedIdealPoint[i]) {
                        updatedIdealPoint[i] = minValue;
                    }
                }
            }
        }
        return updatedIdealPoint;
    }

    public int getRemainingCount(List<List<Individual>> fronts) {
        int individualsCount = 0;
        // Get the last front Fl index
        int lastFrontIndex = -1;
        while (individualsCount < optimizationProblem.getPopulationSize()) {
            lastFrontIndex++;
            individualsCount += fronts.get(lastFrontIndex).size();
        }
        // Determine the number of individuals required to complete the population (REMAINING)
        int remaining;
        if (individualsCount == optimizationProblem.getPopulationSize()) {
            remaining = 0;
        } else {
            individualsCount -= fronts.get(lastFrontIndex).size();
            remaining = optimizationProblem.getPopulationSize() - individualsCount;
        }
        return remaining;
    }

    public int getLimitingFrontIndex(List<List<Individual>> fronts) {
        int individualsCount = 0;
        // Get the last front Fl index
        int lastFrontIndex = -1;
        while (individualsCount <= optimizationProblem.getPopulationSize()) {
            lastFrontIndex++;
            individualsCount += fronts.get(lastFrontIndex).size();
        }
        return lastFrontIndex;
    }

    /**
     * This method fills the destination population with individuals from the
     * source population in two steps: (1)STEP-1: Add any individual in
     * sourcePopulation whose front will be completely accommodated in the
     * destinationPopulation. (2)Step-2: Add all the individuals of
     * lastFrontSubset to the destinationPopulation. (Note: this method is not
     * responsible for selecting individuals from the last front to be
     * accommodated. It is assumed that selection has already been made and that
     * the individuals of lastFrontSubset are the result of that selection).
     * (Note: any previous items in destinationPopulation will be overwritten)
     *
     * @param destinationPopulation population at which selected members are
     * stored
     * @param sourcePopulation population from which top-ranked members are
     * selected
     * @param lastFrontSubset contains selected individuals from the last front
     * to be accommodated
     * @param lastFrontIndex the index of the limiting front (the front before
     * @throws EvaluationException
     */
    public void reFillPopulation(
            Individual[] destinationPopulation,
            Individual[] sourcePopulation,
            Individual[] lastFrontSubset,
            int lastFrontIndex) throws EvaluationException {
        int index = 0;
        for (Individual individual : sourcePopulation) {
            if (individual.getRank() - 1 < lastFrontIndex) { // Remember that ranks starts at 1 not at Zero
                destinationPopulation[index] = individual;
                index++;
            }
        }
        for (int i = 0; i < lastFrontSubset.length; i++) {
            destinationPopulation[index] = /*new Individual(optimizationProblem,*/ lastFrontSubset[i]/*, individualEvaluator)*/;
            index++;
        }
    }

    /**
     * This method fills the destination fronts with individuals from the fronts
     * in order i.e. individuals of smaller-rank fronts are taken before
     * individuals of larger-rank fronts. Individuals are taken from the first
     * front as much as needed and in order of their presence (storage) in that
     * front. It is worth noting that there is no preference information taken
     * into account when selecting individuals from the last front. They are
     * taken in order as mentioned before. So, this method should only be used
     * when selecting individuals from a group of fronts that has a number of
     * feasible solutions less than (case 1) or exactly equal to (case) the
     * number of individuals required in the destination population. Here we use
     * it only for case 1. In case 1, each of the last (k) fronts, where (k) is
     * smaller than the population size (N), contains only one infeasible
     * solution. (Remember: in normal domination sorting, each infeasible
     * solution will lie in a separate front) (Note: any previous items in
     * destinationPopulation will be overwritten)
     *
     * @param destinationPopulation population at which selected members are
     * stored
     * @param fronts fronts from which individuals are selected
     */
    public void reFillPopulation(Individual[] destinationPopulation, List<List<Individual>> fronts) {
        int count = 0;
        int frontIndex = 0;
        int withinFrontIndex = 0;
        while (count < destinationPopulation.length) {
            if (withinFrontIndex < fronts.get(frontIndex).size()) {
                destinationPopulation[count] = fronts.get(frontIndex).get(withinFrontIndex);
                count++;
                withinFrontIndex++;
            } else {
                frontIndex++;
                withinFrontIndex = 0;
            }
        }
    }

    /**
     * Refills the destination population with the top ranked individuals of the
     * source population. The method should be used when the number of
     * individuals in the pre-limiting fronts is exactly equal to the number of
     * individuals required in the destination fronts. (Note: any previous items
     * in destinationPopulation will be overwritten)
     *
     * @param destinationPopulation population at which selected members are
     * stored
     * @param sourcePopulation population from which top-ranked members are
     * selected
     * @param lastFrontIndex the index of the limiting front (the front before
     * which the method should stop).
     * @throws EvaluationException
     */
    public void reFillPopulation(
            Individual[] destinationPopulation,
            Individual[] sourcePopulation,
            int lastFrontIndex) throws EvaluationException {
        int index = 0;
        for (Individual individual : sourcePopulation) {
            if (individual.getRank() - 1 < lastFrontIndex) { // Remember that ranks starts at 1 not at Zero
                destinationPopulation[index] = individual;
                index++;
            }
        }
    }

    protected Individual[] resolveCommonReferences(Individual[] individuals) {
        Individual[] resolvedIndividuals = new Individual[individuals.length];
        for (int i = 0; i < individuals.length; i++) {
            resolvedIndividuals[i] = new Individual(optimizationProblem, individuals[i], individualEvaluator);
        }
        return resolvedIndividuals;
    }

    protected void checkCommonReferences(Individual[] individuals, String message) {
        for (int i = 0; i < individuals.length - 1; i++) {
            for (int j = i + 1; j < individuals.length; j++) {
                if (individuals[i] == individuals[j]) {
                    throw new UnsupportedOperationException(String.format("Common Reference (IND(%d) & IND(%d)): %s", i, j, message));
                }
            }
        }
    }

    protected void checkCommonReferences(List<List<Individual>> fronts, String message) {
        List<Individual> individualsList = new ArrayList<Individual>();
        for (List<Individual> front : fronts) {
            for (Individual individual : front) {
                individualsList.add(individual);
            }
        }
        Individual[] arr = new Individual[individualsList.size()];
        individualsList.toArray(arr);
        checkCommonReferences(arr, message);
    }

    public void reportGenerationWiseInfo(
            Individual[] individuals,
            int generationsCounter,
            int runIndex,
            double epsilon,
            IndividualEvaluator evaluator,
            double[] idealPoint,
            double[] intercepts,
            double[] ultimateIdealPoint,
            double[] ultimateIntercepts,
            double utopianEpsilon,
            List<ReferenceDirection> refDirsList,
            String outputDir) throws IOException {
        if (DUMP_ALL_GENERATIONS_DECISION_SPACE || DUMP_ALL_GENERATIONS_OBJECTIVE_SPACE || DUMP_ALL_GENERATIONS_MATLAB_SCRIPTS || DUMP_ALL_GENERATIONS_NORMALIZED_MATLAB_SCRIPTS) {

            // <editor-fold desc="Create a directory for the current run" defaultstate="collapsed">
            // Note that re-creating the same directory at each generation does
            // not delete its contents. Actually, this line is only effective
            // at the first generation, otherwise it is doing nothing, as it is
            // trying to create an already existing directory.
            File generationWiseDumpDir = new File(
                    outputDir + String.format("generation_wise_run%03d/", runIndex));
            if (!generationWiseDumpDir.exists()) {
                generationWiseDumpDir.mkdir();
            }
            // </editor-fold>

            // <editor-fold desc="Extract non-dominated solutions of the current generation" defaultstate="collapsed">
            // Extract non-dominated solutions of the current generation
//            Individual[] nonDominatedIndividuals
//                    = OptimizationUtilities.getNonDominatedIndividuals(individuals, epsilon);
            // </editor-fold>

            // <editor-fold desc="Collect results only in single-objective problems" defaultstate="collapsed">
            // This step is done here to facilitate statistics collection
            // in single objective problems. Otherwise the objective values will
            // need to be collected post-optimization from files, which is very
            // time consuming especially for NSGA-III, which usually has
            // so many files, because of its small population size and large
            // number of generations, in single-objective problems.
            if (optimizationProblem.objectives.length == 1) {
                TestScript_SingleObjective.appendResults(
                        individuals,
                        runIndex,
                        generationsCounter,
                        optimizationProblem.getGenerationsCount(),
                        optimizationProblem.getPopulationSize());
            }
            // </editor-fold>

            Individual[] individualsToBePrinted = 
                    OptimizationUtilities.getNonDominatedIndividuals(individuals, epsilon);

            if (DUMP_ALL_GENERATIONS_OBJECTIVE_SPACE) {
                // <editor-fold desc="Dumping objective-space data" defaultstate="collapsed">
                String objSpaceFilePath = generationWiseDumpDir
                        + "/"
                        + String.format("gen_%04d_obj.dat", generationsCounter);
                InputOutput.dumpPopulationObjectiveSpace_Undecorated_SpaceSeparated(optimizationProblem, individualsToBePrinted, objSpaceFilePath);
                // </editor-fold>
            }

            // <editor-fold desc="Generate NON-normalized Matlab plots for each generation" defaultstate="collapsed">
            if (DUMP_ALL_GENERATIONS_MATLAB_SCRIPTS) {
                // Dump Matlab plotting script of the current generation
                // (ONLY NON_DOMINATED POINTS ARE CONSIDERED HERE)
                String matlabScriptFilePath
                        = generationWiseDumpDir
                        + "/"
                        + String.format("gen_%04d_obj_matlab.m", generationsCounter);
                // Dump the Matlab plotting script
                if (optimizationProblem.objectives.length == 2) {
                    //InputOutput.dumpMatlabPlottinScriptFor2dPoints(matlabScriptFilePath, evaluator, OptimizationUtilities.getNonDominatedIndividuals(individuals, epsilon));
                } else if(optimizationProblem.objectives.length == 3) {
                    //InputOutput.dumpMatlabPlottinScriptFor3dPoints(matlabScriptFilePath, evaluator, OptimizationUtilities.getNonDominatedIndividuals(individuals, epsilon));
                }
            }
            // </editor-fold>

            // <editor-fold desc="Dumping variable-space data">
            if (DUMP_ALL_GENERATIONS_DECISION_SPACE) {
                // Dump Decision Space
                String decisionSpaceFilePath
                        = generationWiseDumpDir
                        + "/"
                        + String.format("gen_%04d_var.dat", generationsCounter);
                InputOutput.dumpPopulationRealDecisionSpace_Undecorated_SpaceSeparated(optimizationProblem, individualsToBePrinted, decisionSpaceFilePath);
            }
            // </editor-fold>
        }
    }

    public abstract Individual[] niching(List<List<Individual>> fronts, int remainingIndvsCount) throws DoubleAssignmentException;

    public abstract Individual[] start(
            String outputDir,
            int runIndex,
            double epsilon
    )
            throws
            EvaluationException,
            FileNotFoundException,
            DoubleAssignmentException,
            IOException;

    public abstract Individual[] start(
            String outputDir,
            int runIndex,
            double epsilon,
            double hvLimit,
            int funcEvaluationsLimit
    )
            throws
            EvaluationException,
            FileNotFoundException,
            DoubleAssignmentException,
            IOException;

    public abstract Individual[] start(
            String outputDir,
            int runIndex,
            double epsilon,
            double hvLimit,
            int funcEvaluationsLimit,
            double kktMetricLimit
    )
            throws
            EvaluationException,
            FileNotFoundException,
            DoubleAssignmentException,
            IOException;

    public abstract String getAlgorithmName();
}
