/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package engines;

import Utilities.InputOutput;
import Utilities.Mathematics;
import Utilities.PerformanceMetrics;
import Utilities.RandomNumberGenerator;
import emo.DoubleAssignmentException;
import emo.Individual;
import emo.IndividualsSet;
import reference_directions.InvalidReferenceDirectionValue;
import emo.OptimizationProblem;
import emo.OptimizationUtilities;
import static engines.AbstractGeneticEngine.DEBUG_ALL;
import static engines.AbstractGeneticEngine.DEBUG_ASSOCIATION;
import static engines.AbstractGeneticEngine.DEBUG_IDEAL_POINT;
import static engines.AbstractGeneticEngine.DEBUG_INTERCEPTS;
import static engines.AbstractGeneticEngine.DEBUG_POPULATIONS;
import static engines.AbstractGeneticEngine.DEBUG_RANKING;
import static engines.AbstractGeneticEngine.DEBUG_REFERENCE_DIRECTIONS;
import static engines.AbstractGeneticEngine.DEBUG_TRANSLATION;
import static engines.AbstractGeneticEngine.EXTREME_POINTS_DEEP_DEBUG;
import static engines.AbstractGeneticEngine.MIN_DOUBLE_VALUE;
import java.io.BufferedReader;
import reference_directions.ReferenceDirection;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import net.sourceforge.jeval.EvaluationException;
import optimization.TestScript_FixedGenerations;
import parsing.IndividualEvaluator;
import reference_directions.NestedReferenceDirectionsFactory;
import reference_directions.ReferenceDirectionsFactory;

/**
 *
 * @author toshiba
 */
public class UnifiedNSGA3Engine extends AbstractGeneticEngine {

    public static final boolean REMOVE_ADDITIONAL_INDIVIDUALS = false;

    List<ReferenceDirection> referenceDirectionsList;
    protected int currentGenerationsCount = -1;

    public UnifiedNSGA3Engine(OptimizationProblem optimizationProblem, IndividualEvaluator individualEvaluator, int[] divisions) throws EvaluationException {
        super(optimizationProblem, individualEvaluator);
        // Create Reference Diretions
        referenceDirectionsList = new NestedReferenceDirectionsFactory(optimizationProblem.objectives.length).generateDirections(divisions);
        if (DEBUG_ALL || DEBUG_REFERENCE_DIRECTIONS) {
            InputOutput.displayReferenceDirections("Reference Directions", referenceDirectionsList);
        }
    }

    public UnifiedNSGA3Engine(OptimizationProblem optimizationProblem, IndividualEvaluator individualEvaluator) throws EvaluationException {
        super(optimizationProblem, individualEvaluator);
        // Create Reference Diretions
        referenceDirectionsList = new ReferenceDirectionsFactory(optimizationProblem.objectives.length).generateDirections(optimizationProblem.getSteps());
        if (DEBUG_ALL || DEBUG_REFERENCE_DIRECTIONS) {
            InputOutput.displayReferenceDirections("Reference Directions", referenceDirectionsList);
        }
    }

    /**
     *
     * @param optimizationProblem
     * @param individualEvaluator
     * @param directionsFilePath
     */
    public UnifiedNSGA3Engine(OptimizationProblem optimizationProblem, IndividualEvaluator individualEvaluator, String directionsFilePath) throws EvaluationException, IOException {
        super(optimizationProblem, individualEvaluator);
        // Read Directions from File
        BufferedReader reader = null;
        try {
            referenceDirectionsList = new ArrayList<ReferenceDirection>();
            reader = new BufferedReader(new FileReader(directionsFilePath));
            String line;
            while ((line = reader.readLine()) != null) {
                String[] splits = line.trim().split(" ");
                if (!line.isEmpty()) {
                    double[] dirValues = new double[splits.length];
                    for (int i = 0; i < splits.length; i++) {
                        dirValues[i] = Double.parseDouble(splits[i]);
                    }
                    referenceDirectionsList.add(new ReferenceDirection(dirValues));
                }
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
        InputOutput.displayReferenceDirections("Reference Directions", referenceDirectionsList);
    }

    public static int totalSelectionsCount = 0;
    public static int bothFeasibleCount = 0;

    @Override
    protected Individual tournamentSelect(IndividualsSet subset) {
        totalSelectionsCount++;
        Individual individual1 = subset.getIndividual1();
        Individual individual2 = subset.getIndividual2();
        // If only one of the solutions is infeasible, return the feasible solution
        if (optimizationProblem.constraints != null
                && optimizationProblem.constraints.length != 0
                && (individual1.isFeasible() ^ individual2.isFeasible())) {
            // If the problem is constrained and one of the individuals
            // under investigation is feasible while the other is infeasible,
            // return the feasible one (which is normally the dominating
            // individual).
            if (individual1.isFeasible()) {
                return individual1;
            } else {
                return individual2;
            }
        } else if (!individual1.isFeasible() && !individual2.isFeasible()) {
            // If both the two solutions are infeasible, return the less violating.
            if (Mathematics.compare(
                    individual1.getTotalConstraintViolation(),
                    individual2.getTotalConstraintViolation()) == 1) {
                // individual1 is less violating (remember: the more negative
                // the value, the more the violation)
                return individual1;
            } else {
                // individual2 is less violating
                return individual2;
            }
        } //else if (optimizationProblem.objectives.length == 1) {
        // If we have only one objective and both solutions are feasible,
        // return the better in terms of objective value (i.e. the
        // dominating solution). If the two ranks are equal this means that
        // the two individuals are identical so return any of them.
        // (Remember: in the case of single objective, one idividual must
        // dominate the other unless both are identical to each other. This
        // is the only case where they will have the same rank)
        //if (individual1.dominates(individual2)) {
        //return individual1;
        //} else {
        //return individual2;
        //}
        //}
        else if (currentGenerationsCount != 0 && individual1.getReferenceDirection().equals(individual2.getReferenceDirection())) {
            bothFeasibleCount++;
            // If both the two solutions are feasible. You have the
            // following two options:
            // If the two individuals belong to the same reference direction,
            // return the one with lower rank.
            if (individual1.getRank() < individual2.getRank()) {
                return individual1;
            } else if (individual2.getRank() < individual1.getRank()) {
                return individual2;
            } else {
                // If they both belong to the same rank return the one
                // closest to the reference direction.
                if (Mathematics.compare(
                        individual1.getPerpendicularDistance(),
                        individual2.getPerpendicularDistance()) == -1) {
                    return individual1;
                } else {
                    return individual2;
                }
            }
        } else {
            // If the two individuals are associated with two different reference
            // directions, then return one of them randomly
            if (RandomNumberGenerator.randomperc() <= 0.5) {
                return individual1;
            } else {
                return individual2;
            }
        }
    }

    protected void translate(Individual[] individuals, double[] idealPoint) {
        for (Individual individual : individuals) {
            if (individual.translated) {
                throw new UnsupportedOperationException("Translated individuals should NOT be re-translated");
            }
            for (int i = 0; i < idealPoint.length; i++) {
                individual.nonTranslatedObjectives[i] = individual.getObjective(i);
                individual.setObjective(i, individual.getObjective(i) - idealPoint[i]);
            }
            individual.translated = true;
        }
    }

    protected Individual[] getExtremePoints(Individual[] individuals, Individual[] previousExtremePoints, double[] idealPoint, double[] prevIdealPoint) throws EvaluationException {
        if (EXTREME_POINTS_DEEP_DEBUG) {
            System.out.println("---------------------------");
            System.out.println("Extreme Points Calculations");
            System.out.println("---------------------------");
        }
        // Calculate unit directions
        double[][] unitDirections
                = new double[optimizationProblem.objectives.length][optimizationProblem.objectives.length];
        for (int i = 0; i < optimizationProblem.objectives.length; i++) {
            for (int j = 0; j < optimizationProblem.objectives.length; j++) {
                if (i == j) {
                    unitDirections[i][j] = 1;
                } else {
                    unitDirections[i][j] = Math.pow(10, -6);
                }
            }
        }
        if (EXTREME_POINTS_DEEP_DEBUG) {
            System.out.format("* Unit Directions Created Successfully:%n");
            for (int i = 0; i < unitDirections.length; i++) {
                System.out.format("(");
                for (int j = 0; j < unitDirections[i].length; j++) {
                    System.out.format("%7.3f", unitDirections[i][j]);
                    if (j != unitDirections[i].length - 1) {
                        System.out.format(",");
                    }
                }
                System.out.format(")%n");
            }
        }
        // If the previous extreme points are not null then let's take them into
        // consideration because we might not get better extreme points from the
        // current population (remember: we are re-considering them because the
        // individual representing the extreme point might have been lost in
        // transition from the previous generation to this one)
        double[] previousMaxValues = new double[optimizationProblem.objectives.length];
        if (previousExtremePoints != null) {
            // Re-translate the extreme points, because new ideal point may
            // have been found since there previous translation (in the
            // previous generation).
            for (int i = 0; i < optimizationProblem.objectives.length; i++) {
                for (int j = 0; j < optimizationProblem.objectives.length; j++) {
                    double retranslatedObjValue
                            = previousExtremePoints[i].getObjective(j) - (idealPoint[j] - prevIdealPoint[j]);
                    previousExtremePoints[i].setObjective(j, retranslatedObjValue);
                }
            }
            // Re-Calculate the previous MAX values of the previous extreme points
            for (int i = 0; i < optimizationProblem.objectives.length; i++) {
                // Set the unit direction (unit direction j)
                double[] wDirection = unitDirections[i];
                previousMaxValues[i] = previousExtremePoints[i].getObjective(0) / wDirection[0];
                for (int k = 1; k < optimizationProblem.objectives.length; k++) {
                    double nextValue = previousExtremePoints[i].getObjective(k) / wDirection[k];
                    if (nextValue > previousMaxValues[i]) {
                        previousMaxValues[i] = nextValue;
                    }
                }
            }
        }
        if (EXTREME_POINTS_DEEP_DEBUG) {
            System.out.format("* Previous MAX Values:%n");
            System.out.println("[");
            for (int k = 0; k < optimizationProblem.objectives.length; k++) {
                System.out.format("%12.3f", previousMaxValues[k]);
                if (k != optimizationProblem.objectives.length - 1) {
                    System.out.println(",");
                }
            }
            System.out.format("%n]%n");
        }
        // Create an empty 2D array to store the intercepts
        Individual[] extremePoints = new Individual[optimizationProblem.objectives.length];
        // Iterate over all the basic directions (# of directions = # of objectives)
        for (int i = 0; i < optimizationProblem.objectives.length; i++) {
            if (EXTREME_POINTS_DEEP_DEBUG) {
                System.out.format("* OBJ(%d)%n", i);
            }
            // Set the unit direction (unit direction j)
            double[] wDirection = unitDirections[i];
            if (EXTREME_POINTS_DEEP_DEBUG) {
                System.out.format("Direction = (%6.2f, %6.2f, %6.2f)%n", wDirection[0], wDirection[1], wDirection[2]);
            }
            // Each slot in the following array stores the largest value of
            // x.obj(k)/w(k) for 0 <= k <= #objectives-1, for some individual x.
            double[] maxArr = new double[individuals.length];
            // Iterate over all the members of the populations
            for (int j = 0; j < individuals.length; j++) {
                if (EXTREME_POINTS_DEEP_DEBUG) {
                    System.out.format("INDV(%d)%n", j);
                    System.out.format("  max(ind[%d].obj[0]/dir[0], ind[%d].obj[1]/dir[1], ind[%d].obj[2]/dir[2])%n", j, j, j);
                    System.out.format("= max(%10.2f/%-9.2f, %10.2f/%-9.2f, %10.2f/%-9.2f)%n",
                            individuals[j].getObjective(0),
                            wDirection[0],
                            individuals[j].getObjective(1),
                            wDirection[1],
                            individuals[j].getObjective(2),
                            wDirection[2]);
                    System.out.format("= max(%20.2f, %20.2f, %20.2f)%n",
                            individuals[j].getObjective(0) / wDirection[0],
                            individuals[j].getObjective(1) / wDirection[1],
                            individuals[j].getObjective(2) / wDirection[2]);
                }
                double max = individuals[j].getObjective(0) / wDirection[0];
                for (int k = 1; k < optimizationProblem.objectives.length; k++) {
                    double nextValue = individuals[j].getObjective(k) / wDirection[k];
                    if (nextValue > max) {
                        max = nextValue;
                    }
                }
                maxArr[j] = max;
                if (EXTREME_POINTS_DEEP_DEBUG) {
                    System.out.format("= %6.2f%n%n", max);
                }
            }
            // Select the minimum value out of maxArr
            int minIndex = 0;
            for (int j = 1; j < maxArr.length; j++) {
                if (maxArr[j] < maxArr[minIndex]) {
                    minIndex = j;
                }
            }
            if (EXTREME_POINTS_DEEP_DEBUG) {
                System.out.format("Smallest Max Index = %d (value: %6.2f)%n", minIndex, maxArr[minIndex]);
            }
            if (previousExtremePoints != null && previousMaxValues[i] < maxArr[minIndex]) {
                // This means that the previous extreme point was better than
                // the current extreme point and we should retain the previous
                // extreme point instead of replacing it with a new weaker one.
                extremePoints[i] = previousExtremePoints[i];
                if (EXTREME_POINTS_DEEP_DEBUG) {
                    System.out.format(">>> Extreme Point remains unchanged%n");
                }
            } else {
                // Now the individual whose index in minIndex in the population is
                // the one representing the extreme factor in the current directions.
                extremePoints[i] = new Individual(optimizationProblem, individuals[minIndex], individualEvaluator);
                if (EXTREME_POINTS_DEEP_DEBUG) {
                    System.out.format(">>> Extreme Point(%d) = Ind(%d)%n", i, minIndex);
                }
            }
        }
        // Return the extreme points in all basic directions
        return extremePoints;
    }

    protected double[] getIntercepts(Individual[] extremePoints, Individual[] population) {
        // Calculating the vector of maximum objective values & the Nadir point
        // Initialize the structures & set all their initial values to
        // Negative Infinity
        double[] maxObjValues = new double[optimizationProblem.objectives.length];
        double[] nadirPoint = new double[optimizationProblem.objectives.length];
        for (int i = 0; i < nadirPoint.length; i++) {
            maxObjValues[i] = -1 * Double.MIN_VALUE;
            nadirPoint[i] = -1 * Double.MIN_VALUE;
        }
        // Traverse all the individuals of the population and get their maximum
        // value of objective (The simplest way of calculating the nadir point 
        // is to get these maximum values among the first front individuals)
        for (int i = 0; i < population.length; i++) {
            for (int j = 0; j < optimizationProblem.objectives.length; j++) {
                if (maxObjValues[j] < population[i].getObjective(j)) {
                    maxObjValues[j] = population[i].getObjective(j);
                }
                if (population[i].getRank() == 1) {
                    if (nadirPoint[j] < population[i].getObjective(j)) {
                        nadirPoint[j] = population[i].getObjective(j);
                    }
                }
            }
        }

        if (DEBUG_ALL || DEBUG_INTERCEPTS) {
            System.out.println("-----------");
            System.out.println("Nadir Point");
            System.out.println("-----------");
            System.out.print("(");
            for (int i = 0; i < nadirPoint.length; i++) {
                System.out.format("%5.2f", nadirPoint[i]);
                if (i != nadirPoint.length - 1) {
                    System.out.print(",");
                }
            }
            System.out.println(")");

            System.out.println("-----------");
            System.out.println("Max Vector");
            System.out.println("-----------");
            System.out.print("(");
            for (int i = 0; i < maxObjValues.length; i++) {
                System.out.format("%5.2f", maxObjValues[i]);
                if (i != maxObjValues.length - 1) {
                    System.out.print(",");
                }
            }
            System.out.println(")");
        }

        // Caculating the intercepts
        double[] intercepts = new double[optimizationProblem.objectives.length];
        // Create the hyperplane
        // Prepare your arrays for gaussian elimination
        double[][] Z = new double[optimizationProblem.objectives.length][optimizationProblem.objectives.length];
        for (int i = 0; i < Z.length; i++) {
            for (int j = 0; j < Z[i].length; j++) {
                Z[i][j] = extremePoints[i].getObjective(j);
            }
        }
        double[] u = new double[optimizationProblem.objectives.length];
        for (int i = 0; i < u.length; i++) {
            u[i] = 1;
        }
        boolean useNadir = false;
        // Solve the system of equations using gaussian elimination
        try {
            intercepts = Mathematics.gaussianElimination(Z, u);
        } catch (Mathematics.SingularMatrixException ex) {
            useNadir = true;
        }
        // If gaussian elimination resulted in -ve or infinite intercepts,
        // again use the nadir point.
        if (useNadir) {
            for (double intercept : intercepts) {
                if (intercept < MIN_DOUBLE_VALUE) {
                    useNadir = true;
                    break;
                }
            }
        }
        // Now get the intercept = 1/alpha
        if (!useNadir) {
            for (int i = 0; i < intercepts.length; i++) {
                intercepts[i] = 1 / intercepts[i];
                if (intercepts[i] == Double.NaN) {
                    useNadir = true;
                }
            }
        }
        // If the follwing condition is true this means that you have to resort to the nadir point
        if (useNadir) {
            System.arraycopy(nadirPoint, 0, intercepts, 0, intercepts.length);
        }
        // If any of the intercepts is still Zero (which means that one of
        // the nadir values is Zero), then use the maximum value of each
        // objective instead (remember that these values were calculated among
        // all the individuals, not just the first-front individuals)
        for (double intercept : intercepts) {
            if (intercept < MIN_DOUBLE_VALUE) {
                System.arraycopy(maxObjValues, 0, intercepts, 0, intercepts.length);
                break;
            }
        }
        // Return the intercepts
        return intercepts;
    }

    protected double[][] associate(Individual[] individuals, List<ReferenceDirection> referenceDirections, double[] intercepts) {
        // Reset Reference Directions Information
        for (ReferenceDirection referenceDirection : referenceDirections) {
            referenceDirection.surroundingIndividuals.clear();
        }
        double d1, perpendicularDistance, refDirNorm;
        double[][] distanceMatrix = new double[individuals.length][referenceDirections.size()];
        for (Individual individual : individuals) {
            // Initialize the array of distances
            individual.distancesFromDirs = new double[referenceDirections.size()];
        }
        /* Perpendicular distance calculation */
        for (int i = 0; i < referenceDirections.size(); i++) {
            for (int j = 0; j < individuals.length; j++) {
                // Calculate the distance from the current individual to the current direction
                d1 = 0.0; // Eventually will be the DOT product of the individual and the direction
                refDirNorm = 0.0; // Eventually will be the NORM of the direction
                for (int k = 0; k < optimizationProblem.objectives.length; k++) {
                    d1 += individuals[j].getObjective(k) * referenceDirections.get(i).direction[k] / intercepts[k];
                    refDirNorm += Math.pow(referenceDirections.get(i).direction[k], 2);
                }
                refDirNorm = Math.sqrt(refDirNorm); // After this line refDirNorm will be the NORM of the reference direction
                double scalarProjection = d1 / refDirNorm; // Scalar projection of the individual on the reference direction (i.e. the length of the projection of the point on the reference direction)
                perpendicularDistance = 0.0; // Eventaully will be the required perpendicular distance
                for (int k = 0; k < optimizationProblem.objectives.length; k++) {
                    perpendicularDistance += Math.pow(
                            (scalarProjection * referenceDirections.get(i).direction[k] / refDirNorm - individuals[j].getObjective(k) / intercepts[k]),
                            2.0);
                }
                perpendicularDistance = Math.sqrt(perpendicularDistance); // After this line perpendicularDistance will be the required distance
                // Store the distance in the matrix and in the individual object
                distanceMatrix[j][i] = perpendicularDistance;
                individuals[j].distancesFromDirs[i] = distanceMatrix[j][i];
            }
        }
        if (DEBUG_ALL || DEBUG_ASSOCIATION) {
            InputOutput.displayDistanceMatrix(optimizationProblem, referenceDirections, individuals, distanceMatrix);
        }
        // Find the closest reference direction to each individual
        for (int i = 0; i < individuals.length; i++) {
            double minDistance = distanceMatrix[i][0];
            ReferenceDirection refDir = referenceDirections.get(0);
            for (int j = 1; j < referenceDirections.size(); j++) {
                if (distanceMatrix[i][j] < minDistance) {
                    minDistance = distanceMatrix[i][j];
                    refDir = referenceDirections.get(j);
                }
            }
            individuals[i].setReferenceDirection(refDir);
            individuals[i].setPerpendicularDistance(minDistance);
            refDir.surroundingIndividuals.add(individuals[i]);
            individuals[i].validReferenceDirection = true;
        }
        return distanceMatrix;
    }

    @Override
    public Individual[] niching(
            List<List<Individual>> fronts, int remainingIndvsCount) {
        int individualsCount = 0;
        // Get the last front Fl index
        int lastFrontIndex = -1;
        while (individualsCount < optimizationProblem.getPopulationSize()) {
            lastFrontIndex++;
            individualsCount += fronts.get(lastFrontIndex).size();
        }
        // Get the last front Fl (using the index)
        List<Individual> lastFront = fronts.get(lastFrontIndex);
        // Make a copy of the directionsList
        List<ReferenceDirection> referenceDirectionsListCopy
                = new ArrayList<ReferenceDirection>(referenceDirectionsList);
        // Make a copy of the last front Fl (just shallow copies of the objects)
        List<Individual> lastFrontCopy = new ArrayList<Individual>(lastFront);
        // Create the final array of members that will complete the remaining part of the population
        Individual[] remainingIndividuals = new Individual[remainingIndvsCount];
        // Create associationCountArr array (the number of individuals 
        // associated with each Direction in P(t+1)):
        // Remember that the value refDirection.surroundingIndividuals is not
        // what we want here, because it contains the number of individuals
        // associated with the direction out of all the candidates (St), while
        // what we want here is the number of individuals associated with the
        // direction out of the confirmed members only (Pt+1).
        List<Integer> associationCountList = new ArrayList<Integer>();
        for (int i = 0; i < referenceDirectionsListCopy.size(); i++) {
            associationCountList.add(new Integer(0));
            for (List<Individual> front : fronts) {
                if (front.equals(lastFront)) {
                    break;
                }
                for (Individual individual : front) {
                    if (individual.getReferenceDirection().equals(referenceDirectionsListCopy.get(i))) {
                        associationCountList.set(i, associationCountList.get(i) + 1);
                    }
                }
            }
        }
        // The following list of lists contains all the associations of the last front Fl
        List<List<Individual>> lastFrontRefAssociations = new ArrayList<List<Individual>>();
        for (int d = 0; d < referenceDirectionsListCopy.size(); d++) {
            // Find the set of members associated with dir(j) in Fl (nextFrontAssociatedMembers)
            List<Individual> lastFrontAssociatedMembersList = new ArrayList<Individual>();
            for (Individual individual : lastFrontCopy) {
                if (individual.getReferenceDirection().equals(referenceDirectionsListCopy.get(d))) {
                    lastFrontAssociatedMembersList.add(individual);
                }
            }
            lastFrontRefAssociations.add(lastFrontAssociatedMembersList);
        }

        // Repeat until REMAINING members are added to the final list
        while (remainingIndvsCount > 0) {
            // Find the reference directions with the lowest number of associated members
            int minDirClusterSize = associationCountList.get(0);
            for (int i = 1; i < referenceDirectionsListCopy.size(); i++) {
                if (associationCountList.get(i).intValue() < minDirClusterSize) {
                    minDirClusterSize = associationCountList.get(i).intValue();
                }
            }

            // The following part is unique to this method and is not available
            // in the deterministic version of it. Just get all the directions
            // having the same minimum associationCount and pick one of them
            // randomly.
            List<Integer> minAssociationsDirsIndices = new ArrayList<Integer>();
            for (int i = 0; i < referenceDirectionsListCopy.size(); i++) {
                if (associationCountList.get(i).intValue() == minDirClusterSize) {
                    minAssociationsDirsIndices.add(i);
                }
            }

            /**
             * => TWO OPTION: 1-SELECT THE FIRST dirIndex, 2-SELECT A RANDOM
             * dirIndex. (OPTION(2) IS THE ONE CONFORMING WITH THE PAPER)
             */
            int dirIndex = minAssociationsDirsIndices.get(/*new Random().nextInt(minAssociationsDirsIndices.size())*/
                    RandomNumberGenerator.rnd(0, minAssociationsDirsIndices.size() - 1));
            //int dirIndex = minAssociationsDirsIndices.get(0);

            // If nextFrontAssociatedPoints is Empty, Discard reference
            // direction from the set of reference directions (for the current generation).
            if (lastFrontRefAssociations.get(dirIndex).isEmpty()) {
                referenceDirectionsListCopy.remove(dirIndex);
                associationCountList.remove(dirIndex);
                lastFrontRefAssociations.remove(dirIndex);
            } else {
                int newMemberIndex;
                if (associationCountList.get(dirIndex).intValue() == 0) {
                    // If the number of members(of Pt+1) associated with direction Pj is Zero,
                    // Select the member of nextFrontAssociatedMembers with the shortest 
                    // perpendicular distance from Pj (to be added to the final list)
                    int minDistanceIndividualIndex;

                    minDistanceIndividualIndex = 0;
                    double minDistance = lastFrontRefAssociations.get(dirIndex).get(0).getPerpendicularDistance();
                    for (int i = 1; i < lastFrontRefAssociations.get(dirIndex).size(); i++) {
                        if (lastFrontRefAssociations.get(dirIndex).get(i).getPerpendicularDistance() < minDistance) {
                            minDistanceIndividualIndex = i;
                            minDistance = lastFrontRefAssociations.get(dirIndex).get(i).getPerpendicularDistance();
                        }
                    }

                    newMemberIndex = minDistanceIndividualIndex;
                } else {
                    // Otherwise, select the first individual of nextFrontAssociatedMembers (to be added to the final list)
                    newMemberIndex = 0;
                }
                // Add the new Member to the final list
                remainingIndividuals[remainingIndividuals.length - remainingIndvsCount] = lastFrontRefAssociations.get(dirIndex).get(newMemberIndex);
                // Add the new member to the list of members associated with Pj
                associationCountList.set(dirIndex, associationCountList.get(dirIndex).intValue() + 1);
                // Exclude the new member from the list Fl
                lastFrontRefAssociations.get(dirIndex).remove(newMemberIndex);
                // Decrement the counter of remaining members
                remainingIndvsCount--;
            }
        }
        return remainingIndividuals;
    }

    protected void resetObjectiveValues(Individual[] individuals) {
        for (int i = 0; i < individuals.length; i++) {
            if (!individuals[i].translated) {
                System.out.println("STOP");
            }
        }
        for (int i = 0; i < individuals.length; i++) {
            if (!individuals[i].translated) {
                throw new UnsupportedOperationException("Un-translated individuals connot be reset (they are already un-translated i.e. reset)");
            }
            for (int j = 0; j < optimizationProblem.objectives.length; j++) {
                individuals[i].setObjective(j, individuals[i].nonTranslatedObjectives[j]);
            }
            individuals[i].translated = false;
        }
    }

    @Override
    public Individual[] start(
            String outputDir,
            int runIndex,
            double epsilon
    ) throws EvaluationException, FileNotFoundException, IOException, DoubleAssignmentException {
        return start(outputDir, runIndex, epsilon, Double.MAX_VALUE, Integer.MAX_VALUE);
    }

    @Override
    public Individual[] start(
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
            IOException {
        return start(outputDir,
                runIndex,
                epsilon,
                hvLimit,
                funcEvaluationsLimit,
                -1);
    }

    @Override
    public Individual[] start(
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
            IOException {
        double[] idealPoint;
        double[] prevIdealPoint;
        double[] intercepts = null;
        Individual[] extremePoints = null;
        // Generate initial population
        Individual[] parentPopulation = generateInitialPopulation();
        if (DEBUG_ALL || DEBUG_POPULATIONS) {
            InputOutput.displayPopulation(optimizationProblem, "0", parentPopulation);
        }
        // Perform initial ranking
        assign_rank(parentPopulation, epsilon);
        if (DEBUG_ALL || DEBUG_RANKING) {
            InputOutput.displayRanks("0", parentPopulation);
        }
        // Create the initial ideal point
        idealPoint = getInitialIdealPoint(parentPopulation);
        if (DEBUG_ALL || DEBUG_IDEAL_POINT) {
            InputOutput.displayIdealPoint("0", idealPoint);
        }
        // Translate all the objective values
        translate(parentPopulation, idealPoint);
        // Save the current ideal point for future use (we need to improve
        // on the current ideal points as generations go. Creating a completely
        // new ideal point from each generation independent from previous
        // generations, may cause degradation in the ideal point values)
        prevIdealPoint = new double[optimizationProblem.objectives.length];
        if (OptimizationUtilities.getFeasibleIndividuals(parentPopulation).length > 0) {
            // Normalization
            extremePoints = getExtremePoints(OptimizationUtilities.getFeasibleIndividuals(parentPopulation), extremePoints, idealPoint, prevIdealPoint);
            intercepts = getIntercepts(extremePoints, OptimizationUtilities.getFeasibleIndividuals(parentPopulation));
            //extremePoints = getExtremePoints(parentPopulation, extremePoints, idealPoint, prevIdealPoint);
            //intercepts = getIntercepts(extremePoints, parentPopulation);
            if (DEBUG_ALL || DEBUG_INTERCEPTS) {
                InputOutput.displayIntercepts("Initial Population Intercepts", intercepts);
            }
        } else {

        }
        // Reset obj. values to the pre-translation values
        // (in order to start the loop with valid objective values)
        resetObjectiveValues(parentPopulation);
        // Get the actual reference point (to calculate hypervolume)
        double[] hvReferencePoint = null;
        double[] hvIdealPoint = null;
        if (hvLimit != Double.MAX_VALUE) {
            hvReferencePoint = individualEvaluator.getReferencePoint();
            // Adjust reference point so that end points will contribute to the
            // hypervolume
            for (int i = 0; i < hvReferencePoint.length; i++) {
                hvReferencePoint[i] *= 1.01;
            }
            hvIdealPoint = individualEvaluator.getIdealPoint();
        }
        // Initialize currentHV
        double currentHV = -1;
        // Initialize currentMedianKKTMetric
        double currentKktMetric = Double.MAX_VALUE;
        // Generations counter
        int generationsCounter = 0;
        // Each loop represents one generation
        for (; generationsCounter < optimizationProblem.getGenerationsCount()
                && currentHV < hvLimit
                && (currentKktMetric != currentKktMetric || currentKktMetric > kktMetricLimit)
                && individualEvaluator.getIndividualEvaluationsCount() < funcEvaluationsLimit;
                generationsCounter++) {
            currentGenerationsCount = generationsCounter;
            InputOutput.displayGenerationCount(generationsCounter);
            // Create the offspring (tournament selection & crossover)
            Individual[] offspringPopulation = getOffspringPopulation(parentPopulation);
            // Mutation (binary & real)
            mutate(offspringPopulation);
            // Update the objective values & constraints violation of these offspring
            for (Individual individual : offspringPopulation) {
                individualEvaluator.updateIndividualObjectivesAndConstraints(optimizationProblem, individual);
            }
            if (DEBUG_ALL || DEBUG_POPULATIONS) {
                InputOutput.displayPopulation(optimizationProblem, generationsCounter + "-offspring", offspringPopulation);
            }
            // Merge the parents and the offspring in a single population
            Individual[] mergedPopulation
                    = merge(parentPopulation, offspringPopulation);
            // Resolve common references in the merged population. Notice that
            // selection may cause copying some individuals to be reference
            // copied from the parent population to the offspring population.
            mergedPopulation = resolveCommonReferences(mergedPopulation);
            if (DEBUG_ALL || DEBUG_POPULATIONS) {
                InputOutput.displayPopulation(optimizationProblem, String.format("%d-population+%d-offspring", generationsCounter, generationsCounter), mergedPopulation);
            }
            // Non-dominated sorting
            List<List<Individual>> fronts = assign_rank(mergedPopulation, epsilon);
            if (DEBUG_ALL || DEBUG_RANKING) {
                InputOutput.displayRanks(String.format("%d-population+%d-offspring", generationsCounter, generationsCounter), mergedPopulation);
            }
            // Update ideal point
            prevIdealPoint = idealPoint;
            idealPoint = getUpdatedIdealPoint(mergedPopulation, idealPoint);
            if (DEBUG_ALL || DEBUG_IDEAL_POINT) {
                InputOutput.displayIdealPoint(String.format("Ideal Point(merged population)"), idealPoint);
            }
            // Translate objective values
            translate(mergedPopulation, idealPoint);
            if (DEBUG_ALL || DEBUG_TRANSLATION) {
                InputOutput.displayPopulationObjectiveSpace(optimizationProblem, String.format("After Translation (megred population)"), mergedPopulation);
            }
            // Count feasible solutions
            int feasibleCount = 0;
            for (Individual individual : mergedPopulation) {
                if (individual.isFeasible()) {
                    feasibleCount++;
                }
            }
            // Get reamining individuals
            int remainingIndividualsCount = getRemainingCount(fronts);
            int limitingFrontIndex = getLimitingFrontIndex(fronts);
            // Retrieve all the individuals in the fronts to be accomodated,
            // even those individuals of the last front that might be partially
            // accommodated.
            Individual[] candidates = OptimizationUtilities.getCandidates(mergedPopulation, fronts, optimizationProblem.getPopulationSize());
            // Retrieve the feasible solutions only from all the candidates
            Individual[] feasibleCandidates = OptimizationUtilities.getFeasibleIndividuals(candidates);

            if (feasibleCandidates.length != 0) {
                // Normalization (infeasible solutions should not go through
                // normalization, association and niching. This destroys
                // normalization and degrades results significantly)
                extremePoints = getExtremePoints(feasibleCandidates, extremePoints, idealPoint, prevIdealPoint);
                if (DEBUG_ALL || DEBUG_INTERCEPTS) {
                    InputOutput.displayExtremePoints(optimizationProblem, "Merged Population Extreme Points", extremePoints);
                }
                intercepts = getIntercepts(extremePoints, /*candidates*/ feasibleCandidates);
                if (DEBUG_ALL || DEBUG_INTERCEPTS) {
                    InputOutput.displayIntercepts("Merged Population Intercepts", intercepts);
                }
                // Association
                double[][] distanceMatrix = associate(feasibleCandidates, referenceDirectionsList, intercepts);
                if (DEBUG_ALL || DEBUG_ASSOCIATION) {
                    InputOutput.displayAssociationResluts(optimizationProblem, "Merged Population", /*candidates*/ feasibleCandidates);
                }
            }
            // Niching is required only if the number of feasible solutions in
            // the merged population is greater the population size, and the
            // number of remaining individuals to complete the population is
            // greater than Zero. (Notice that the having a number of feasible
            // solutions greater than the population size does not mean that
            // niching is necessary. Assume that population size is 8 and the 
            // first front contains 5 individuals while the second front 
            // contains 3 individuals. In this case even if there exist more
            // feasible individuals in subsequent fronts, they are useless,
            // because the new population is now full and niching is not
            // required.
            if (feasibleCount > optimizationProblem.getPopulationSize() && remainingIndividualsCount > 0) {
                // Niching
                Individual[] lastFrontSubset = niching(fronts, remainingIndividualsCount);
                if (DEBUG_ALL) {
                    System.out.println("---------------");
                    System.out.println("Niching Results");
                    System.out.println("---------------");
                    for (Individual individual : lastFrontSubset) {
                        System.out.println(individual.getShortVariableSpace());
                    }
                }
                // Re-fill the new population after niching
                reFillPopulation(parentPopulation, mergedPopulation, lastFrontSubset, limitingFrontIndex);
            } else {
                // Re-fill the new population in the case when niching is not
                // required.
                reFillPopulation(parentPopulation, fronts);
            }
            // Reset the values of the objective values for the next iteration
            resetObjectiveValues(parentPopulation);
            if (DEBUG_ALL || DEBUG_POPULATIONS) {
                System.out.format("--------------------------------------%n");
                System.out.format("Final Population After Generation(%d)%n", generationsCounter);
                System.out.format("--------------------------------------%n");
                InputOutput.displayPopulation(optimizationProblem, String.valueOf(generationsCounter + 1), parentPopulation);
            }
            // Report generation-wise information
            if (REMOVE_ADDITIONAL_INDIVIDUALS) {
                reportGenerationWiseInfo(
                        OptimizationUtilities.selectOnlyOneIndividualForEachDirection(
                                parentPopulation,
                                referenceDirectionsList),
                        generationsCounter,
                        runIndex,
                        epsilon,
                        individualEvaluator,
                        idealPoint,
                        intercepts,
                        (optimizationProblem.objectives.length == 1) ? null : individualEvaluator.getIdealPoint(),
                        (optimizationProblem.objectives.length == 1) ? null : individualEvaluator.getReferencePoint(),
                        0,
                        referenceDirectionsList,
                        outputDir);
            } else {
                reportGenerationWiseInfo(
                        parentPopulation,
                        generationsCounter,
                        runIndex,
                        epsilon,
                        individualEvaluator,
                        idealPoint,
                        intercepts,
                        (optimizationProblem.objectives.length == 1) ? null : individualEvaluator.getIdealPoint(),
                        (optimizationProblem.objectives.length == 1) ? null : individualEvaluator.getReferencePoint(),
                        0,
                        referenceDirectionsList,
                        outputDir);
            }
            // Update the current hypervolume (if hypervolume is used as a stopping criterion)
            if (hvLimit != Double.MAX_VALUE) {
                // Calculate current hypervolume
                if (optimizationProblem.objectives.length == 2) {
                    if (REMOVE_ADDITIONAL_INDIVIDUALS) {
                        currentHV = PerformanceMetrics.calculateHyperVolumeForTwoObjectivesOnly(
                                this,
                                parentPopulation,
                                referenceDirectionsList,
                                hvReferencePoint,
                                hvIdealPoint,
                                epsilon);
                    } else {
                        currentHV = PerformanceMetrics.calculateHyperVolumeForTwoObjectivesOnly(
                                this,
                                parentPopulation,
                                hvReferencePoint,
                                hvIdealPoint,
                                epsilon);
                    }
                } else {
                    throw new UnsupportedOperationException("Hypervolume can only be calculated for 2 objectives");
                }
            }
        } // End of generations loop (main loop)
        // Report last generation info
        Individual[] finalPopulation;
        if (REMOVE_ADDITIONAL_INDIVIDUALS) {
            finalPopulation = OptimizationUtilities.getNonDominatedIndividuals(
                    OptimizationUtilities.selectOnlyOneIndividualForEachDirection(
                            parentPopulation,
                            referenceDirectionsList),
                    epsilon
            );
        } else {
            finalPopulation = OptimizationUtilities.getNonDominatedIndividuals(parentPopulation, epsilon);
        }
        // Display some quick info.
        System.out.format("----------------%n");
        System.out.format("Final Population%n", generationsCounter);
        System.out.format("----------------%n");
        InputOutput.displayPopulation(optimizationProblem, String.valueOf(generationsCounter - 1), parentPopulation);
        InputOutput.displayPopulationUndecoratedObjectiveSpace(optimizationProblem, "Final Population Obj. Space", parentPopulation);
        // Display some simple statistics only for the 2 objectives case
        if (optimizationProblem.objectives.length == 2) {
            double minObj1 = OptimizationUtilities.getMinObjectiveValue(parentPopulation, 0);
            double maxObj1 = OptimizationUtilities.getMaxObjectiveValue(parentPopulation, 0);
            double minObj2 = OptimizationUtilities.getMinObjectiveValue(parentPopulation, 1);
            double maxObj2 = OptimizationUtilities.getMaxObjectiveValue(parentPopulation, 1);
            System.out.format("Objective-1 (MIN = %-15.3f, MIN = %-15.3f)%n", minObj1, maxObj1);
            System.out.format("Objective-2 (MIN = %-15.3f, MIN = %-15.3f)%n", minObj2, maxObj2);
        }
        // Return the final population
        return finalPopulation;
    }

    @Override
    public String getAlgorithmName() {
        return "unified_nsga3";
    }
}
