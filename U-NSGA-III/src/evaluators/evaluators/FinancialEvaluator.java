package evaluators;

import emo.Individual;
import emo.OptimizationProblem;
import emo.VirtualIndividual;

import javax.xml.stream.XMLStreamException;
import net.sourceforge.jeval.EvaluationException;
import parsing.IndividualEvaluator;
import parsing.InvalidOptimizationProblemException;

/*
*
* Farhan Hormasji, Bonnie Reiff
* CSE 848: Survey of Evolutionary Computing
* 
*/
public class FinancialEvaluator extends IndividualEvaluator 
{
	@Override
	/* Check that this exception is correct! */
	public double[] getReferencePoint() 
    {
    	throw new UnsupportedOperationException("Nadir point is not defined for this problem.");
    }

	@Override
	/* Check that this exception is correct! */
	public double[] getIdealPoint() 
    {
    	throw new UnsupportedOperationException("Ideal point is not defined for this problem.");
    }

    @Override
	public VirtualIndividual[] getParetoFront(int objectivesCount, int n) throws InvalidOptimizationProblemException, XMLStreamException, EvaluationException 
    {
        throw new UnsupportedOperationException("getParetoFront() is not defined for this problem.");
    }
    
    public void updateIndividualObjectivesAndConstraints(
            OptimizationProblem problem,
            Individual individual)
            throws EvaluationException 
    {
        double[] x = individual.real;
        double[] objs = getObjectives(calculateindicators(x));
        individual.setObjective(0, objs[0]);
        individual.setObjective(1, objs[1]);
        // Announce that objective function values are valid
        individual.validObjectiveFunctionsValues = true;
        // Update constraint violations if constraints exist
        if (problem.constraints != null) {
            // Three non-boundary constraints   
            double[] g = new double[3];
            g[0] = x[1] - x[0]; // EMA_long >= EMA =_short (DMAC)
            g[1] = x[3] - x[2]; // EMA_long >= EMA =_short (MACD)
            g[2] = x[2] - x[4]; // Signal <= EMA_short (MACD)
            // Set constraints violations
            for (int i = 0; i < g.length; i++) {
                if (g[i] < 0) {
                    individual.setConstraintViolation(i, g[i]);
                } else {
                    individual.setConstraintViolation(i, 0);
                }
            }
            // Announce that objective function values are valid
            individual.validConstraintsViolationValues = true;
        }
        // Increase Evaluations Count by One (counted per individual)
        individualEvaluationsCount++;
    }
    
    // TODO: calculate the value of  the four indicator functions using the values in the chromosome
    double[] calculateindicators(double[] x) {
    	throw new UnsupportedOperationException();
    	
    	/* DEMAC */
    	
    	/* MADC */
    	
    	/* RSI */
    	
    	/* MARSI */
    }

    // TODO: run the simulation using the new indicator values and output the value of the objective functions
    double[] getObjectives(double[] indicators) {
    	throw new UnsupportedOperationException();
    	
    	/* Annual Return */
    	
    	/* Sharpe Ratio */
    }
}
