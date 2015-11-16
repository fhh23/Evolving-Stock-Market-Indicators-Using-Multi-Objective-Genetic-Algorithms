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
    /* TODO: override this appropriately */
	public VirtualIndividual[] getParetoFront(int objectivesCount, int n) throws InvalidOptimizationProblemException, XMLStreamException, EvaluationException 
    {
        throw new UnsupportedOperationException("getParetoFront() method must be overrided with appropriate logic.");
    }
    
    public void updateIndividualObjectivesAndConstraints(
            OptimizationProblem problem,
            Individual individual)
            throws EvaluationException 
    {
    	// TODO: see example code in ZDT1Evaluator.java
    }

}
