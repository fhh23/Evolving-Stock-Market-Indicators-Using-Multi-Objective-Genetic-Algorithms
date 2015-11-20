package evaluators;

import emo.Individual;
import emo.OptimizationProblem;
import emo.VirtualIndividual;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;

import java.util.List;
import java.util.ArrayList;

import javax.xml.stream.XMLStreamException;
import net.sourceforge.jeval.EvaluationException;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

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
	private List<Double> stockData_highs;
	private List<Double> stockData_lows;
	private List<Double> stockData_opens;
	private List<Double> stockData_closes;
	
	public FinancialEvaluator() throws IOException {
		// Initialize the ArrayLists
		stockData_highs = new ArrayList<Double>();
		stockData_lows = new ArrayList<Double>();
		stockData_opens = new ArrayList<Double>();
		stockData_closes = new ArrayList<Double>();
		
		// TODO: Determine the relative path for this file
		// Change this line for the location on your computer!
		File djiaDataFile = new File("C:/Users/breif/Documents/MSU/CSE848/cse_848_project/DJI_Data/AAPL.csv");
		
		// Read the data into the ArrayLists
		// Note: the CSV file should be formatted with: (a) headers as the first line,
		// (b) data starting at the second line, (c) no extra lines at the end of the file,
		// (d) the data is ordered from earliest date to most current date!
		CSVParser parser = CSVParser.parse(djiaDataFile, StandardCharsets.UTF_8, CSVFormat.RFC4180);
		for (CSVRecord record : parser)
		{
			// Ignore the first record because it contains the data headers
			if(parser.getRecordNumber() != 1)
			{
				stockData_highs.add(Double.parseDouble(record.get(2)));
				stockData_lows.add(Double.parseDouble(record.get(3)));
				stockData_opens.add(Double.parseDouble(record.get(1)));
				stockData_closes.add(Double.parseDouble(record.get(4)));
			}
		}
	}
	
	@Override
	public double[] getReferencePoint() 
    {
    	throw new UnsupportedOperationException("Nadir point is not defined for this problem.");
    }

	@Override
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
        double[] objs = getObjectives(evaluateIndicators(x));
        individual.setObjective(0, objs[0]);
        individual.setObjective(1, objs[1]);
        // Announce that objective function values are valid
        individual.validObjectiveFunctionsValues = true;
        // Update constraint violations if constraints exist
        if (problem.constraints != null) {
            // Three non-boundary constraints   
            double[] g = new double[3];
//            g[0] = x[1] - x[0]; // EMA_long >= EMA =_short (DMAC)
//            g[1] = x[3] - x[2]; // EMA_long >= EMA =_short (MACD)
//            g[2] = x[2] - x[4]; // Signal <= EMA_short (MACD)
            // Temporary substitute for the three constraints to use in the problem
            // TODO: remove these and uncomment the above
            g[0] = 0; g[1] = 0; g[2] = 0;
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
    
    /*
     * Calculate the Simple Moving Average (SMA) of the stock for the specified number of days
     * Input: n, the number of days for look-back period
     *        averageType, the values to use when calculating the average
     */
    double calculateSMA( int n, String averageType )
    {
    	double sum = 0;
    	for( int i = 1; i <= n; i++ ) {
    		if (averageType.equals("price"))
    		{
    			// TODO: read the closing price from the given data
    			sum += 0;
    		} else if (averageType.equals("marsi"))
    		{
    			// TODO: implement this option
    			sum += 0;
    		} 
    	}
    	return sum / ((double) n);
    }
    
    /*
     *  Calculate the Exponential Moving Average (EMA) of the stock for the specified number of days
     *  Input: n, the number of days for look-back period
     *         averageType, the values to use when calculating the average
     */
    double calculateEMA( int n, String averageType ) {
    	// TODO: implement!!
    	return 0;
    }
    
    /*
     * Calculate the Relative Strength Index (RSI) of the stock for the specified number of days
     * Input: currDay, the current data for which we are determining the buy/sell signal
     *        n, the number of days for look-back period
     */
    double calculateRSI( int currDay, int n )
    {
    	double upDays = 0;
    	double downDays = 0;

    	for ( int observeDay = currDay; observeDay > currDay - n; observeDay-- ) {
    		if (observeDay < 0) {
    			// Decrement n so that the relative strength calculation is not affected by the non-existent value
    			n -= 1;
    		} else if (stockData_closes.get(observeDay) < stockData_opens.get(observeDay)) {
    			downDays++;
    		} else if (stockData_closes.get(observeDay) > stockData_opens.get(observeDay)){
    			upDays++;
    		}
    		// Do nothing if there was no change during the day
    	}
    	double rs = upDays/downDays;
    	double rsi = 100 - (100.0 / (1.0 + rs));

    	return rsi;
    }
    
    /*
     * Evaluate buy/sell signals from the technical indicators based on the chromosome 
     * Input: x, the chromosome containing the parameters for all technical indicators
     * Output: buy/sell signals of the four technical indicators
     *         -1 = sell, 0 = no action, 1 = buy
     */ 
    int[] evaluateIndicators( double[] x ) {
    	// throw new UnsupportedOperationException();
		// TODO: change this variable type so that it can hold buy/sell signals for all four indicators
		int[] signals = new int[stockData_closes.size()];
    	
    	for ( int currentDay = 0; currentDay < stockData_closes.size(); currentDay++ ) {

    		/* DEMAC */
//    		// Calculate EMA_short
//    		double ema_short_demac = calculateEMA((int) Math.floor(x[0]), "price");
//    		// Calculate EMA_long
//    		double ema_long_demac = calculateEMA((int) Math.floor(x[1]), "price");
//    		if (ema_short_demac > ema_long_demac) {
//    			signals[0] = 1;
//    		} else if (ema_short_demac < ema_long_demac){
//    			signals[0] = -1;
//    		} else {
//    			signals[0] = 0;
//    		}

    		/* MACD */
//    		double ema_short_macd = calculateEMA((int) Math.floor(x[2]), "price");
//    		double ema_long_macd = calculateEMA((int) Math.floor(x[3]), "price");
//    		double macd_line = ema_short_macd - ema_long_macd;
//    		double signal_line = calculateEMA((int) Math.floor(x[4]), "macd");
//    		// TODO: set the buy/sell signal value
//    		signals[1] = 0;

    		/* RSI */
    		// Change to x[0] to x[5] for the full set of indicators
    		double rsi_rsi = calculateRSI(currentDay, (int) Math.floor(x[0])); 
    		if (rsi_rsi < x[1]) { 				 // Change to x[6] for the full set of indicators
    			signals[2] = 1; // Over-sold condition: buy
    		} else if (rsi_rsi > x[2]) {		 // Change to x[7] for the full set of indicators
    			signals[2] = -1; // Over-bought condition: sell
    		} else {
    			signals[2] = 0;
    		}

    		/* MARSI */
//    		double rs_marsi = 0; // TODO: needs to read in data from the text file
//    		double rsi_marsi = 100 - (100.0 / (1.0 + rs_marsi));
//    		// TODO: determine how the RSI indicator relates to the MARSI calculation
//    		double marsi = calculateSMA((int) Math.floor(x[11]), "marsi");
//    		if (marsi < x[9]) {
//    			signals[3] = 1; // Oversold condition: buy
//    		} else if (rsi_rsi > x[10]) {
//    			signals[3] = -1; // Overbought condition: sell
//    		} else {
//    			signals[3] = 0;
//    		}

    	}
    	
		return signals;
    }

    // TODO: run the simulation using the new indicator values and output the value of the objective functions
    double[] getObjectives( int[] indicators ) {
    	throw new UnsupportedOperationException();
    	
    	/* Annual Return */
    	
    	/* Sharpe Ratio */
    }
}
