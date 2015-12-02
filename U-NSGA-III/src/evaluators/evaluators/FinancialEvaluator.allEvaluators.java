package evaluators;

import emo.Individual;
import emo.OptimizationProblem;
import emo.VirtualIndividual;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

import javax.xml.stream.XMLStreamException;

import net.sourceforge.jeval.EvaluationException;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

import parsing.IndividualEvaluator;
import parsing.InvalidOptimizationProblemException;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

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
	private List<Double> returns;
	private int[][] signals;
	
	public FinancialEvaluator() throws IOException {
		// Initialize the ArrayLists to contain the stock data
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
		signals = new int[stockData_closes.size()][4];
	}
	
	@Override
	public double[] getReferencePoint() 
    {
		//throw new UnsupportedOperationException("Nadir point is not defined for this problem.");
		return null;
    }

	@Override
	public double[] getIdealPoint() 
    {
		//throw new UnsupportedOperationException("Ideal point is not defined for this problem.");
		return null;
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
        evaluateIndicators(x);
        double[] objs = getObjectives(x);
        individual.setObjective(0, objs[0]);
        individual.setObjective(1, objs[1]);
        // Announce that objective function values are valid
        individual.validObjectiveFunctionsValues = true;
        // Update constraint violations if constraints exist
        if (problem.constraints != null) {
            // Three non-boundary constraints   
            double[] g = new double[5];
            g[0] = x[1] - x[0]; // EMA_long >= EMA =_short (DMAC)
            g[1] = x[3] - x[2]; // EMA_long >= EMA =_short (MACD)
            g[2] = x[2] - x[4]; // Signal <= EMA_short (MACD)
            g[3] = x[12] + x[13] + x[14] + x[15] - 1.0; // Sum of the indicator weights >= 1
            g[4] = 1.0 - x[12] - x[13] - x[14] - x[15]; // 1 >= Sum of the indicator weights
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
     * Calculate the Simple Moving Average (SMA) of the data provided as input for the specified number of days
     * Input: n, the number of days for look-back period (n >= 0)
     *        dataValues, the values to use for the calculation of the EMA
     */
    double calculateSMA( int currDay, int n, List<Double> dataValues )
    {
    	int observeDay; double sum = 0; 
    	for ( observeDay = currDay; observeDay > currDay - n; observeDay-- ) {
    		if (observeDay < 0)
    			break;
    		sum += dataValues.get(observeDay);
    	}
    	// Note: currDay - observeDay is used instead of n for a correct calculation in the case that
    	// the loop is broken early
    	return sum / ((double) (currDay - observeDay));
    }
    
    /*
     *  Calculate the Exponential Moving Average (EMA) of the data provided as input for the specified number of days
     *  Input: n, the number of days for look-back period
     *         dataValues, the values to use for the calculation of the EMA
     */
    List<Double> calculateEMA( int n, List<Double> dataValues ) {
    	double currentMultiplier = 2.0 / (1.0 + n);
    	double ema; double previousEMA = 0;
    	
    	List<Double> stockData_EMAs = new ArrayList<Double>();
    	
		// Calculate the Simple Moving Average over the first n days of data
    	for ( int observeDay = 0; observeDay < n; observeDay++ )
    	{
    		stockData_EMAs.add(calculateSMA(observeDay, n, dataValues));
    	}
    	previousEMA = stockData_EMAs.get(n - 1);
    	for ( int observeDay = n; observeDay < dataValues.size(); observeDay++ )
    	{
    		ema = currentMultiplier * dataValues.get(observeDay) + (1.0 - currentMultiplier) * previousEMA;
    		stockData_EMAs.add(ema);
    		previousEMA = ema;
    	}
    	
    	return stockData_EMAs;
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
    			break;
    		} else if (stockData_closes.get(observeDay) < stockData_opens.get(observeDay)) {
    			downDays++;
    		} else if (stockData_closes.get(observeDay) > stockData_opens.get(observeDay)){
    			upDays++;
    		}
    		// Do nothing if there was no change during the day
    	}
    	double rsi = 100; // Set to maximum initially
    	// If downDays = 0, then RS is infinite, meaning that RSI = 100
    	if (downDays != 0) {
    		double rs = upDays/downDays;
    		rsi = 100 - (100.0 / (1.0 + rs));
    	}
    	
    	return rsi;
    }
    
    /*
     * Evaluate buy/sell signals from the technical indicators based on the chromosome 
     * Input: x, the chromosome containing the parameters for all technical indicators
     * Output: buy/sell signals of the four technical indicators
     *         -1 = sell, 0 = no action, 1 = buy
     */ 
    void evaluateIndicators( double[] x ) {
    	
    	for ( int currentDay = 0; currentDay < stockData_closes.size(); currentDay++ ) {

    		/* DEMAC */
    		// Calculate EMA_short
    		double ema_short_demac = (calculateEMA((int) Math.round(x[0]), stockData_closes)).get(currentDay);
    		// Calculate EMA_long
    		double ema_long_demac = (calculateEMA((int) Math.round(x[1]), stockData_closes)).get(currentDay);
    		if (ema_short_demac > ema_long_demac) {
    			signals[currentDay][0] = 1;
    		} else if (ema_short_demac < ema_long_demac){
    			signals[currentDay][0] = -1;
    		} else {
    			signals[currentDay][0] = 0;
    		}

    		/* MACD */
    		List<Double> ema_short_macd = calculateEMA((int) Math.round(x[2]), stockData_closes);
    		List<Double> ema_long_macd = calculateEMA((int) Math.round(x[3]), stockData_closes);
    		List<Double> macd_line = new ArrayList<Double>();
    		if (ema_short_macd.size() == ema_long_macd.size())
    		{
    			for (int idx = 0; idx < ema_short_macd.size(); idx ++)
    			{
    				macd_line.add(ema_short_macd.get(idx) - ema_long_macd.get(idx));
    			}
    		}
    		List<Double> signal_line = calculateEMA((int) Math.round(x[4]), macd_line);
    		if (currentDay != 0)
    		{
    			if (((macd_line.get(currentDay) > signal_line.get(currentDay)) && (macd_line.get(currentDay - 1) <= signal_line.get(currentDay - 1))) || 
    					((macd_line.get(currentDay) > 0) && (macd_line.get(currentDay - 1) <= 0)))
    			{
    				signals[currentDay][1] = 1;
    			}
    			else if (((macd_line.get(currentDay) < signal_line.get(currentDay)) && (macd_line.get(currentDay - 1) >= signal_line.get(currentDay - 1))) || 
    					((macd_line.get(currentDay) < 0) && (macd_line.get(currentDay - 1) >= 0)))
    			{
    				signals[currentDay][1] = -1;
    			}
    			else
    				signals[currentDay][1] = 0;
    		}
    		else
				signals[currentDay][1] = 0;
    		

    		/* RSI */
    		double rsi = calculateRSI(currentDay, (int) Math.round(x[5])); 
    		if (rsi < x[6]) {
    			signals[currentDay][2] = 1; // Over-sold condition: buy
    		} else if (rsi > x[7]) {
    			signals[currentDay][2] = -1; // Over-bought condition: sell
    		} else {
    			signals[currentDay][2] = 0;
    		}

    		/* MARSI */
    		int sma_days = (int) Math.round(x[11]);
    		double rsi_sum = 0;
    		// TODO: account for observeDays less than 0 here
    		for (int observeDay = currentDay - sma_days; observeDay < currentDay; observeDay++)
    		{
    			rsi_sum  += calculateRSI(observeDay, (int) Math.round(x[9]));
    		}
    		double marsi = rsi_sum / sma_days;
    		if (marsi < x[9]) {
    			signals[currentDay][3] = 1; // Oversold condition: buy
    		} else if (marsi > x[10]) {
    			signals[currentDay][3] = -1; // Overbought condition: sell
    		} else {
    			signals[currentDay][3] = 0;
    		}

    	}
    }

    // TODO: run the simulation using the new indicator values and output the value of the objective functions
    double[] getObjectives( double[] x ) {
    	
    	int n = stockData_closes.size(); //number of trading days
    	// TODO: base this value on the values in the data set
    	double capital = 20000; // amount initially invested
    	double wallet = capital;
    	returns = new ArrayList<Double>();
    	
    	// Initialize the signals to indicate that a buy must happen before a sell of the stock
    	int buy = 0;
    	int sell = 1;
    	double buyvalue = 0;
    	double sellvalue = 0;
    	
    	// Find all returns using indicator voting to determine whether to buy or sell
    	for ( int currentDay = 0; currentDay < stockData_closes.size(); currentDay++ ) {
    		double determinedSignal = (x[12] * signals[currentDay][0] + x[13] * signals[currentDay][1] + x[14] * signals[currentDay][2] + x[15] * signals[currentDay][3]) / 4.0;
    		if ((int) Math.round(determinedSignal) == 1) {
    			if ((stockData_opens.get(currentDay) < wallet) && (buy == 0)) {
    				buyvalue = stockData_opens.get(currentDay);
    				wallet -= buyvalue;
    				sell = 0;
    				buy = 1;
    			}
    		}
    		else if ((int) Math.round(determinedSignal) == -1) { 	
    			if (sell == 0) {
    				sellvalue = stockData_opens.get(currentDay);
    				wallet += sellvalue;
    				buy = 0;
    				sell = 1;
    				returns.add(sellvalue-buyvalue);
    			}
    		}
    				
    	}
    	
    	// Sum all returns
    	double sum = 0;
    	for ( int index = 0; index < returns.size(); index++ ) {
    		sum += returns.get(index);
    	}
    	
    	// Average of returns
    	double average = sum/returns.size();
    	
    	// Standard deviation of returns
    	double stddev_sum = 0;
    	for ( int index = 0; index < returns.size(); index++ ) {
    		stddev_sum += pow(returns.get(index) - average, 2);
    	}
    	double stddev = sqrt(stddev_sum/returns.size());

    	/* Annual Return */
    	double base = wallet / capital;
    	int exponent_denom = n / 250;
    	double annual_return = (pow(base, 1.0/exponent_denom)-1) * 100;
    	
    	/* Sharpe Ratio */
    	double sharpe_ratio = average/stddev;
    	
    	double[] objs = {(-1) * annual_return, (-1) * sharpe_ratio};
		return objs;
    }
}
