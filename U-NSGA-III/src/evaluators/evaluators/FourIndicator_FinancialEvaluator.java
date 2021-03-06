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
* File Description: Evaluator for the full stock market predictor system with
*     four technical indicators (DEMAC, MACD, RSI, MARSI)
* Associated XML file: fourIndicator_financial.xml 
* 
*/
public class FourIndicator_FinancialEvaluator extends IndividualEvaluator 
{
	private List<Double> stockData_highs;
	private List<Double> stockData_lows;
	private List<Double> stockData_opens;
	private List<Double> stockData_closes;
	private List<Double> returns;
	private int[][] signals;
	
	public FourIndicator_FinancialEvaluator() throws IOException {
		// Initialize the ArrayLists to contain the stock data
		stockData_highs = new ArrayList<Double>();
		stockData_lows = new ArrayList<Double>();
		stockData_opens = new ArrayList<Double>();
		stockData_closes = new ArrayList<Double>();
		
		// TODO: Determine the relative path for this file
		// Change this line for the location on your computer!
		File djiaDataFile = new File("C:/Users/breif/Documents/MSU/CSE848/cse_848_project/DJI_Data/AAPL_Short.csv");
		
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
        individual.setObjective(2, objs[2]);
        // Announce that objective function values are valid
        individual.validObjectiveFunctionsValues = true;
        // Update constraint violations if constraints exist
        if (problem.constraints != null) {
            // Three non-boundary constraints   
            double[] g = new double[3];
            g[0] = Math.round(x[1]) - Math.round(x[0]); // EMA_long >= EMA =_short (DMAC)
            g[1] = Math.round(x[3]) - Math.round(x[2]); // EMA_long >= EMA =_short (MACD)
            g[2] = Math.round(x[2]) - Math.round(x[4]); // Signal <= EMA_short (MACD)
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
     * Input: n, the number of days for look-back period (n >= 1)
     *        dataValues, the values to use for the calculation of the SMA
     * Output: Simple Moving Average of the given data over all valid days in the look-back period
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
     *  Input: n, the number of days for look-back period (n >= 1)
     *         dataValues, the values to use for the calculation of the EMA
     *  Output: n-day Exponential Moving Average of the given data for every trading day
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
     *        n, the number of days for look-back period (n >= 1)
     */
    double calculateRSI( int currDay, int n )
    {
    	double ups = 0; double downs = 0;
    	double gain = 0; 

    	for ( int observeDay = currDay; observeDay > currDay - n; observeDay-- ) {
    		if (observeDay <= 0)
    			break;
    		else
    			gain = stockData_closes.get(observeDay) - stockData_closes.get(observeDay - 1);
    		
    		if (gain > 0) // closing price of observeDay is higher than that of the previous day
    		{
    			ups += gain;
    		}
    		else if (gain < 0) // closing price of observeDay is lower than that of the previous day
    		{
    			downs += Math.abs(gain);
    		}
    		// Do nothing if there was no change during the day
    	}
    	double rsi = 100; // Set to maximum initially
    	// If downs = 0, then RS is infinite, meaning that RSI = 100
    	if (downs != 0) {
    		double rs = ups/downs;
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
    	
    	/* Perform all calculations  and initializations that are not 
    	 * based on the value of currentDay and that only need to be done once */
    	
    	/* DEMAC (initializations only) */
    	double previousEMA_short = 0; double previousEMA_long = 0;
    	double ema_short_demac; double ema_long_demac;
		int short_amt_days = (int) Math.round(x[0]);
		double short_multiplier = 2.0 / (1.0 + short_amt_days);
		int long_amt_days = (int) Math.round(x[1]);
		double long_multiplier = 2.0 / (1.0 + long_amt_days);
    	
    	/* MACD */
    	List<Double> ema_short_macd = calculateEMA((int) Math.round(x[2]), stockData_closes);
    	List<Double> ema_long_macd = calculateEMA((int) Math.round(x[3]), stockData_closes);
    	List<Double> macd_line = new ArrayList<Double>();
    	if (ema_short_macd.size() == ema_long_macd.size()) // Error Checking: this should never be false
    	{
    		for (int idx = 0; idx < ema_short_macd.size(); idx ++)
    		{
    			macd_line.add(ema_short_macd.get(idx) - ema_long_macd.get(idx));
    		}
    	}
    	List<Double> signal_line = calculateEMA((int) Math.round(x[4]), macd_line);
		
		/* MARSI */
		// Precalculation of the RSI values provides significant speedup in SMA calculation
    	List<Double> marsi_rsi_values = new ArrayList<Double>();
    	// For loop skips day 0 of the data because n for RSI calculations is >= 1
    	marsi_rsi_values.add(0.0);
    	for ( int currentDay = 1; currentDay < stockData_closes.size(); currentDay++ )
    	{
    		marsi_rsi_values.add(calculateRSI(currentDay, (int) Math.round(x[8])));
    	}
    	
    	/* Perform all evaluator calculations for the current day and set all 4 buy/sell signals */
    	for ( int currentDay = 0; currentDay < stockData_closes.size(); currentDay++ ) {

    		/* DEMAC (Double EMA Crossovers) */
    		
    		// Calculate the EMA values based on the closing stock price data for the current day
    		// EMAs are calculated inline to save storage space (no need to save every value)
    		if (short_amt_days > currentDay)
    		{
    			ema_short_demac = calculateSMA(currentDay, short_amt_days, stockData_closes);
    		}
    		else // currentDay >= short_amt_days
    		{
    			ema_short_demac = short_multiplier * stockData_closes.get(currentDay) + (1.0 - short_multiplier) * previousEMA_short;
    		}
			
    		if (long_amt_days > currentDay)
    		{
    			ema_long_demac = calculateSMA(currentDay, long_amt_days, stockData_closes);
    		}
    		else // currentDay >= long_amt_days
    		{
    			ema_long_demac = long_multiplier * stockData_closes.get(currentDay) + (1.0 - long_multiplier) * previousEMA_long;
    		}
    		
    		// Determine the buy/sell signal based on the two calculated values
    		if ((ema_short_demac > ema_long_demac) && (previousEMA_short <= previousEMA_long)) {
    			signals[currentDay][0] = 1;
    		} else if ((ema_short_demac < ema_long_demac) && (previousEMA_short >= previousEMA_long)) {
    			signals[currentDay][0] = -1;
    		} else {
    			signals[currentDay][0] = 0;
    		}
    		
    		previousEMA_short = ema_short_demac;
    		previousEMA_long = ema_long_demac;

    		/* MACD (Moving Average Convergence/Divergence) */
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
    		

    		/* RSI (Relative Strength Index) */
    		if ( currentDay != 0 )
			{
				double rsi = calculateRSI(currentDay, (int) Math.round(x[5])); 
				if (rsi < x[6]) {
					signals[currentDay][2] = 1; // Over-sold condition: buy
				} else if (rsi > x[7]) {
					signals[currentDay][2] = -1; // Over-bought condition: sell
				} else {
					signals[currentDay][2] = 0;
				}
			}
			else
			{
				signals[currentDay][2] = 0;
			}

    		/* MARSI */
			if ( currentDay != 0 )
			{
				int sma_days = (int) Math.round(x[11]);
				int observeDay; double rsi_sum = 0;
				for (observeDay = currentDay; observeDay > currentDay - sma_days; observeDay--)
				{
					if (observeDay <= 0)
						break;
					rsi_sum  += marsi_rsi_values.get(observeDay);
				}
				double marsi = rsi_sum / ((double) (currentDay - observeDay));
				if (marsi < x[9]) {
					signals[currentDay][3] = 1; // Over-sold condition: buy
				} else if (marsi > x[10]) {
					signals[currentDay][3] = -1; // Over-bought condition: sell
				} else {
					signals[currentDay][3] = 0;
				}
			}
			else
			{
				signals[currentDay][3] = 0;
			}

    	}
    }

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
    		// Majority voting implementation
    		int determinedSignal = signals[currentDay][0] + signals[currentDay][1] + signals[currentDay][2] + signals[currentDay][3];
    		if ( determinedSignal >= 1 ) {
    			if ((stockData_opens.get(currentDay) < wallet) && (buy == 0)) {
    				buyvalue = stockData_opens.get(currentDay);
    				wallet -= buyvalue;
    				sell = 0;
    				buy = 1;
    			}
    		}
    		else if ( determinedSignal <= -1 ) {
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
    	double sharpe_ratio = 0;
    	if (returns.size() > 4)
    	{
    		sharpe_ratio = average/stddev;
    	}
    	// else: sharpe_ratio remains equal to 0
    	
    	double[] objs = {(-1) * annual_return, (-1) * sharpe_ratio, returns.size()};
		return objs;
    }
}
