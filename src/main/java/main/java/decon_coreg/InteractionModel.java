package main.java.decon_coreg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

/*
 * InteractionModel contains the information (observed values, p-value, regression parameter values etc)
 * of one single interaction model. All the InteractionModels that are used in the deconvolution are collected 
 * in InteractionModelCollection
 */
public class InteractionModel {
	/* independentVariables = the names of the independent variables, e.g. neut%, mono%, neut%:GT */
	private List<String> independentVariableNames = new ArrayList<String>();
	private List<int[]> celltypeVariablesIndex = new ArrayList <int[]>();
	private double[][] observedValues;
	private String modelName;
	// Initialise with number so that we can test if it has been set
	private Double sumOfSquares;
	private Double pvalue;
	private Integer degreesOfFreedom;
	private Integer modelLength;
	private Integer numberOfTerms;
	private double[] estimatedRegressionParameters;
	private double[] residuals;
	@SuppressWarnings("unused")
	private double[] estimatedBetaStandardErrors;
	private double estimatedStandardError;
	private String celltypeName;
	private int sampleSize;
	
	/**
	 * Initialise object by setting the observed values size. Per QTL for each sample the observed values are each term of the 
	 * linear model. E.g. if the model is y = mono% + neut% + mono%:gene, the observedValues are
	 * [mono%, neut%, mono% * gene]
	 * 
	 * @param sampleSize 	number of samples
	 * @param numberOfTerms 	number of terms that the interaction model has
	 */
	public InteractionModel( int sampleSize, int numberOfTerms){
	    this.observedValues = new double[sampleSize][numberOfTerms];
	    this.numberOfTerms = numberOfTerms;
	    this.sampleSize = sampleSize;
	  }
	
	public void addObservedValue( double observedValue, int sampleIndex, int termIndex){
	    this.observedValues[sampleIndex][termIndex] = observedValue;
	  }
	
	/**
	 * Get the observed values. Per QTL for each sample the observed values are each term of the 
	 * linear model. E.g. if the model is y = mono% + neut% + mono%:GT, the observedValues are
	 * [mono%, neut%, mono% * GT]
	 * 
	 * @return Matrix of observed values
	 */
	public double[][] getObservedValues(){
	    return this.observedValues;
	  }
	
	/**
	 * Add the index of the celltype variables of the linear model e.g.
	 * the index of the celltype% and celltype%:GT of the model. If 
	 *    y = neut% + mono% + eos% + neut%:GT + eos%:GT
	 *    celltypeTerms = [[0,3],[1],[2,4]
	 * This can be used to sum up the Beta * variable per cell type  
	 * 
	 * @param values	The index of the celltype variables of the linear model
	 */
	public void addCelltypeVariablesIndex(int[] values){
		celltypeVariablesIndex.add(values);
	  }
	
	/** 
	 * Add the name of the independent variable name at the end of the existing list
	 * of independent variables.
	 * 
	 * @param independentVariables	Variable name to add to independentVariableNames list
	 */
	public void addIndependentVariableName(String independentVariables){
	    this.independentVariableNames.add(independentVariables);
	  }
	
	/**
	 * Add the name of the independent variable name at index of the existing list
	 * of independent variables. 
	 * 
	 * @param index	Index of list at which to add the variable name
	 * 
	 * @param independentVariables	Variable name to add
	 */
	public void addIndependentVariableName(int index, String independentVariables){
	    this.independentVariableNames.add(index, independentVariables);
	  }
	
	/** 
	 * Get a list of the independent variables of the interaction model e.g.
	 * 		[neut%, mono%, neut%:GT]
	 * 
	 * @return List of independent variable names
	 */
	public List<String> getIndependentVariableNames(){
	    return this.independentVariableNames;
	  }
	
	public void setModelName(String modelName){
	    this.modelName = modelName;
	  }
	
	public String getModelName(){
	    return this.modelName;
	}	

	public void setPvalue(double pvalue){
		this.pvalue = pvalue;
	}
	public double getPvalue(){
		return this.pvalue;
	}
	public void setAlltIndependentVariableNames(List<String> list){
		for (int celltypeIndex = 0; celltypeIndex < list.size(); ++celltypeIndex) {
			int[] index = new int[] {celltypeIndex, list.size() + celltypeIndex};
			addCelltypeVariablesIndex(index);
			addIndependentVariableName(celltypeIndex, list.get(celltypeIndex));
			addIndependentVariableName(list.get(celltypeIndex)+":GT");
		}
	}

	public void setModelLength(){
		this.modelLength = this.observedValues.length;
	}
	
	public int getModelLength() throws IllegalAccessException {
		return this.modelLength;
	}

	public void setSumOfSquares(double sumOfSquares) {
		this.sumOfSquares = sumOfSquares;
	}
	
	public double getSumOfSquares() throws IllegalAccessException {
		return(this.sumOfSquares);
	}

	public void setDegreesOfFreedom(int degreesOfFreedom) {
		this.degreesOfFreedom  = degreesOfFreedom;
	}
	public int getDegreesOfFreedom() throws IllegalAccessException {
		return(this.degreesOfFreedom);
	}
	
	public int getNumberOfTerms(){
		return this.numberOfTerms;
	}
	
	
	public void setEstimatedRegressionParameters(double[] estimatedRegressionParamters){
		this.estimatedRegressionParameters = estimatedRegressionParamters;
	}
	
	public double[] getEstimateRegressionParameters() throws IllegalAccessException{
		return(this.estimatedRegressionParameters);
	}
	
	private void setResiduals(double[] residuals) {
		this.residuals = residuals;
	}
	public double[] getResiduals() {
		return(this.residuals);
	}
	
	/**
	 * Calculate the sum of squares, using Non-Negative Linear Regression, given a y expression vector with y ~
	 * model.
	 * Remove of the intercept (equivalent to y ~ model -1 in R) is hard-code in
	 * 
	 * This uses NNLS from Rochester Institute of technology. Documentation here: https://www.cs.rit.edu/~ark/pj/doc/edu/rit/numeric/NonNegativeLeastSquares.html
	 * Can download JAR from: https://www.cs.rit.edu/~ark/pj.shtml#installed
	 * 
	 * @param expressionValues	Vector of expression values
	 * 
	 * @throws IllegalAccessException	Exception thrown when observed values can't be retrieved	
	 */
	public void calculateSumOfSquaresNNLS(double[] expressionValues) throws IllegalAccessException {
		NonNegativeLeastSquares nnls = new NonNegativeLeastSquares();
		
		nnls.newSampleData(expressionValues, this.getObservedValues());

		
		// results contain:
		// normsqr: sqroot of the norm error vector
		// x: the parameters
		// For more, check out the Class documentation
		double[] estimatedRegressionParameters = nnls.estimateRegressionParameters();
		setEstimatedRegressionParameters(estimatedRegressionParameters);

		setSumOfSquares(nnls.calculateResidualSumOfSquares());
		setDegreesOfFreedom(expressionValues.length - (getNumberOfTerms() + 1));
		
		double[] residuals = nnls.estimateResiduals();
		setEstimatedBetaStandardErrors(nnls.estimateBetaStandardErrors(this.sampleSize, this.numberOfTerms,
																	   residuals, estimatedRegressionParameters,
																	   this.getObservedValues()));
		setResiduals(residuals);
	}
	public void setEstimatedBetaStandardErrors(double[] estimatedBetaStandardErrors){
		this.estimatedBetaStandardErrors = estimatedBetaStandardErrors;

	}

	public void setEstimatedStandardError(double estimatedStandardError){
		this.estimatedStandardError = estimatedStandardError;
	}
	public double getEstimatedStandardError(){
		return this.estimatedStandardError;
	}

	public void cleanUp() {
		this.observedValues = null;
		this.residuals = null;
	}

	public void setCelltypeName(String celltypeName) {
		this.celltypeName = celltypeName;
	}
	
	public String getCelltypeName(){
		return this.celltypeName;
	}
	
	/* Calculate the sum of squares, using Ordinary Linear Regression, given a y expression vector with y ~
	 * model. Remove of the intercept (equivalent to y ~ model -1 in R) is hard-code in
	 * 
	 * @param model An InteractionModel object including the y vector expression values and ObservedValues (model)
	 * Such that
	 * test_trait ~ geno_A + lymph% + geno_A:geno_B it can be for one QTL
	 * [[2, 43.4, 86.8], [2, 40.3, 80.6]], for another QTL [[0, 46.7, 0],
	 * [0, 51.5, 0] [0, 48.7, 0]]
	 * @param expressionValues Vector of expression values to use
	 * @return A regression object
	 */
	public void calculateSumOfSquaresOLS(double[] expressionValues) throws IOException, IllegalAccessException {
		// OLS = Ordinary Least Squares
		OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
		// if GetIntercept is false, remove the intercept (Beta1) from the linear model
		regression.setNoIntercept(true);
		regression.newSampleData(expressionValues, this.getObservedValues());
		/*for(int i = 0; i < expressionValues.length; i++) {
		 * 
			System.out.println(expressionValues[i]+"\t"+ this.getObservedValues()[i][0]+"\t"+this.getObservedValues()[i][1]+"\t"+
							   this.getObservedValues()[i][2]+"\t"+ this.getObservedValues()[i][3]+"\t"+this.getObservedValues()[i][4]+
							   "\t"+this.getObservedValues()[i][5]+"\t"+this.getObservedValues()[i][6]);
		}*/
		this.setSumOfSquares(regression.calculateResidualSumOfSquares());
		this.setDegreesOfFreedom(expressionValues.length - (this.getNumberOfTerms() + 1));
		setResiduals(regression.estimateResiduals());
		setEstimatedRegressionParameters(regression.estimateRegressionParameters());
	}
}

