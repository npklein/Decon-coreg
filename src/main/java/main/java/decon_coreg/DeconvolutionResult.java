package main.java.decon_coreg;

import java.util.ArrayList;
import java.util.List;

public class DeconvolutionResult {
	private List<String> celltypes;
	private String genePairName;
	private List<Double> pvalues = new ArrayList<Double>();
	private InteractionModelCollection interactionModelCollection;
	private double wholeBloodCorrelation;
	private double wholeBloodCorrelationpvalue;
	
	public DeconvolutionResult(){};
	
	/**
	 * Set the deconvolutionResult with the InteractionModels.
	 * 
	 * @param interactionModelCollection Collection of the interaction models used to get the deconvolution result
	 * 
	 * @param wholeBloodCorrelation Spearman correlation of expression of the genes in whole blood
	 * 
	 * @param wholeBloodCorrelationpvalue p-value of the spearman correlation of whole blood gene-gene correlation
	 * 
	 * @throws IllegalAccessException	QTL name or p-value can not be retrieved from interactionModelCollection
	 */
	public DeconvolutionResult( InteractionModelCollection interactionModelCollection, double wholeBloodCorrelation, 
								double wholeBloodCorrelationpvalue) throws IllegalAccessException{

		celltypes = interactionModelCollection.getAllCelltypes();
		this.genePairName = interactionModelCollection.getQtlName();
		for (int i = 0; i < celltypes.size(); i++){
			String modelName = celltypes.get(i);
			Double pvalue = interactionModelCollection.getPvalue(modelName);
			this.pvalues.add(pvalue);
		}
		this.interactionModelCollection = interactionModelCollection;
		this.wholeBloodCorrelation = wholeBloodCorrelation;
		this.wholeBloodCorrelationpvalue = wholeBloodCorrelationpvalue;
	}

	public String getGenePairName() throws IllegalAccessException{
		if(this.genePairName == null){
			throw new IllegalAccessException("QTL name not set for this model");
		}
		return this.genePairName;
	}
	
	public List<Double> getPvalues() throws IllegalAccessException{
		if(this.pvalues == null){
			throw new IllegalAccessException("pvalues not set for this model");
		}
		return this.pvalues;
	}
	
	public InteractionModelCollection getInteractionModelCollection() throws IllegalAccessException{
		if(this.interactionModelCollection == null){
			throw new IllegalAccessException("interactionModelCollection not set");
		}
		return this.interactionModelCollection;
	}

	public double  getWholeBloodCorrelation() throws IllegalAccessException{
		return this.wholeBloodCorrelation;
	}
	public double  getWholeBloodCorrelationpvalue() throws IllegalAccessException{
		return this.wholeBloodCorrelationpvalue;
	}
}
