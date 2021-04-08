package main.java.decon_coreg;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.lang3.time.DurationFormatUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import main.java.decon_coreg.CellCount;

public class Deconvolution {
	private String outputFolder;
	public CellCount cellCounts;
	public ExpressionData expressionData;
	public HashMap<String,ArrayList<String>> geneGenePairs;
	private CommandLineOptions commandLineOptions;

	public Deconvolution(CommandLineOptions commandLineOptions) {
		this.commandLineOptions = commandLineOptions;
	}

	/*
	 * Read all the input data
	 * 
	 * @throw IllegalAccessException	If file can't be accessed
	 * @throw IOException	If file can not be read
	 */
	public void readInputData() throws IllegalAccessException, IOException {
		outputFolder = commandLineOptions.getOutfolder();
		cellCounts = new CellCount();
		cellCounts.readCellCountData(commandLineOptions.getCellcountFile());
		
		geneGenePairs = Utils.parseSnpPerGeneFile(commandLineOptions.getGenesToTestFile());
		String expressionFile = commandLineOptions.getExpressionFile();
		DeconvolutionLogger.log.info(String.format("Parse expression data from %s",expressionFile));
		expressionData = new ExpressionData(expressionFile);
		DeconvolutionLogger.log.info("Done");

	}

	/**
	 * For each of the gene-SNP pair in the SnpsToTestFile run deconvolution
	 * 
	 * @param commandLineOptions	commandLineOptions to run it with 
	 * 
	 * @return Deconvolution results
	 * 
	 * @throws RuntimeException
	 * @throws IllegalAccessException 
	 * @throws IOException 
	 */
	public List<DeconvolutionResult> runDeconPerGeneGenePair() throws RuntimeException, IllegalAccessException, IOException{
		int whileIndex = 0;
		long time = System.currentTimeMillis();
		List<DeconvolutionResult> deconvolutionResults = new ArrayList<DeconvolutionResult>();
		int genePairsTotal = 0;
		HashMap<String, double[]> geneExpressionLevels = expressionData.getGeneExpression();
		// model: geneY ~ cc1 + cc2 + cc1*geneX + cc2*geneX
		for(String geneY : geneGenePairs.keySet()){
			for(String geneX : geneGenePairs.get(geneY)){
				if(commandLineOptions.getTestRun() && whileIndex == 100){
					break;
				}
				if (whileIndex % 500 == 0) {
					long completedIn = System.currentTimeMillis() - time;
					DeconvolutionLogger.log.info(String.format("Processed %d gene-gene pairs - %s", whileIndex, DurationFormatUtils.formatDuration(completedIn, "HH:mm:ss:SS")));
				}
				++whileIndex;
				String genePairName = geneY+'_'+geneX;
				++genePairsTotal;

				double[] expressionLevelsGeneY = geneExpressionLevels.get(geneY);
				double[] expressionLevelsGeneX = geneExpressionLevels.get(geneX);
				
				if(expressionLevelsGeneY == null || expressionLevelsGeneX == null){
					DeconvolutionLogger.log.info(String.format("Error: Gene %s and %s are included in gene/gene combinations to test, but at least one is not available in the expression file!",geneY, geneX));
					throw new RuntimeException(String.format("Gene %s and %s included in gene/gene combinations to test, but not available in the expression file!",geneY, geneX));
				}
				DeconvolutionResult deconResult = deconvolution(expressionLevelsGeneY, expressionLevelsGeneX, genePairName);
				deconvolutionResults.add(deconResult);
			}
		}
		DeconvolutionLogger.log.info(String.format("Total: %d",genePairsTotal));
		return deconvolutionResults;
	}

	/**
	 * Write the deconvolution results
	 * 
	 * @param deconvolutionResult The deconvolution result
	 */
	public void writeDeconvolutionResults(List<DeconvolutionResult> deconvolutionResults) throws IllegalAccessException, IOException{
		List<String> celltypes = cellCounts.getAllCelltypes();
		String header = "\t"+Utils.listToTabSeparatedString(celltypes, "_pvalue");

		DeconvolutionLogger.log.info("Getting decon result with full model info for writing the header");
		// celltypes.size()*2 because there are twice as many betas as celltypes (CC% & CC%:GT)
		InteractionModelCollection firstInteractionModelCollection = deconvolutionResults.get(0).getInteractionModelCollection();
		InteractionModel bestFullModelForHeaderOnly = firstInteractionModelCollection.getBestFullModel();

		for(int i = 1; i < cellCounts.getNumberOfCelltypes()*2 + 1; ++i){
			header += "\tBeta" + Integer.toString(i) +"_"+bestFullModelForHeaderOnly.getIndependentVariableNames().get(i-1);
		}

		if(commandLineOptions.getWholeBloodCorrelation()){
			header += "\tSpearman correlation expression~gene\tSpearman correlation p-value";
		}

		
		//header += "\tStandardError";
		List<String> output = new ArrayList<String>();
		output.add(header);
		for(DeconvolutionResult deconvolutionResult : deconvolutionResults){
			InteractionModelCollection interactionModelCollection = deconvolutionResult.getInteractionModelCollection();

			String results = "";
			results += deconvolutionResult.getGenePairName()+"\t"+Utils.listToTabSeparatedString(deconvolutionResult.getPvalues());
			InteractionModel bestFullModel = null;

			bestFullModel = interactionModelCollection.getBestFullModel();


			double[] estimateRegressionParameters = bestFullModel.getEstimateRegressionParameters();

			// check what the genotype configuration is and the beta of the interaction term. 
			// If genotype configuration == 0 and beta == positive, dosage2 effect = positive
			// If genotype configuration == 1 and beta == negative, dosage2 effect = positive
			// else is negative
			int numberOfCelltypes = cellCounts.getNumberOfCelltypes();
			// first write out the beta of the cell proportion term
			for(int i = 0; i < numberOfCelltypes; ++i){
				results += "\t"+estimateRegressionParameters[i];
			}
			
			// then write out cell proportion-genotype interaction term with correct sign
			for(int i = 0; i < numberOfCelltypes; ++i){
				double interactionTermCurrentCelltype = estimateRegressionParameters[i+numberOfCelltypes];
				results += "\t"+interactionTermCurrentCelltype;
			}

			//results += "\t"+bestFullModel.getGenotypeConfiguration();
			//for(String celltype : cellCounts.getAllCelltypes()){
			//	InteractionModel bestCtModel = deconvolutionResult.getInteractionModelCollection().getBestCtModel(celltype); 
			//	results += "\t"+bestCtModel.getGenotypeConfiguration();
			//}
			if(commandLineOptions.getWholeBloodCorrelation()){
				results += "\t"+deconvolutionResult.getWholeBloodCorrelation();
				results += "\t"+deconvolutionResult.getWholeBloodCorrelationpvalue();
			}

			//results += "\t"+bestFullModel.getEstimatedStandardError();
			output.add(results);	
		}

		Path file = Paths.get(outputFolder+"/"+commandLineOptions.getOutfile());
		Files.write(file, output, Charset.forName("UTF-8"));

		DeconvolutionLogger.log.info(String.format("Deconvolution output written to %s", file.toAbsolutePath()));
		DeconvolutionLogger.log.info(String.format("Files with additional info in  %s", outputFolder));
	}

	/**
	 * Make the linear regression models and then do an Anova of the sum of
	 * squares
	 * 
	 * Full model: geneY ~ celltype_1 + celltype_2 + ... + celltype_n + celltype_1:geneX +
	 *             		   celltype_2:geneX + ... + celltype_n:gene X <- without intercept
	 * 
	 * Compare with anova to: 
	 * geneY ~ celltype_1 + celltype_2 + celtype_n + celltype_1: + celltype_2:geneX + .. + 
	 * 		    celltype_n-1:geneX
	 * 
	 * @param expressionGeneY A vector with the expression value per sample for geneY in geneY ~ cc1 + cc1:geneX
	 * @param exp5essionGeneX A vector with the expression value per sample for geneX in geneY ~ cc1 + cc1:geneX
	 * 
	 * @param genotypes A vector with the expression levels of all
	 * samples for *one* eQTL-gene pair. This should include qtl names as in first column, and sample names in first row
	 * 
	 * @param genePairName Name of the gene pair
	 * 
	 * @return A list with for each celltype a p-value for the celltype
	 * specific eQTL for one eQTL
	 */
	private DeconvolutionResult deconvolution(double[] expressionGeneY, double[] expressionGeneX, String genePairName) 
			throws RuntimeException, IllegalAccessException, IOException {

		InteractionModelCollection interactionModelCollection = new InteractionModelCollection(cellCounts,
																							   commandLineOptions.getUseOLS());
		interactionModelCollection.setGenePairName(genePairName);
		interactionModelCollection.setExpressionValues(expressionGeneY, expressionGeneX);
		
		/**
		 * For each cell type model, e.g. ctModel 1 -> y = neut% + mono% + neut%:GT; ctModel 2 -> y = neut% + mono% + mono%:GT, one for each cell type, 
		 * where the interaction term (e.g mono%:GT) of the celltype:genotype to test is removed, calculate and save the observations in an observation vector
		 * where the observation vector for the example ctModel 1 is
		 *  
		 * 		celltypeModel = [[sample1_neut%, sample1_mono%, sample1_neut%*sample1_genotype], [sample2_neut%, sample2_mono%, sample2_neut%*sample2_genotype]]
		 *  
		 * with for each sample a cellcount percentage for each cell type and the genotype of the QTL that is being testetd. 
		 * 
		 * Using this observation vector calculate the sum of squares and test with Anova if it is significantly different from the sum of squares of the full model. 
		 * Here the full model includes all interaction terms of the cell type models, e.g. fullModel -> y = neut% + mono% + neut%:GT + mono%:GT so the observation vector
		 * 
		 * 		fullModel = [[sample1_neut%, sample1_mono%, sample1_neut%*sample1_genotype, sample1_mono%*sample1_genotype], [sample2_neut%, ..., etc]]
		 * 
		 */
		interactionModelCollection.createObservedValueMatricesFullModel();
		interactionModelCollection.findBestFullModel();		
		interactionModelCollection.createObservedValueMatricesCtModels();
		interactionModelCollection.findBestCtModel();
		calculateDeconvolutionPvalue(interactionModelCollection);

		double wholeBloodCorrelation = 0;
		double wholeBloodCorrelationPvalue = 0;
		if(commandLineOptions.getWholeBloodCorrelation()){
			// if true calculate spearman correlation between gene Y and gene X
			wholeBloodCorrelation = new SpearmansCorrelation().correlation(expressionGeneY, expressionGeneX);
			wholeBloodCorrelationPvalue = Statistics.calculateSpearmanTwoTailedPvalue(wholeBloodCorrelation, cellCounts.getNumberOfSamples());
		}
		DeconvolutionResult deconResult =  new DeconvolutionResult();

		interactionModelCollection.cleanUp();
		deconResult = new DeconvolutionResult(interactionModelCollection, wholeBloodCorrelation, wholeBloodCorrelationPvalue);
		return deconResult;
	}


	/**
	 * get pvalue for each ctmodel
	 * 
	 * @param interactionModelCollection InteractionModelCollection object that has fullModel and ctModels for ANOVA comparison
	 */
	private void calculateDeconvolutionPvalue(InteractionModelCollection interactionModelCollection) 
			throws IllegalAccessException, IOException {
		for (int modelIndex = 0; modelIndex < cellCounts.getNumberOfCelltypes(); ++modelIndex) {
			String celltypeName = cellCounts.getCelltype(modelIndex);
			InteractionModel fullModel;

			fullModel = interactionModelCollection.getBestFullModel();

			InteractionModel ctModel = interactionModelCollection.getBestCtModel(celltypeName);
			double pval = Statistics.anova(fullModel.getSumOfSquares(), ctModel.getSumOfSquares(), 
					fullModel.getDegreesOfFreedom(),ctModel.getDegreesOfFreedom(), 
					true);

			ctModel.setPvalue(pval);
			interactionModelCollection.setPvalue(pval, ctModel.getCelltypeName());
			// TODO: why is this method called twice?
			interactionModelCollection.setPvalue(pval,interactionModelCollection
					.getBestCtModel(cellCounts.getCelltype(modelIndex)).getModelName());

		}

	}
}
