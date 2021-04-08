package main.java.decon_coreg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/*
 *  Collection of all the interaction models and their shared data (expression etc)
 *  There are n + 1 interaction models, where n = the number of cell types. One full model
 *  with all cell types, and for each cell type one model with the interaction term for that
 *  model removed
 */
public class InteractionModelCollection {
    private double[] expressionValuesGeneY;
    private double[] expressionValuesGeneX;
    private String genePairName;
    private HashMap<String, InteractionModel> interactionModels = new HashMap<String, InteractionModel>();
    private HashMap<String, Double> pvalues = new HashMap<String, Double>();
    private ArrayList<String> fullModelNames = new ArrayList<String>();
    private HashMap<String, ArrayList<String>> ctModelNames = new HashMap<String, ArrayList<String>>();
    private HashMap<String, ArrayList<String>> fullModelNamesByCelltype = new HashMap<String, ArrayList<String>>();
    private String bestFullModelName;
    private HashMap<String, String> bestCtModel = new HashMap<String, String>();
    private CellCount cellCount;
    private HashMap<String, String> modelCelltype = new HashMap<String, String>();
    private List<String> celltypes = new ArrayList<String>();
    private List<String> sampleNames = new ArrayList<String>();
    private Boolean useOLS;


    /*
     * Have to initialize instance with if NNLS or OLS will be used, and for that we need cellCounts
     */
    public InteractionModelCollection(CellCount cellCount, Boolean useOLS) throws IllegalAccessException {
        setCellCount(cellCount);
        this.useOLS = useOLS;

    }

    public List<String> getAllCelltypes() {
        return celltypes;
    }

    /*
     * Get interaction model with modelName
     */
    public InteractionModel getInteractionModel(String modelName) throws IllegalAccessException {
        InteractionModel interactionModel = this.interactionModels.get(modelName);
        return interactionModel;
    }

    /*
     * Remove interaction model with modelName
     */
    public void removeInteractionModel(String modelName) throws IllegalAccessException {
        this.interactionModels.remove(modelName);
    }

    /*
     * Set CellCount
     */
    private void setCellCount(CellCount cellCount) {
        this.cellCount = cellCount;
        celltypes = cellCount.getAllCelltypes();
        sampleNames = cellCount.getSampleNames();
    }

    public CellCount getCellCount() throws IllegalAccessException {
        return this.cellCount;
    }

    /**
     * Set the expression values (y) for all the interaction models.
     *
     * @param expression Expression vector
     */
    public void setExpressionValues(double[] expressionGeneY, double[] expressionGeneX) {
        this.expressionValuesGeneY = expressionGeneY;
        this.expressionValuesGeneX = expressionGeneX;
    }

    /**
     * Get the expression values (y) of all the interaction models.
     *
     * @return Expression vector
     */
    public double[] getExpessionValuesGeneY() {
        return this.expressionValuesGeneY;
    }
    
    /**
     * Get the expression values (y) of all the interaction models.
     *
     * @return Expression vector
     */
    public double[] getExpessionValuesGeneX() {
        return this.expressionValuesGeneX;
    }


    public void setGenePairName(String genePairName) {
        this.genePairName = genePairName;
    }

    public String getQtlName() throws IllegalAccessException {
        return this.genePairName;
    }

    /*
     * Each ctModel will have a p-value from ANOVA test with fullmodel, save it per ctModel
     */
    public void setPvalue(Double pvalue, String modelName) {
        this.pvalues.put(modelName, pvalue);
    }

    public Double getPvalue(String modelName) throws IllegalAccessException {
        Double pvalue = this.pvalues.get(modelName);
        return pvalue;
    }

    public ArrayList<String> getFullModelNames() throws IllegalAccessException {
        return this.fullModelNames;
    }

    public ArrayList<String> getCtModelNames(String celltype) throws IllegalAccessException {
        return this.ctModelNames.get(celltype);
    }

    private void setBestFullModelName(String modelName) {
        this.bestFullModelName = modelName;
    }

    public InteractionModel getBestFullModel() throws IllegalAccessException {
        return this.getInteractionModel(this.bestFullModelName);
    }

    private void setBestCtModel(String celltype, String modelName) {
        this.bestCtModel.put(celltype, modelName);
    }

    // per cell type there is one best Ct model
    public InteractionModel getBestCtModel(String celltype) throws IllegalAccessException {
        return this.getInteractionModel(this.bestCtModel.get(celltype));
    }

    /*
     * Add interaction model to the collections
     */
    private void addInteractionModel(InteractionModel interactionModel, String modelName, Boolean isFullModel) {
        this.interactionModels.put(modelName, interactionModel);
        String cellType = interactionModel.getCelltypeName();
        if (isFullModel) {
            fullModelNames.add(modelName);
            fullModelNamesByCelltype.putIfAbsent(cellType, new ArrayList<String>());
            fullModelNamesByCelltype.get(cellType).add(modelName);
        } else {
            ctModelNames.putIfAbsent(cellType, new ArrayList<String>());
            ctModelNames.get(cellType).add(modelName);
        }
    }

    /*
     * Go through all full models, calculate the regression statistics and
     * select the model with the highest R2 as the new full model
     */
    public void findBestFullModel() throws IllegalAccessException, IOException {
        // set to -1 so that first loop can be initialised
        double sumOfSquares = -1;
        for (String modelName : getFullModelNames()) {
            InteractionModel fullModel = getInteractionModel(modelName);

            if (useOLS) {
                fullModel.calculateSumOfSquaresOLS(getExpessionValuesGeneY());
            } else {
                fullModel.calculateSumOfSquaresNNLS(getExpessionValuesGeneY());
            }

            if (sumOfSquares == -1) {
                sumOfSquares = fullModel.getSumOfSquares();
            }
            if (fullModel.getSumOfSquares() <= sumOfSquares) {
                sumOfSquares = fullModel.getSumOfSquares();
                setBestFullModelName(fullModel.getModelName());
            } else {
                removeInteractionModel(fullModel.getModelName());
            }
        }

    }

    /*
     * Go through all cell type models, calculate the regression statistics and select the model with the highest R2 as the new cell type model
     * TODO: merge with findBestFullModel()
     */
    public void findBestCtModel() throws IllegalAccessException, IOException {
        // set to -1 so that first loop can be initialised
        for (String celltype : celltypes) {
            double sumOfSquares = -1;
            for (String modelName : getCtModelNames(celltype)) {
                InteractionModel ctModel = getInteractionModel(modelName);
                modelCelltype.put(modelName, celltype);

                if (useOLS) {
                    ctModel.calculateSumOfSquaresOLS(getExpessionValuesGeneY());
                } else {
                    ctModel.calculateSumOfSquaresNNLS(getExpessionValuesGeneY());
                }

                if (sumOfSquares == -1) {
                    sumOfSquares = ctModel.getSumOfSquares();
                    setBestCtModel(ctModel.getCelltypeName(), ctModel.getModelName());
                }

                double ctSumOfSquares = ctModel.getSumOfSquares();
                if (ctSumOfSquares <= sumOfSquares) {
                    sumOfSquares = ctSumOfSquares;
                    setBestCtModel(ctModel.getCelltypeName(), ctModel.getModelName());
                }
            }
        }
    }


    /**
     * Construct the observed value matrices that are used for calculating the regression for the full model.
     *
     * @throws IllegalAccessException If cell counts file can not be read
     */
    public void createObservedValueMatricesFullModel()
            throws IllegalAccessException {
        CellCount cellCount = getCellCount();
        int numberOfCelltypes = cellCount.getNumberOfCelltypes();
        int numberOfSamples = cellCount.getNumberOfSamples();
        int numberOfTerms = numberOfCelltypes * 2;

        // things neded for fullModel defined outside of loop because every celltype model (ctModel) has to be compared to it
        InteractionModel fullModel = new InteractionModel(numberOfSamples,
                numberOfTerms);

        String modelName = String.format("fullModel");
        fullModel.setModelName(modelName);
        addInteractionModel(fullModel, modelName, true);

        for (int celltypeIndex = 0; celltypeIndex < numberOfCelltypes; ++celltypeIndex) {
            /** save the index of the variables related to current celltype so that this can be used later to calculate
            * Beta1 celltype% + Beta2 * celltype%:gene. For fullModel not so necesarry as it's always <numberOfCelltypes> away,
            * but for ctModel this is easiest method
            */
            int[] index = new int[]{celltypeIndex, numberOfCelltypes + celltypeIndex};
            fullModel.addCelltypeVariablesIndex(index);
            // add the celltype name at position i so that it gets in front of the celltype:GT
            fullModel.addIndependentVariableName(celltypeIndex, cellCount.getCelltype(celltypeIndex));
            fullModel.addIndependentVariableName(cellCount.getCelltype(celltypeIndex) + ":GT");
        }

        // number of terms + 1 because for full model all cell types are included
        for (int sampleIndex = 0; sampleIndex <= numberOfSamples - 1; ++sampleIndex) {
            for (int celltypeIndex = 0; celltypeIndex < numberOfCelltypes; ++celltypeIndex) {

                double celltypePerc = cellCount.getCellCountPercentages()[sampleIndex][celltypeIndex];
                // if i (cell type index) is the same as m (model index), don't add the interaction term of celltype:GT
                fullModel.addObservedValue(celltypePerc, sampleIndex, celltypeIndex);
                fullModel.addObservedValue(celltypePerc * this.expressionValuesGeneX[sampleIndex],
                							sampleIndex, numberOfCelltypes + celltypeIndex);


            }
            fullModel.setModelLength();
        }
    }

    /**
     * Construct the observed value matrices that are used for calculating the regression
     * <p>
     *
     * @throws IllegalAccessException If cell counts file can not be read
     */
    public void createObservedValueMatricesCtModels()
            throws IllegalAccessException {
        CellCount cellCount = getCellCount();
        int numberOfCelltypes = cellCount.getNumberOfCelltypes();
        int numberOfSamples = cellCount.getNumberOfSamples();
        int interactionTermCount = cellCount.getNumberOfCelltypes();
        // -1 because one interaction term is removed
        int numberOfTerms = (numberOfCelltypes * 2) - 1;
        //System.out.println("Number of cell types: "+numberOfCelltypes);
        //System.out.println("Number of terms: "+numberOfTerms);
        // m = model, there are equally many models as celltypes
        for (int modelIndex = 0; modelIndex < numberOfCelltypes; modelIndex++) {
        	//System.out.println("Model nr: "+modelIndex);
            InteractionModel ctModel = new InteractionModel(numberOfSamples, numberOfTerms);
            // calculate p-value and save it, with other information, in a ctModel object.
            // Then, add it to a list of these models to return as decon results
            String celltypeName = cellCount.getCelltype(modelIndex);
            String modelName = String.format("ctModel_%s", celltypeName);
            ctModel.setModelName(modelName);
            ctModel.setCelltypeName(celltypeName);
            addInteractionModel(ctModel, ctModel.getModelName(), false);
            for (int sampleIndex = 0; sampleIndex <= numberOfSamples - 1; sampleIndex++) {
            	//System.out.println("Sample index: "+sampleIndex);
                for (int celltypeIndex = 0; celltypeIndex < numberOfCelltypes; celltypeIndex++) {
                    // There is one fullModel including all celltypes add values for celltypePerc and interaction term of
                    // celltypePerc * gene so that you get [[0.3, 0.6], [0.4, 0.8], [0.2, 0.4], [0.1, 0.2]]
                    // where numberOfSamples = 1 and numberOfCellTypes = 4 with celltypePerc = 0.3, 0.4, 0.2, and 0.1 and gene = 2
                    // for each cell type is 1 model, celltype% *gene without 1 celltype.
                    // j+1 because j==0 is header
                    double celltype_perc = cellCount.getCellCountPercentages()[sampleIndex][celltypeIndex];
                    ctModel.addObservedValue(celltype_perc, sampleIndex, celltypeIndex);
                    if (sampleIndex == 0) {
                        // add the celltype name at position i so that it gets in front of the celltype:GT, but once
                        try {
                            ctModel.addIndependentVariableName(celltypeIndex, celltypeName);
                        } catch (NullPointerException e) {
                            DeconvolutionLogger.log.info(String.format("Nullpoint exception with celltype %s", celltypeIndex));
                            throw e;
                        }
                    }

                    // if celltypeIndex is the same as m modelIndex, don't add the interaction term of celltype:gene
                    if (celltypeIndex != modelIndex) {
                        // Only add IndependentVariableName once per gene-gene pair (j==0)
                        if (sampleIndex == 0) {

                            // Add the interaction term of celltype:gene
                            ctModel.addIndependentVariableName(cellCount.getCelltype(celltypeIndex) + ":gene");
                            // save the index of the variables related to current celltype so that this can be used later to calculate
                            // Beta1 celltype% + Beta2 * celltype%:gene. For fullModel not so necessary as it's always <numberOfCelltypes> away,
                            // but for ctModel this is easiest method
                            int[] index = new int[]{celltypeIndex, numberOfCelltypes - 1 + celltypeIndex};
                            ctModel.addCelltypeVariablesIndex(index);
                         }
                        //System.out.println("interaction count: "+interactionTermCount);
                        ctModel.addObservedValue(celltype_perc * this.expressionValuesGeneX[sampleIndex], sampleIndex,interactionTermCount);
                        interactionTermCount++;
                        //System.out.println("debug");
                        //System.out.println(modelIndex+"\t"+celltypeIndex+"\t"+celltype_perc+"\t"+this.expressionValuesGeneX[sampleIndex]+"\t"+celltype_perc * this.expressionValuesGeneX[sampleIndex]);
                    }
                    // if i==m there is not celltype:GT interaction term so only one index added to CelltypeVariables
                    else if (sampleIndex == 0) {
                        int[] index = new int[]{celltypeIndex};
                        ctModel.addCelltypeVariablesIndex(index);
                    }
                }
                ctModel.setModelLength();

                interactionTermCount = cellCount.getNumberOfCelltypes();
            }
        }
    }


    public void cleanUp() throws IllegalAccessException {
        this.expressionValuesGeneX = null;
        this.expressionValuesGeneY = null;
        for (InteractionModel interactionModel : this.interactionModels.values()) {
            interactionModel.cleanUp();
        }
        this.cellCount = null;
    }

    public List<String> getSampleNames() {
        return sampleNames;
    }
}
