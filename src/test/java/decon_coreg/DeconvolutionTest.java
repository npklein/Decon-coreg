package decon_coreg;

import static org.junit.Assert.*;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.hamcrest.CoreMatchers;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import main.java.decon_coreg.CommandLineOptions;
import main.java.decon_coreg.Deconvolution;
import main.java.decon_coreg.DeconvolutionResult;
import main.java.decon_coreg.Main;

public class DeconvolutionTest {
	String outputDir = "src/test/resources/deconvolutionTestResults/";
	String counts;
	String genotypes;
	String geneSnpList;
	String expression;
	CommandLineOptions commandLineOptions;
	
	@Before
	public void init() {
		commandLineOptions = new CommandLineOptions();
		File countsFile = new File("src/test/resources/cellcount_files/cellcounts.txt");
		File genotypesFile = new File("src/test/resources/genotype_files/genotype_dosages.txt");
		File geneSnpListFile = new File("src/test/resources/gene_snp_list_files/gene_snp_list.txt");
		File expressionFile = new File("src/test/resources/expression_files/expression_levels.txt");
		counts = countsFile.getAbsolutePath();
		genotypes = genotypesFile.getAbsolutePath();
		geneSnpList = geneSnpListFile.getAbsolutePath();
		expression = expressionFile.getAbsolutePath();
	}
	
	@After
	public void tearDown() throws Exception {
		//deleteDir(new File(outputDir));
	}	
	
	@Test
	public void mainTest() throws Exception {
		/*
		String[] args = {"-o",outputDir+"deconvolutionTestResultsTestRun","-c",counts,
						 "-e",expression, "-sn", geneSnpList};

		Main.main(args);*/
	}
	
	/*
	 * Give error when expression and genotype file have different names
	 */
	@Test
	public void readInputDataTest() throws Exception {
		/*File cellCountsSmall = new File("src/test/resources/cellcount_files/cellcounts_small.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts,cellCountsSmall.getAbsolutePath(),
						 "-e",expression, 
						 "-sn", geneSnpList};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);

		deconvolution.readInputData();*/
	}

	
	@Test
	public void runDeconPerGeneSnpPairTestRunTest() throws Exception {
		/*
		File geneSnpListFile = new File("src/test/resources/gene_snp_list_files/gene_snp_list_long.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResultsTestRun","-c",counts,
						 "-e",expression, "-g",genotypes, "-sn", geneSnpListFile.getAbsolutePath(),
						 "-t"};

		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		List<DeconvolutionResult> deconvolutionResults = deconvolution.runDeconPerGeneGenePair();
		deconvolution.writeDeconvolutionResults(deconvolutionResults);
		Path path = Paths.get(outputDir+"deconvolutionTestResultsTestRun/deconvolutionResults.csv");
		long lineCount = Files.lines(path).count();
		assertEquals("100 example lines written", lineCount, 101);*/
	}
		
	@Test
	public void runDeconPerGenePairNotInExpressionFileTest() throws Exception {
		/*
		File geneSnpList = new File("src/test/resources/gene_snp_list_files/gene_snp_list_non_existing_gene.txt");
		String[] args = {"-o",outputDir+"deconvolutionTestResults","-c",counts,
						 "-e",expression,
						 "-sn", geneSnpList.getAbsolutePath()};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		try {
			deconvolution.runDeconPerGeneGenePair();
			fail( "My method didn't throw when I expected it to" );
		} catch (RuntimeException expectedException) {
			assertThat(expectedException.getMessage(), CoreMatchers.containsString("included in gene/gene combinations to test, but not available in the expression file"));
		}*/
	}

	@Test
	public void writeDeconvolutionResultWithSpearmanCorrelationTest() throws Exception {
		/*
		String[] args = {"-o",outputDir+"deconvolutionSpearmanResults","-c",counts,
						 "-e",expression,
						 "-sn", geneSnpList, "-w"};
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		List<DeconvolutionResult> deconvolutionResults = deconvolution.runDeconPerGeneGenePair();
		deconvolution.writeDeconvolutionResults(deconvolutionResults);

		LineIterator deconResults = FileUtils.lineIterator(new File(outputDir+"deconvolutionSpearmanResults/deconvolutionResults.csv"), "UTF-8");
		LineIterator deconExpected = FileUtils.lineIterator(new File("src/test/resources/expected_results/deconSpearmanExpected.txt"), "UTF-8");
		//test if header is same
		assertEquals("File header the same",deconExpected.next(),deconResults.next());
		while (deconResults.hasNext() && deconExpected.hasNext()){
			ArrayList<String> deconResultsStringVector = new ArrayList<String>(Arrays.asList(deconResults.next().split("\t")));
			ArrayList<String> deconExpectedStringVector = new ArrayList<String>(Arrays.asList(deconExpected.next().split("\t")));
			assertEquals("Deconresult same as expected", deconExpectedStringVector, deconResultsStringVector);
			assertEquals("QTL name the same", deconExpectedStringVector.remove(0), deconResultsStringVector.remove(0));
		}*/
	}
	
	
	void deleteDir(File file) {
	    File[] contents = file.listFiles();
	    if (contents != null) {
	        for (File f : contents) {
	            deleteDir(f);
	        }
	    }
	    file.delete();
	}
}