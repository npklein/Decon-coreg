package main.java.decon_coreg;

import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.ParseException;

public class Main {
	/**
	 * Deconvolute gene-gene co-regulation given the expression levels and
	 * and cell counts. Calculates the p-values for the deconvoluted coregulation
	 * and writes them to an out file
	 * 
	 * @param args List of command line arguments
	 * 
	 * @throws ParseException	If cell count file is not in right format to be parsed correctly
	 * @throws IllegalAccessException	If out folder can not be retrieved from commandLineOptions
	 * @throws IOException	If cell counts file can not be found or read
	 */
	public static void main(String[] args) throws ParseException, IllegalAccessException, IOException {
		CommandLineOptions commandLineOptions = new CommandLineOptions();
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		List<DeconvolutionResult> deconvolutionResults = deconvolution.runDeconPerGeneGenePair();
		deconvolution.writeDeconvolutionResults(deconvolutionResults);

	}

}
