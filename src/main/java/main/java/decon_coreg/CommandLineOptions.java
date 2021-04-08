package main.java.decon_coreg;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Properties;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


public class CommandLineOptions {
	private String expressionFile;
	private String cellCountFile;
	private String genesToTestFile;
	private String outfile = "deconvolutionResults.csv";
	private String outfolder;
	private Boolean testRun = false;
	private Boolean wholeBloodCorrelation = false;
	private Boolean noConsole = false;
	private String programVersion;
	private Boolean useOLS = false;

	/**
	 * Standard command line parsing.
	 * 
	 * @param args A string vector of all arguments given to the command line, e.g. for `java -jar Deconvolution.jar -o` args = ["-o"]
	 * 
	 * @throws ParseException	If command line options can not be parsed
	 * @throws IOException 
	 */
	public void parseCommandLine(String[] args) throws ParseException, IOException {
		// load the properties file so that version can be printed
		Properties properties = new Properties();
		properties.load(this.getClass(). getClassLoader().getResourceAsStream("project.properties"));
		programVersion = properties.getProperty("artifactId")+"_"+properties.getProperty("version");
		Options options = new Options();
		Option help = new Option("help", "print this message");
		Option cellcount = Option.builder("c").required(true).hasArg().longOpt("cellcount").desc("Cellcount file name")
				.argName("file").build();
		Option expression = Option.builder("e").required(true).hasArg().longOpt("expression")
				.desc("Expression file name").argName("file").build();
		Option noConsoleOption = Option.builder("no").required(false).longOpt("no_console")
				.desc("Do not output logging info to the console").build();
		Option outfolder = Option.builder("o").required(true).hasArg().longOpt("outfolder").desc("Path to folder to write output to")
				.argName("path").build();
		Option outfile = Option.builder("of").required(false).hasArg().longOpt("outfile").desc("Outfile name of deconvolution results (will be written in outfolder)")
				.argName("file").build();
		Option genesToTestOption = Option.builder("sn").required(true).hasArg().longOpt("genesToTest").argName("file")
				.desc("Tab delimited file with first column gene name, second column gene name.").build();
		Option doTestRun = Option.builder("t").required(false).longOpt("test_run")
				.desc("Only run deconvolution for 100 QTLs for quick test run").build();
		Option version = Option.builder("v").required(false).longOpt("version")
				.desc("Print the version of the program").build();
		Option useOlsOption = Option.builder("uo").required(false).longOpt("use_OLS")
				.desc("use OLS option").build();
		Option wholeBloodCorrelationOption = Option.builder("wc").required(false).longOpt("whole_blood_correlation")
				.desc("Calculate the whole blood correlation").build();
		
		options.addOption(help);
		options.addOption(outfile);
		options.addOption(expression);
		options.addOption(cellcount);
		options.addOption(outfolder);
		options.addOption(doTestRun);
		options.addOption(genesToTestOption);
		options.addOption(noConsoleOption);
		options.addOption(version);
		options.addOption(useOlsOption);
		options.addOption(wholeBloodCorrelationOption);

		CommandLineParser cmdLineParser = new DefaultParser();
		try{
			CommandLine cmdLine = cmdLineParser.parse(options, args);
			if (cmdLine.hasOption("help")) {
				// automatically generate the help statement
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp("deconvolution", options, true);
			}
			parseOptions (cmdLine);
			printArgumentValues(cmdLine);
		}
			catch(MissingOptionException e){
				HelpFormatter formatter = new HelpFormatter();
				DeconvolutionLogger.log.info(e.toString());
				formatter.printHelp("deconvolution", options, true);
				System.exit(0);
		}
	}
	
	private void parseOptions(CommandLine cmdLine) throws IOException{

		expressionFile = cmdLine.getOptionValue("expression");
		cellCountFile = cmdLine.getOptionValue("cellcount");
		genesToTestFile = cmdLine.getOptionValue("genesToTest");
		// check if all input files exist before starting the program to return error as early as possible
		if(!new File(expressionFile).exists() || new File(expressionFile).isDirectory()) { 
		    throw new FileNotFoundException(expressionFile+" does not exist");
		}
		if(!new File(cellCountFile).exists() || new File(cellCountFile).isDirectory()) { 
		    throw new FileNotFoundException(cellCountFile+" does not exist");
		}
		if(!new File(genesToTestFile).exists() || new File(genesToTestFile).isDirectory()) { 
		    throw new FileNotFoundException(genesToTestFile+" does not exist");
		}
				
		if (cmdLine.hasOption("outfile")) {
			outfile = cmdLine.getOptionValue("outfile");
		}
		
		outfolder = cmdLine.getOptionValue("outfolder")+"/";
		
		if (cmdLine.hasOption("no_console")) {
			noConsole = !noConsole;
		}


		if (cmdLine.hasOption("test_run")) {
			testRun = !testRun;
		}
		
		if (cmdLine.hasOption("whole_blood_correlation")){
			wholeBloodCorrelation = !wholeBloodCorrelation;
		}
		
		if (cmdLine.hasOption("version")) {
			DeconvolutionLogger.log.info("Version: "+programVersion);
		}
		
		
		if (cmdLine.hasOption("use_OLS")){
			useOLS = !useOLS;
		}
	}
	

	private void printArgumentValues(CommandLine cmdLine){
	    try {
	    	File outfolderDir = new File(outfolder);
	    	// if the directory does not exist, create it
	    	Boolean dirDidNotExist = false;
	    	if (!outfolderDir.exists()) {
	    		outfolderDir.mkdir();
	    		dirDidNotExist = true;
	    	}
	    	DeconvolutionLogger.setup(outfolder, noConsole);
	    	if(dirDidNotExist){
	    		DeconvolutionLogger.log.info("Created directory "+outfolder);
	    	}
	    	DeconvolutionLogger.log.info("Writing output and logfile to "+outfolder);
	    } catch (IOException e) {
	      e.printStackTrace();
	      throw new RuntimeException("Problems with creating the log files");
	    }
	    DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	    Date date = new Date();
	    DeconvolutionLogger.log.info("Starting deconvolution");
	    DeconvolutionLogger.log.info(dateFormat.format(date));
	    DeconvolutionLogger.log.info(String.format("Running deconvolution version: %s", programVersion));
	    DeconvolutionLogger.log.info("======= DECONVOLUTION paramater settings =======");
		DeconvolutionLogger.log.info(String.format("Expression file (-e): %s", expressionFile));
		DeconvolutionLogger.log.info(String.format("Cellcount file (-c): %s", cellCountFile));
		DeconvolutionLogger.log.info(String.format("Genes to test file (-sn): %s", genesToTestFile));
		DeconvolutionLogger.log.info(String.format("Outfolder (-o): %s", outfolder));
		DeconvolutionLogger.log.info(String.format("Outfile (-of): %s", outfile));
		DeconvolutionLogger.log.info(String.format("test run doing only 100 QTL (-t): %s", testRun));
		DeconvolutionLogger.log.info(String.format("Do not ouput logging info to console (-no): %s", noConsole));
		DeconvolutionLogger.log.info(String.format("Use OLS(-uo): %s", useOLS));
		DeconvolutionLogger.log.info(String.format("Calculate whole blood correlation (-wc): %s", wholeBloodCorrelation));
		DeconvolutionLogger.log.info("=================================================");
	}
	public String getExpressionFile(){
		return (expressionFile);
	}
	public String getGenesToTestFile(){
		return (genesToTestFile);
	}	

	public String getCellcountFile(){
		return cellCountFile;
	}
	public String getOutfile(){
		return outfile;
	}

	public String getOutfolder() throws IllegalAccessException{
		if(this.outfolder == null){
			throw new IllegalAccessException("Outfolder has not been set");
		}
		return outfolder;
	}

	public void setOutfolder(String newOutfolder){
		outfolder = newOutfolder;
	}

	public Boolean getTestRun(){
		return testRun;
	}

	public Boolean getWholeBloodCorrelation(){
		return wholeBloodCorrelation;
	}
	
	public Boolean getUseOLS(){
		return(useOLS);
	}
}





