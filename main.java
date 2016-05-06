import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Scanner;
import java.util.Vector;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.imageio.ImageIO;
import javax.media.opengl.GLOffscreenAutoDrawable;
import javax.media.opengl.GLProfile;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JTextField;

import org.apache.log4j.BasicConfigurator;
import org.eclipse.swt.widgets.Display;

import com.algosome.eutils.blast.BlastParser;
import com.algosome.eutils.blast.GetCommand;
import com.algosome.eutils.blast.PutCommand;
import com.metsci.glimpse.canvas.FBOGlimpseCanvas;
import com.metsci.glimpse.examples.Example;
import com.metsci.glimpse.gl.util.GLUtils;
import com.metsci.glimpse.plot.ColorAxisPlot2D;
import com.metsci.glimpse.support.settings.SwingLookAndFeel;

import com.metsci.glimpse.examples.Example;

public class main extends Frame {
	//Variable initialization
	public static int nmer = 3;
	public static int numShifts = 0;
	public static int numShiftsMinus = 0;
	public static int windowSize = 5000;
	public static String file1 = "sequence1.fa";
	public static String file2 = "sequence2.fa";
	public static float[][] regressionValues = null;
	public static double numberIterationsX = 0;
	public static double numberIterationsY = 0;
	public static String seqXname = "sequence1";
	public static String seqYname = "sequence2";
	public static String y = "";
	public static String x = "";
	public static double [] maxValues_X;
	public static double [] maxValues_Y;
	public static double [] averageValues_X;
	public static double [] averageValues_Y;
	public static Vector<Integer> interestingX = new Vector<Integer>();
	public static Vector<Integer> interestingY = new Vector<Integer>();
	
	//Code starts here
	public static void main(String[] args) throws Exception {		
		//GUI runs, inputting from user. After input, goes to RunOurBaby() 
		FrontGUI window = new FrontGUI();
		window.runGUI();
		//RunOurBaby();
		
	}
	
	//Runner method to get kmers from window: ignores kmers with n's
	public static float[] runGetKmers (String sequence) {
		float[] kmerComp = new float[(int) Math.pow(4, nmer)];
		String[] toRun = sequence.split("N");
		for (int i =0; i<toRun.length; i++) {
			if (toRun[i].length()>=nmer) {
				kmerComp = getKmers(toRun[i], kmerComp);
				}
			}
		return kmerComp;
	}

	/*
	 * Automated BLAST method
	 * Looks through the array of correlation values generated for each window in X and Y
	 * If the averageValue of that window and maxValue of that window is below a certain threshold, sends a 
	 * BLAST query to NCBI via runBLAST method
	 */
	public static void automatedBLAST () {
		int intIterationsX = (int)numberIterationsX;
		int intIterationsY = (int)numberIterationsY;
		System.out.println();
		System.out.println("averageValues_X");
		for(int xv = 0;xv<intIterationsX;xv++){
			averageValues_X[xv] = averageValues_X[xv]/intIterationsX;
			//if thresholds met, add the indice to interesting X for later blasting
			if (averageValues_X[xv]<.5 && maxValues_X[xv]<.8) {
				interestingX.add(xv);
				System.out.print(averageValues_X[xv]);
			}
		}
		System.out.println();
		System.out.println("averageValues_Y");
		for(int xv = 0;xv<intIterationsY;xv++){
			averageValues_Y[xv] = averageValues_Y[xv]/intIterationsY;
			//if thresholds met, add the indice to interesting Y for later blasting
			if (averageValues_Y[xv]<.5 && maxValues_Y[xv]<.8) {
				interestingY.add(xv);
				System.out.print(averageValues_Y[xv]);
			}
		}
		System.out.println();
		
		//Blast interestingYs
		for (int i = 0; i<interestingY.size(); i++) {
			int temp = interestingY.get(i);
			boolean stillContig = true;
			while (stillContig){
				if(interestingY.contains(temp+1)){
					interestingY.removeElement(temp+1);
					temp++;
				}
				else{
					stillContig = false;
				}
			}
			String seq = y.substring(interestingY.get(i) * windowSize, (temp+1) * windowSize);
			seq = seq.toUpperCase();
			String filename = seqYname + "_" + interestingY.get(i)+".txt";
			try {
				runBlast(seq, filename);
			} catch (FileNotFoundException e1) {
				e1.printStackTrace();
			}
		}
		
		//Blast interestingXs
		for (int i = 0; i<interestingX.size(); i++) {
			int temp = interestingX.get(i);
			boolean stillContig = true;
			while (stillContig){
				if(interestingX.contains(temp+1)){
					interestingX.removeElement(temp+1);
					temp++;
				}
				else{
					stillContig = false;
				}
			}
			String seq = x.substring(interestingX.get(i) * windowSize, (temp + 1) * windowSize);
			seq = seq.toUpperCase();
			String filename = seqXname + "_" + interestingX.get(i)+".txt";
			try {
				runBlast(seq, filename);
			} catch (FileNotFoundException e1) {
				e1.printStackTrace();
			}
		}
	}
	
	/*
	 * Similar to autoBLAST,
	 * but instead, this method writes all interesting sequences to file.
	 */
	public static void automatedFile () throws FileNotFoundException {
		int intIterationsX = (int)numberIterationsX;
		int intIterationsY = (int)numberIterationsY;
		System.out.println();
		System.out.println("averageValues_X");
		for(int xv = 0;xv<intIterationsX;xv++){
			averageValues_X[xv] = averageValues_X[xv]/intIterationsX;
			if (averageValues_X[xv]<.5 && maxValues_X[xv]<.8) {
				interestingX.add(xv);
				System.out.print(averageValues_X[xv]);
			}
		}
		System.out.println();
		System.out.println("averageValues_Y");
		for(int xv = 0;xv<intIterationsY;xv++){
			averageValues_Y[xv] = averageValues_Y[xv]/intIterationsY;
			if (averageValues_Y[xv]<.5 && maxValues_Y[xv]<.8) {
				interestingY.add(xv);
				System.out.print(averageValues_Y[xv]);
			}
		}
		System.out.println();
		
		String filename = seqYname+".txt";
		PrintStream outt = new PrintStream(new File(filename));
		for (int i = 0; i<interestingY.size(); i++) {
			int temp = interestingY.get(i);
			boolean stillContig = true;
			while (stillContig){
				if(interestingY.contains(temp+1)){
					interestingY.removeElement(temp+1);
					temp++;
				}
				else{
					stillContig = false;
				}
			}
			String seq = y.substring(interestingY.get(i) * windowSize, (temp+1) * windowSize);
			seq = seq.toUpperCase();
			outt.println(">Window:"+interestingY.get(i)+","+(temp+1));
			int numIter = seq.length()/50;
			int j = 0;
			for (j = 0; j<numIter; j++) {
				outt.println(seq.substring(j*50, (j+1)*50));
			}
			if (j*50<seq.length()) {
				outt.println(seq.substring((j)*50, seq.length()));
			}
		}
		outt.close();
		
		
		filename = seqXname+".txt";
		PrintStream out = new PrintStream(new File(filename));
		for (int i = 0; i<interestingX.size(); i++) {
			int temp = interestingX.get(i);
			boolean stillContig = true;
			while (stillContig){
				if(interestingX.contains(temp+1)){
					interestingX.removeElement(temp+1);
					temp++;
				}
				else{
					stillContig = false;
				}
			}
			String seq = x.substring(interestingX.get(i) * windowSize, (temp + 1) * windowSize);
			seq = seq.toUpperCase();
			out.println(">Window:"+interestingX.get(i)+","+(temp+1));
			int numIter = seq.length()/50;
			int j = 0;
			for (j = 0; j<numIter; j++) {
				out.println(seq.substring(j*50, (j+1)*50));
			}
			if (j*50<seq.length()) {
				out.println(seq.substring((j)*50, seq.length()));
			}
		}
		out.close();
		}
	
	/*
	 * Uses JeUtil package to send query to NCBI's Servers and gets 
	 * files back with BLAST results
	 */
	public static void runBlast(String seq, String filename) throws FileNotFoundException {
    	BasicConfigurator.configure();
    	//logger.info("Blast utility test");
    	PutCommand put = new PutCommand();
    	PrintStream out = new PrintStream(new File(filename));
    	put.setQuery(seq);
    	put.setProgram("blastn");
    	put.setDatabase("nr");
    	
    	GetCommand get = new GetCommand(new BlastParser(){

			@Override
			public void parseBlastOutput(String output) {
				out.println(output);
			}
    		
    	});
    	get.setFormatType("Text");
    	//logger.info("Blasting");
    	Blast blast = new Blast(put, get);
    	blast.run();
	}
	
	/*
	 * Replaces those funcky nucleotides with Ns
	 */
	public static String replaceNucs(String sequence) {
		sequence = sequence.replaceAll("M", "N");
		sequence = sequence.replaceAll("K", "N");
		sequence = sequence.replaceAll("R", "N");
		sequence = sequence.replaceAll("S", "N");
		sequence = sequence.replaceAll("W", "N");
		sequence = sequence.replaceAll("Y", "N");
		sequence = sequence.replaceAll("B", "N");
		sequence = sequence.replaceAll("D", "N");
		sequence = sequence.replaceAll("H", "N");
		sequence = sequence.replaceAll("V", "N");
		sequence = sequence.toUpperCase();
		return sequence;
	}
	
	/*
	 * This is where stuff happens
	 */
	public static void RunOurBaby() throws Exception{
		//This chunk imports files and replaces weird nucleotides
		System.out.println("Importing file");
		File file;
		file = new File(file1);
		InputStream fis = null;
		try {
			fis = new FileInputStream(file);
		} catch (FileNotFoundException e1) {
			JOptionPane.showMessageDialog(null, "Sequence1 not found. Please rerun SPlot");
			System.exit(0);
		}
		BufferedReader reader = new BufferedReader(new InputStreamReader(fis));
		StringBuilder out = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            out.append(line.trim());
        }
        x = out.toString();   //Prints the string content read from input stream
        reader.close();
		x = replaceNucs(x);
		
		File fileTwo;
		fileTwo = new File(file2);
		InputStream is = null;
		try {
			is = new FileInputStream(fileTwo);
		} catch (FileNotFoundException e1) {
			JOptionPane.showMessageDialog(null, "Sequence2 not found. Please rerun SPlot");
			System.exit(0);
		}
		BufferedReader reader2 = new BufferedReader(new InputStreamReader(is));
		StringBuilder out2 = new StringBuilder();
		String line2;
		while ((line = reader2.readLine()) != null) {
			out2.append(line);
		}
		y = out2.toString();
		reader2.close();
		y = replaceNucs(y);
		
		//Some prints to know how long are the things we are dealing with
		System.out.println(y.length());
		System.out.println(y.length()/windowSize);
		
		//calculate number of iterations for loops
		numberIterationsX = (x.length()/windowSize);
		int intIterationsX = (int)numberIterationsX;
		numberIterationsY = (y.length()/windowSize);
		int intIterationsY = (int)numberIterationsY;
		//calculate bitshifting backwards/forwards numbers
		numShifts = 64 - (2*(nmer-1));
		numShiftsMinus = numShifts - 2;
		
		//Vectors and such to be used for autoID of RUCPS
		maxValues_X = new double[intIterationsX];
		maxValues_Y = new double[intIterationsY];
		averageValues_X = new double[intIterationsX];
		averageValues_Y = new double[intIterationsY];
		Vector<float[]> xStorage = new Vector<float[]>();
		Vector<float[]> yStorage = new Vector<float[]>();

		/*
		 * We start filtering things here
		 */
		int badX = 0;
		int badY = 0;
		int sum = 0;
		
		//adds windows that do not have many Rs (less than 20% Rs, to be processed later)
		for (int iii = 0; iii<intIterationsX; iii++) {
			sum = 0;
			float[] toAdd = runGetKmers(x.substring(windowSize*iii, windowSize*(iii+1)));
			for (int i = 0; i<toAdd.length; i++) {
				sum += (int)toAdd[i];
			}
			if (sum < windowSize * .8) {
				badX++;
			}
			else {
				xStorage.add(toAdd);
			}
		}
		
		//adds windows that do not have many Rs (less than 20% Rs, to be processed later)
		for (int iii = 0; iii<intIterationsY; iii++) {
			sum = 0;
			float[] toAdd = runGetKmers(y.substring(windowSize*iii, windowSize*(iii+1)));
			for (int i = 0; i<toAdd.length; i++) {
				sum += (int)toAdd[i];
			}
			if (sum < windowSize * .8) {
				badY++;
			}
			else {
				yStorage.add(toAdd);
			}
		}
		System.out.println(badX);
		System.out.println(badY);
		//adjust number of iterations after throwing out bad windows
		intIterationsX = intIterationsX - badX;
		intIterationsY = intIterationsY - badY;
		numberIterationsX = (double)intIterationsX;
		numberIterationsY = (double)intIterationsY;
		regressionValues = new float[intIterationsX][intIterationsY];
		
		//Begin the bit shifting magic
		System.out.println("Begin Bitting");
		for (int ii = 0; ii<intIterationsX; ii++) {
			for (int jj = 0; jj<intIterationsY; jj++) {
				//for each y and x, get a regression value and store it in our 2D array
				float rvalue = getR(xStorage.get(ii), yStorage.get(jj));
				regressionValues[ii][jj] = rvalue;
				//tracks lows/highs/averages for automatic analysis
				if(rvalue > maxValues_X[ii]){
					maxValues_X[ii] = rvalue;
				}
				if(rvalue > maxValues_Y[jj]){
					maxValues_Y[jj] = rvalue;
				}
				averageValues_X[ii] += rvalue;
				averageValues_Y[jj] += rvalue;
			}
			System.out.println("I=" + ii);
		}		
		System.out.println("Regression Completed");
		

		//initiate graphing
		System.out.println("Graphing");
    	HeatMapExample a = new HeatMapExample();
        Example b = null;
        Example.showWithSwing(a);
        
        //run GUI in seperate thread
        Thread z = new Thread (new Runnable() {
			@Override
			public void run() {
				try {
					System.out.println("0");
			        BackGUI processing = new BackGUI();
					processing.runGUI();
					System.out.println("1");
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		});
		z.start();        
		//Take picture of our beautiful SPlot
        takePicture(a);
	}
	
	//Glimpse code to take a picture of our pretty SPlot
	public static void takePicture (HeatMapExample a) throws IOException {
       // a.cursorPainter.setVisible(false);
        GLProfile glProfile = GLUtils.getDefaultGLProfile( );
        
        // generate a GLContext by constructing a small offscreen framebuffer
        final GLOffscreenAutoDrawable glDrawable = GLUtils.newOffscreenDrawable( glProfile );

        // create an offscreen GlimpseCanvas which shares an OpenGL context with the above drawable
        // (its size is 1000 by 1000 pixels)
        final FBOGlimpseCanvas canvas = new FBOGlimpseCanvas( glDrawable.getContext( ), 1000, 1000);

        // set the Glimpse look and feed of the canvas just like we would for an onscreen canvas
        canvas.setLookAndFeel( new SwingLookAndFeel( ) );
        
        // use one of the previous examples to build a simple plot to draw
        ColorAxisPlot2D layout = a.getLayout( );
        // add the layout to the offscreen canvas
        layout.getCrosshairPainter().setVisible(false);
        
        canvas.addLayout( layout );
        // draw the canvas to a BufferedImage and write the image to a file
        BufferedImage image = canvas.toBufferedImage( );
        ImageIO.write( image, "PNG", new File( "Chr21HumanVsChimp.png" ) );
        //a.cursorPainter.setVisible(true);

	}
	
	//Returns the Pearson Correlation coefficient
	public static float getR (float[] x, float[] y) {
		float xAverage = 0;
		float yAverage = 0;
		for (int i =0; i<x.length; i++) {
			xAverage += x[i];
			yAverage += y[i];
		}
		xAverage = xAverage/x.length;
		yAverage = yAverage/y.length;
		float xy = 0;
		float xSq = 0;
		float ySq = 0;
		for (int i = 0; i<x.length; i++) {
			x[i] = (x[i] - xAverage);
			y[i] = (y[i] - yAverage);
			xy += (x[i]*y[i]);
			xSq += x[i]*x[i];
			ySq += y[i]*y[i];
		}
		return (float) (xy/Math.sqrt(xSq*ySq));
	}

//This method is the runner for each windows' bitshifts	
	public static float[] getKmers(String sequence, float[] kmerComp) {
		//initialize first set of kmers
		Long temp = null;
		Long full = Long.parseUnsignedLong("0");
		int i = 0;
		for (i=0; i<nmer; i++) {
			temp  = nucToNum(sequence.charAt(i));
			full = full + temp;
			if (i<nmer-1) {
				full = full<<2;
			}
		}
		//add it and its reverse kmer to count array
		kmerComp[full.intValue()] += 1;
		kmerComp[reverser(full)] += 1;

		//delete first nucleotide and add to the end of it
		//add it and its reverse complement to count array
		while (i<sequence.length()) {
			temp = nucToNum(sequence.charAt(i));
			full = fancyShift(full);
			full = full + temp;
			kmerComp[full.intValue()] += 1;
			kmerComp[reverser(full)] += 1;
			i++;
		}
		//return kmerComposition array
		return kmerComp;
	}
	
	//shifts kmer off the cliff and back
	public static Long fancyShift (Long a) {
		a = a << numShifts;
		a = a >>> numShiftsMinus;
		return a;
	}
	
	//reverses kmer and makes the reverse complement
	public static int reverser(Long xx) {
		xx = ~xx;
		xx = Long.reverse(xx);
		xx = xx >>> numShiftsMinus;
		return xx.intValue();
	}
	
	//Silly initializations for optimization
	public static Long aLong = Long.parseUnsignedLong("0");
	public static Long cLong = Long.parseUnsignedLong("1");
	public static Long gLong = Long.parseUnsignedLong("2");
	public static Long tLong = Long.parseUnsignedLong("3");

	//turns a nucleotide character to a binary representation
	public static Long nucToNum (char a) {
		switch (a) {
			case 'A':
				return aLong;
			case 'C':
				return cLong;
			case 'G':
				return gLong;
			case 'T':
				return tLong;
			default:
				return null;
		}
	}

}
