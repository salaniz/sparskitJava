package edu.tuberlin.sparskit;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;
import java.util.Random;
import java.util.Scanner;

import edu.tuberlin.sparskit.formats.AdaptedCompRowMatrix;
import edu.tuberlin.sparskit.formats.JaggedDiagonalMatrix;
import edu.tuberlin.sparskit.formats.ModifiedCompRowMatrix;
import no.uib.cipr.matrix.DenseVector;


public class DatasetTestsPlain {
	
	private static AdaptedCompRowMatrix readToCsr(Path file, int offset) throws Exception {
		
//		DataInputStream in = new DataInputStream(new FileInputStream(file.toString()));
		
		Scanner in = new Scanner(new File(file.toString()));
		in.useLocale(Locale.US);
		
		String line = "";
		
		while(in.hasNextLine() && (line = in.nextLine()).startsWith("%")){
			//skip comment lines
		}
				
		String[] lineItems = line.split(" ");
		
        int numColumns = Integer.parseInt(lineItems[0]);
        int numRows = Integer.parseInt(lineItems[1]);
        int nnz = Integer.parseInt(lineItems[2]);
                
//        List<Double> data = new ArrayList<Double>();
//        List<Integer> columnIndex = new ArrayList<Integer>();
//        List<Integer> rowPointer = new ArrayList<Integer>();
        double[] data = new double[nnz];
        int[] columnIndex = new int[nnz];
        int[] rowPointer = new int[numRows+1];
        
        int prevRow = 0, currRow, i = 0;
//        while(in.hasNext()){
//        	columnIndex.add(in.nextInt());
//        	currRow = in.nextInt();
//        	if(currRow > prevRow) {
//        		for(int j = 0; j < currRow - prevRow; j++) {
//        			rowPointer.add(i);
//        		}
//        	}
//        	if(in.hasNextDouble()) { // does not work since it also detects int
//        		data.add(in.nextDouble());
//        	} else {
//        		data.add(1D);
//        	}
//        	in.nextLine();
//        	prevRow = currRow;
//        	i++;
//		}
        while(in.hasNextLine()){
        	lineItems = in.nextLine().split(" ");
        	
        	columnIndex[i] = Integer.parseInt(lineItems[0]) - (1-offset);
        	currRow = Integer.parseInt(lineItems[1]);
        	if(currRow > prevRow) {
        		for(int j = prevRow; j < currRow; j++) {
        			rowPointer[j] = i + offset;
        		}
        	} else if (currRow < prevRow) {
        		in.close();
        		throw new Exception("Matrix Market Format not sorted, aborting");
        	}
        	if(lineItems.length > 2) {
        		data[i] = Double.parseDouble(lineItems[2]);
        	} else {
        		data[i] = 1D;
        	}
        	prevRow = currRow;
        	i++;
		}
        rowPointer[prevRow] = i + offset;
        
		in.close();
		return new AdaptedCompRowMatrix(numRows, numColumns, data, columnIndex, rowPointer);
	}
	
	private static double[] createRandomDoubleArray(int size) {
		Random rand = new Random();
		double[] array = new double[size];
		for(int i = 0; i < size; i++) {
			array[i] = rand.nextDouble();
		}
		return array;
	}
	
	public static void main(String[] args) throws IOException {
		
//		double[][] dense = new double[][]{{2,0,0},{0,0,0},{0,0,2}};
//		
//		double[] data = new double[2];
//		int[] coli = new int[2];
//		int[] rowp = new int[4];
//		
//		Sparskit.dnscsr(3, 3, 2, dense, 3, data, coli, rowp);
//		System.out.println(Arrays.toString(data));
//		System.out.println(Arrays.toString(coli));
//		System.out.println(Arrays.toString(rowp));
		
		
		Files.walk(Paths.get("data")).forEach(filePath -> {
		    if (Files.isRegularFile(filePath)) {
		    	AdaptedCompRowMatrix csrMatrix = null;
		    	try {
		    		csrMatrix = readToCsr(filePath, 1);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		    	
		    	
		    	
		    	if(csrMatrix != null) {
		    	
			    	DenseVector xVector = new DenseVector(createRandomDoubleArray(csrMatrix.numColumns()));
			    	DenseVector yVector = new DenseVector(csrMatrix.numRows());
			    	
			    	ModifiedCompRowMatrix modCsrMatrix = csrMatrix.toModifiedCompRowMatrix();
			    	JaggedDiagonalMatrix jaggedMatrix = csrMatrix.toJaggedDiagonalMatrix();
//			    	EllpackMatrix ellpackMatrix = csrMatrix.toEllpackMatrix();
			    	System.out.println("Transformation Done");
			    	csrMatrix.multSparskit(xVector, yVector);
			    	modCsrMatrix.mult(xVector, yVector);
			    	jaggedMatrix.mult(xVector, yVector);
//			    	ellpackMatrix.mult(xVector, yVector);
			    	System.out.println("Multiplication Done");

		    	}
		    	
		    }
		});
	}
}
