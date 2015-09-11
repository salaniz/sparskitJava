package edu.tuberlin.sparskit.formats;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;

import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExternalResource;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import edu.tuberlin.sparskit.Sparskit;
import edu.tuberlin.sparskit.formats.AdaptedCompRowMatrix;
import edu.tuberlin.sparskit.formats.EllpackMatrix;
import edu.tuberlin.sparskit.formats.JaggedDiagonalMatrix;
import edu.tuberlin.sparskit.formats.ModifiedCompRowMatrix;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@RunWith(Parameterized.class)
public class FormatsDatasetTests {
	
	private static final long seed = new Random().nextLong();

    // Keep Track of global stats
    static Map<Path, Map<Experiment, Long>> convertRuntime;
    static Map<Path, Map<Experiment, Long>> multRuntime;
    static Map<Path, Map<Experiment, Long>> formatMemory;
    static Map<Path, Map<String, Double>> mtxStatistics;
    
    public AdaptedCompRowMatrix csrMatrix;
    public AdaptedCompRowMatrix transCsrMatrix;
    public boolean symmetric = false; //all matrices are not symmetric currently, can be optimized when reading the matrix into CSR format
    public DenseVector xVector;
	public DenseVector yVector;
	public DenseVector solutionVector;
	public static final String[] statsNames = new String[]
			{"rc_dimension",
			"lower_bandwidth", "upper_bandwidth", "max_bandwidth", "average_bandwidth", 
			"max_column_length", "min_column_length", "max_row_length", "min_row_length", "zero_column_number", "zero_row_number",
			"diag_domi_column_percentage", "diag_domi_row_percentage",
			"Frobenius_norm",
			"sym_part_Frobenius_norm", "nonsym_part_Frobenius_norm", "matching_elements_number", "relative_sym_match",
			"average_dist_of_a", "standard_deviation",
			"nonzero_number_in_skyline",
			//"element_number_in_eachdiag",
			"nonzero_number_of_lower_part", "nonzero_number_of_upper_part", "nonzero_number_of_maindiag"};

    @Parameterized.Parameter
    public Path testMatrixFile;

    @Parameterized.Parameters(name = "{index}: Dataset {0}")
    public static List<Object> data() throws IOException {
    	
    	List<Object> parameters = new ArrayList<Object>();
    	Files.walk(Paths.get("data")).forEach(filePath -> {
		    if (Files.isRegularFile(filePath)) {
		    	parameters.add(filePath);
		    }
    	});

        return parameters;
    }

    @BeforeClass
    public static void setup() throws Exception {
    	convertRuntime = new HashMap<Path, Map<Experiment, Long>>();
    	multRuntime = new HashMap<Path, Map<Experiment, Long>>();
    	formatMemory = new HashMap<Path, Map<Experiment, Long>>();
    	mtxStatistics = new HashMap<Path, Map<String, Double>>();

        for (Object file : data()) {
        	convertRuntime.put((Path) file, new HashMap<Experiment, Long>());
        	multRuntime.put((Path) file, new HashMap<Experiment, Long>());
        	formatMemory.put((Path) file, new HashMap<Experiment, Long>());
        	mtxStatistics.put((Path) file, new HashMap<String, Double>());
        }
        
        System.out.println("Test Seed: " + Math.abs(seed) % 1000000);
        
    }

    @AfterClass
    public static void tearDown() throws IOException {
    	
		Pattern p = Pattern.compile("/(.*?)\\.");
    	
    	//write statistics to file
    	FileWriter multStatsWriter = new FileWriter("stats/" + "MatVecMultStats_" + Math.abs(seed) % 1000000 + ".txt");
    	FileWriter convertStatsWriter = new FileWriter("stats/" + "MatVecConvertStats_" + Math.abs(seed) % 1000000 + ".txt");
    	FileWriter formatMemoryWriter = new FileWriter("stats/" + "MatVecFormatMemoryStats_" + Math.abs(seed) % 1000000 + ".txt");
    	FileWriter mxtStatsWriter = new FileWriter("stats/" + "MtxStats_" + Math.abs(seed) % 1000000 + ".txt");
    	multStatsWriter.write("Matrix"); // CSR ModCSR Jagged Ellpack\n
    	convertStatsWriter.write("Matrix");
    	formatMemoryWriter.write("Matrix");
    	for(Experiment exp : Experiment.values()) {
        	multStatsWriter.write(" " + exp.toString());
        	convertStatsWriter.write(" " + exp.toString());
    		formatMemoryWriter.write(" " + exp.toString());
        }
    	mxtStatsWriter.write("Matrix");
    	for(String statsName : statsNames) {
    		mxtStatsWriter.write(" " + statsName);
    	}
    	multStatsWriter.write("\n");
    	convertStatsWriter.write("\n");
		formatMemoryWriter.write("\n");
    	mxtStatsWriter.write("\n");
    	for(Object f : data()) {
    		Path file = (Path)f;
    		
    		Matcher m = p.matcher(file.toString());
    		String matrix;
    		if(m.find()) {
    			matrix = m.group(1);
    		} else {
    			matrix = file.toString();
    		}
    		
    		multStatsWriter.write(matrix);
    		convertStatsWriter.write(matrix);
    		formatMemoryWriter.write(matrix);
    		mxtStatsWriter.write(matrix);
    		
    		Map<Experiment,Long> mtxMultRuntime = multRuntime.get(file);
    		Map<Experiment,Long> mtxConvertRuntime = convertRuntime.get(file);
    		Map<Experiment,Long> mtxFormatMemory = formatMemory.get(file);
    		Map<String,Double> mtxStats = mtxStatistics.get(file);
    		for(Experiment exp : Experiment.values()) {
    			if(mtxMultRuntime.containsKey(exp)) {
    				multStatsWriter.write(" " + mtxMultRuntime.get(exp));
    			}
    			if(mtxConvertRuntime.containsKey(exp)) {
    				convertStatsWriter.write(" " + mtxConvertRuntime.get(exp));
    			}
    			if(mtxFormatMemory.containsKey(exp)) {
    				formatMemoryWriter.write(" " + mtxFormatMemory.get(exp));
    			}
    		}
    		
    		for(String statsName : statsNames) {
    			if(mtxStats.containsKey(statsName)) {
    				mxtStatsWriter.write(" " + mtxStats.get(statsName));
    			}
    		}
    		
    		multStatsWriter.write("\n");
    		convertStatsWriter.write("\n");
    		formatMemoryWriter.write("\n");
    		mxtStatsWriter.write("\n");
    	}
    	multStatsWriter.flush();
    	convertStatsWriter.flush();
    	formatMemoryWriter.flush();
    	mxtStatsWriter.flush();
    	multStatsWriter.close();
    	convertStatsWriter.close();
    	formatMemoryWriter.close();
    	mxtStatsWriter.close();
    }
    
    @Rule
    public ExternalResource externalResource = new ExternalResource() {
    	//TODO: find solution that this is done only once per parameter! currently it's executed before every test
        protected void before() throws Throwable {
        	// reading first 60x60 chunk of each matrix (or less if matrix is smaller)
        	csrMatrix = readToCsr(testMatrixFile, 1, new int[]{60,60});
        	
        	transCsrMatrix = csrMatrix.getTranspose();
        	xVector = new DenseVector(createRandomDoubleArray(csrMatrix.numColumns()));
        	
        	solutionVector = new DenseVector(csrMatrix.numRows());
        	csrMatrix.mult(xVector, solutionVector);
        }
        protected void after() {}
    };

    private void testMatVecMult(Matrix matrix, Experiment experiment) throws Exception {
    	yVector = new DenseVector(csrMatrix.numRows());
    	long start = System.nanoTime();
    	matrix.mult(xVector, yVector);
    	long runtimeNanoSec = System.nanoTime() - start;
    	
    	
    	multRuntime.get(testMatrixFile).put(experiment, runtimeNanoSec);
    	
    	Assert.assertArrayEquals("Matrix vector multiplication result differs from expected output!", solutionVector.getData(), yVector.getData(), 0.1D);
    }

    @Test
    public void testMultiplyCompRowMatrix() throws Exception {
    	Experiment exp = Experiment.COMPRESSEDROW;

    	testMatVecMult(csrMatrix, exp);
    	formatMemory.get(testMatrixFile).put(exp, csrMatrix.getMemorySizeInBytes());
    }
    
    @Test
    public void testMultiplyModifiedCompRowMatrix() throws Exception {

    	Experiment exp = Experiment.MODIFIEDCOMPRESSEDROW;
    	
    	long start = System.nanoTime();
    	ModifiedCompRowMatrix modCsrMatrix = csrMatrix.toModifiedCompRowMatrix();
    	long runtimeNanoSec = System.nanoTime() - start;
    	
    	convertRuntime.get(testMatrixFile).put(exp, runtimeNanoSec);
    	
    	testMatVecMult(modCsrMatrix, exp);
    	formatMemory.get(testMatrixFile).put(exp, modCsrMatrix.getMemorySizeInBytes());
    }
    
    @Test
    public void testMultiplyJaggedCompDiagMatrix() throws Exception {
    	Experiment exp = Experiment.JAGGED_DIAGONAL;
    	
    	long start = System.nanoTime();
    	JaggedDiagonalMatrix jaggedMatrix = csrMatrix.toJaggedDiagonalMatrix();
    	long runtimeNanoSec = System.nanoTime() - start;

    	convertRuntime.get(testMatrixFile).put(exp, runtimeNanoSec);
    	
    	testMatVecMult(jaggedMatrix, exp);
    	formatMemory.get(testMatrixFile).put(exp, jaggedMatrix.getMemorySizeInBytes());
    }
    
    @Test
    public void testMultiplyEllpackMatrix() throws Exception {
    	Experiment exp = Experiment.ELLPACK;
    	
    	long start = System.nanoTime();
    	EllpackMatrix ellpackMatrix = csrMatrix.toEllpackMatrix();
    	long runtimeNanoSec = System.nanoTime() - start;

    	convertRuntime.get(testMatrixFile).put(exp, runtimeNanoSec);
    	
    	testMatVecMult(ellpackMatrix, exp);
    	formatMemory.get(testMatrixFile).put(exp, ellpackMatrix.getMemorySizeInBytes());
    }
    
    //TODO: maybe move this to "externalResource" Rule -> but first fix that it is performed once per parameter
    @Test
    public void calcMatrixStatistics() {
    	if(csrMatrix.numRows() != csrMatrix.numColumns()) {
    		return;
    	}
    	
    	//index for names of stats, see statsNames variable
    	int statsNameIdx = 0;
    	
    	//row and column dimension
    	mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)csrMatrix.numRows());
    	
    	//bandwidth
    	int[] ml = new int[1];
		int[] mu = new int[1];
		int[] iband = new int[1];
		double[] bndav = new double[1];
		Sparskit.bandwidth(csrMatrix.numColumns(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), ml, mu, iband, bndav);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)ml[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)mu[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)iband[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)bndav[0]);

		//non zeros
		int[] nzmaxc = new int[1];
		int[] nzminc = new int[1];
		int[] nzmaxr = new int[1];
		int[] nzminr = new int[1];
		int[] nzcol = new int[1];
		int[] nzrow = new int[1];
		Sparskit.nonz(csrMatrix.numColumns(), symmetric ? 1 : 0, transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), csrMatrix.columnIndex.array(), nzmaxc, nzminc, nzmaxr, nzminr, nzcol, nzrow);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nzmaxc[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nzminc[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nzmaxr[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nzminr[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nzcol[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nzrow[0]);

		//dominant diagonals
		double[] ddomc = new double[1];
		double[] ddomr = new double[1];
		Sparskit.diagdomi(csrMatrix.numColumns(), symmetric ? 1 : 0, 1, transCsrMatrix.data.array(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), csrMatrix.data.array(), csrMatrix.columnIndex.array(), csrMatrix.rowPointer.array(), ddomc, ddomr);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)ddomc[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)ddomr[0]);
		
		//Frobenius norm
		double frobnorm = Sparskit.frobnorm(csrMatrix.numColumns(), symmetric ? 1 : 0, transCsrMatrix.data.array(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array());
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)frobnorm);

		//Frobenius norm, symmetric and asymmetric
		int[] imatch = new int[1];
		double[] av = new double[1];
		double[] fas = new double[1];
		double[] fan = new double[1];
		Sparskit.ansym(csrMatrix.numColumns(), symmetric ? 1 : 0, transCsrMatrix.data.array(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), new double[csrMatrix.data.array().length], new int[csrMatrix.columnIndex.array().length], new int[csrMatrix.rowPointer.array().length], imatch, av, fas, fan);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)fas[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)fan[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)imatch[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)av[0]);

		//distance of elements
		double[] dist = new double[1];
		double[] std = new double[1];
		Sparskit.distaij(csrMatrix.numColumns(), transCsrMatrix.data.array().length, symmetric ? 1 : 0, transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), dist, std);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)dist[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)std[0]);

		//non zero elements in skyline format
		int nsky = Sparskit.skyline(csrMatrix.numColumns(), symmetric ? 1 : 0, transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), csrMatrix.columnIndex.array(), csrMatrix.rowPointer.array());
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nsky);
		
		/*
		//number of elements in each diagonal
		int[] dists = new int[csrMatrix.numRows() + csrMatrix.numColumns() - 1];
		Sparskit.distdiag(csrMatrix.numRows(), csrMatrix.numColumns(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), dists);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nsky[0]);
		*/
		
		//number of elements in main diagonal, lower and upper part
		int[] nlower = new int[1];
		int[] nupper = new int[1];
		int[] ndiag = new int[1];
		Sparskit.nonzlud(csrMatrix.numColumns(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), nlower, nupper, ndiag);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nlower[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)nupper[0]);
		mtxStatistics.get(testMatrixFile).put(statsNames[statsNameIdx++], (double)ndiag[0]);
    }
    
    //TODO: how to use MatrixVectorReader correctly. this does not seem to work. own implementation is used.
//    private AdaptedCompRowMatrix readToCsr(Path file, int offset) throws Exception {
//    	AdaptedCompRowMatrix matrix = new AdaptedCompRowMatrix(new MatrixVectorReader(new FileReader(file.toString())));
//    	matrix.incrAllIndices();
//    	return matrix;
//    }
    
    //read full matrix
    private AdaptedCompRowMatrix readToCsr(Path file, int offset) throws Exception {
    	return readToCsr(file, offset, new int[]{-1,-1});
    }
    
    //read matrix with dimension limits
    //own implementation of matrix reader: second column in data files are sorted -> this is used as the row dimension for easy CSR
    private AdaptedCompRowMatrix readToCsr(Path file, int offset, int[] dimensionLimits) throws Exception {
				
		Scanner in = new Scanner(new File(file.toString()));
		in.useLocale(Locale.US);
		
		String line = "";
		
		while(in.hasNextLine() && (line = in.nextLine()).startsWith("%")){
			//skip comment lines
		}
				
		String[] lineItems = line.split(" ");
		
        int numColumns = Integer.parseInt(lineItems[0]);
        int numRows = Integer.parseInt(lineItems[1]);
//        int nnz = Integer.parseInt(lineItems[2]);
        
        if(dimensionLimits[0] != -1 && dimensionLimits[1] != -1) {
        	if(numRows > dimensionLimits[0]) {
        		numRows = dimensionLimits[0];
        	}
        	if(numColumns > dimensionLimits[1]) {
        		numColumns = dimensionLimits[1];
        	}
        }
        
                
        List<Double> data = new ArrayList<Double>();
        List<Integer> columnIndex = new ArrayList<Integer>();
        List<Integer> rowPointer = new ArrayList<Integer>();
        
        int prevRow = 0, currRow, i = 0;
        while(in.hasNextLine()){
        	lineItems = in.nextLine().split(" ");
        	
        	int columnID = Integer.parseInt(lineItems[0]);
        	int rowID = Integer.parseInt(lineItems[1]);
        	
        	double value = 0;
        	if(lineItems.length > 2) {
        		value = Double.parseDouble(lineItems[2]);
        	}
        	
        	if(columnID <= numColumns && rowID <= numRows && (lineItems.length <= 2 || value != 0)) {
	        	columnIndex.add(columnID - (1-offset));
	        	currRow = rowID;
	        	if(currRow > prevRow) {
	        		for(int j = prevRow; j < currRow; j++) {
	        			rowPointer.add(i + offset);
	        		}
	        	} else if (currRow < prevRow) {
	        		in.close();
	        		throw new Exception("Matrix Market Format not sorted, aborting");
	        	}
	        	if(lineItems.length > 2) {
	        		data.add(value);
	        	} else {
	        		data.add(1D);
	        	}
	        	prevRow = currRow;
	        	i++;
        	}
		}
        for(int j = rowPointer.size(); j < numRows+1; j++) {
        	rowPointer.add(i + offset);
        }
        
		in.close();
		
		double[] dataArray = new double[data.size()];
		int[] columnIndexArray = new int[columnIndex.size()];
		for(int k = 0; k < data.size(); k++) {
			dataArray[k] = data.get(k).doubleValue();
			columnIndexArray[k] = columnIndex.get(k).intValue();
		}
		
		int[] rowPointerArray = new int[rowPointer.size()];
		for(int k = 0; k < rowPointer.size(); k++) {
			rowPointerArray[k] = rowPointer.get(k).intValue();
		}
		
		return new AdaptedCompRowMatrix(numRows, numColumns, dataArray, columnIndexArray, rowPointerArray);
	}

	private double[] createRandomDoubleArray(int size) {
		Random rand = new Random(seed);
		double[] array = new double[size];
		for(int i = 0; i < size; i++) {
			array[i] = rand.nextDouble();
		}
		return array;
	}
	
    enum Experiment {
        COMPRESSEDROW,
        MODIFIEDCOMPRESSEDROW,
        JAGGED_DIAGONAL,
        ELLPACK
    }
}
