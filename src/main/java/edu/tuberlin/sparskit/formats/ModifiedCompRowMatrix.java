package edu.tuberlin.sparskit.formats;
import edu.tuberlin.sparskit.Sparskit;
import edu.tuberlin.sparskit.formats.utils.Utils;
import no.uib.cipr.matrix.*;

public class ModifiedCompRowMatrix extends AbstractMatrix{
    
	/**
     * Matrix data
     */
    double[] data;

    /**
     * Column indices. These are kept sorted within each row.
     */
    int[] columnIndex;
    
    
	public ModifiedCompRowMatrix(int numRows, int numColumns, double[] data, int[] columnIndex) {
		super(numRows, numColumns);
		this.data = data;
		this.columnIndex = columnIndex;
	}
	
    public long getMemorySizeInBytes() {
    	return data.length * Double.BYTES + columnIndex.length * Integer.BYTES;
    }
	
	@Override
	public Vector mult(Vector x, Vector y) {
		checkMultAdd(x, y);
		DenseVector yDense = Utils.getDenseVector(y);
		Sparskit.amuxms(numRows, Utils.getDenseVector(x).getData(), yDense.getData(), data, columnIndex);
		return yDense;
	}
	
	public AdaptedCompRowMatrix toAdaptedCompRowMatrix() {
    	double[] csrData = new double[data.length]; // TODO: true size has to be determined
    	int[] csrColumnIndex = new int[data.length]; // TODO: true size has to be determined
    	int[] csrRowPointer = new int[numRows + 1];
    	int[] intWorkArray;
    	double[] doubleWorkArray;
    	// use in place computations for work arrays if possible
    	if(csrData.length >= numRows) {
    		doubleWorkArray = csrData;
    	} else {
    		doubleWorkArray = new double[numRows];
    	}
    	if(csrRowPointer.length >= numRows + 1) {
    		intWorkArray = csrRowPointer;
    	} else {
    		intWorkArray = new int[numRows + 1];
    	}
    	Sparskit.msrcsr(numRows, data, columnIndex, csrData, csrColumnIndex, csrRowPointer, doubleWorkArray, intWorkArray);
    	return new AdaptedCompRowMatrix(numRows, numColumns, csrData, csrColumnIndex, csrRowPointer);
	}
}
