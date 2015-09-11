package edu.tuberlin.sparskit.formats;
import edu.tuberlin.sparskit.Sparskit;
import edu.tuberlin.sparskit.formats.utils.Utils;
import no.uib.cipr.matrix.*;

public class EllpackMatrix extends AbstractMatrix{

	/**
     * Matrix data in ellpack format
     */
    double[][] ell;

    /**
     * Column indices of the ellpack format.
     */
    int[][] columnIndex;
    
    /**
     * Number of active 'diagonals'
     */
    int numDiagonals;
    
	public EllpackMatrix(int numRows, int numColumns, double[][] ell, int[][] columnIndex, int numDiagonals) {
		super(numRows, numColumns);
		this.ell = ell;
		this.columnIndex = columnIndex;
		this.numDiagonals = numDiagonals;
	}
	
    public long getMemorySizeInBytes() {
    	return ell.length * ell[0].length * Double.BYTES + columnIndex.length * columnIndex[0].length * Integer.BYTES + Integer.BYTES;
    }
	
	@Override
	public Vector mult(Vector x, Vector y) {
		checkMultAdd(x, y);
		DenseVector yDense = Utils.getDenseVector(y);
		Sparskit.amuxe(numRows, Utils.getDenseVector(x).getData(), yDense.getData(), numRows, numDiagonals, ell, columnIndex);
		return yDense;
	}
	
	public AdaptedCompRowMatrix toAdaptedCompRowMatrix() throws Exception {
		double[] csrData = new double[ell.length * ell[0].length];
        int[] csrColumnIndex = new int[csrData.length];
        int[] csrRowPointer = new int[numRows + 1];
        int[] ierr = new int[1];
        int nzmax = Sparskit.ellcsr(numRows, ell, columnIndex, numRows, numDiagonals, csrData, csrColumnIndex, csrRowPointer, ierr);
    	// TODO: use nzmax to crop csrData and csrColumnIndex to size nzmax (space optimization)
        if(ierr[0] == 1) {
    		throw new Exception("Sparskit Ellpack to CSR transformation failed. There is not enough space in CSR data structures to store output matrix.");
        }
        return new AdaptedCompRowMatrix(numRows, numColumns, csrData, csrColumnIndex, csrRowPointer);
	}
}
