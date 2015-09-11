package edu.tuberlin.sparskit.formats;
import edu.tuberlin.sparskit.Sparskit;
import edu.tuberlin.sparskit.formats.utils.Utils;
import no.uib.cipr.matrix.*;

public class JaggedDiagonalMatrix extends AbstractMatrix{
    
	/**
     * Matrix data in jagged diagonal format
     */
    double[] jdiag;

    /**
     * Column indices. These are kept sorted within each row.
     */
    int[] columnIndex;

    /**
     * Row indices of permutation
     */
    int[] iperm;
    
    /**
     * Pointer to beginnings of each jagged diagonal
     */
    int[] jdPointer;
    
    /**
     * Number of jagged diagonals
     */
    int numJaggedDiagonals;
    
    
	public JaggedDiagonalMatrix(int numRows, int numColumns, double[] jdiag, int[] columnIndex, int[] iperm, int[] jdPointer, int numJaggedDiagonals) {
		super(numRows, numColumns);
		this.jdiag = jdiag;
		this.columnIndex = columnIndex;
		this.iperm = iperm;
		this.jdPointer = jdPointer;
		this.numJaggedDiagonals = numJaggedDiagonals;
	}
	
    public long getMemorySizeInBytes() {
    	return jdiag.length * Double.BYTES + columnIndex.length * Integer.BYTES + iperm.length * Integer.BYTES + jdPointer.length * Integer.BYTES + Integer.BYTES;
    }
	
	@Override
	public Vector mult(Vector x, Vector y) {
		checkMultAdd(x, y);
		DenseVector yDense = Utils.getDenseVector(y);
		Sparskit.amuxj(numRows, Utils.getDenseVector(x).getData(), yDense.getData(), numJaggedDiagonals, jdiag, columnIndex, jdPointer);
		Sparskit.dvperm(numRows, yDense.getData(), iperm);
		return yDense;
	}
	
	public AdaptedCompRowMatrix toAdaptedCompRowMatrix() {
		double[] csrData = new double[jdiag.length];
        int[] csrColumnIndex = new int[jdiag.length];
        int[] csrRowPointer = new int[numRows + 1];
        Sparskit.jadcsr(numRows, numJaggedDiagonals, jdiag, columnIndex, jdPointer, iperm, csrData, csrColumnIndex, csrRowPointer);
    	return new AdaptedCompRowMatrix(numRows, numColumns, csrData, csrColumnIndex, csrRowPointer);
	}
}
