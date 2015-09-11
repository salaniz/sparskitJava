package edu.tuberlin.sparskit.formats;
import edu.tuberlin.sparskit.Sparskit;
import edu.tuberlin.sparskit.formats.utils.Utils;
import no.uib.cipr.matrix.*;

public class DiagonalMatrix extends AbstractMatrix{

	/**
     * Matrix data in diagonal format
     */
    double[][] diag;

    /**
     * Column indices. These are kept sorted within each row.
     */
    int[] diagonalOffsets;
    
    /**
     * Number of diagonals
     */
    int numDiagonals;

    /**
     * Remainder of the matrix in CSR Format.
     */
    double[] remainderData;
    
    /**
     * Column indices of remainder data.
     */
    int[] remainderColumnIndex;
    
    /**
     * Row pointer of remainder data.
     */
    int[] remainderRowPointer;
    
	public DiagonalMatrix(int numRows, int numColumns, double[][] diag, int[] diagonalOffsets, double[] remainderData, int[] remainderColumnIndex, int[] remainderRowPointer) {
		super(numRows, numColumns);
		this.diag = diag;
		this.diagonalOffsets = diagonalOffsets;
		this.remainderData = remainderData;
		this.remainderColumnIndex = remainderColumnIndex;
		this.remainderRowPointer = remainderRowPointer;
	}
	
	@Override
	public Vector mult(Vector x, Vector y) {
		checkMultAdd(x, y);
		DenseVector yDense = Utils.getDenseVector(y);
		Sparskit.amuxd(numRows, Utils.getDenseVector(x).getData(), yDense.getData(), diag, diag[0].length, numDiagonals, diagonalOffsets);
		return yDense;
	}
	
	public AdaptedCompRowMatrix toAdaptedCompRowMatrix() {
		double[] csrData = new double[diag.length + diag[0].length + remainderData.length];
        int[] csrColumnIndex = new int[csrData.length];
        int[] csrRowPointer = new int[numRows + 1];
        Sparskit.diacsr(numRows, Utils.DiaCsrJob.SKIP_ZEROS_IN_DIAGONAL.getValue(), numDiagonals, diag, diag[0].length, diagonalOffsets, csrData, csrColumnIndex, csrRowPointer);
    	return new AdaptedCompRowMatrix(numRows, numColumns, csrData, csrColumnIndex, csrRowPointer);
	}
}
