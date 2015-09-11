package edu.tuberlin.sparskit.formats;
import edu.tuberlin.sparskit.Sparskit;
import edu.tuberlin.sparskit.formats.utils.Utils;
import no.uib.cipr.matrix.*;

public class VbrMatrix extends AbstractMatrix{

	/**
     * Matrix data in vbr format
     */
    double[] vbr;

    int[] index1;
    
    int[] index2;
    
    int[] index3;
    
    /**
     * First row number for each block
     */
    int[] kvstr;
    		
    /**
     * First column number for each block
     */
    int[] kvstc;
    
    /**
     * Matrix block row dimension
     */
    int nr;
    
    /**
     *  Matrix block column dimension
     */
    int nc;
    
	public VbrMatrix(int numRows, int numColumns, double[] vbrData, int[] index1, int[] index2, int[] index3, int[] kvstr, int[] kvstc, int nr, int nc) {
		super(numRows, numColumns);
		vbr = vbrData;
		this.index1 = index1;
		this.index2 = index2;
		this.index3 = index3;
		this.kvstr = kvstr;
		this.kvstc = kvstc;
		this.nr = nr;
		this.nc = nc;
	}
	
	@Override
	public Vector mult(Vector x, Vector y) {
		checkMultAdd(x, y);
		DenseVector yDense = Utils.getDenseVector(y);
		Sparskit.vbrmv(nr, nc, index1, index2, index3, vbr, kvstr, kvstc, Utils.getDenseVector(x).getData(), yDense.getData());;
		return yDense;
	}
	
	public AdaptedCompRowMatrix toAdaptedCompRowMatrix() throws Exception {
		double[] csrData = new double[vbr.length];
        int[] csrColumnIndex = new int[csrData.length];
        int[] csrRowPointer = new int[numRows + 1];
        int ierr = Sparskit.vbrcsr(csrColumnIndex, csrRowPointer, csrData, nr, kvstr, kvstc, index1, index2, index3, vbr, csrData.length);
        if(ierr != 0) {
    		throw new Exception("Sparskit VBR to CSR transformation failed. out of space in CSR data structures (found in row " + -ierr + ")");
        }
        return new AdaptedCompRowMatrix(numRows, numColumns, csrData, csrColumnIndex, csrRowPointer);
	}
}
