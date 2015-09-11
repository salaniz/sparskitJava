package edu.tuberlin.sparskit.formats.utils;

import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseVector;

public class Utils {
	public static DenseVector getDenseVector(Vector x) {
		if(x instanceof DenseVector) {
			return (DenseVector)x;
		} else {
			return new DenseVector(x);
		}
	}
	
    public enum CsrDiaJob {
		FIXED_DIAGONALS_WITHOUT_REMAINDER_DATA(0), 
        SELECT_DIAGONALS_WITHOUT_REMAINDER_DATA(10),
        FIXED_DIAGONALS_WITH_REMAINDER_DATA(1),
        SELECT_DIAGONALS_WITH_REMAINDER_DATA(11);

        private final int val;
        CsrDiaJob(int val) { this.val = val; }
        public int getValue() { return val; }
    }
    
    public enum DiaCsrJob {
		SKIP_ZEROS_IN_DIAGONAL(0), 
        CONVERT_ALL_DIAGONAL_VALUES(1);

        private final int val;
        DiaCsrJob(int val) { this.val = val; }
        public int getValue() { return val; }
    }
    
    public enum CsrEllJob {
		SUPPLY_BLOCK_STRUCTURE(0), 
        CALCULATE_BLOCK_PARTITIONING(1),
        ALT_CALCULATE_BLOCK_CONFORMAL_PARTITIONING(2);

        private final int val;
        CsrEllJob(int val) { this.val = val; }
        public int getValue() { return val; }
    }
}
