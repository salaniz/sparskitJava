package edu.tuberlin.sparskit.formats.utils;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

import edu.tuberlin.sparskit.Sparskit;
import edu.tuberlin.sparskit.formats.AdaptedCompRowMatrix;
import edu.tuberlin.sparskit.formats.EllpackMatrix;
import edu.tuberlin.sparskit.formats.JaggedDiagonalMatrix;
import edu.tuberlin.sparskit.formats.ModifiedCompRowMatrix;
import no.uib.cipr.matrix.DenseMatrix;
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
	
//	public static void main(String[] args) {
//		// invoke the native method
//		
//		double[] a = {1, 2, 3, 4, 5, 6, 7};
//		int[] ja = {1, 3, 3, 1, 2, 3, 3};
//		int[] ia = {1, 3, 4, 7, 8};
//		
////		double[] a = {2, 3, 4, 6};
////		int[] ja = {2, 3, 1, 3};
////		int[] ia = {1, 3, 5};
//		
////		double[] a = {2, 3, 4, 6, 7, 8};
////		int[] ja = {2, 3, 1, 3, 1, 2};
////		int[] ia = {1, 3, 5, 7};
//		
//		int size = 100;
//		double[][] dense = new double[size][size];
//		Scanner in = null;
//		try {
//			in = new Scanner(new FileReader("data/mhd1280b.mtx"));
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
//		int nnz = 0;
//		while(in.hasNext()) {
//			String line = in.nextLine();
//			if(! (line.startsWith("%") || line.equals("1280 1280 12029")) ) {
//				String[] items = line.split(" ");
//				if(Integer.parseInt(items[0])-1 < size && Integer.parseInt(items[1])-1 < size)
//				dense[Integer.parseInt(items[0])-1][Integer.parseInt(items[1])-1] = Double.parseDouble(items[2]);
//				nnz++;
//			}
//		}
//		
//		Random rand = new Random();
//		double[] xArray = new double[size];
//		for (int i = 0; i < size; i++) {
//			xArray[i] = (rand.nextDouble() - 0.5D) * 2D;
//		}
//		
//		DenseVector x = new DenseVector(xArray);
//		DenseVector y = new DenseVector(size);
//		
////		DenseVector x = new DenseVector(new double[]{5,1,2});
////		DenseVector y = new DenseVector(ia.length-1);
//		
////		double[] diag = new double[3];
////		double[] al = new double[5];
////		int[] jal = new int[5];
////		int[] ial = new int[5];
////		double[] au = new double[5];
////		double[][] dns = new double[3][4];
////		boolean sorted = false;
//		//new Sparskit().amux(4,x,y,a,ja,ia);
//		//new Sparskit().csrsss(n, a, ja, ia, sorted ? 1 : 0, diag, al, jal, ial, au);		
//		
//		
//		// CSR TEST
//		DenseVector result;
//		AdaptedCompRowMatrix csrMatrix = new AdaptedCompRowMatrix(ia.length-1,3,a,ja,ia);
//		AdaptedCompRowMatrix transCsrMatrix = csrMatrix.getTranspose();
//		AdaptedCompRowMatrix test = csrMatrix.getTranspose().getTranspose();
//		System.out.println(Arrays.toString(test.data.array()));
//		int[] ml = new int[1];
//		int[] mu = new int[1];
//		int[] iband = new int[1];
//		double[] bndav = new double[1];
//		Sparskit.bandwidth(csrMatrix.numColumns(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), ml, mu, iband, bndav);
//		System.out.println(ml[0] + " " + mu[0] + " " + iband[0] + " " + bndav[0]);
//		int[] nzmaxc = new int[1];
//		int[] nzminc = new int[1];
//		int[] nzmaxr = new int[1];
//		int[] nzminr = new int[1];
//		int[] nzcol = new int[1];
//		int[] nzrow = new int[1];
//		Sparskit.nonz(csrMatrix.numColumns(), false ? 1 : 0, transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), csrMatrix.columnIndex.array(), nzmaxc, nzminc, nzmaxr, nzminr, nzcol, nzrow);
//		System.out.println(nzmaxc[0] + " " + nzminc[0] + " " + nzmaxr[0] + " " + nzminr[0] + " " + nzcol[0] + " " + nzrow[0]);
//		double[] ddomc = new double[1];
//		double[] ddomr = new double[1];
//		Sparskit.diagdomi(csrMatrix.numColumns(), false ? 1 : 0, true ? 1 : 0, transCsrMatrix.data.array(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), csrMatrix.data.array(), csrMatrix.columnIndex.array(), csrMatrix.rowPointer.array(), ddomc, ddomr);
//		System.out.println(ddomc[0]  + " " + ddomr[0]);
//		double frobnorm = Sparskit.frobnorm(csrMatrix.numColumns(), false ? 1 : 0, transCsrMatrix.data.array(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array());
//		System.out.println(frobnorm);
//		int[] imatch = new int[1];
//		double[] av = new double[1];
//		double[] fas = new double[1];
//		double[] fan = new double[1];
//		Sparskit.ansym(csrMatrix.numColumns(), false ? 1 : 0, transCsrMatrix.data.array(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), new double[csrMatrix.data.array().length], new int[csrMatrix.columnIndex.array().length], new int[csrMatrix.rowPointer.array().length], imatch, av, fas, fan);
//		System.out.println(imatch[0] + " " + av[0] + " " + fas[0] + " " + fan[0]);
//		double[] dist = new double[1];
//		double[] std = new double[1];
//		Sparskit.distaij(csrMatrix.numColumns(), transCsrMatrix.data.array().length, false ? 1 : 0, transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), dist, std);
//		System.out.println(dist[0] + " " + std[0]);
//		int nsky = Sparskit.skyline(csrMatrix.numColumns(), false ? 1 : 0, transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), csrMatrix.columnIndex.array(), csrMatrix.rowPointer.array());
//		System.out.println(nsky);
//		int[] dists = new int[csrMatrix.numRows() + csrMatrix.numColumns() - 1];
//		Sparskit.distdiag(csrMatrix.numRows(), csrMatrix.numColumns(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), dists);
//		System.out.println(Arrays.toString(dists));
//		int[] nlower = new int[1];
//		int[] nupper = new int[1];
//		int[] ndiag = new int[1];
//		Sparskit.nonzlud(csrMatrix.numColumns(), transCsrMatrix.columnIndex.array(), transCsrMatrix.rowPointer.array(), nlower, nupper, ndiag);
//		System.out.println(nlower[0] + " " + nupper[0] + " " + ndiag[0]);
//		
//		System.out.println(Arrays.toString(csrMatrix.data.array()));
//		System.out.println(Arrays.toString(csrMatrix.columnIndex.array()));
//		System.out.println(Arrays.toString(csrMatrix.rowPointer.array()));
////		result = (DenseVector) csrMatrix.multSparskit(x, y);
////		System.out.println(Arrays.toString(result.getData()));
//		
//		// Dense TEST
//		double[][] dnsMatrix;
//		try {
//			dnsMatrix = csrMatrix.getTranspose().toDense();
//			for(double[] row : dnsMatrix) {
//				System.out.println(Arrays.toString(row));
//			}
//		} catch (Exception e1) {
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		}
////		csrMatrix = new AdaptedCompRowMatrix(dnsMatrix, a.length);
//		System.out.println(nnz);
//		csrMatrix = new AdaptedCompRowMatrix(new DenseMatrix(dense));
//		result = (DenseVector) csrMatrix.multSparskit(x, y);
////		System.out.println(Arrays.toString(result.getData()));
//
//		
//		// Modified CSR TEST
//		System.out.println("Modified CSR");
//		ModifiedCompRowMatrix modMatrix = csrMatrix.toModifiedCompRowMatrix();
////		y = new DenseVector(ia.length-1);
//		y = new DenseVector(size);
//		result = (DenseVector) modMatrix.mult(x, y);
//		System.out.println(Arrays.toString(result.getData()));
////		y = new DenseVector(ia.length-1);
//		y = new DenseVector(size);
//		modMatrix.toAdaptedCompRowMatrix().multSparskit(x, y);
//		System.out.println(Arrays.toString(result.getData()));	
//
//		
//		// Jagged Diagonal TEST
//		System.out.println("Jagged");
//		JaggedDiagonalMatrix jadMatrix = csrMatrix.toJaggedDiagonalMatrix();
////		y = new DenseVector(ia.length-1);
//		y = new DenseVector(size);
//		result = (DenseVector) jadMatrix.mult(x, y);
//		System.out.println(Arrays.toString(result.getData()));
////		y = new DenseVector(ia.length-1);
//		y = new DenseVector(size);
//		jadMatrix.toAdaptedCompRowMatrix().multSparskit(x, y);
//		System.out.println(Arrays.toString(result.getData()));
//		
////		// Diagonal TEST
////		System.out.println("DIAGONAL");
////		DiagonalMatrix diaMatrix = csrMatrix.toDiagonalMatrix(Utils.CsrDiaJob.SELECT_DIAGONALS_WITH_REMAINDER_DATA, size, 14400);
//////		y = new DenseVector(ia.length-1);
////		y = new DenseVector(size);
////		result = (DenseVector) diaMatrix.mult(x, y);
////		System.out.println("Result: " + Arrays.toString(result.getData()));
////		System.out.println("DIAG:");
////		for(double[] row : diaMatrix.diag)
////			System.out.println(Arrays.toString(row));
////		System.out.println("Diag Offsets: " + Arrays.toString(diaMatrix.diagonalOffsets));
//////		y = new DenseVector(ia.length-1);
////		y = new DenseVector(size);
////		diaMatrix.toAdaptedCompRowMatrix().multSparskit(x, y);
////		csrMatrix = diaMatrix.toAdaptedCompRowMatrix();
////		System.out.println("Diag to CSR");
////		System.out.println(Arrays.toString(csrMatrix.data.array()));
////		System.out.println(Arrays.toString(result.getData()));
//		
//		
//		//Ellpack TEST
//		System.out.println("Ellpack");
//		EllpackMatrix ellMatrix;
//		try {
//			ellMatrix = csrMatrix.toEllpackMatrix();
////			y = new DenseVector(ia.length-1);
//			y = new DenseVector(size);
//			result = (DenseVector) ellMatrix.mult(x, y);
//			System.out.println(Arrays.toString(result.getData()));
////			y = new DenseVector(ia.length-1);
//			y = new DenseVector(size);
//			ellMatrix.toAdaptedCompRowMatrix().multSparskit(x, y);
//			System.out.println(Arrays.toString(result.getData()));
//		} catch (Exception e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		
//		
//////		System.out.println("VBR");
//////		//VBR TEST
////		VbrMatrix vbrMatrix = csrMatrix.toVbrMatrix(Utils.CsrEllJob.CALCULATE_BLOCK_PARTITIONING);
//////		y = new DenseVector(ia.length-1);
////		y = new DenseVector(size);
////		result = (DenseVector) vbrMatrix.mult(x, y);
////		System.out.println(Arrays.toString(result.getData()));
//////		y = new DenseVector(ia.length-1);
////		y = new DenseVector(size);
////		vbrMatrix.toAdaptedCompRowMatrix().multSparskit(x, y);
////		System.out.println(Arrays.toString(result.getData()));
//	}
}
