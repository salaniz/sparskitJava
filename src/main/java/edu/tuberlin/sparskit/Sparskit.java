package edu.tuberlin.sparskit;
public class Sparskit {
	static {
		// sparskit.dll on Windows or libsparskit.so on Linux
	   System.loadLibrary("sparskit");
	}
	
	//private native void readunf(int nmax, int nzmax, int n, int nnz, int a, int[] ja, int[] ia, String iounit);
	
	// ##### Matrix Vector Multiplication #####
	public static native void amux(int n, double[] x, double[] y, double[] a, int[] ja, int[] ia); //A times a vector. Compressed Sparse Row (CSR) format.
	public static native void amuxms(int n, double[] x, double[] y, double[] a, int[] ja); //A times a vector. Modified Compress Sparse Row format
	public static native void atmux(int n, double[] x, double[] y, double[] a, int[] ja, int[] ia); //Transp(A) times a vector. CSR format
	public static native void atmuxr(int m, int n, double[] x, double[] y, double[] a, int[] ja, int[] ia); //Transp(A) times a vector. CSR format. A rectangular
	public static native void amuxe(int n, double[] x, double[] y, int na, int ncol, double[][] a, int[][] ja); //A times a vector. Ellpack/Itpack (ELL) format
	public static native void amuxd(int n, double[] x, double[] y, double[][] diag, int ndiag, int idiag, int[] ioff); //A times a vector. Diagonal (DIA) format
	public static native void amuxj(int n, double[] x, double[] y, int jdiag, double[] a, int[] ja, int[] ia); //A times a vector. Jagged Diagonal (JAD) format
	public static native void vbrmv(int nr, int nc, int[] ia, int[] ja, int[] ka, double[] a, int[] kvstr, int[] kvstc, double[] x, double[] b); //Sparse matrix-full vector product, in VBR format

	// ##### Matrix Format Conversion #####
	public static native int csrdns(int nrow, int ncol, double[] a, int[] ja, int[] ia, double[][] dns, int ndns); //converts a row-stored sparse matrix into the dense format
	public static native int dnscsr(int n, int ncol, int nzmax, double[][] dns, int ndns, double[] a, int[] ja, int[] ia); //converts a dense matrix to a sparse storage format
	public static native void csrmsr(int n, double[] a, int[] ja, int[] ia, double[] ao, int[] jao, double[] wk, int[] iwk); //converts compressed sparse row format to modified sparse row format
	public static native void msrcsr(int n, double[] a, int[] ja, double[] ao, int[] jao, int[] iao, double[] wk, int[]iwk); //converts modified sparse row format to compressed sparse row format
	public static native void csrcsc(int n, int job, int ipos, double[] a, int[] ja,int[] ia, double[] ao, int[] jao, int[] iao); //converts compressed sparse row format to compressed sparse column format (transposition)
	public static native void csrcsc2(int n, int n2, int job, int ipos, double[] a, int[] ja,int[] ia, double[] ao, int[] jao, int[] iao); //rectangular version of csrcsc
	public static native int csrell(int nrow, double[] a, int[] ja, int[] ia, int maxcol, double[][] coef, int[][] jcoef, int ncoef, int[] ierr); //converts compressed sparse row to ellpack format
	public static native int ellcsr(int n, double[][] coef, int[][] jcoef, int ncoef, int ndiag, double[] a, int[] ja, int[] ia, int[] ierr); //converts ellpack format to compressed sparse row format
	public static native void csrdia(int n, int idiag, int job, double[] a, int[] ja, int[] ia, int ndiag, double[][] diag, int[] ioff, double[] ao, int[] jao, int[] iao, int[] ind); //converts a compressed sparse row format into a diagonal format
	public static native void diacsr(int n, int job, int idiag, double[][] diag, int ndiag, int[] ioff, double[] a, int[] ja, int[] ia); //converts a diagonal format into a compressed sparse row format
	public static native int csrjad(int nrow, double[] a, int[] ja, int[] ia, int[] iperm, double[] ao, int[] jao, int[] iao); //converts the csr format into the jagged diagonal format
	public static native void jadcsr(int nrow, int idiag, double[] a, int[] ja, int[] ia, int[] iperm, double[] ao, int[] jao, int[] iao); //converts the jagged-diagonal format into the csr format
	public static native int csrvbr(int n, int[] ia, int[] ja, double[] a, int[] nr, int[] nc, int[] kvstr, int[] kvstc, int[] ib, int[] jb, int[] kb, double[] b, int job, int[] iwk, int nkmax, int nzmax); //converts compressed sparse row to var block row format
	public static native int vbrcsr(int[] ia, int[] ja, double[] a, int nr, int[] kvstr, int[] kvstc, int[] ib, int[] jb, int[] kb, double[] b, int nzmax); //converts var block row to compressed sparse row format
	
	// ##### Statistics Calculations #####
	public static native void bandwidth(int n, int[] ja, int[] ia, int[] ml, int[] mu, int[] iband, double[] bndav);	//Computes  ml     = lower_bandwidth(A)
																														//          mu     = upper_bandwidth(A)
																														//          iband  = max_bandwidth(A)
																														//          bndav  = average_bandwidth(A)
	public static native void nonz(int n, int sym, int[] ja, int[] ia, int[] iao, int[] nzmaxc, int[] nzminc, int[] nzmaxr, int[] nzminr, int[] nzcol, int[] nzrow); //Computes  nzmaxc = max_column_length(A)
																																							             //          nzminc = min_column_length(A)
																																							             //          nzmaxr = max_row_length(A)
																																							             //          nzminr = min_row_length(A)
																																							             //          nzcol  = zero_column_number(A)
																																							             //          nzrow  = zero_row_number(A)
	public static native void diagdomi(int n, int sym, int valued, double[] a, int[] ja, int[] ia, double[] ao, int[] jao, int[] iao, double[] ddomc, double[] ddomr);	//Computes  ddomc  = diag_domi_column_percentage(A)
	             																																									//          ddomr  = diag_domi_row_percentage(A)
	public static native double frobnorm(int n, int sym, double[] a, int[] ja, int[] ia); //Computes  Fnorm  = Frobenius_norm(A)
	public static native void ansym(int n, int sym, double[] a, int[] ja, int[] ia, double[] ao, int[] jao, int[] iao, int[] imatch, double[] av, double[] fas, double[] fan);  //Computes  fas    = sym_part_Frobenius_norm(A)
																																													//          fan    = nonsym_part_Frobenius_norm(A)
																																													//          imatch = matching_elements_number(A)
																																													//          av     = relative_sym_match(A)
	public static native void distaij(int n, int nnz, int sym, int[] ja, int[] ia, double[] dist, double[] std); //Computes  dist   = average_dist_of_a(i,j)(A)
																										//          std    = standard_deviation(A)
	public static native int skyline(int n, int sym, int[] ja, int[] ia, int[] jao, int[] iao); //Computes  nsky   = nonzero_number_in_skyline(A)
	public static native void distdiag(int nrow, int ncol, int[] ja, int[] ia, int[] dist); //Computes  dist   = element_number_in_eachdiag(A)
	public static native int bandpart(int n, int[] ja, int[] ia, int[] dist, int nper); //Computes  band   = bandwidth_width(A)
	public static native int nimpdiag(int n, int nnz, int[] dist, int ipar1, int[] ioff, double[] dcount); //Computes  ndiag  = important_diag_number(A)
	public static native void nonzlud(int n, int[] ja, int[] ia, int[] nlower, int[] nupper, int[] ndiag); //Computes  nlower = nonzero_number_of_lower_part(A)
																											//          nupper = nonzero_number_of_upper_part(A)
																											//          ndiag  = nonzero_number_of_maindiag(A)
	public static native void avnzcol(int n, int[] ja, int[] ia, int[] iao, int ndiag, double[] av, double[] st); 	//Computes  av     = average_nonzero_number_in_column(A)


	// ##### Utilization Functions #####
	public static native void dvperm(int n, double[] x, int[] perm); //permutes a real vector (in-place)
	
	public static native void csrsss(int nrow, double[] a, int[] ja, int[] ia, int sorted, double[] diag, double[] al, int[] jal, int[] ial, double[] au);
	public static native void exper(double[] a);
}
