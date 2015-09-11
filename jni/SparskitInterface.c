#include <jni.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "edu_tuberlin_sparskit_Sparskit.h"

// ##### matrix vector multiplication #####
extern void amux_(int *, double[], double[], double[], int[], int[]); //A times a vector. Compressed Sparse Row (CSR) format.
extern void amuxms_(int *, double[], double[], double[], int[]); //A times a vector. Modified Compress Sparse Row format
extern void atmux_(int *, double[], double[], double[], int[], int[]); //Transp(A) times a vector. CSR format
extern void atmuxr_(int *, int *, double[], double[], double[], int[], int[]); //Transp(A) times a vector. CSR format. A rectangular
extern void amuxe_(int *, double[], double[], int *, int *, double *, int *); //A times a vector. Ellpack/Itpack (ELL) format
extern void amuxd_(int *, double[], double[], double *, int *, int *, int[]); //A times a vector. Diagonal (DIA) format
extern void amuxj_(int *, double[], double[], int *, double[], int[], int[]); //A times a vector. Jagged Diagonal (JAD) format
extern void vbrmv_(int *, int *, int[], int[], int[], double[], int[], int[], double[], double[]); //Sparse matrix-full vector product, in VBR format

// ##### format conversion #####
extern void csrdns_(int *, int *, double[], int[], int[], double *, int *, int *); //converts a row-stored sparse matrix into the dense format
extern void dnscsr_(int *, int *, int *, double *, int *, double[], int[], int[], int *); //converts a dense matrix to a sparse storage format
extern void csrmsr_(int *, double[], int[], int[], double[], int[], double[], int[]); //converts compressed sparse row format to modified sparse row format
extern void msrcsr_(int *, double[], int[], double[], int[], int[], double[], int[]); //converts modified sparse row format to compressed sparse row format
extern void csrcsc_(int *, int *, int *, double[], int[], int[], double[], int[], int[]); //converts compressed sparse row format to compressed sparse column format (transposition)
extern void csrcsc2_(int *, int *, int *, int *, double[], int[], int[], double[], int[], int[]); //rectangular version of csrcsc
extern void csrell_(int *, double[], int[], int[], int *, double *, int *, int *, int *, int *); //converts compressed sparse row to ellpack format
extern void ellcsr_(int *, double *, int *, int *, int *, double[], int[], int[], int *, int *); //converts ellpack format to compressed sparse row format
extern void csrdia_(int *, int *, int *, double[], int[], int[], int *, double *, int[], double[], int[], int[], int[]); //converts a compressed sparse row format into a diagonal format
extern void diacsr_(int *, int *, int *, double *, int *, int[], double[], int[], int[]); //converts a diagonal format into a compressed sparse row format
extern void csrjad_(int *, double[], int[], int[], int *, int[], double[], int[], int[]); //converts the csr format into the jagged diagonal format
extern void jadcsr_(int *, int *, double[], int[], int[], int[], double[], int[], int[]); //converts the jagged-diagonal format into the csr format
extern void csrvbr_(int *, int[], int[], double[], int *, int *, int[], int[], int[], int[], int[], double[], int *, int[], int *, int *, int *); //converts compressed sparse row to var block row format
extern void vbrcsr_(int[], int[], double[], int *, int[], int[], int[], int[], int[], double[], int *, int *); //converts var block row to compressed sparse row format

// ##### statistics #####
extern void bandwidth_(int *, int[], int[], int *, int *, int *, double *);	//Computes  ml     = lower_bandwidth(A)
											//          mu     = upper_bandwidth(A)
											//          iband  = max_bandwidth(A)
											//          bndav  = average_bandwidth(A)
extern void nonz_(int *, int *, int[], int[], int[], int *, int *, int *, int *, int *, int *); //Computes  nzmaxc = max_column_length(A)
             //          nzminc = min_column_length(A)
             //          nzmaxr = max_row_length(A)
             //          nzminr = min_row_length(A)
             //          nzcol  = zero_column_number(A)
             //          nzrow  = zero_row_number(A)
extern void diag_domi_(int *, int *, int *, double[], int[], int[], double[], int[], int[], double *, double *);  //Computes  ddomc  = diag_domi_column_percentage(A)
             //          ddomr  = diag_domi_row_percentage(A)
extern void frobnorm_(int *, int *, double[], int[], int[], double *); //Computes  Fnorm  = Frobenius_norm(A)
extern void ansym_(int *, int *, double[], int[], int[], double[], int[], int[], int *, double *, double *, double *); //Computes  fas    = sym_part_Frobenius_norm(A)
             //          fan    = nonsym_part_Frobenius_norm(A)
             //          imatch = matching_elements_number(A)
             //          av     = relative_sym_match(A)
extern void distaij_(int *, int *, int *, int[], int[], double *, double *); //Computes  dist   = average_dist_of_a(i,j)(A)
             //          std    = standard_deviation(A)
extern void skyline_(int *, int *, int[], int[], int[], int[], int *); //Computes  nsky   = nonzero_number_in_skyline(A)
extern void distdiag_(int *, int *, int[], int[], int[]); //Computes  dist   = element_number_in_eachdiag(A)
extern void bandpart_(int *, int[], int[], int[], int *, int *); //Computes  band   = bandwidth_width(A)
extern void n_imp_diag_(int *, int *, int[], int *, int *, int[], double[]); //Computes  ndiag  = important_diag_number(A)
extern void nonz_lud_(int *, int[], int[], int *, int *, int *); //Computes  nlower = nonzero_number_of_lower_part(A)
             //          nupper = nonzero_number_of_upper_part(A)
             //          ndiag  = nonzero_number_of_maindiag(A)
extern void avnz_col_(int *, int[], int[], int[], int *, double *, double *); //Computes  av     = average_nonzero_number_in_column(A)
             //          st     = standard_deviation(A)

// ##### utils #####
extern void dvperm_(int *, double[], int[]); //permutes a real vector (in-place)

//not needed for now
extern void csrsss_(int *, double[], int[], int[], int *, double[], double[], int[], int[], double[]);

void dense_transpose(double *matrix, double *destination_matrix, int nrows, int ncols) {
	for(int i = 0; i < nrows; i++) {
		for(int j = 0; j < ncols; j++) {
			destination_matrix[j*ncols+i] = matrix[i*nrows+j];
		}
	}
}

// ##### matrix vector multiplication #####
JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_amux(JNIEnv *env, jobject thisObj, jint n, jdoubleArray x, jdoubleArray y, jdoubleArray a, jintArray ja, jintArray ia) {

	jdouble *c_x = (*env)->GetPrimitiveArrayCritical(env, x, NULL);
	jdouble *c_y = (*env)->GetPrimitiveArrayCritical(env, y, NULL);
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);

	amux_(&n, c_x, c_y, c_a, c_ja, c_ia);

	(*env)->ReleasePrimitiveArrayCritical(env, x, c_x, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, y, c_y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_amuxms(JNIEnv *env, jobject thisObj, jint n, jdoubleArray x, jdoubleArray y, jdoubleArray a, jintArray ja) {
	jdouble *c_x = (*env)->GetPrimitiveArrayCritical(env, x, NULL);
	jdouble *c_y = (*env)->GetPrimitiveArrayCritical(env, y, NULL);
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);

	amuxms_(&n, c_x, c_y, c_a, c_ja);

	(*env)->ReleasePrimitiveArrayCritical(env, x, c_x, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, y, c_y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_atmux(JNIEnv *env, jobject thisObj, jint n, jdoubleArray x, jdoubleArray y, jdoubleArray a, jintArray ja, jintArray ia) {
	jdouble *c_x = (*env)->GetPrimitiveArrayCritical(env, x, NULL);
	jdouble *c_y = (*env)->GetPrimitiveArrayCritical(env, y, NULL);
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);

	atmux_(&n, c_x, c_y, c_a, c_ja, c_ia);

	(*env)->ReleasePrimitiveArrayCritical(env, x, c_x, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, y, c_y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_atmuxr(JNIEnv *env, jobject thisObj, jint m, jint n, jdoubleArray x, jdoubleArray y, jdoubleArray a, jintArray ja, jintArray ia) {
	jdouble *c_x = (*env)->GetPrimitiveArrayCritical(env, x, NULL);
	jdouble *c_y = (*env)->GetPrimitiveArrayCritical(env, y, NULL);
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);

	atmuxr_(&m, &n, c_x, c_y, c_a, c_ja, c_ia);

	(*env)->ReleasePrimitiveArrayCritical(env, x, c_x, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, y, c_y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_amuxe(JNIEnv *env, jobject thisObj, jint n, jdoubleArray x, jdoubleArray y, jint na, jint ncol, jobjectArray a, jobjectArray ja) {
	jdouble *c_x = (*env)->GetPrimitiveArrayCritical(env, x, NULL);
	jdouble *c_y = (*env)->GetPrimitiveArrayCritical(env, y, NULL);

	jsize c_na = na;
	jsize c_ncol = ncol;
	jsize i, j;

	double c_a[ncol][na];
	for(i=0; i<c_ncol; i++) {
		jdoubleArray a1_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, a, i));
		jdouble *a2_tmp = (*env)->GetPrimitiveArrayCritical(env, a1_tmp, NULL);
		for(j=0; j<c_na; j++) {
			c_a[i][j] = a2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, a1_tmp, a2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, a1_tmp);
	}

	int c_ja[ncol][na];
	for(i=0; i<c_ncol; i++) {
		jintArray ja1_tmp = (jintArray)((*env)->GetObjectArrayElement(env, ja, i));
		jint *ja2_tmp = (*env)->GetPrimitiveArrayCritical(env, ja1_tmp, NULL);
		for(j=0; j<c_na; j++) {
			c_ja[i][j] = ja2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, ja1_tmp, ja2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, ja1_tmp);
	}


	amuxe_(&n, c_x, c_y, &na, &ncol, (double*)c_a, (int*)c_ja);

	(*env)->ReleasePrimitiveArrayCritical(env, x, c_x, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, y, c_y, 0);

	for(i=0; i<c_ncol; i++) {
		jdoubleArray a_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, a, i));
		(*env)->SetDoubleArrayRegion(env, a_tmp, (jsize)0, c_na, (jdouble*)c_a[i]);
		(*env)->SetObjectArrayElement(env, a, i, a_tmp);
		(*env)->DeleteLocalRef(env, a_tmp);
	}

	for(i=0; i<c_ncol; i++) {
		jintArray ja_tmp = (jintArray)((*env)->GetObjectArrayElement(env, ja, i));
		(*env)->SetIntArrayRegion(env, ja_tmp, (jsize)0, c_na, (jint*)c_ja[i]);
		(*env)->SetObjectArrayElement(env, ja, i, ja_tmp);
		(*env)->DeleteLocalRef(env, ja_tmp);
	}
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_amuxd(JNIEnv *env, jobject thisObj, jint n, jdoubleArray x, jdoubleArray y, jobjectArray diag, jint ndiag, jint idiag, jintArray ioff) {
	jdouble *c_x = (*env)->GetPrimitiveArrayCritical(env, x, NULL);
	jdouble *c_y = (*env)->GetPrimitiveArrayCritical(env, y, NULL);
	jint *c_ioff = (*env)->GetPrimitiveArrayCritical(env, ioff, NULL);

	jsize c_ndiag = ndiag;
	jsize c_idiag = idiag;
	jsize i, j;

	double c_diag[idiag][ndiag];
	for(i=0; i<c_idiag; i++) {
		jdoubleArray diag1_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, diag, i));
		jdouble *diag2_tmp = (*env)->GetPrimitiveArrayCritical(env, diag1_tmp, NULL);
		for(j=0; j<c_ndiag; j++) {
			c_diag[i][j] = diag2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, diag1_tmp, diag2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, diag1_tmp);
	}

	amuxd_(&n, c_x, c_y, (double*)c_diag, &ndiag, &idiag, c_ioff);

	(*env)->ReleasePrimitiveArrayCritical(env, x, c_x, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, y, c_y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ioff, c_ioff, JNI_ABORT);

	for(i=0; i<c_idiag; i++) {
		jintArray diag_tmp = (jintArray)((*env)->GetObjectArrayElement(env, diag, i));
		(*env)->SetIntArrayRegion(env, diag_tmp, (jsize)0, c_ndiag, (jint*)c_diag[i]);
		(*env)->SetObjectArrayElement(env, diag, i, diag_tmp);
		(*env)->DeleteLocalRef(env, diag_tmp);
	}
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_amuxj(JNIEnv *env, jobject thisObj, jint n, jdoubleArray x, jdoubleArray y, jint jdiag, jdoubleArray a, jintArray ja, jintArray ia) {
	jdouble *c_x = (*env)->GetPrimitiveArrayCritical(env, x, NULL);
	jdouble *c_y = (*env)->GetPrimitiveArrayCritical(env, y, NULL);
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);

	amuxj_(&n, c_x, c_y, &jdiag, c_a, c_ja, c_ia);

	(*env)->ReleasePrimitiveArrayCritical(env, x, c_x, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, y, c_y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);

}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_vbrmv(JNIEnv *env, jobject thisObj, jint nr, jint nc, jintArray ia, jintArray ja, jintArray ka, jdoubleArray a, jintArray kvstr, jintArray kvstc, jdoubleArray x, jdoubleArray b) {
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_ka = (*env)->GetPrimitiveArrayCritical(env, ka, NULL);
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_kvstr = (*env)->GetPrimitiveArrayCritical(env, kvstr, NULL);
	jint *c_kvstc = (*env)->GetPrimitiveArrayCritical(env, kvstc, NULL);
	jdouble *c_x = (*env)->GetPrimitiveArrayCritical(env, x, NULL);
	jdouble *c_b = (*env)->GetPrimitiveArrayCritical(env, b, NULL);

	vbrmv_(&nr, &nc, c_ia, c_ia, c_ka, c_a, c_kvstr, c_kvstc, c_x, c_b);

	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ka, c_ka, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, kvstr, c_kvstr, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, kvstc, c_kvstc, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, x, c_x, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, b, c_b, 0);
}















// ##### format conversion #####
JNIEXPORT jint JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrdns(JNIEnv *env, jobject thisObj, jint nrow, jint ncol, jdoubleArray a, jintArray ja, jintArray ia, jobjectArray dns, jint ndns) {

	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);

	jsize c_ndns = ndns;
	jsize c_ncol = ncol;
	jsize i, j;
	double c_dns[ncol][ndns];
	for(i=0; i<c_ncol; i++) {
		jdoubleArray dns1_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, dns, i));
		jdouble *dns2_tmp = (*env)->GetPrimitiveArrayCritical(env, dns1_tmp, NULL);
		for(j=0; j<c_ndns; j++) {
			c_dns[i][j] = dns2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, dns1_tmp, dns2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, dns1_tmp);
	}

	jint ierr;

	csrdns_(&nrow, &ncol, c_a, c_ja, c_ia, (double*)c_dns, &ndns, &ierr);


	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);

	for(i=0; i<c_ncol; i++) {
		jdoubleArray dns_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, dns, i));
		(*env)->SetDoubleArrayRegion(env, dns_tmp, (jsize)0, c_ndns, (jdouble*)c_dns[i]);
		(*env)->SetObjectArrayElement(env, dns, i, dns_tmp);
		(*env)->DeleteLocalRef(env, dns_tmp);
	}

	return ierr;
}


JNIEXPORT jint JNICALL Java_edu_tuberlin_sparskit_Sparskit_dnscsr(JNIEnv *env, jobject thisObj, jint nrow, jint ncol, jint nzmax, jobjectArray dns, jint ndns, jdoubleArray a, jintArray ja, jintArray ia) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);

	jsize c_ndns = ndns;
	jsize c_ncol = ncol;
	jsize i, j;
	double c_dns[ncol][ndns];
	for(i=0; i<c_ncol; i++) {
		jdoubleArray dns1_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, dns, i));
		jdouble *dns2_tmp = (*env)->GetPrimitiveArrayCritical(env, dns1_tmp, NULL);
		for(j=0; j<c_ndns; j++) {
			c_dns[i][j] = dns2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, dns1_tmp, dns2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, dns1_tmp);
	}

	jint ierr;

	dnscsr_(&nrow, &ncol, &nzmax, (double*)c_dns, &ndns, c_a, c_ja, c_ia, &ierr);

	//TODO: release c_dns?

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, 0);

	return ierr;
}


JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrmsr(JNIEnv *env, jobject thisObj, jint n, jdoubleArray a, jintArray ja, jintArray ia, jdoubleArray ao, jintArray jao, jdoubleArray wk, jintArray iwk) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jdouble *c_ao = (*env)->GetPrimitiveArrayCritical(env, ao, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jdouble *c_wk = (*env)->GetPrimitiveArrayCritical(env, wk, NULL);
	jint *c_iwk = (*env)->GetPrimitiveArrayCritical(env, iwk, NULL);

	csrmsr_(&n, c_a, c_ja, c_ia, c_ao, c_jao, c_wk, c_iwk);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ao, c_ao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, wk, c_wk, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, iwk, c_iwk, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_msrcsr(JNIEnv *env, jobject thisObj, jint n, jdoubleArray a, jintArray ja, jdoubleArray ao, jintArray jao, jintArray iao, jdoubleArray wk, jintArray iwk) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jdouble *c_ao = (*env)->GetPrimitiveArrayCritical(env, ao, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);
	jdouble *c_wk = (*env)->GetPrimitiveArrayCritical(env, wk, NULL);
	jint *c_iwk = (*env)->GetPrimitiveArrayCritical(env, iwk, NULL);

	msrcsr_(&n, c_a, c_ja, c_ao, c_jao, c_iao, c_wk, c_iwk);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ao, c_ao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, wk, c_wk, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, iwk, c_iwk, JNI_ABORT);
}


JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrcsc(JNIEnv *env, jobject thisObj, jint n, jint job, jint ipos, jdoubleArray a, jintArray ja, jintArray ia, jdoubleArray ao, jintArray jao, jintArray iao) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jdouble *c_ao = (*env)->GetPrimitiveArrayCritical(env, ao, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);

	csrcsc_(&n, &job, &ipos, c_a, c_ja, c_ia, c_ao, c_jao, c_iao);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ao, c_ao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, 0);

}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrcsc2(JNIEnv *env, jobject thisObj, jint n, jint n2, jint job, jint ipos, jdoubleArray a, jintArray ja, jintArray ia, jdoubleArray ao, jintArray jao, jintArray iao) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jdouble *c_ao = (*env)->GetPrimitiveArrayCritical(env, ao, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);

	csrcsc2_(&n, &n2, &job, &ipos, c_a, c_ja, c_ia, c_ao, c_jao, c_iao);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ao, c_ao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, 0);
}

JNIEXPORT jint JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrell(JNIEnv *env, jobject thisObj, jint nrow, jdoubleArray a, jintArray ja, jintArray ia, jint maxcol, jobjectArray coef, jobjectArray jcoef, jint ncoef, jintArray ierr) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_ierr = (*env)->GetPrimitiveArrayCritical(env, ierr, NULL);

	jsize c_ncoef = ncoef;
	jsize c_maxcol = maxcol;
	jsize i, j;
	double c_coef[maxcol][ncoef];
	for(i=0; i<c_maxcol; i++) {
		jdoubleArray coef1_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, coef, i));
		jdouble *coef2_tmp = (*env)->GetPrimitiveArrayCritical(env, coef1_tmp, NULL);
		for(j=0; j<c_ncoef; j++) {
			c_coef[i][j] = coef2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, coef1_tmp, coef2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, coef1_tmp);
	}

	int c_jcoef[maxcol][ncoef];
	for(i=0; i<c_maxcol; i++) {
		jintArray jcoef1_tmp = (jintArray)((*env)->GetObjectArrayElement(env, jcoef, i));
		jint *jcoef2_tmp = (*env)->GetPrimitiveArrayCritical(env, jcoef1_tmp, NULL);
		for(j=0; j<c_ncoef; j++) {
			c_jcoef[i][j] = jcoef2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, jcoef1_tmp, jcoef2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, jcoef1_tmp);
	}

	jint ndiag;

	csrell_(&nrow, c_a, c_ja, c_ia, &maxcol, (double*)c_coef, (int*)c_jcoef, &ncoef, &ndiag, (int*)c_ierr);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ierr, c_ierr, 0);

	for(i=0; i<c_maxcol; i++) {
		jdoubleArray coef_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, coef, i));
		(*env)->SetDoubleArrayRegion(env, coef_tmp, (jsize)0, c_ncoef, (jdouble*)c_coef[i]);
		(*env)->SetObjectArrayElement(env, coef, i, coef_tmp);
		(*env)->DeleteLocalRef(env, coef_tmp);
	}

	for(i=0; i<c_maxcol; i++) {
		jintArray jcoef_tmp = (jintArray)((*env)->GetObjectArrayElement(env, jcoef, i));
		(*env)->SetIntArrayRegion(env, jcoef_tmp, (jsize)0, c_ncoef, (jint*)c_jcoef[i]);
		(*env)->SetObjectArrayElement(env, jcoef, i, jcoef_tmp);
		(*env)->DeleteLocalRef(env, jcoef_tmp);
	}

	return ndiag;
}

JNIEXPORT jint JNICALL Java_edu_tuberlin_sparskit_Sparskit_ellcsr(JNIEnv *env, jobject thisObj, jint nrow, jobjectArray coef, jobjectArray jcoef, jint ncoef, jint ndiag, jdoubleArray a, jintArray ja, jintArray ia, jintArray ierr) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_ierr = (*env)->GetPrimitiveArrayCritical(env, ierr, NULL);

	jsize c_ncoef = ncoef;
	jsize c_ndiag = ndiag;
	jsize i, j;
	double c_coef[ndiag][ncoef];
	for(i=0; i<c_ndiag; i++) {
		jdoubleArray coef1_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, coef, i));
		jdouble *coef2_tmp = (*env)->GetPrimitiveArrayCritical(env, coef1_tmp, NULL);
		for(j=0; j<c_ncoef; j++) {
			c_coef[i][j] = coef2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, coef1_tmp, coef2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, coef1_tmp);
	}

	int c_jcoef[ndiag][ncoef];
	for(i=0; i<c_ndiag; i++) {
		jintArray jcoef1_tmp = (jintArray)((*env)->GetObjectArrayElement(env, jcoef, i));
		jint *jcoef2_tmp = (*env)->GetPrimitiveArrayCritical(env, jcoef1_tmp, NULL);
		for(j=0; j<c_ncoef; j++) {
			c_jcoef[i][j] = jcoef2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, jcoef1_tmp, jcoef2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, jcoef1_tmp);
	}

	jint nzmax;

	ellcsr_(&nrow, (double*)c_coef, (int*)c_jcoef, &ncoef, &ndiag, c_a, c_ja, c_ia, &nzmax, (int*)c_ierr);

	//TODO: release c_coef, c_jcoef?

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ierr, c_ierr, 0);

	return nzmax;
}


JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrdia(JNIEnv *env, jobject thisObj, jint n, jint idiag, jint job, jdoubleArray a, jintArray ja, jintArray ia, jint ndiag, jobjectArray diag, jintArray ioff, jdoubleArray ao, jintArray jao, jintArray iao, jintArray ind) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jdouble *c_ao = (*env)->GetPrimitiveArrayCritical(env, ao, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);
	jint *c_ioff = (*env)->GetPrimitiveArrayCritical(env, ioff, NULL);
	jint *c_ind = (*env)->GetPrimitiveArrayCritical(env, ind, NULL);

	jsize c_ndiag = ndiag;
	jsize c_idiag = idiag;
	jsize i, j;
	double c_diag[idiag][ndiag];
	for(i=0; i<c_idiag; i++) {
		jdoubleArray diag1_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, diag, i));
		jdouble *diag2_tmp = (*env)->GetPrimitiveArrayCritical(env, diag1_tmp, NULL);
		for(j=0; j<c_ndiag; j++) {
			c_diag[i][j] = diag2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, diag1_tmp, diag2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, diag1_tmp);
	}

	csrdia_(&n, &idiag, &job, c_a, c_ja, c_ia, &ndiag, (double*)c_diag, c_ioff, c_ao, c_jao, c_iao, c_ind);

	for(i=0; i<c_idiag; i++) {
		jdoubleArray diag_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, diag, i));
		(*env)->SetDoubleArrayRegion(env, diag_tmp, (jsize)0, c_ndiag, (jdouble*)c_diag[i]);
		(*env)->SetObjectArrayElement(env, diag, i, diag_tmp);
		(*env)->DeleteLocalRef(env, diag_tmp);
	}


	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ao, c_ao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ioff, c_ioff, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ind, c_ind, JNI_ABORT);
}


JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_diacsr(JNIEnv *env, jobject thisObj, jint n, jint job, jint idiag, jobjectArray diag, jint ndiag, jintArray ioff, jdoubleArray a, jintArray ja, jintArray ia) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_ioff = (*env)->GetPrimitiveArrayCritical(env, ioff, NULL);

	jsize c_ndiag = ndiag;
	jsize c_idiag = idiag;
	jsize i, j;
	double c_diag[idiag][ndiag];
	for(i=0; i<c_idiag; i++) {
		jdoubleArray diag1_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, diag, i));
		jdouble *diag2_tmp = (*env)->GetPrimitiveArrayCritical(env, diag1_tmp, NULL);
		for(j=0; j<c_ndiag; j++) {
			c_diag[i][j] = diag2_tmp[j];
		}
		(*env)->ReleasePrimitiveArrayCritical(env, diag1_tmp, diag2_tmp, JNI_ABORT);
		(*env)->DeleteLocalRef(env, diag1_tmp);
	}

	diacsr_(&n, &job, &idiag, (double*)diag, &ndiag, c_ioff, c_a, c_ja, c_ia);

	//TODO: release c_diag?

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ioff, c_ioff, JNI_ABORT);
}

JNIEXPORT jint JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrjad(JNIEnv *env, jobject thisObj, jint nrow, jdoubleArray a, jintArray ja, jintArray ia, jintArray iperm, jdoubleArray ao, jintArray jao, jintArray iao) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jdouble *c_ao = (*env)->GetPrimitiveArrayCritical(env, ao, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);
	jint *c_iperm = (*env)->GetPrimitiveArrayCritical(env, iperm, NULL);

	jint idiag;

	csrjad_(&nrow, c_a, c_ja, c_ia, &idiag, c_iperm, c_ao, c_jao, c_iao);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ao, c_ao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iperm, c_iperm, 0);

	return idiag;
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_jadcsr(JNIEnv *env, jobject thisObj, jint nrow, jint idiag, jdoubleArray a, jintArray ja, jintArray ia, jintArray iperm, jdoubleArray ao, jintArray jao, jintArray iao) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jdouble *c_ao = (*env)->GetPrimitiveArrayCritical(env, ao, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);
	jint *c_iperm = (*env)->GetPrimitiveArrayCritical(env, iperm, NULL);

	jadcsr_(&nrow, &idiag, c_a, c_ja, c_ia, c_iperm, c_ao, c_jao, c_iao);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ao, c_ao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iperm, c_iperm, JNI_ABORT);
}

JNIEXPORT jint JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrvbr(JNIEnv *env, jobject thisObj, jint n, jintArray ia, jintArray ja, jdoubleArray a, jintArray nr, jintArray nc, jintArray kvstr, jintArray kvstc, jintArray ib, jintArray jb, jintArray kb, jdoubleArray b, jint job, jintArray iwk, jint nkmax, jint nzmax) {

	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_kvstr = (*env)->GetPrimitiveArrayCritical(env, kvstr, NULL);
	jint *c_kvstc = (*env)->GetPrimitiveArrayCritical(env, kvstc, NULL);
	jint *c_ib = (*env)->GetPrimitiveArrayCritical(env, ib, NULL);
	jint *c_jb = (*env)->GetPrimitiveArrayCritical(env, jb, NULL);
	jint *c_kb = (*env)->GetPrimitiveArrayCritical(env, kb, NULL);
	jdouble *c_b = (*env)->GetPrimitiveArrayCritical(env, b, NULL);
	jint *c_iwk = (*env)->GetPrimitiveArrayCritical(env, iwk, NULL);
	jint *c_nr = (*env)->GetPrimitiveArrayCritical(env, nr, NULL);
	jint *c_nc = (*env)->GetPrimitiveArrayCritical(env, nc, NULL);

	jint ierr;

	csrvbr_(&n, c_ia, c_ja, c_a, (int*)c_nr, (int*)c_nc, c_kvstr, c_kvstc, c_ib, c_jb, c_kb, c_b, &job, c_iwk, &nkmax, &nzmax, &ierr);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, kvstr, c_kvstr, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, kvstc, c_kvstc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ib, c_ib, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jb, c_jb, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, kb, c_kb, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, b, c_b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iwk, c_iwk, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, nr, c_nr, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, nc, c_nc, 0);

	return ierr;
}

JNIEXPORT jint JNICALL Java_edu_tuberlin_sparskit_Sparskit_vbrcsr(JNIEnv *env, jobject thisObj, jintArray ia, jintArray ja, jdoubleArray a, jint nr, jintArray kvstr, jintArray kvstc, jintArray ib, jintArray jb, jintArray kb, jdoubleArray b, jint nzmax) {

	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_kvstr = (*env)->GetPrimitiveArrayCritical(env, kvstr, NULL);
	jint *c_kvstc = (*env)->GetPrimitiveArrayCritical(env, kvstc, NULL);
	jint *c_ib = (*env)->GetPrimitiveArrayCritical(env, ib, NULL);
	jint *c_jb = (*env)->GetPrimitiveArrayCritical(env, jb, NULL);
	jint *c_kb = (*env)->GetPrimitiveArrayCritical(env, kb, NULL);
	jdouble *c_b = (*env)->GetPrimitiveArrayCritical(env, b, NULL);

	jint ierr;

	vbrcsr_(c_ia, c_ja, c_a, &nr, c_kvstr, c_kvstc, c_ib, c_jb, c_kb, c_b, &nzmax, &ierr);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, kvstr, c_kvstr, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, kvstc, c_kvstc, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ib, c_ib, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, jb, c_jb, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, kb, c_kb, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, b, c_b, JNI_ABORT);

	return ierr;
}


// ##### statistics #####
JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_bandwidth(JNIEnv *env, jobject thisObj, jint n, jintArray ja, jintArray ia, jintArray ml, jintArray mu, jintArray iband, jdoubleArray bndav) {
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_ml = (*env)->GetPrimitiveArrayCritical(env, ml, NULL);
	jint *c_mu = (*env)->GetPrimitiveArrayCritical(env, mu, NULL);
	jint *c_iband = (*env)->GetPrimitiveArrayCritical(env, iband, NULL);
	jdouble *c_bndav = (*env)->GetPrimitiveArrayCritical(env, bndav, NULL);

	bandwidth_(&n, c_ja, c_ia, (int*)c_ml, (int*)c_mu, (int*)c_iband, (double*)c_bndav);

	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ml, c_ml, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, mu, c_mu, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, iband, c_iband, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, bndav, c_bndav, 0);
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_nonz(JNIEnv *env, jobject thisObj, jint n, jint sym, jintArray ja, jintArray ia, jintArray iao, jintArray nzmaxc, jintArray nzminc, jintArray nzmaxr, jintArray nzminr, jintArray nzcol, jintArray nzrow) {
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);
	jint *c_nzmaxc = (*env)->GetPrimitiveArrayCritical(env, nzmaxc, NULL);
	jint *c_nzminc = (*env)->GetPrimitiveArrayCritical(env, nzminc, NULL);
	jint *c_nzmaxr = (*env)->GetPrimitiveArrayCritical(env, nzmaxr, NULL);
	jint *c_nzminr = (*env)->GetPrimitiveArrayCritical(env, nzminr, NULL);
	jint *c_nzcol = (*env)->GetPrimitiveArrayCritical(env, nzcol, NULL);
	jint *c_nzrow = (*env)->GetPrimitiveArrayCritical(env, nzrow, NULL);

	nonz_(&n, &sym, c_ja, c_ia, c_iao, (int*)c_nzmaxc, (int*)c_nzminc, (int*)c_nzmaxr, (int*)c_nzminr, (int*)c_nzcol, (int*)c_nzrow);

	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, nzmaxc, c_nzmaxc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, nzminc, c_nzminc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, nzmaxr, c_nzmaxr, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, nzminr, c_nzminr, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, nzcol, c_nzcol, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, nzrow, c_nzrow, 0);
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_diagdomi(JNIEnv *env, jobject thisObj, jint n, jint sym, jint valued, jdoubleArray a, jintArray ja, jintArray ia, jdoubleArray ao, jintArray jao, jintArray iao, jdoubleArray ddomc, jdoubleArray ddomr) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jdouble *c_ao = (*env)->GetPrimitiveArrayCritical(env, ao, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);
	jdouble *c_ddomc = (*env)->GetPrimitiveArrayCritical(env, ddomc, NULL);
	jdouble *c_ddomr = (*env)->GetPrimitiveArrayCritical(env, ddomr, NULL);

	diag_domi_(&n, &sym, &valued, c_a, c_ja, c_ia, c_ao, c_jao, c_iao, (double*)c_ddomc, (double*)c_ddomr);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ao, c_ao, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ddomc, c_ddomc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ddomr, c_ddomr, 0);
}

JNIEXPORT jdouble JNICALL Java_edu_tuberlin_sparskit_Sparskit_frobnorm(JNIEnv *env, jobject thisObj, jint n, jint sym, jdoubleArray a, jintArray ja, jintArray ia) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);

	jdouble Fnorm;

	frobnorm_(&n, &sym, c_a, c_ja, c_ia, &Fnorm);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);

	return Fnorm;
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_ansym(JNIEnv *env, jobject thisObj, jint n, jint sym, jdoubleArray a, jintArray ja, jintArray ia, jdoubleArray ao, jintArray jao, jintArray iao, jintArray imatch, jdoubleArray av, jdoubleArray fas, jdoubleArray fan) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jdouble *c_ao = (*env)->GetPrimitiveArrayCritical(env, ao, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);
	jint *c_imatch = (*env)->GetPrimitiveArrayCritical(env, imatch, NULL);
	jdouble *c_av = (*env)->GetPrimitiveArrayCritical(env, av, NULL);
	jdouble *c_fas = (*env)->GetPrimitiveArrayCritical(env, fas, NULL);
	jdouble *c_fan = (*env)->GetPrimitiveArrayCritical(env, fan, NULL);

	ansym_(&n, &sym, c_a, c_ja, c_ia, c_ao, c_jao, c_iao, (int*)c_imatch, (double*)c_av, (double*)c_fas, (double*)c_fan);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ao, c_ao, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, imatch, c_imatch, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, av, c_av, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, fas, c_fas, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, fan, c_fan, 0);
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_distaij(JNIEnv *env, jobject thisObj, jint n, jint nnz, jint sym, jintArray ja, jintArray ia, jdoubleArray dist, jdoubleArray std) {
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jdouble *c_dist = (*env)->GetPrimitiveArrayCritical(env, dist, NULL);
	jdouble *c_std = (*env)->GetPrimitiveArrayCritical(env, std, NULL);

	distaij_(&n, &nnz, &sym, c_ja, c_ia, (double*)c_dist, (double*)c_std);

	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, dist, c_dist, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, std, c_std, 0);
}

JNIEXPORT int JNICALL Java_edu_tuberlin_sparskit_Sparskit_skyline(JNIEnv *env, jobject thisObj, jint n, jint sym, jintArray ja, jintArray ia, jintArray jao, jintArray iao) {
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_jao = (*env)->GetPrimitiveArrayCritical(env, jao, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);

	jint nsky;

	skyline_(&n, &sym, c_ja, c_ia, c_jao, c_iao, &nsky);

	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, jao, c_jao, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, JNI_ABORT);

	return nsky;
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_distdiag(JNIEnv *env, jobject thisObj, jint nrow, jint ncol, jintArray ja, jintArray ia, jintArray dist) {
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_dist = (*env)->GetPrimitiveArrayCritical(env, dist, NULL);

	distdiag_(&nrow, &ncol, c_ja, c_ia, c_dist);

	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, dist, c_dist, 0);
}

JNIEXPORT int JNICALL Java_edu_tuberlin_sparskit_Sparskit_bandpart(JNIEnv *env, jobject thisObj, jint n, jintArray ja, jintArray ia, jintArray dist, jint nper) {
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_dist = (*env)->GetPrimitiveArrayCritical(env, dist, NULL);

	jint band;

	bandpart_(&n, c_ja, c_ia, c_dist, &nper, &band);

	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, dist, c_dist, JNI_ABORT);

	return band;
}

JNIEXPORT int JNICALL Java_edu_tuberlin_sparskit_Sparskit_nimpdiag(JNIEnv *env, jobject thisObj, jint n, jint nnz, jintArray dist, jint ipar1, jintArray ioff, jdoubleArray dcount) {
	jint *c_dist = (*env)->GetPrimitiveArrayCritical(env, dist, NULL);
	jint *c_ioff = (*env)->GetPrimitiveArrayCritical(env, ioff, NULL);
	jdouble *c_dcount = (*env)->GetPrimitiveArrayCritical(env, dcount, NULL);

	jint ndiag;

	n_imp_diag_(&n, &nnz, c_dist, &ipar1, &ndiag, c_ioff, c_dcount);

	(*env)->ReleasePrimitiveArrayCritical(env, dist, c_dist, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ioff, c_ioff, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, dcount, c_dcount, 0);

	return ndiag;
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_nonzlud(JNIEnv *env, jobject thisObj, jint n, jintArray ja, jintArray ia, jintArray nlower, jintArray nupper, jintArray ndiag) {
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_nlower = (*env)->GetPrimitiveArrayCritical(env, nlower, NULL);
	jint *c_nupper = (*env)->GetPrimitiveArrayCritical(env, nupper, NULL);
	jint *c_ndiag = (*env)->GetPrimitiveArrayCritical(env, ndiag, NULL);

	nonz_lud_(&n, c_ja, c_ia, (int*)c_nlower, (int*)c_nupper, (int*)c_ndiag);

	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, nlower, c_nlower, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, nupper, c_nupper, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ndiag, c_ndiag, 0);
}

JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_avnzcol(JNIEnv *env, jobject thisObj, jint n, jintArray ja, jintArray ia, jintArray iao, jint ndiag, jdoubleArray av, jdoubleArray st) {
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);
	jint *c_iao = (*env)->GetPrimitiveArrayCritical(env, iao, NULL);
	jdouble *c_av = (*env)->GetPrimitiveArrayCritical(env, av, NULL);
	jdouble *c_st = (*env)->GetPrimitiveArrayCritical(env, st, NULL);

	avnz_col_(&n, c_ja, c_ia, c_iao, &ndiag, (double*)c_av, (double*)c_st);

	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, iao, c_iao, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, av, c_av, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, st, c_st, 0);
}

// ##### utils #####
JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_dvperm(JNIEnv *env, jobject thisObj, jint n, jdoubleArray x, jintArray perm) {
	jdouble *c_x = (*env)->GetPrimitiveArrayCritical(env, x, NULL);
	jint *c_perm = (*env)->GetPrimitiveArrayCritical(env, perm, NULL);

	dvperm_(&n, c_x, c_perm);

	(*env)->ReleasePrimitiveArrayCritical(env, x, c_x, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, perm, c_perm, JNI_ABORT);
}





JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrsss(JNIEnv *env, jobject thisObj, jint nrow, jdoubleArray a, jintArray ja, jintArray ia, jint sorted, jdoubleArray diag, jdoubleArray al, jintArray jal, jintArray ial, jdoubleArray au) {
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, NULL);
	jint *c_ja = (*env)->GetPrimitiveArrayCritical(env, ja, NULL);
	jint *c_ia = (*env)->GetPrimitiveArrayCritical(env, ia, NULL);

	jdouble *c_diag = (*env)->GetPrimitiveArrayCritical(env, diag, NULL);
	jdouble *c_al = (*env)->GetPrimitiveArrayCritical(env, al, NULL);
	jint *c_jal = (*env)->GetPrimitiveArrayCritical(env, jal, NULL);
	jint *c_ial = (*env)->GetPrimitiveArrayCritical(env, ial, NULL);
	jdouble *c_au = (*env)->GetPrimitiveArrayCritical(env, au, NULL);

	csrsss_(&nrow, c_a, c_ja, c_ia, &sorted, c_diag, c_al, c_jal, c_ial, c_au);

	(*env)->ReleasePrimitiveArrayCritical(env, a, c_a, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, c_ja, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ia, c_ia, JNI_ABORT);

	(*env)->ReleasePrimitiveArrayCritical(env, diag, c_diag, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, al, c_al, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jal, c_jal, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ial, c_ial, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, au, c_au, 0);
}

//JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_amux(JNIEnv *env, jobject thisObj, jint n, jdoubleArray x, jdoubleArray y, jdoubleArray a, jintArray ja, jintArray ia) {
//
//	jdouble *c_x = (*env)->GetDoubleArrayElements(env, x, NULL);
//	jdouble *c_y = (*env)->GetDoubleArrayElements(env, y, NULL);
//	jdouble *c_a = (*env)->GetDoubleArrayElements(env, a, NULL);
//	jint *c_ja = (*env)->GetIntArrayElements(env, ja, NULL);
//	jint *c_ia = (*env)->GetIntArrayElements(env, ia, NULL);
//
//	amux_(&n, c_x, c_y, c_a, c_ja, c_ia);
//
//	(*env)->ReleaseDoubleArrayElements(env, x, c_x, JNI_ABORT);
//	(*env)->ReleaseDoubleArrayElements(env, y, c_y, 0);
//	(*env)->ReleaseDoubleArrayElements(env, a, c_a, JNI_ABORT);
//	(*env)->ReleaseIntArrayElements(env, ja, c_ja, JNI_ABORT);
//	(*env)->ReleaseIntArrayElements(env, ia, c_ia, JNI_ABORT);
//}

//JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrsss(JNIEnv *env, jobject thisObj, jint nrow, jdoubleArray a, jintArray ja, jintArray ia, jint sorted, jdoubleArray diag, jdoubleArray al, jintArray jal, jintArray ial, jdoubleArray au) {
//	jdouble *c_a = (*env)->GetDoubleArrayElements(env, a, NULL);
//	jint *c_ja = (*env)->GetIntArrayElements(env, ja, NULL);
//	jint *c_ia = (*env)->GetIntArrayElements(env, ia, NULL);
//
//	jdouble *c_diag = (*env)->GetDoubleArrayElements(env, diag, NULL);
//	jdouble *c_al = (*env)->GetDoubleArrayElements(env, al, NULL);
//	jint *c_jal = (*env)->GetIntArrayElements(env, jal, NULL);
//	jint *c_ial = (*env)->GetIntArrayElements(env, ial, NULL);
//	jdouble *c_au = (*env)->GetDoubleArrayElements(env, au, NULL);
//
//	csrsss_(&nrow, c_a, c_ja, c_ia, &sorted, c_diag, c_al, c_jal, c_ial, c_au);
//
//	(*env)->ReleaseDoubleArrayElements(env, a, c_a, JNI_ABORT);
//	(*env)->ReleaseIntArrayElements(env, ja, c_ja, JNI_ABORT);
//	(*env)->ReleaseIntArrayElements(env, ia, c_ia, JNI_ABORT);
//
//	(*env)->ReleaseDoubleArrayElements(env, diag, c_diag, 0);
//	(*env)->ReleaseDoubleArrayElements(env, al, c_al, 0);
//	(*env)->ReleaseIntArrayElements(env, jal, c_jal, 0);
//	(*env)->ReleaseIntArrayElements(env, ial, c_ial, 0);
//	(*env)->ReleaseDoubleArrayElements(env, au, c_au, 0);
//}

//JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_csrdns(JNIEnv *env, jobject thisObj, jint nrow, jint ncol, jdoubleArray a, jintArray ja, jintArray ia, jobjectArray dns, jint nrowdns, jint ncoldns, jint ierr) {
//	jdouble *c_a = (*env)->GetDoubleArrayElements(env, a, NULL);
//	jint *c_ja = (*env)->GetIntArrayElements(env, ja, NULL);
//	jint *c_ia = (*env)->GetIntArrayElements(env, ia, NULL);
//
//	jsize c_nrowdns = nrowdns;
//	jsize c_ncoldns = ncoldns;
//	jsize i, j;
//	double c_dns[ncoldns][nrowdns];
//	for(i=0; i<c_ncoldns; i++) {
//		jdoubleArray dns1_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, dns, i));
//		jdouble *dns2_tmp = (*env)->GetDoubleArrayElements(env, dns1_tmp, NULL);
//		for(j=0; j<c_nrowdns; j++) {
//			c_dns[i][j] = dns2_tmp[j];
//		}
//	}
//
//	//transpose because of column-major order in Fortran compared to row-major order in C
////	double c_dnsT[ncoldns][nrowdns];
////	dense_transpose(*c_dns, *c_dnsT, nrowdns, ncoldns);
////	for(int k = 0; k < nrowdns; k++) {
////		for(int l = 0; l < ncoldns; l++) {
////			printf("%d, %d, %f \n", k, l, c_dns[k][l]);
////			printf("%d, %d, %f \n", k, l, c_dnsT[k][l]);
////		}
////	}
//
//	csrdns_(&nrow, &ncol, c_a, c_ja, c_ia, (double*)c_dns, &nrowdns, &ierr);
//
//	//transpose back
////	dense_transpose(*c_dnsT, *c_dns, ncoldns, nrowdns);
//
//
//	(*env)->ReleaseDoubleArrayElements(env, a, c_a, JNI_ABORT);
//	(*env)->ReleaseIntArrayElements(env, ja, c_ja, JNI_ABORT);
//	(*env)->ReleaseIntArrayElements(env, ia, c_ia, JNI_ABORT);
//
//	for(i=0; i<c_ncoldns; i++) {
//		jdoubleArray dns_tmp = (jdoubleArray)((*env)->GetObjectArrayElement(env, dns, i));
//		(*env)->SetDoubleArrayRegion(env, dns_tmp, (jsize)0, c_nrowdns, (jdouble*)c_dns[i]);
//		//(*env)->ReleaseDoubleArrayElements(env, dns_tmp, (jdouble*)c_dns[i], 0);
//		(*env)->SetObjectArrayElement(env, dns, i, dns_tmp);
//	}
//}

//experiment
JNIEXPORT void JNICALL Java_edu_tuberlin_sparskit_Sparskit_exper(JNIEnv *env, jobject thisObj, jdoubleArray a) {
	clock_t begin, end;
	double time_spent;

	begin = clock();
	// do something
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("%f seconds \n", time_spent);

	jboolean isCopy;
	jdouble *c_a = (*env)->GetPrimitiveArrayCritical(env, a, &isCopy);
	if(isCopy) {
		printf("TRUE");
	} else {
		printf("FALSE");
	}
}
