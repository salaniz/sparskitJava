package edu.tuberlin.sparskit.formats;
/*
 * Copyright (C) 2003-2006 Bjørn-Ove Heimsund
 * 
 * This file is part of MTJ.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */



import no.uib.cipr.matrix.*;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.io.MatrixInfo;
import no.uib.cipr.matrix.io.MatrixSize;
import no.uib.cipr.matrix.io.MatrixVectorReader;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.*;

import edu.tuberlin.sparskit.Sparskit;
import edu.tuberlin.sparskit.formats.utils.Utils;
import edu.tuberlin.sparskit.formats.utils.Utils.CsrDiaJob;
import edu.tuberlin.sparskit.formats.utils.Utils.CsrEllJob;

/**
 * Compressed row storage (CRS) matrix
 */
public class AdaptedCompRowMatrix extends AbstractMatrix {

    ByteBuffer byteBuffer;

    /**
     * Matrix data
     */
    DoubleBuffer data;

    /**
     * Column indices. These are kept sorted within each row.
     */
    IntBuffer columnIndex;

    /**
     * Indices to the start of each row
     */
    IntBuffer rowPointer;

    public AdaptedCompRowMatrix() {
        super(0,0);
    }

    /**
     * Constructor for CompRowMatrix
     * 
     * @param r
     *            Reader to get sparse matrix from
     */
    public AdaptedCompRowMatrix(MatrixVectorReader r) throws IOException {
        // Start with a zero-sized matrix
        super(0, 0);

        // Get matrix information. Use the header if present, else just assume
        // that the matrix stores real numbers without any symmetry
        MatrixInfo info = null;
        if (r.hasInfo())
            info = r.readMatrixInfo();
        else
            info = new MatrixInfo(true, MatrixInfo.MatrixField.Real,
                    MatrixInfo.MatrixSymmetry.General);

        // Check that the matrix is in an acceptable format
        if (info.isPattern())
            throw new UnsupportedOperationException(
                    "Pattern matrices are not supported");
        if (info.isDense())
            throw new UnsupportedOperationException(
                    "Dense matrices are not supported");
        if (info.isComplex())
            throw new UnsupportedOperationException(
                    "Complex matrices are not supported");

        // Resize the matrix to correct size
        MatrixSize size = r.readMatrixSize(info);
        numRows = size.numRows();
        numColumns = size.numColumns();

        // Start reading entries
        int numEntries = size.numEntries();
        int[] row = new int[numEntries];
        int[] column = new int[numEntries];
        double[] entry = new double[numEntries];
        r.readCoordinate(row, column, entry);

        // Shift the indices from 1 based to 0 based
        r.add(-1, row);
        r.add(-1, column);

        // Find the number of entries on each row
        List<Set<Integer>> rnz = new ArrayList<Set<Integer>>(numRows);
        for (int i = 0; i < numRows; ++i)
            rnz.add(new HashSet<Integer>());

        for (int i = 0; i < numEntries; ++i)
            rnz.get(row[i]).add(column[i]);

        // Allocate some more in case of symmetry
        if (info.isSymmetric() || info.isSkewSymmetric())
            for (int i = 0; i < numEntries; ++i)
                if (row[i] != column[i])
                    rnz.get(column[i]).add(row[i]);

        int[][] nz = new int[numRows][];
        for (int i = 0; i < numRows; ++i) {
            nz[i] = new int[rnz.get(i).size()];
            int j = 0;
            for (Integer colind : rnz.get(i))
                nz[i][j++] = colind;
        }

        // Create the sparse matrix structure
        construct(nz);

        // Insert the entries
        for (int i = 0; i < size.numEntries(); ++i)
            set(row[i], column[i], entry[i]);

        // Put in extra entries from symmetry or skew symmetry
        if (info.isSymmetric())
            for (int i = 0; i < numEntries; ++i) {
                if (row[i] != column[i])
                    set(column[i], row[i], entry[i]);
            }
        else if (info.isSkewSymmetric())
            for (int i = 0; i < numEntries; ++i) {
                if (row[i] != column[i])
                    set(column[i], row[i], -entry[i]);
            }
    }

    public AdaptedCompRowMatrix(int numRows, int numColumns, double[] data, int[] columnIndex, int[] rowPointer) {
        super(numRows, numColumns);
        this.data = DoubleBuffer.wrap(data);
        this.columnIndex = IntBuffer.wrap(columnIndex);
        this.rowPointer = IntBuffer.wrap(rowPointer);
    }

    /**
     * Constructor for CompRowMatrix
     *
     * @param numRows
     *            Number of rows
     * @param numColumns
     *            Number of columns
     * @param nz
     *            The nonzero column indices on each row
     */
    public AdaptedCompRowMatrix(int numRows, int numColumns, int[][] nz) {
        super(numRows, numColumns);
        construct(nz);
    }

    /**
     * Constructor for CompRowMatrix
     *
     * @param numRows
     *            Number of rows
     * @param numColumns
     *            Number of columns
     * @param nz
     *            The nonzero column indices on each row
     */
    public AdaptedCompRowMatrix(int numRows, int numColumns, int[][] nz, double[] data) {
        super(numRows, numColumns);
        construct(nz, data);
    }

    /**
     * Constructor for CompRowMatrix
     *
     * @param A
     *            Copies from this matrix
     * @param deep
     *            True if the copy is to be deep. If it is a shallow copy,
     *            <code>A</code> must be a <code>CompRowMatrix</code>
     */
    public AdaptedCompRowMatrix(Matrix A, boolean deep) {
        super(A);
        construct(A, deep);
    }

    /**
     * Constructor for CompRowMatrix
     *
     * @param A
     *            Copies from this matrix. The copy will be deep
     */
    public AdaptedCompRowMatrix(Matrix A) {
        this(A, true);
    }

    private void construct(int[][] nz) {
        int nnz = 0;
        for (int i = 0; i < nz.length; ++i)
            nnz += nz[i].length;

        int[] rowPointer = new int[numRows + 1];
        int[] columnIndex = new int[nnz];
        data = DoubleBuffer.wrap(new double[nnz]);

        if (nz.length != numRows)
            throw new IllegalArgumentException("nz.length != numRows");

        for (int i = 1; i <= numRows; ++i) {
            rowPointer[i] = rowPointer[i - 1] + nz[i - 1].length;

            for (int j = rowPointer[i - 1], k = 0; j < rowPointer[i]; ++j, ++k) {
                columnIndex[j] = nz[i - 1][k];
                if (nz[i - 1][k] < 0 || nz[i - 1][k] >= numColumns)
                    throw new IllegalArgumentException("nz[" + (i - 1) + "]["
                            + k + "]=" + nz[i - 1][k]
                            + ", which is not a valid column index");
            }

            Arrays.sort(columnIndex, rowPointer[i - 1], rowPointer[i]);
        }

        this.rowPointer = IntBuffer.wrap(rowPointer);
        this.columnIndex = IntBuffer.wrap(columnIndex);
    }

    private void construct(int[][] nz, double[] data) {
        int nnz = 0;
        for (int i = 0; i < nz.length; ++i)
            nnz += nz[i].length;

        int[] rowPointer = new int[numRows + 1];
        int[] columnIndex = new int[nnz];
        this.data = DoubleBuffer.wrap(data);

        if (nz.length != numRows)
            throw new IllegalArgumentException("nz.length != numRows");

        for (int i = 1; i <= numRows; ++i) {
            rowPointer[i] = rowPointer[i - 1] + nz[i - 1].length;

            for (int j = rowPointer[i - 1], k = 0; j < rowPointer[i]; ++j, ++k) {
                columnIndex[j] = nz[i - 1][k];
                if (nz[i - 1][k] < 0 || nz[i - 1][k] >= numColumns)
                    throw new IllegalArgumentException("nz[" + (i - 1) + "]["
                            + k + "]=" + nz[i - 1][k]
                            + ", which is not a valid column index");
            }

            Arrays.sort(columnIndex, rowPointer[i - 1], rowPointer[i]);
        }
        this.rowPointer = IntBuffer.wrap(rowPointer);
        this.columnIndex = IntBuffer.wrap(columnIndex);
    }

    private void construct(Matrix A, boolean deep) {
        if (deep) {
            if (A instanceof AdaptedCompRowMatrix) {
                AdaptedCompRowMatrix Ac = (AdaptedCompRowMatrix) A;
                data = Ac.data.duplicate();//new double[Ac.data.capacity()];

                this.columnIndex = Ac.columnIndex.duplicate();
                this.rowPointer = Ac.rowPointer.duplicate();;
            } else {

                List<Set<Integer>> rnz = new ArrayList<>(numRows);
                for (int i = 0; i < numRows; ++i)
                    rnz.add(new HashSet<Integer>());

                for (MatrixEntry e : A)
                    rnz.get(e.row()).add(e.column());

                int[][] nz = new int[numRows][];
                for (int i = 0; i < numRows; ++i) {
                    nz[i] = new int[rnz.get(i).size()];
                    int j = 0;
                    for (Integer colind : rnz.get(i))
                        nz[i][j++] = colind;
                }

                construct(nz);
                set(A);

            }
        } else {
            AdaptedCompRowMatrix Ac = (AdaptedCompRowMatrix) A;
            columnIndex = Ac.getColumnIndices();
            rowPointer = Ac.getRowPointers();
            data = Ac.data;//getData();
        }
    }

    /**
     * Returns the column indices
     */
    public IntBuffer getColumnIndices() {
        return columnIndex;
    }

    /**
     * Returns the row pointers
     */
    public IntBuffer getRowPointers() {
        return rowPointer;
    }

    /**
     * Returns the internal data storage
     */
   /* public double[] getData() {
        return data;
    }*/

    @Override
    public Matrix mult(Matrix B, Matrix C) {
        checkMultAdd(B, C);
        C.zero();

        // optimised a little bit to avoid zeros in rows, but not to
        // exploit sparsity of matrix B
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < C.numColumns(); ++j) {
                double dot = 0;
                for (int k = rowPointer.array()[i]; k < rowPointer.array()[i + 1]; ++k) {
                    dot += data.get(k) * B.get(columnIndex.array()[k], j);
                }
                if (dot != 0) {
                    C.set(i, j, dot);
                }
            }
        }
        return C;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.multAdd(alpha, x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData();
        double[] yd = ((DenseVector) y).getData();

        for (int i = 0; i < numRows; ++i) {
            double dot = 0;
            for (int j = rowPointer.array()[i]; j < rowPointer.array()[i + 1]; ++j)
                dot += data.get(j) * xd[columnIndex.array()[j]];
            yd[i] += alpha * dot;
        }

        return y;
    }

    @Override
    public Vector transMult(Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.transMult(x, y);

        checkTransMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData();
        double[] yd = ((DenseVector) y).getData();

        y.zero();

        for (int i = 0; i < numRows; ++i)
            for (int j = rowPointer.array()[i]; j < rowPointer.array()[i + 1]; ++j)
                yd[columnIndex.array()[j]] += data.get(j) * xd[i];

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.transMultAdd(alpha, x, y);

        checkTransMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData();
        double[] yd = ((DenseVector) y).getData();

        // y = 1/alpha * y
        y.scale(1. / alpha);

        // y = A'x + y
        for (int i = 0; i < numRows; ++i)
            for (int j = rowPointer.array()[i]; j < rowPointer.array()[i + 1]; ++j)
                yd[columnIndex.array()[j]] += data.get(j) * xd[i];

        // y = alpha*y = alpha*A'x + y
        return y.scale(alpha);
    }

    @Override
    public void set(int row, int column, double value) {
        check(row, column);

        int index = getIndex(row, column);

        data.put(index,value);
    }

    @Override
    public void add(int row, int column, double value) {
        check(row, column);

        int index = getIndex(row, column);
        throw new UnsupportedOperationException("Sorry");

        //data[index] += value;
    }

    @Override
    public double get(int row, int column) {
        check(row, column);

        int index = Arrays.binarySearch(columnIndex.array(),
                rowPointer.array()[row], rowPointer.array()[row + 1], column);

        if (index >= 0)
            return data.get(index);
        else
            return 0;
    }

    /**
     * Finds the insertion index
     */
    private int getIndex(int row, int column) {
        int i = Arrays.binarySearch(columnIndex.array(), rowPointer.array()[row], rowPointer.array()[row + 1], column);

        if (i != -1 && columnIndex.array()[i] == column)
            return i;
        else
            throw new IndexOutOfBoundsException("Entry (" + (row + 1) + ", "
                    + (column + 1) + ") is not in the matrix structure");
    }

    @Override
    public AdaptedCompRowMatrix copy() {

        throw new UnsupportedOperationException("Sorry");

        //return new AdaptedCompRowMatrix(this);
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new CompRowMatrixIterator();
    }

    @Override
    public AdaptedCompRowMatrix zero() {
        data = DoubleBuffer.wrap(new double[data.limit()]);
        return this;
    }

    @Override
    public Matrix set(Matrix B) {

        if (!(B instanceof AdaptedCompRowMatrix))
            return super.set(B);

        checkSize(B);

        AdaptedCompRowMatrix Bc = (AdaptedCompRowMatrix) B;

        // Reallocate matrix structure, if necessary
        if (Bc.columnIndex.limit() != columnIndex.limit()
                || Bc.rowPointer.limit() != rowPointer.limit()) {
            data = DoubleBuffer.wrap(new double[Bc.data.capacity()]);
            columnIndex = IntBuffer.wrap(new int[Bc.columnIndex.limit()]);
            rowPointer = IntBuffer.wrap(new int[Bc.rowPointer.limit()]);
        }

        System.arraycopy(Bc.data.array(), 0, data.array(), 0, data.capacity());
        System.arraycopy(Bc.columnIndex.array(), 0, columnIndex.array(), 0, columnIndex.limit());
        System.arraycopy(Bc.rowPointer.array(), 0, rowPointer.array(), 0, rowPointer.limit());

        return this;
    }

    public IntBuffer getRowPointer() {
        return rowPointer;
    }

    /**
     * Iterator over a compressed row matrix
     */
    private class CompRowMatrixIterator implements Iterator<MatrixEntry> {

        private int row, cursor;

        private CompRowMatrixEntry entry = new CompRowMatrixEntry();

        public CompRowMatrixIterator() {
            // Find first non-empty row
            nextNonEmptyRow();
        }

        /**
         * Locates the first non-empty row, starting at the current. After the
         * new row has been found, the cursor is also updated
         */
        private void nextNonEmptyRow() {
            while (row < numRows() && rowPointer.get(row) == rowPointer.get(row + 1))
                row++;
            cursor = rowPointer.get(row);
        }

        public boolean hasNext() {
            return cursor < data.capacity();
        }

        public MatrixEntry next() {
            entry.update(row, cursor);

            // Next position is in the same row
            if (cursor < rowPointer.get(row + 1) - 1)
                cursor++;

            // Next position is at the following (non-empty) row
            else {
                row++;
                nextNonEmptyRow();
            }

            return entry;
        }

        public void remove() {
            entry.set(0);
        }

    }

    /**
     * Entry of a compressed row matrix
     */
    private class CompRowMatrixEntry implements MatrixEntry {

        private int row, cursor;

        /**
         * Updates the entry
         */
        public void update(int row, int cursor) {
            this.row = row;
            this.cursor = cursor;
        }

        public int row() {
            return row;
        }

        public int column() {
            return columnIndex.get(cursor);
        }

        public double get() {
            return data.get(cursor);
        }

        public void set(double value) {
            throw new UnsupportedOperationException("Sorry");

//            data[cursor] = value;
        }
    }
    
    // STEPHAN: NEW CODE FOR SPARSKIT INTEGRATION
    
    @Override
    public Vector mult(Vector x, Vector y) {
    	//CHANGED TO SPARSKIT IMPLEMENTATION (old implementation in comments)
    	
//        // check dimensions
//        checkMultAdd(x, y);
//        // can't assume this, unfortunately
//        y.zero();
//
//        if (x instanceof DenseVector) {
//            // DenseVector optimisations
//            double[] xd = ((DenseVector) x).getData();
//            for (int i = 0; i < numRows; ++i) {
//                double dot = 0;
//                for (int j = rowPointer.array()[i]; j < rowPointer.array()[i + 1]; j++) {
//                    dot += data.get(j) * xd[columnIndex.array()[j]];
//                }
//                if (dot != 0) {
//                    y.set(i, dot);
//                }
//            }
//            return y;
//        }
//        // use sparsity of matrix (not vector), as get(,) is slow
//        // TODO: additional optimisations for mult(ISparseVector, Vector)
//        // note that this would require Sparse BLAS, e.g. BLAS_DUSDOT(,,,,)
//        // @see http://www.netlib.org/blas/blast-forum/chapter3.pdf
//        for (int i = 0; i < numRows; ++i) {
//            double dot = 0;
//            for (int j = rowPointer.get(i); j < rowPointer.get(i + 1); j++) {
//                dot += data.get(j) * x.get(columnIndex.get(j));
//            }
//            y.set(i, dot);
//        }
//        return y;
    	
    	//assuming the indices start from 1
    	//otherwise incrAllIndices() has to be called - TODO: make option
		checkMultAdd(x, y);
    	DenseVector yDense = Utils.getDenseVector(y);
    	Sparskit.amux(numRows, Utils.getDenseVector(x).getData(), yDense.getData(), data.array(), columnIndex.array(), rowPointer.array());
    	return yDense;
    }
    
    
    public AdaptedCompRowMatrix(double[][] denseMatrix, int numberNonZeros) throws Exception {
    	super(denseMatrix[0].length, denseMatrix.length);
    	double[] data = new double[numberNonZeros];
    	int[] columnIndex = new int[numberNonZeros];
    	int[] rowPointer = new int[numRows + 1];
    	int error = Sparskit.dnscsr(numRows, numColumns, numberNonZeros, denseMatrix, numRows, data, columnIndex, rowPointer);
    	if(error != 0) {
    		throw new Exception("Sparskit Dense to CSR transformation failed. Code stopped while processing row number " + error + ", because there was no space left in CSR data structures");
    	}
    	this.data = DoubleBuffer.wrap(data);
        this.columnIndex = IntBuffer.wrap(columnIndex);
        this.rowPointer = IntBuffer.wrap(rowPointer);
    }
    
    public void incrAllIndices() {
    	int[] csrColumnIndex = columnIndex.array();
    	int[] csrRowPointer = rowPointer.array();
        for(int i = 0; i < csrColumnIndex.length; i++) {
        	if(i < csrRowPointer.length) {
        		csrRowPointer[i]++;
        	}
        	csrColumnIndex[i]++;
        }
        if(csrRowPointer.length > csrColumnIndex.length) {
        	for(int i = csrColumnIndex.length; i < csrRowPointer.length; i++) {
        		csrRowPointer[i]++;
        	}
        }
    }
    
    public long getMemorySizeInBytes() {
    	return data.capacity() * Double.BYTES + columnIndex.capacity() * Integer.BYTES + rowPointer.capacity() * Integer.BYTES;
    }
    
    //use this when matrix is in 0 base and you want to use Sparskit multiplication
    public Vector multSparskit(Vector x, Vector y) {
    	//TODO: make option for "add 1 to each index":
    	incrAllIndices(); //add 1 to each index
    	
		checkMultAdd(x, y);
    	DenseVector yDense = Utils.getDenseVector(y);
    	Sparskit.amux(numRows, Utils.getDenseVector(x).getData(), yDense.getData(), data.array(), columnIndex.array(), rowPointer.array());
    	return yDense;
    }
    
    public double[][] toDense() throws Exception {
    	double[][] denseMatrix = new double[numColumns][numRows];
    	int error = Sparskit.csrdns(numRows, numColumns, data.array(), columnIndex.array(), rowPointer.array(), denseMatrix, numRows);
    	if(error != 0) {
    		throw new Exception("Sparskit CSR to Dense transformation failed. Code has stopped when processing row number " + error + ", because it found a column number greater than number of columns (" + numColumns + ")");
    	}
    	return denseMatrix;
    	
    }
    
    public AdaptedCompRowMatrix getTranspose() {
    	double[] data = new double[this.data.capacity()];
    	int[] rowIndex = new int[this.data.capacity()];
    	int[] columnPointer = new int[numColumns + 1];
    	Sparskit.csrcsc2(numRows, numColumns, 1, 1, this.data.array(), columnIndex.array(), rowPointer.array(), data, rowIndex, columnPointer);
    	return new AdaptedCompRowMatrix(numColumns, numRows, data, rowIndex, columnPointer);
    }
    
    public ModifiedCompRowMatrix toModifiedCompRowMatrix() {
    	double[] mData = new double[numRows + data.capacity() + 1];
    	int[] mColumnIndex = new int[numRows + data.capacity() + 1];
    	Sparskit.csrmsr(numRows, data.array(), columnIndex.array(), rowPointer.array(), mData, mColumnIndex, mData, mColumnIndex);
    	return new ModifiedCompRowMatrix(numRows, numColumns, mData, mColumnIndex);
    }
    
    public DiagonalMatrix toDiagonalMatrix(Utils.CsrDiaJob job, int numDiagonals, int sizeRemainder) {
        int n = numRows;
    	double[][] diag = new double[numDiagonals][n];
        int[] diagonalOffsets = new int[numDiagonals];
        double[] remainderData = new double[sizeRemainder];
        int[] remainderColumnIndex = new int[sizeRemainder];
        int[] remainderRowPointer = new int[numRows + 1];
        int[] workArray;
        if(job == Utils.CsrDiaJob.SELECT_DIAGONALS_WITH_REMAINDER_DATA || job == Utils.CsrDiaJob.SELECT_DIAGONALS_WITHOUT_REMAINDER_DATA) {
        	workArray = new int[(n)*2-1];
        } else {
        	workArray = new int[1];
        }
    	Sparskit.csrdia(n, numDiagonals, job.getValue(), data.array(), columnIndex.array(), rowPointer.array(), n, diag, diagonalOffsets, remainderData, remainderColumnIndex, remainderRowPointer, workArray);
    	return new DiagonalMatrix(numRows, numColumns, diag, diagonalOffsets, remainderData, remainderColumnIndex, remainderRowPointer);
    }
    
    public JaggedDiagonalMatrix toJaggedDiagonalMatrix() {
        double[] jdiag = new double[data.capacity()];
        int[] jColumnIndex = new int[data.capacity()];
        int[] iperm = new int[numRows];
        int[] jdPointer = new int[numColumns];
        int numJaggedDiagonals = Sparskit.csrjad(numRows, data.array(), columnIndex.array(), rowPointer.array(), iperm, jdiag, jColumnIndex, jdPointer);
    	return new JaggedDiagonalMatrix(numRows, numColumns, jdiag, jColumnIndex, iperm, jdPointer, numJaggedDiagonals);
    }
    
    // not working when data gets bigger, reason unresolved. Sparskit.csrell gets called and nothing happens, not a single line in the C interface is executed
    public EllpackMatrix toEllpackMatrix() throws Exception {
    	int ellColDimension = (int)(data.capacity()/numRows * 15);
    	double[][] ellData = new double[ellColDimension][numRows];
    	int[][] ellColumnIndex = new int[ellColDimension][numRows];
    	int[] ierr = new int[1];
    	int numDiags = Sparskit.csrell(numRows, data.array(), columnIndex.array(), rowPointer.array(), ellColDimension, ellData, ellColumnIndex, numRows, ierr);
    	if(ierr[0] == 1) {
    		throw new Exception("Sparskit CSR to Ellpack transformation failed. Init arrays not big enough for ellpack convertion");
    	}
    	
    	double[][] ellData2 = null;
    	int[][] ellColumnIndex2 = null;
    	for(int i = 0; i < ellData.length; i++) {
        	double sum = 0D;
        	for(int j = 0; j < numRows; j++) {
        		sum += ellData[i][j];
        	}
        	if(sum == 0D) {
        		ellData2 = new double[i][numRows];
            	ellColumnIndex2 = new int[i][numRows];
        		for(int k = 0; k < i; k++) {
        			for(int l = 0; l < numRows; l++) {
        				ellData2[k][l] = ellData[k][l];
        				ellColumnIndex2[k][l] = ellColumnIndex[k][l];
                	}
        		}
        		break;
        	}
    	}
    	if(ellData2 == null){
    		ellData2 = ellData;
    		ellColumnIndex2 = ellColumnIndex;
    	}
    	return new EllpackMatrix(numRows, numColumns, ellData2, ellColumnIndex2, numDiags);
    }
    
    public VbrMatrix toVbrMatrix(Utils.CsrEllJob job) throws Exception {
        double[] vbrData = new double[numRows * numColumns];
        int[] vbrIndex1 = new int[vbrData.length];
        int[] vbrIndex2 = new int[vbrData.length];
        int[] vbrIndex3 = new int[vbrData.length];
        // nr,nc   = matrix block row and block column dimension
        int[] nr = new int[1];
        int[] nc = new int[1];
        int[] kvstr = new int[numRows];
        int[] kvstc = new int[numColumns];
        int ierr = Sparskit.csrvbr(numRows, columnIndex.array(), rowPointer.array(), data.array(), nr, nc, kvstr, kvstc, vbrIndex1, vbrIndex2, vbrIndex3, vbrData, job.getValue(), new int[vbrData.length], vbrIndex1.length, vbrData.length);
        if(ierr == 1) {
    		throw new Exception("Sparskit CSR to VBR transformation failed. Out of space in jb and/or kb arrays");
        } else if (ierr == 2) {
    		throw new Exception("Sparskit CSR to VBR transformation failed. Out of space in b array");
        } else if (ierr == 3) {
    		throw new Exception("Sparskit CSR to VBR transformation failed. Nonsquare matrix used with job=2");
    	}
        return new VbrMatrix(numRows, numColumns, vbrData, vbrIndex1, vbrIndex2, vbrIndex3, kvstr, kvstc, nr[0], nc[0]);
    }
}