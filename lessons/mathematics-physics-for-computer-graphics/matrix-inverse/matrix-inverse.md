In this lesson, we will show how the inverse of a matrix can be derived using a technique known as the Gauss-Jordan (or reduced row) elimination method. Finding the inverse of a matrix implies a couple of things, starting with the fact that the matrix is invertible in the first place (a matrix is not necessarily invertible). First, to be invertible, a matrix has to be a **square matrix** (it has as many rows as it has columns, for instance, 2x2, 3x3, 4x4, etc.), and also the determinant of the matrix has to be different than zero (to learn about the determinant of a matrix refer to [Geometry: Matrix Operations](/lessons/mathematics-physics-for-computer-graphics/geometry/matrix-operations)). Any matrix with a zero determinant is said to be singular (meaning it is not invertible).

Before we get started, remember that a matrix can be seen as a system of linear equations that can be solved using what we call [**row elementary operations**](https://en.wikipedia.org/wiki/Elementary_matrix#Operations). Row elementary operations have the property to preserve the solution set by the matrix. There are three such operations: 

- **Operation 1**: swapping rows of a matrix.
- **Operation 2**: multiplying each coefficient of a row by a non-zero constant \(k\). Note that we say multiplying here, which implies that this can be a division too. We can multiply coefficients by \(\frac{1}{2}\), which is the same as dividing the coefficients by 2. Similarly, this constant can be positive or negative, 2 or -2.
- **Operation 3**: replacing a row with the sum of itself and the multiple of another row. For instance, we can add the coefficients of row 4 to those of row 1. But we can also, in the process, multiply the coefficients of row 4 by some constant \(k\) as follows:

$$
\begin{array}{l}
M_{11} += M_{41} * k,\\
M_{12} += M_{42} * k,\\
M_{13} += M_{43} * k,\\
M_{14} += M_{44} * k
\end{array}
$$

These are illustrated with the following examples (row switching, row multiplication, and row addition):

![](/images/matrix-inverse/mi-rowop.png)

In the second example, we multiply all the coefficients of row 1 by 1/2 (or divide them by 2) and the coefficients of row 3 by 2. Finally, in the third example, we add the coefficients from row 1 to the coefficients of row 2. 

The idea behind the Gauss-Jordan elimination is to use these elementary row operations to transform the 4x4 matrix on the left inside of the augmented matrix into the identity matrix (we say that M is row-reduced), as shown below. We say that we [augment M by identity](https://en.wikipedia.org/wiki/Augmented_matrix). We obtain the inverse matrix by performing the same row operations to the 4x4 identity matrix on the right inside of the augmented matrix.

The Gauss-Jordan elimination method works as follows:

- We will start with two matrices. We want to invert the matrix \(M\) alongside another matrix that we will set at the start of the process to the identity matrix. An identity matrix is a square matrix with ones on the diagonal and zeros elsewhere (the Matrix class constructor sets the matrix to the identity matrix by default). Our matrix \(M\) is augmented. 

$$
[MI]=
\begin{bmatrix}
M_{11} & M_{12} & M_{13} & M_{14}\\
M_{21} & M_{22} & M_{23} & M_{24}\\
M_{31} & M_{32} & M_{33} & M_{34}\\
M_{41} & M_{42} & M_{43} & M_{44}
\end{bmatrix}
\begin{bmatrix}
1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
$$

- We then have to use row elementary operations to transform the original \(M\) matrix into the identity matrix. This process is called **row reduction**. While we do so, we apply the same operations to the second matrix in parallel. 

- Once the process is finished (if we have been successful in reducing \(M\) to the identity matrix), the second matrix will end up being \(M^{-1}\) our inverted matrix. At this point, you can decide to either return \(M^{-1}\) (the method named `inverse` below) or copy the inverted matrix into \(M\) (the method called `invertIt` below).

$$
[IM^{-1}]=
\begin{bmatrix}
1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
M'_{11} & M'_{12} & M'_{13} & M'_{14}\\
M'_{21} & M'_{22} & M'_{23} & M'_{24}\\
M'_{31} & M'_{32} & M'_{33} & M'_{34}\\
M'_{41} & M'_{42} & M'_{43} & M'_{44}
\end{bmatrix}
$$

```
Matrix44& Matrix44::invertIt()
{
    Matrix44 inv; // start from identity matrix
    ...
    *this = inv;
    return *this;
}

Matrix44 Matrix44::inverse(Matrix44& M) const
{ Matrix44 temp(*this); temp.invertIt(); return temp; }

```

Let's see how this is done.

## Step 1: Set The Row So That The Pivot Is Different Than Zero

The coefficients making the diagonal of the matrix are called the **pivots** of the matrix.

![](/images/matrix-inverse/mi-pivots.png?)

As we mentioned earlier, the goal of the matrix inversion process is to use the row elementary operations to set the pivot of each column to 1 and all the other coefficients to 0 (at the end of this process, we will get the identity matrix). To achieve this, the best is to row-reduce each column one after the other starting from the left. Proceeding from left to right guarantees that as we move from column to column, we do not change the values of the columns which have already been row-reduced.

![](/images/matrix-inverse/mi-columns.gif)

We will first examine the value of the pivot coefficient for the current column (the column being processed). This value has to be different than zero. If it is other than zero, we can directly move to step two. If it is equal to zero, we will need to find another row in the matrix for this column whose value is different than zero and swap these two rows (the first of the permitted row elementary operations). However, more than one row may have values other than zero. So which one should we pick in this case? We will select the row whose absolute value is the highest.

![](/images/matrix-inverse/mi-swap.gif)

In the example on the right, we are looking at the pivot coefficient of the second column. Because the value of that coefficient is zero, we swap row 2 with the row of the second column whose absolute value is the greatest (in this example, row 3).

If we can't find another value for the pivot coefficient (that is, if all the other values in the column are zero), then the matrix can't be inverted and is singular (line 11).

If the current column we are processing has a valid pivot coefficient, we can proceed to the next step.

```
Matrix4<T>& invertIt() 
{ 
    Matrix4<T> inv; 
    for (uint32_t column = 0; column < 4; ++column) { 
        // Swap row in case our pivot point is not working
        if (m[column][column] == 0) { 
            uint32_t big = column; 
            for (uint32_t row = 0; row < 4; ++row) 
                if (fabs(m[row][column]) > fabs(m[big][column])) big = row; 
            // Print this as a singular matrix, return identity?
            if (big == column) fprintf(stderr, "Singular matrix\n"); 
            // Swap rows                               
            else for (uint32_t j = 0; j < 4; ++j) { 
                std::swap(m[column][j], m[big][j]); 
                std::swap(inv[column][j], inv[big][j]); 
            } 
        } 
        ... 
    } 
} 
```

## Step 2: Eliminate All Numbers Under the Diagonal Element

The Gauss-Jordan elimination method at this step processes the column one by one from left to right.

In this step, we will first row-reduce the coefficients (in other words, eliminate the numbers) of each column apart under the column pivot coefficient, which we won't touch for now. In order, \(M_{21}\) (we skip \(M_{11}\) since it's a pivot coefficient), \(M_{31}\), \(M_{41}\), then \(M_{32}\) (we skip \(M_{22}\) since it's a pivot coefficient), \(M_{42}\), and finally  \(M_{43}\). There is no coefficient under \(M_{44}\) so note that we can skip the last column. We are about to eliminate the numbers in green in the below matrix.

$$
\begin{bmatrix}
pivot  & M_{12} & M_{13} & M_{14}\\
\color{green}{M_{21}} & pivot  & M_{23} & M_{24}\\
\color{green}{M_{31}} & \color{green}{M_{32}} & pivot  & M_{34}\\
\color{green}{M_{41}} & \color{green}{M_{42}} & \color{green}{M_{43}} & pivot
\end{bmatrix}
$$

Now, remember that we can only use **row elementary operations**. That means you can only modify one single coefficient through a series of operations on that coefficient by changing all the other coefficients from the same row with the same set of operations. If you don't do so, you do not preserve the solution set by the matrix. That's the condition by which **row-reduction** works. Do whatever you want to your coefficient, as long you process all the coefficients in the row that this coefficient is on with the same set of row elementary operations. You can multiply all the coefficients in a row by the same constant and swap rows, and add the coefficients of a given row to those of another row.

Let's see how this works in practice. Let's imagine we are about to process the first column. Let's begin with the first coefficient \(\textcolor{red}{M_{21}}\) (remember we skip the pivot coefficient, so we ignore \(M_{11}\) for now). Let's call it @@\rJoe@@. Let's imagine that @@\rJoe@@'s value is 2. Now, remember that the goal is to somehow make @@\rJoe@@ 0 using one or a combination of the three possible row elementary operations. To make @@\rJoe@@ 0 (in other words, "eliminate" @@\rJoe@@), we should subtract the value 2. That's pretty elementary, but how? Imagine that the coefficient \(\textcolor{blue}{M_{11}}\), the coefficient above @@\rJoe@@, is equal to 4. Note that coincidentally, \(\textcolor{blue}{M_{11}}\) is also @@\rJoe@@'s column pivot coefficient. @@\rJoe@@ is in column 1, and the pivot coefficient of column 1 is \(\textcolor{blue}{M_{11}}\). 

$$
\begin{bmatrix}
\textcolor{blue}{M_{11}} & M_{12} & M_{13} & M_{14}\\
\textcolor{red}{M_{21}} & M_{22} & M_{23} & M_{24}\\
M_{31} & M_{32} & M_{33} & M_{34}\\
M_{41} & M_{42} & M_{34} & M_{44}
\end{bmatrix}
\;=\;
\begin{bmatrix}
\textcolor{blue}{4} & x & x & x\\
\textcolor{red}{2} & x & x & x\\
x & x & x & x\\
x & x & x & x
\end{bmatrix}
\;=\;
\begin{bmatrix}
\textcolor{blue}{pivot}  & x & x & x\\
\textcolor{red}{Joe} & x & x & x\\
x & x & x & x\\
x & x & x & x
\end{bmatrix}
$$

Operation 3 says that we can add coefficients of any of the matrix rows to our current row coefficients. If we use the coefficient of the row that is above @@\rJoe@@, \(\textcolor{blue}{M_{11}}\) (which we know is equal to 4), we get 6:

$$2 + 4 = 6$$

That's different from what we are looking for. However, remember that operation 2 says that you can also multiply the coefficients of a row by a constant. Let's imagine that this constant, for now, is -1. If we multiply \(\textcolor{blue}{M_{11}}\) by this constant, we get:

$$2 + -1 * 4 = -2$$

Not yet 0, but we are getting a step closer. What if our constant was equal to \(-\frac{1}{2}\) or \(-\frac{2}{4}\)

$$2 - \frac{2}{4}4 = 0$$

Now we canceled out the value of coefficient \(\textcolor{red}{M_{21}}\) and made it 0. Maybe you have noticed that 2 is the value of the coefficient \(\textcolor{red}{M_{21}}\), Joe itself, and 4 is the value of the coefficient \(\textcolor{blue}{M_{11}}\), which happens to be the pivot coefficient of the first column, the column than Joe belongs to.

Remember that the operations we applied to Joe also need to be applied to all the coefficients of the row Joe belongs to. That's the rule if you are to use row-elementary operations. So every coefficient on the row that the coefficient being canceled out belongs to needs to receive the same treatment as the one given to that coefficient (Joe). Like so for our Joe example:

$$
\begin{array}{l}
\textcolor{red}{M_{21}} &-=& k * M_{11},\\
M_{22} &-=& k * M_{12},\\
M_{23} &-=& k * M_{13},\\
M_{24} &-=& k * M_{14}.
\end{array}
$$

With:

$$k = -\frac{\textcolor{red}{M_{21}}}{\textcolor{blue}{M_{11}}}$$

Let's now generalize:

!!!
- For every coefficient \(M_{ij}\) in the matrix under the pivot \(M_{ii}\), where \(i\) is the coefficient's row and \(j\) its column, proceeding in order through the columns first, and then for each column from top to bottom.
- Create a constant \(k\) that's made up of \(M_{ij}\) divided by the pivot coefficient of the row above \(M_{ij}\), \(i-1\): \(M_{(i-1)(i-1)}\). If you process \(M_{21}\) (row 2, column 1), then use the pivot coefficient of the first row \(M_{11}\). If you process \(M_{32}\) (row 3, column 2), then use the pivot coefficient \(M_{22}\) (row 2). Etc.
- Negate that number. Let's call this number \(k\).
- Add the coefficients from row \(j\) multiplied by \(k\) to the coefficients of row \(i\) (\(M_{ij}\)'s row). That's allowed by operation 3.
- \(M_{ij}\) is now equal to 0.
!!!

Keep in mind that we loop through the columns first (we loop through \(j\)) and that for each column, we process the coefficients under the pivot for that column from top to bottom (for each given \(j\), we loop through \(i\)), and skip the column's pivot coefficient (when \(i == j\)).

And the end of this process, all of these coefficients **under the diagonal** should be equal to 0 (except the pivot coefficients, of course).

In C++ code, we get the following:

```
// Go left to right. No need to process the last column since
// there's no coefficient under pivot m[3][3].
for (uint32_t column = 0; column < 3; ++column) {
...
    // Make each row in column = 0. Proceed top to bottom
    for (uint32_t row = column + 1; row < 4; ++row) {
			// that's our constant k
            T constant = m[row][column] / m[column][column]; 
			// process each coefficient on the current row
            for (uint32_t j = 0; j < 4; ++j) {
				m[row][j] -= constant * m[column][j]; 
                inv[row][j] -= constant * inv[column][j]; 
            } 
            // Set the element to 0 for safety
            m[row][column] = 0; 
        } 
    }
}
```

Note how we apply the same operations to the right-hand side of the augmented matrix. This step is called the **forward substitution** step. Our matrix should now look like this:

$$
\begin{bmatrix}
pivot & x & x & x\\
0 & pivot & x & x\\
0 & 0 & pivot & x & x\\
0 & 0 & 0 & pivot
\end{bmatrix}
$$

## Step 3: Scale All The Pivots Coefficients To 1

We can now set the pivot coefficients to 1. To do so, we need to scale the coefficients by their value. Similarly to what we have done in step 2, we will loop over each row, then each column, and divide the current coefficient by the pivot coefficient value for the current row. Let's remember to apply the same operation to the inverted matrix.

```
for (uint32_t row = 0; row < 4; ++row) {
	float divisor = m[row][row];
    for (uint32_t column = 0; column < 4; ++column) { 
        m[row][column] = m[row][column] / divisor;
		inv[row][column] = inv[row][column]  / divisor;
    } 
} 
```

Our matrix should now look like this:

$$
\begin{bmatrix}
1 & x & x & x\\
0 & 1 & x & x\\
0 & 0 & 1 & x\\
0 & 0 & 0 & 1
\end{bmatrix}
$$

## Step 4 (and Final): Eliminate All Numbers Above the Diagonal

We are left with eliminating the numbers above the diagonal. The process here is simpler than in step 2. We will subtract the coefficient itself from the coefficient we wish to eliminate (reduce to 0). As with step 2, we can only do so by using one of the three possible row operations. We will use operation 3: replacing the row with the sum of its elements and another row's elements (multiplied by some constant). For the constant part, we will choose the coefficient we are eliminating (the one we called @@\rJoe@@ before). If this coefficient (@@\rJoe@@) is \(M_{ij}\) where as usual \(i\) is the row index and \(j\) the column index, we will add the coefficients of row \(j\) multiplied by \(-M_{ij}\) to the coefficients of row \(i\). In code, this gives the following:

```
// for each row
for (uint32_t row = 0; row < 4; ++row) {
	// for each coefficient above the diagonal for this row
	for (uint32_t column = row + 1; column < 4; ++column) {
		T constant = m[row][column];
		for (uint32_t k = 0; k < 4; ++k) {
			m[row][k] -= m[column][k] * constant;
			inv[row][k] -= inv[column][k] * constant;
		}
		// in case of a round-off error
		m[row][column] = 0.f;
	}
}

``` 

This step is sometimes called **backward substitution**. Note that we step through the rows and then through the columns this time. This is because we only process the elements of the matrix that are above the diagonal. Either way, you can process the column from left to right (our example) or right to left. The maths work here because we have already set the matrix elements below the diagonal to 0 (which is why we can't apply this method in step 2).

Finally, we set the result of the current matrix with the result of our right inside matrix (`inv` in the code). We have inverted the matrix using the Gauss-Jordan elimination technique.

Here is the complete code for the method:

```
Matrix4 Matrix4::Inverse() const
{
    Matrix4 s;
    Matrix4 t(*this);

    // Forward elimination
    for (uint32_t i = 0; i < 3; i++) {
		
		// Step 1: Choose a pivot
        uint32_t pivot = i;

        float pivotsize = t[i][i];

        if (pivotsize < 0) pivotsize = -pivotsize;

        for (uint32_t j = i + 1; j < 4; j++) {
            float tmp = t[j][i];

            if (tmp < 0) tmp = -tmp;

            if (tmp > pivotsize) {
                pivot = j;
                pivotsize = tmp;
            }
        }

        if (pivotsize == 0) { return Matrix4(); }

        if (pivot != i) {
            for (uint32_t j = 0; j < 4; j++) {
                float tmp;

                tmp = t[i][j];
                t[i][j] = t[pivot][j];
                t[pivot][j] = tmp;

                tmp = s[i][j];
                s[i][j] = s[pivot][j];
                s[pivot][j] = tmp;
            }
        }

		// Step 2: eliminate all the numbers below the diagonal
        for (uint32_t j = i + 1; j < 4; j++) {
            float f = t[j][i] / t[i][i];

            for (uint32_t k = 0; k < 4; k++) {
                t[j][k] -= f * t[i][k];
                s[j][k] -= f * s[i][k];
            }
			// Set the column value to exactly 0 in case
			// numeric roundoff left it a very tiny number
			t[j][i] = 0.f;
        }
    }
	
	// Step 3: set elements along the diagonal to 1.0
	for (uint32_t i = 0; i < 4; i++) {
		float divisor = t[i][i];
		for (uint32_t j = 0; j < 4; j++) {
			t[i][j] = t[i][j] / divisor;
			s[i][j] = s[i][j] / divisor;
		}
		// set the diagonal to 1.0 exactly to avoid
		// possible round-off error
		t[i][i] = 1.f;
	}

	// Step 4: eliminate all the numbers above the diagonal
	for (uint32_t i = 0; i < 3; i++) {
		for (uint32_t j = i + 1; j < 4; j++) {
			float constant = t[i][j];
			for (uint32_t k = 0; k < 4; k++) {
				t[i][k] -= t[j][k] * constant;
				s[i][k] -= s[j][k] * constant;
			}
			t[i][j] = 0.f;	// in case of round-off error
		}
	}

    return s;
}
```

## What's Next?

There is one other popular method to compute the inverse of a matrix. It is based on computing determinants and will be explained in a future revision of this lesson.

Check the lesson in the Intermediate Rendering Section called "Warm Up Lap: Reviewing Projections" in which we look into how various libraries, such as Imath, implement matrix inversion.
