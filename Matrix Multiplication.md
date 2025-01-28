# Matrix Multiplication

## Basics and Notation

The dot product of two $n$-vectors involves $n$ multiplications and #n# additions. The dot product operation is an "$O(n)$" operation meaning that the amount of work scales linearly with the dimension. 

If $A \in \mathbb{R}^{m \times n}$, $x \in \mathbb{R}^n$ and $y \in \mathbb{R}^m$ the the alogirthm to compute $y = y + Ax$ is a generalized $\textit{saxpy}$ operation referred to as $\textit{gaxpy}$, which often proceeds by updating the components one at a time in a row or column oriented way. This makes these operations $O(mn)$ work, where if each dimension of $A$ is doubled, the amount of arithmetic increases by a factor of $4$. A column oriented gaxpy arrives from if we regard $Ax$ as a linear combination of $A$'s columns. 

The vector space of complex $m$-by-$n$ matrices is designated by $\mathbb{C}^{m \times n}$. While the scaling, addition and multiplication of complex matrices correspond exactly to the real case, the transposition becomes $\textit{conjugate transposition}$

$$ C = A^{H} \implies c_{ij} = \bar{a}_{ji} $$

and the dot product of complex $n$-vectors $x$ and $y$ is prescribed by 

$$s = x^{H}y = \sum_{i=1}^n \bar{x}_i y_i $$

## Structure and Efficiency

A $\textit{band matrix}$ is when say $A \in \mathbb{R}^{m \times n}$ has $\textit{lower bandwidth} \: p$ if $a_{ij} = 0$ whenever $i < j + p$ and $\textit{upper bandwidth} \: q$ if $j > i + q$ implies $a_{ij} = 0$. So, diagonal matrices and a lower and upper bandwidth of 0, upper triangular matrices have a lower bandwidth of 0 but upper bandwidth of $n-1$. Lower triangular matrices have a lower bandwidth of $m-1$ and an upper bandwidth of 0. 

Defining $\textit{flop}$ as a floating point add, subtract, multiply or divide, triangular matrix multiplication requires one-sixth the number of flops as full matrix multiplication:

$$ \sum_{i=1}^n \sum_{j=1}^n 2(j - i + 1) = \sum_{i=1}^n \sum_{j=1}^{n-i+1} 2j \approx \sum_{i=1}^n \frac{2(n-i+1)^2}{2} = \sum_{i=1}^n i^2 \approx \frac{n^3}{3} $$

Suppose $A \in \mathbb{R}^{n \times n}$ has lower bandwidth $p$ and upper bandwidth $q$. Such a matrix can be stored in a $(p + q + 1)$-by-$n$ array $A.\textit{band}$ with the convection that

$$ a_{ij} = A.\textit{band}(i - j + q + 1, j) $$ 

for all $(i,j)$ that fall inside the band

## Diagonal, Symmetry and Permutation Matrices

If $A \in \mathbb{R}^{m \times n}$ then pre-multiplication by $D = \text{diag}(d_1,\cdots, d_m) \in \mathbb{R}^{m \times m}$ scales rows, while post multiplication by $D = \text{diag}(d_1, \cdots, d_m) \in \mathbb{R}^{n \times n}$ scales columns. 

A matrix $A \in \mathbb{R}^{n \times n}$ is $\textit{symmetric}$ if $A^T = A$ and $\textit{skew-symmetric}$ if $A^T = -A$. Similarly a matrix $A \in \mathbb{C}^{n \times n}$ is $\textit{Hermitian}$ if $A^H = A$ and $\textit{skew-Hermitian}$ if $A^H = -A$. 

If the rows of an $n$-by-$n$ identity matrix are reordered, then the resulting matrix is said to be a $\textit{permutation}$ matrix.

## Block Matrices

A Block matrix is a matrix whose entries are themselves matrices. Block matrices can be scaled and transposed

$$ \begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \\ A_{31} & A_{32}\end{bmatrix}^T = \begin{bmatrix} A_{11}^T & A_{21}^T & A_{31}^T \\ \\ A_{12}^T & A_{22}^T & A_{32}^T\end{bmatrix} $$

In addition, block matrices can be added by summing the corresponding blocks and matrix multiplication between block matrices follows a similar scheme with nested matrices instead of values, given that the row and column dimensions of the blocks satisfy all the necessary constraints. 

## Submatrices

Suppose $A \in \mathbb{R}^{m \times n}$. If $\alpha = [\alpha_1, \cdots, \alpha_s]$ and $\beta = [\beta_1, \cdots, \beta_t]$ are integer vectors with distinct components that satisfy $1 \leq \alpha_i \leq m$ and $1 \leq B_i \leq n$, then

$$ A(\alpha, \beta) = \begin{bmatrix} a_{\alpha_1, \beta_1} & \cdots & a_{\alpha_1, \beta_t} \\ \vdots & \ddots & \vdots \\ a_{\alpha_s, \beta_1} & \cdots & a_{\alpha_s, \beta_t} \end{bmatrix} $$

is an $s$-by-$t$ submatrix of A. For example, if $A \in \mathbb{R}^{8 \times 6}$, $\alpha = [2\: 4\: 6\: 8]$ and $\beta = [4 \: 5 \: 6]$ then

$$ A(\alpha, \beta) = \begin{bmatrix} a_{24} & a_{25} & a_{26} \\ a_{44} & a_{45} & a_{46} \\ a_{64} & a_{65} & a_{66} \\ a_{84} & a_{85} & a_{86} \end{bmatrix} $$

If $\alpha = \beta$ then $A(\alpha, \beta)$ is a $\textit{principal submatrix}$.

## The Kronecker product

$A$ is a $\textit{Kronecker product}$ if the entries in a block matrix $A$ are all scalar multiples of the same matrix. So if $B \in \mathbb{R}^{m_1 \times n_1}$ and $C \in \mathbb{R}^{m_2 \times n_2}$, then their Kronecker product $B \otimes C$ is an $m_1$-by-$n_1$ block matrix whose $(i, j)$ block is the $m_2$-by-$n_2$ matrix $b_{ij}C$. Thus if, 

$$ A = \begin{bmatrix} b_{11} & b_{12} \\ b_{21} & b_{22} \\ b_{31} & b_{32} \\ \end{bmatrix} \otimes \begin{bmatrix} c_{11} & c_{12} & c_{13} \\ c_{21} & c_{22} & c_{23} \\ c_{31} & c_{32} & c_{33} \\ \end{bmatrix} $$

Then, 

$$ A = \begin{bmatrix} b_{11}C & b_{12}C \\ b_{21}C & b_{22}C \\ b_{31}C & b_{32}C \end{bmatrix} $$

So if $B$ has a band structure it inherits that structure in $B \otimes C$ in block form.

Other important properties of the Kronecker product are

$$ (B \otimes C)^T = B^T \otimes C^T, \\
(B \otimes C)(D \otimes F) = BD \otimes CF, \\
(B \otimes C)^{-1} = B^{-1} \otimes C^{-1}, \\
B \otimes (C \otimes D) = (B \otimes C) \otimes D $$

In general, $B \otimes C \neq C \otimes B$, but we can use the perfect shuffle permutation matrix to reoder the rows and columns and expressive commutativity. The perfect shuffle permuation matrix for $B \in \mathbb{R}^{m_1 \times n_1}$ and $C \in \mathbb{R}^{m_2 \times n_2}$, gives us

$$ P_{m_1, m_2} (B \otimes C)P_{n_1, n_2}^T = C \otimes B $$

where $P_{m_1, m_2}$ can be defined explicitly using indices

1. Let $I = \{0, 1, \cdots, m_1-1\}$ and $J = \{0, 1, \cdots, m_2-1\}$
2. For each pair $(i,j)$ in $B \otimes C$ maps to the position $(j,i)$ in $C \otimes B$

Then we can define 

$$ P_{m,n}[k, l] = \begin{cases} 1, \: \text{if} \: k=m_2(i) + j \: \text{and} \: l = m_1(j) + i, \\ 0, \: \text{otherwise}  \end{cases} $$

We can also reshape Kroncker products to matrix-matrix-matrix product expressions. If $X \in \mathbb{R}^{m \times n}$ then $\text{vec}(X)$ is an $nm$-by-$1$ vector obtained by "stacking" $X$'s columns i.e. $\text{vec}(X) = \begin{bmatrix} X(:, 1) \\ \vdots \\ X(:,n) \end{bmatrix}$

If $B \in \mathbb{R}^{m_1 \times n_1}$, $C \in \mathbb{R}^{m_2 \times n_2}$ and $X \in \mathbb{R}^{n_1 \times m_2}$, then

$$ Y = CXB^T \iff \text{vec}(Y) = (B \otimes C)\text{vec}(X) $$

where $CXB^T$ costs $O(n^3)$ to evaluate while the Kronecker structure product would give an $O(n^4)$ calculation. 

## Hamiltonian and Symplectic Matrices

A matrix $M \in \mathbb{R}^{2n \times 2n}$ is a $\textit{Hamiltonian matrix}$ if it has the form

$$ M = \begin{bmatrix} A & G \\ F & -A^T \end{bmatrix} $$

where $A, F, G \in \mathbb{R}^{n \times n}$ and $F$ and $G$ are symmetric. An equivalent definition can be given in terms of the permutation matrix

$$ J = \begin{bmatrix} 0 & I_n \\ -I_n & 0 \end{bmatrix} $$

where a Hamiltonian matrix $M$ satisfies $JMJ^T = -M^T$ or equivalently it saitsfies

$$ M^TJ + JM = 0 $$

We can show this by computing $M^T$

$$ M^T = \begin{bmatrix}A^T & F^T \\ G^T & -A \end{bmatrix} $$

and since $F = F^T$ and $G = G^T$ we get 

$$ M^T = \begin{bmatrix} A^T & F \\ G & -A \end{bmatrix} $$

thus,

$$ M^TJ = \begin{bmatrix} A^T & F \\ G & -A \end{bmatrix} \begin{bmatrix} 0 & I_n \\ -I_n & 0 \end{bmatrix} = \begin{bmatrix} -F & A^T \\ A & G \end{bmatrix}$$

We also get that

$$ JM = \begin{bmatrix} 0 & I_n \\ -I_n & 0 \end{bmatrix} \begin{bmatrix} A & G \\ F & -A^T \end{bmatrix} = \begin{bmatrix} F & -A^T \\ -A & -G\end{bmatrix} $$

which proves that $M^TJ + JM = 0$

a related class of matrices are the symplectic matrices where a matrix $S \in \mathbb{R}^{2n \times 2n}$ is symplectic if

$$ S^TJS = J $$

If 

$$S = \begin{bmatrix} S_{11} & S_{12} \\ S_{21} & S_{22} \end{bmatrix} $$

where the blocks are $n$-by-$n$, then it follows that both $S_{11}^TS_{21}$ and $S_{22}^TS_{21}$ are symmetric and $S_{11}^TS_22 = I_n + S_{21}^TS_{12}$. 

## Strassen Matrix Multiplication

Let $A$ and $B$ be $2 \times 2$ block matrices:
$$
A = \begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}, \quad
B = \begin{bmatrix} B_{11} & B_{12} \\ B_{21} & B_{22} \end{bmatrix},
$$
where each block (e.g., $A_{11}, B_{12}$) is itself an $n \times n$ matrix.

The product $C = AB$ is given as:
$$
C = \begin{bmatrix} C_{11} & C_{12} \\ C_{21} & C_{22} \end{bmatrix},
$$
where:
$$
C_{11} = A_{11}B_{11} + A_{12}B_{21}, \quad
C_{12} = A_{11}B_{12} + A_{12}B_{22},
$$
$$
C_{21} = A_{21}B_{11} + A_{22}B_{21}, \quad
C_{22} = A_{21}B_{12} + A_{22}B_{22}.
$$

### Strassen's Algorithm
Strassen's algorithm reduces the number of multiplications from 8 to 7 by introducing 7 intermediate products, $P_1, P_2, \dots, P_7$, as follows:
$$
P_1 = (A_{11} + A_{22})(B_{11} + B_{22}),
$$
$$
P_2 = (A_{21} + A_{22})B_{11},
$$
$$
P_3 = A_{11}(B_{12} - B_{22}),
$$
$$
P_4 = A_{22}(B_{21} - B_{11}),
$$
$$
P_5 = (A_{11} + A_{12})B_{22},
$$
$$
P_6 = (A_{21} - A_{11})(B_{11} + B_{12}),
$$
$$
P_7 = (A_{12} - A_{22})(B_{21} + B_{22}).
$$

The blocks of $C$ are then computed as:
$$
C_{11} = P_1 + P_4 - P_5 + P_7,
$$
$$
C_{12} = P_3 + P_5,
$$
$$
C_{21} = P_2 + P_4,
$$
$$
C_{22} = P_1 - P_2 + P_3 + P_6.
$$

### Verification of Strassen's Algorithm
We verify that the computed $C_{11}, C_{12}, C_{21}, C_{22}$ match the standard definitions.

#### $C_{11}$:
$$
C_{11} = P_1 + P_4 - P_5 + P_7.
$$
Substituting the definitions of $P_1, P_4, P_5, P_7$:
$$
P_1 = (A_{11} + A_{22})(B_{11} + B_{22}), \quad
P_4 = A_{22}(B_{21} - B_{11}),
$$
$$
P_5 = (A_{11} + A_{12})B_{22}, \quad
P_7 = (A_{12} - A_{22})(B_{21} + B_{22}).
$$
Expanding:
$$
C_{11} = (A_{11}B_{11} + A_{11}B_{22} + A_{22}B_{11} + A_{22}B_{22})
+ (A_{22}B_{21} - A_{22}B_{11})
- (A_{11}B_{22} + A_{12}B_{22})
+ (A_{12}B_{21} + A_{12}B_{22} - A_{22}B_{21} - A_{22}B_{22}).
$$
Combine like terms:
$$
C_{11} = A_{11}B_{11} + A_{12}B_{21},
$$
which matches the standard definition.

#### $C_{12}$:
$$
C_{12} = P_3 + P_5.
$$
Substituting $P_3$ and $P_5$:
$$
P_3 = A_{11}(B_{12} - B_{22}), \quad
P_5 = (A_{11} + A_{12})B_{22}.
$$
Expanding:
$$
C_{12} = A_{11}B_{12} - A_{11}B_{22} + A_{11}B_{22} + A_{12}B_{22}.
$$
Simplify:
$$
C_{12} = A_{11}B_{12} + A_{12}B_{22},
$$
which matches the standard definition.

#### $C_{21}$:
$$
C_{21} = P_2 + P_4.
$$
Substituting $P_2$ and $P_4$:
$$
P_2 = (A_{21} + A_{22})B_{11}, \quad
P_4 = A_{22}(B_{21} - B_{11}).
$$
Expanding:
$$
C_{21} = A_{21}B_{11} + A_{22}B_{11} + A_{22}B_{21} - A_{22}B_{11}.
$$
Simplify:
$$
C_{21} = A_{21}B_{11} + A_{22}B_{21},
$$
which matches the standard definition.

#### $C_{22}$:
$$
C_{22} = P_1 - P_2 + P_3 + P_6.
$$
Substituting $P_1, P_2, P_3, P_6$:
$$
P_1 = (A_{11} + A_{22})(B_{11} + B_{22}), \quad
P_2 = (A_{21} + A_{22})B_{11},
$$
$$
P_3 = A_{11}(B_{12} - B_{22}), \quad
P_6 = (A_{21} - A_{11})(B_{11} + B_{12}).
$$
Expanding:
$$
C_{22} = (A_{11}B_{11} + A_{11}B_{22} + A_{22}B_{11} + A_{22}B_{22})
- (A_{21}B_{11} + A_{22}B_{11})
+ (A_{11}B_{12} - A_{11}B_{22})
+ (A_{21}B_{11} + A_{21}B_{12} - A_{11}B_{11} - A_{11}B_{12}).
$$
Combine like terms:
$$
C_{22} = A_{21}B_{12} + A_{22}B_{22},
$$
which matches the standard definition.

Suppose $n = 2m$, counting adds and multiplies with conventional matrix multiplication involves $(2m)^3$ multiplies and $(2m)^3-(2m)^2$ adds. Strassen's algorithm multiplies $7m^3$ and $7m^3 + 11m^2$ adds. If $m >> 1$ then the Strassen method involves about $7/8$ the arithmetic which works much better

Additionally, we can recur the Strassen idea, i.e. we can apply the Strassen algorithm to each of the half-sized block multiplications associated with the $P_i$