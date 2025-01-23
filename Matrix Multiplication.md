# Chapter 1
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