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