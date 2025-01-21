# Chapter 1
# Matrix Multiplication

The dot product of two $n$-vectors involves $n$ multiplications and #n# additions. The dot product operation is an "$O(n)$" operation meaning that the amount of work scales linearly with the dimension. 

If $A \in \mathbb{R}^{m \times n}$, $x \in \mathbb{R}^n$ and $y \in \mathbb{R}^m$ the the alogirthm to compute $y = y + Ax$ is a generalized $\textit{saxpy}$ operation referred to as $\textit{gaxpy}$, which often proceeds by updating the components one at a time in a row or column oriented way. This makes these operations $O(mn)$ work, where if each dimension of $A$ is doubled, the amount of arithmetic increases by a factor of $4$. A column oriented gaxpy arrives from if we regard $Ax$ as a linear combination of $A$'s columns. 