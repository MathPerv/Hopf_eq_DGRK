# Hopf_eq_DGRK
Solving Hopf Equation using DGRK

# Requirements

* Polynomial.jl
* LegendrePolynomials.jl

# Our problem

$$
\begin{cases}
    \frac{\partial U}{\partial t} + U \frac{\partial U}{\partial x} = 0, T\in [0, 2], x \in [-3, 3];
    \\
    U(x, 0) = \frac{1}{1 + x^2}
\end{cases}
$$