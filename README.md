In this project, we consider the problem of recovering a sparse vector in the context of a linear regression model. 
We develop an algorithm of primal-dual active set type for a class of nonconvex sparsity-promoting penalties, 
which cover $\ell^0$, bridge, smoothly clipped absolute deviation, capped $\ell^1$ and minimax concavity penalty. 
The solutions to the optimality system are coordinate-wise
minimizers, and under minor conditions, they are also local minimizers. Upon introducing
the dual variable, the active set can be determined from the primal and dual variables. This
relation lends itself to an iterative algorithm of active set type which at each step involves
updating the primal variable only on the active set and then updating the dual variable
explicitly. When combined with a continuation strategy on the regularization parameter,
the primal dual active set method has a global convergence property under the restricted
isometry property. Extensive numerical experiments demonstrate its superior performance in efficiency and accuracy compared with the existing methods.
