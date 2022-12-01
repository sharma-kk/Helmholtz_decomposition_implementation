# Helmholtz_decomposition_implementation
According to the Helmholtz decompostion theorem we can split a vector field into its divergence-free and curl-free components. 
We do this numerically here. 
We take a vector field and calculate its div-free and curl-free components. 

Observations: 
1. If i take 'CG' 1 finite element space instead 2 then the solution is not very accurate. at present the divergence is of the order 10^-6 
comparison to theoretical value of the order 10^-11
Ques: Why normal mesh yields accurate results for only some of the test cases.
Ans: This is because I took inappropriate functions as test cases. u_sol should be tangential at the boundary according to the Helmholtz theorem which in my test cases was not the followed.
Ques: Surprisingly barycentric mesh resolves this issue and solving test cases on this mesh gives very accurate results. Why ? 
Ans: This is because barycentric meshes resolve the div(u) things more accurately compared to other meshes.
