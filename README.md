# Helmholtz_decomposition_implementation
According to the Helmholtz decompostion theorem we can split a vector field into its divergence-free and curl-free components. 
We do this numerically here. 
We take a vector field and calculate its div-free and curl-free components. 

Observations: 
1. If i take 'CG' 1 finite element space instead 2 then the solution is not very accurate. at present the divergence is of the order 10^-6 
comparison to theoretical value of the order 10^-11
2. Normal mesh yields accurate results for only some of the test cases.
3. Surprisingly barycentric mesh resolves this issue and solving test cases on this mesh gives very accurate results. Why ? 
4. Apparantly not ! I used peridic mesh that's why it was giving accurate results. The improvement from normal mesh is not that much. I should think of some other approach or maybe
the example that I am solving has some problems. 