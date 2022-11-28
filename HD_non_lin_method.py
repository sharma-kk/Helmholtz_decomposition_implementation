# solving like a nonlinear problem
# Using Helmholtz decomposition to separate divergence free term from the atmospheric velocity
# u = u_sol + \nabla q
# \Delta q = nabla.u
# n. \nabla q = n.u
from firedrake import *

mesh = UnitSquareMesh(100,100)
V = FunctionSpace(mesh, 'CG', 2)
W = VectorFunctionSpace(mesh, 'CG',2)
phi = TestFunction(V)


x,y = SpatialCoordinate(mesh)

u = Function(W)
u.assign(project(as_vector([-sin(pi*x)*cos(pi*y) + pi*cos(pi*(x+y)), cos(pi*x)*sin(pi*y) + pi*cos(pi*(x+y))]), W))
n = FacetNormal(mesh)

q =  Function(V)
F = -(inner(grad(q), grad(phi)))*dx  - div(u)*phi*dx + dot(n, u)*phi*ds

nullspace = VectorSpaceBasis(constant = True)

solve (F == 0, q, nullspace = nullspace) 

U_ns = Function(W)
U_ns.interpolate(grad(q))
U_ns.rename("U_ns_calculated")

u.rename("u_exact")

u_ns = Function(W)
u_ns.assign(project(as_vector([pi*cos(pi*(x+y)), pi*cos(pi*(x+y))]),W))
u_ns.rename("u_ns_exact")

U_s = Function(W)
U_s.assign(u - U_ns)
U_s.rename("u_s_calculated")

u_s = Function(W)
u_s.assign(project(as_vector([-sin(pi*x)*cos(pi*y), cos(pi*x)*sin(pi*y)]),W))
u_s.rename("u_s_exact")

print("divergence of U_s calculated numerically:", assemble(div(U_s)*dx))
print("divergence of u_s calculated theoretically:", assemble(div(u_s)*dx))
outfile = File("./results/HD_nonlin.pvd")
outfile.write(u, U_ns,u_ns, U_s, u_s)
#Observations: If i take 'CG' 1 finite element space instead 2 then the solution is not very accurate. at present the divergence is of the order 10^-6 
# comparison to theoretical value of the order 10^-11
