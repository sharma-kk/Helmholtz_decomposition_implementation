# solving like a nonlinear problem
# Using Helmholtz decomposition to separate divergence free term from the atmospheric velocity
# u = u_sol + \nabla q
# \Delta q = nabla.u
# n. \nabla q = n.u
from firedrake import *
import math 

mesh = UnitSquareMesh(100,100)
V = FunctionSpace(mesh, 'CG', 2)
W = VectorFunctionSpace(mesh, 'CG',2)
phi = TestFunction(V)


x,y = SpatialCoordinate(mesh)

sq = conditional(And(And(x> 0.25, x < 0.75), And(y > 0.25, y < 0.75)), 1.0, 0.0)

u = Function(W)
u.assign(project(as_vector([sq, 0]), W))
n = FacetNormal(mesh)

q =  Function(V)
F = -(inner(grad(q), grad(phi)))*dx  - div(u)*phi*dx 

nullspace = VectorSpaceBasis(constant = True)

solve (F == 0, q, nullspace = nullspace) 

U_ns = Function(W)
U_ns.interpolate(grad(q))
U_ns.rename("U_ns_calculated")

u.rename("u_exact")

U_s = Function(W)
U_s.assign(u - U_ns)
U_s.rename("u_s_calculated")


print("divergence of U_s calculated numerically:", assemble(div(U_s)*dx))

outfile = File("./results/test_ex_2.pvd")
outfile.write(u, U_ns, U_s)

