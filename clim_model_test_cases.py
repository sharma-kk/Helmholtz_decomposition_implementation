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

print("option 1: Wind blowing over the ocean from left to right in square patch of the domain \n")
print("option 2: Circular velocity field of the atmosphere like a cyclone rotating anti-clockwise \n")
print("option 3: explosion, the air is moving outwards")
opt = input("Pleace choose from the above options (type 1, 2, or 3) for atm vel field: \n")

circ = conditional(sqrt(pow(x-0.5, 2) + pow(y-0.5,2)) < 0.25, 1.0, 0.0)
sq = conditional(And(And(x> 0.25, x < 0.75), And(y > 0.25, y < 0.75)), 1.0, 0.0)

u = Function(W)
if opt=="1":
    u.assign(project(as_vector([sq, 0]), W))
    o_file = "test_case_wind.pvd"
elif opt=="2":
    u.assign(project(as_vector([-circ*(y-0.5), circ*(x-0.5)]), W))
    o_file = "test_case_cyclone.pvd"
elif opt == "3":
    u.assign(project(as_vector([circ*(x-0.5), circ*(y-0.5)]), W))
    o_file = "test_case_explosion.pvd"
else:
    print("wrong choice ! program won't run successfully !")

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

outfile = File("./results/" + o_file)
outfile.write(u, U_ns, U_s)

