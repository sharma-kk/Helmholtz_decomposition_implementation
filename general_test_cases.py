# solving like a nonlinear problem
# Using Helmholtz decomposition to separate divergence free term from the atmospheric velocity
# u = u_sol + \nabla q
# \Delta q = nabla.u
# n. \nabla q = n.u
from firedrake import *
from functions_barycentric_mesh import*

N = 50
M = UnitSquareMesh(N,N)
bmh = BaryMeshHierarchy(M, 0)
mesh = bmh[-1]

V = FunctionSpace(mesh, 'CG', 2)
W = VectorFunctionSpace(mesh, 'CG',2)
phi = TestFunction(V)


x,y = SpatialCoordinate(mesh)

print("option 1: Test example from literature: u = [-sin(pi*x)*cos(pi*y) + pi*cos(pi*(x+y)), cos(pi*x)*sin(pi*y) + pi*cos(pi*(x+y))] \n")
print("option 2: I created this example for testing: u = [x-y, x+y-1] \n")

opt = input("Please choose a vector field whose div-free part is sought (type 1 or 2):")

u = Function(W)

if opt=="1":
    u.assign(project(as_vector([-sin(pi*x)*cos(pi*y) + pi*cos(pi*(x+y)), cos(pi*x)*sin(pi*y) + pi*cos(pi*(x+y))]), W))
    o_file = "general_test_from_literature.pvd"
elif opt=="2":
    u.assign(project(as_vector([x-y, x+y-1]), W))
    o_file = "general_test_own_example.pvd"
else:
    print("wrong choice ! program won't run successfully !")

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
u_s = Function(W)

if opt=="1":
    u_ns.assign(project(as_vector([pi*cos(pi*(x+y)), pi*cos(pi*(x+y))]),W))
    u_s.assign(project(as_vector([-sin(pi*x)*cos(pi*y), cos(pi*x)*sin(pi*y)]),W))
elif opt=="2":
    u_ns.assign(project(as_vector([x-.5,y-.5]),W))
    u_s.assign(project(as_vector([-y+.5,x-.5]),W))
else:
    print("wrong choice ! program won't run successfully !")


u_ns.rename("u_ns_exact")
u_s.rename("u_s_exact")

U_s = Function(W)
U_s.assign(u - U_ns)
U_s.rename("u_s_calculated")

print("divergence of U_s calculated numerically:", assemble(div(U_s)*dx))
print("divergence of u_s calculated theoretically:", assemble(div(u_s)*dx))

# print("Curl of U_ns calculated numerically:", assemble(norm(curl(U_ns))*dx))
# print("Curl of u_ns calculated theoretically:", assemble(norm(curl(u_ns))*dx))

outfile = File("./results/"+ o_file)
outfile.write(u, U_ns,u_ns, U_s, u_s)

