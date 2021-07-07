from firedrake import *
import numpy
import random
try:
  import matplotlib.pyplot as plt
except:
  warning("Matplotlib not imported")

#== Create mesh == 
nx, ny = 100, 40
Lx, Ly = 1.0, 0.4
mesh = RectangleMesh(nx,ny,Lx,Ly)

#== Function spaces == 
#---Double porosity/permeability flow problem---
velSpace = VectorFunctionSpace(mesh,"DG",2)
pSpace = FunctionSpace(mesh,"DG",2)
wSpace = MixedFunctionSpace([velSpace,pSpace,velSpace,pSpace])

#---Advection-diffusion problem---
uSpace = FunctionSpace(mesh,"CG",1)

#---Permeability---
kSpace = FunctionSpace(mesh, "DG", 0)

#== Material properties and parameters ==
mu0, Rc, D = Constant(1e-3), Constant(3.0), Constant(2e-6)
k1_0 = 1.1
k1_1 = 0.9
tol = 1E-14

class myk1(Expression): #Macro-permeability
    def eval(self, values, x):
        if x[1] < Ly/2 + tol:
            values[0] = k1_0
        else: 
            values[0] = k1_1

k1 = interpolate(myk1(),kSpace)

k2_0 = 0.01 * 1.1
k2_1 = 0.01 * 0.9

class myk2(Expression): #Micro-permeability
    def eval(self, values, x):
        if x[1] < Ly/2 + tol:
            values[0] = k2_0
        else: 
            values[0] = k2_1

k2 = interpolate(myk2(),kSpace)

#---Drag coefficients---

def alpha1(c): 
    return mu0 * exp(Rc * (1.0 - c))/k1

def invalpha1(c): 
    return 1/alpha1(c)

def alpha2(c): 
    return mu0 * exp(Rc * (1.0 - c))/k2

def invalpha2(c): 
    return 1/alpha2(c)

#== Boundary and initial conditions ==
v_topbottom = Constant(0.0)
p_L = Constant(10.0) 
p_R = Constant(1.0) 
c_inj = Constant(1.0)

#== Perturbation function for initial concentration ==
#---Needed to trigger the instability---
class c_0(Expression):
    def eval(self, values, x):
        if x[0] < 0.010*Lx:
            values[0] = abs(.10*exp(-x[0]*x[0]) * random.random())
        else: 
            values[0] = 0.0

#== Define trial and test functions ==  
#---DPP flow problem---
(v1,p1,v2,p2) = TrialFunctions(wSpace)
(w1,q1,w2,q2) = TestFunctions(wSpace)
DPP_solution = Function(wSpace) 

#---AD problem---
c1 = TrialFunction(uSpace)
u = TestFunction(uSpace)
conc = Function(uSpace) 
conc_k = interpolate(c_0(),uSpace)

#== Time parameters ==
T = 0.0015	# Total simulation time
dt = 0.00005	# Time step

#== Boundary conditions ==
#---DPP velocity BCs---
bcDPP = []

#---AD concentration BCs---  
bcleft_c = DirichletBC(uSpace,c_inj,1,method = "geometric")

bcAD = [bcleft_c]

#== Define source terms ==
#---DPP model---
rhob1, rhob2 = Constant((0.0,0.0)), Constant((0.0,0.0))

#---AD problem---
f = Constant(0.0)

#== Normal vectors and mesh size ==
n = FacetNormal(mesh)
h = CellSize(mesh)
h_avg = (h('+') + h('-'))/2

#== Penalty parameters ==
eta_p, eta_u = Constant(0.0), Constant(0.0)

#== Define variational forms ==  

#---DPP stabilized mixed DG formulation---
aDPP = dot(w1, alpha1(conc_k) * v1) * dx +\
      dot(w2, alpha2(conc_k) * v2) * dx -\
      div(w1) * p1 * dx -\
      div(w2) * p2 * dx +\
      q1 * div(v1) * dx +\
      q2 * div(v2) * dx +\
      q1 * (p1 - p2) * dx -\
      q2 * (p1 - p2) * dx +\
      jump(w1,n) * avg(p1) * dS +\
      jump(w2,n) * avg(p2) * dS -\
      avg(q1) * jump(v1,n) * dS -\
      avg(q2) * jump(v2,n) * dS +\
      dot(w1,n) * p1 * ds(3) +\
      dot(w2,n) * p2 * ds(3) -\
      q1 * dot(v1,n) * ds(3) -\
      q2 * dot(v2,n) * ds(3) +\
      dot(w1,n) * p1 * ds(4) +\
      dot(w2,n) * p2 * ds(4) -\
      q1 * dot(v1,n) * ds(4) -\
      q2 * dot(v2,n) * ds(4) -\
      0.5 * dot( alpha1(conc_k) * w1 - grad(q1), \
                 invalpha1(conc_k) * (alpha1(conc_k) * v1 + grad(p1)) ) * dx -\
      0.5 * dot( alpha2(conc_k) * w2 - grad(q2), \
                 invalpha2(conc_k) * (alpha2(conc_k) * v2 + grad(p2)) ) * dx +\
      (eta_u * h_avg)  * avg(alpha1(conc_k)) * (jump(v1,n) * jump(w1,n)) * dS +\
      (eta_u * h_avg)  * avg(alpha2(conc_k)) * (jump(v2,n) * jump(w2,n)) * dS +\
      (eta_p / h_avg)  * avg(1 / alpha1(conc_k)) * dot(jump(q1,n),jump(p1,n)) * dS +\
      (eta_p / h_avg)  * avg(1 / alpha2(conc_k)) * dot(jump(q2,n),jump(p2,n)) * dS   

LDPP = dot(w1,rhob1) * dx +\
      dot(w2,rhob2) * dx -\
      dot(w1,n) * p_L * ds(1) -\
      dot(w2,n) * p_L * ds(1) -\
      dot(w1,n) * p_R * ds(2) -\
      dot(w2,n) * p_R * ds(2) -\
      0.5 * dot( alpha1(conc_k) * w1 - grad(q1), \
                 invalpha1(conc_k) * rhob1 ) * dx -\
      0.5 * dot( alpha2(conc_k) * w2 - grad(q2), \
                 invalpha2(conc_k) * rhob2 ) * dx 


#---AD formulation with SUPG Stabilization---
vnorm = sqrt(dot((DPP_solution.sub(0)+DPP_solution.sub(2)),\
                 (DPP_solution.sub(0)+DPP_solution.sub(2))))

taw = h/(2*vnorm)*dot((DPP_solution.sub(0)+DPP_solution.sub(2)),\
                      grad(u))

a_r = taw*(c1 + dt*(dot((DPP_solution.sub(0)+DPP_solution.sub(2)),\
                        grad(c1)) - div(D*grad(c1))))*dx

L_r = taw*(conc_k + dt*f)*dx

#---Weak form (GL + SUPG)---
aAD = a_r + u*c1*dx + dt*(u*dot((DPP_solution.sub(0)+DPP_solution.sub(2)),\
                                grad(c1))*dx + dot(grad(u),D*grad(c1))*dx)

LAD = L_r + u*conc_k*dx + dt*u*f*dx

#---Create files for storing solution---
cfile = File("Concentration.pvd")
v1file = File("Macro_Velocity.pvd")
p1file = File("Macro_Pressure.pvd")
v2file = File("Micro_Velocity.pvd")
p2file = File("Micro_Pressure.pvd")

#== Solver for flow problem ==
solver_parameters = { # Default solver -- medium sized problems
  'ksp_type': 'gmres',
  'pc_type': 'bjacobi',
  'mat_type': 'aij',
  'ksp_rtol': 1e-7,
  'ksp_monitor': True
}

problem_flow = LinearVariationalProblem(aDPP, LDPP, DPP_solution, bcs=bcDPP,
  constant_jacobian=False)
solver_flow = LinearVariationalSolver(problem_flow, options_prefix="flow_",
  solver_parameters=solver_parameters)

#== March the solution over time  ==
t = dt
while t <= T:
    print '=============================='
    print '          time =', t
    print '=============================='
    c_0.t = t
    
    #---Compute DPP model---
    solver_flow.solve()

    #---Compute AD problem---
    solve(aAD == LAD,conc,bcs=bcAD)	
    conc_k.assign(conc)   # update for next iteration
    
    #---Dump solutions for each time step---
    cfile.write(conc, time = t)
    v1file.write(DPP_solution.sub(0), time = t)
    p1file.write(DPP_solution.sub(1), time = t)
    v2file.write(DPP_solution.sub(2), time = t)
    p2file.write(DPP_solution.sub(3), time = t)
    t += dt

print "total time = ", t

v1sol, p1sol, v2sol, p2sol = DPP_solution.split()

#==  Dump solution fields to file in VTK format ==
file = File("Concentration.pvd")
file.write(conc)

file = File('Macro_Velocity.pvd')
file.write(v1sol)

file = File('Macro_Pressure.pvd')
file.write(p1sol)

file = File('Micro_Velocity.pvd')
file.write(v2sol)

file = File('Micro_Pressure.pvd')
file.write(p2soll)

