import rebound
print(rebound.__build__)
import numpy as np
def getJ(sim):
    J = np.zeros(Nsmall)
    n = sim.particles[1].n
    for i in range(2,sim.N):
        p = sim.particles[i]
        p0 = p - sim.particles[0]
        p1 = p - sim.particles[1]
        p0d = np.sqrt(p0.x**2+p0.y**2)
        p1d = np.sqrt(p1.x**2+p1.y**2)
        J[p.hash.value] = (p.vx**2 + p.vy**2)/2. - sim.particles[0].m/p0d - sim.particles[1].m/p1d - n*(p.x*p.vy - p.y*p.vx)
    return J
sim = rebound.Simulation()
sim.add(m=1,r=0.0046524726)
sim.add(m=1e-3,a=1, r=0.00046732617)

sim.collision = "direct"
sim.dt = 1e-2*2.*np.pi
sim.ri_mercurana.kappa = 4./3.
sim.ri_mercurana.epsilon = 0.1
sim.ri_mercurana.Nmaxshells = 30
sim.ri_mercurana.N_dominant = 2
sim.integrator = "mercurana"
sim.collision_resolve = "merge"
sim.collision_resolve_keep_sorted = 1
sim.move_to_com()
sim.N_active = sim.N

np.random.seed(9)

Nsmall = 100
for i in range(Nsmall):
    #    m=1e-6,
    sim.add(m=0, a=1,f="uniform",omega="uniform",e=0.249,primary=sim.particles[0],r=4.2587571e-05,hash=i)

J0 = getJ(sim)
sim.integrate(100)
J1 = getJ(sim)
