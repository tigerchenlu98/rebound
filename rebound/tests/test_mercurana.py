import rebound
import unittest
import warnings

class TestMercurana(unittest.TestCase):
    def test_encounter(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=1e-3,a=1.14,f=-0.2)
        sim.dt = 0.1
        sim.move_to_com()
        E0 = sim.calculate_energy()
        sim2 = sim.copy()
        sim3 = sim.copy()

        sim.integrator = "whfast"
        sim.integrate(50.)
        E1 = sim.calculate_energy()
        self.assertGreater((abs((E1-E0)/E0)),.1)
    
        sim2.integrator = "mercurana"
        sim2.ri_mercurana.Nmaxshells = 30
        sim2.ri_mercurana.kappa = 0.001
        sim2.integrate(50.)
        E1 = sim2.calculate_energy()
        self.assertGreater(sim2.ri_mercurana.Nmaxshellsused,10)
        self.assertLess((abs((E1-E0)/E0)),1e-4)
   
        sim3.integrator = "mercurana"
        sim3.ri_mercurana.Nmaxshells = 30
        sim3.ri_mercurana.N_dominant = 1
        sim3.ri_mercurana.n0=30
        sim3.ri_mercurana.kappa = 0.001
        sim3.integrate(50.)
        E1 = sim3.calculate_energy()
        self.assertGreater(sim3.ri_mercurana.Nmaxshellsused,4)
        self.assertLess((abs((E1-E0)/E0)),1e-4)
   
    def test_merge(self):
        sim = rebound.Simulation()
        sim.add(m=1,x=-1,r=0.1)
        sim.add(m=1,x=1,r=0.1)
        sim.dt = 0.1

        sim.integrator = "mercurana"
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        sim.ri_mercurana.Nmaxshells = 30
        sim.ri_mercurana.kappa = 0.001
        sim.integrate(10.)
        self.assertEqual(sim.N,1)
        self.assertEqual(sim.ri_mercurana.Nmaxshellsused,5)
        self.assertEqual(sim.particles[0].x,0)
        self.assertEqual(sim.particles[0].y,0)
        self.assertEqual(sim.particles[0].vx,0)
        self.assertEqual(sim.particles[0].vy,0)
        self.assertEqual(sim.particles[0].m,2)
    

    def test_merge_three(self):
        sim = rebound.Simulation()
        sim.add(m=1,x=-1,r=0.1)
        sim.add(m=1,x=1,r=0.1)
        sim.add(m=1,a=100,r=0.1)
        sim.move_to_com()
        sim.track_energy_offset = 1
        E0 = sim.calculate_energy()
        sim.dt = 0.1

        sim.integrator = "mercurana"
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        sim.ri_mercurana.Nmaxshells = 30
        sim.ri_mercurana.kappa = 0.001
        sim.integrate(10.)
        E1 = sim.calculate_energy()
        self.assertLess(abs((E1-E0)/E0),1e-4)
        self.assertEqual(sim.N,2)
        self.assertGreater(sim.ri_mercurana.Nmaxshellsused,4)
        self.assertEqual(sim.particles[0].m,2)
    
    def test_merge_N_dominant(self):
        sim = rebound.Simulation()
        sim.add(m=1,x=-1,r=0.1)
        sim.add(m=1,x=1,r=0.1)
        sim.add(m=1,a=100,r=0.1)
        sim.ri_mercurana.N_dominant = 2
        sim.move_to_com()
        sim.track_energy_offset = 1
        sim.collision_resolve_keep_sorted = 1
        E0 = sim.calculate_energy()
        sim.dt = 0.1

        sim.integrator = "mercurana"
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        sim.ri_mercurana.Nmaxshells = 30
        sim.ri_mercurana.kappa = 0.001
        sim.integrate(10.)
        E1 = sim.calculate_energy()
        self.assertLess(abs((E1-E0)/E0),1e-4)
        self.assertEqual(sim.N,2)
        self.assertEqual(sim.ri_mercurana.N_dominant,1)
        self.assertGreater(sim.ri_mercurana.Nmaxshellsused,4)
        self.assertEqual(sim.particles[0].m,2)
    
    def test_merge_N_active(self):
        sim = rebound.Simulation()
        sim.add(m=1,x=-1,r=0.1)
        sim.add(m=1,x=1,r=0.1)
        sim.add(m=0,a=100,r=0.1)
        sim.N_active = 2
        sim.move_to_com()
        sim.track_energy_offset = 1
        sim.collision_resolve_keep_sorted = 1
        E0 = sim.calculate_energy()
        sim.dt = 0.1

        sim.integrator = "mercurana"
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        sim.ri_mercurana.Nmaxshells = 30
        sim.ri_mercurana.kappa = 0.001
        sim.integrate(10.)
        E1 = sim.calculate_energy()
        self.assertLess(abs((E1-E0)/E0),1e-4)
        self.assertEqual(sim.N,2)
        self.assertEqual(sim.N_active,1)
        self.assertGreater(sim.ri_mercurana.Nmaxshellsused,4)
        self.assertEqual(sim.particles[0].m,2)
    
    def test_restart(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=1e-3,a=1.14,f=-0.2)
        sim.dt = 0.1
        sim.move_to_com()
        sim.integrator = "mercurana"
        sim.ri_mercurana.Nmaxshells = 30
        sim.ri_mercurana.N_dominant = 1
        sim.ri_mercurana.n0=30
        sim.ri_mercurana.kappa = 0.001
        sim.integrate(50.)
        sim.save("test.bin")
        sim.integrate(100.)
        sim2 = rebound.Simulation("test.bin")
        sim2.integrate(100.)
  

        self.assertEqual(sim.ri_mercurana._dcrit[0][0],sim2.ri_mercurana._dcrit[0][0])
        self.assertEqual(sim.ri_mercurana._dcrit[1][1],sim2.ri_mercurana._dcrit[1][1])
        self.assertEqual(sim.t,sim2.t)
        self.assertEqual(sim.particles[0].m,sim2.particles[0].m)
        self.assertEqual(sim.particles[1].vx,sim2.particles[1].vx)
        self.assertEqual(sim.particles[2].x,sim2.particles[2].x)

if __name__ == "__main__":
    unittest.main()
