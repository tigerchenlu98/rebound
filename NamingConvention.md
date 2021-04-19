Naming Convention for REBOUND Version 4
=======================================

|                                    | C                                  | Python                 |
| ---                                | ---                                | ---                    |
| module prefix                      | `reb_`                             |                        |
| Macros, constants                  | `REB_PI`                           | `REB_PI`               | 
| Structs                            | `RebSimulation`, `RebParticle`     |                        | 
| Variables                          | `particle_number`                  | `particle_number`      |
| Classes                            |                                    | `TestParticle`         | 
| Functions operating on structs     | `reb_simulation_free()`            |                        | 
| Instance methods                   |                                    | `Simulation.free()`    | 
| Private Functions                  | `_reb_simulation_cache()` (static if possible, then prefix is not needed)  | `Simulation._cache()`  |


Other naming rules
------------------

- If a function operates on a struct, the struct's name appears first: e.g. `reb_simulation_move_to_com()`, **not** `reb_move_simulation_to_com().
- Following GLib/GTK: `_new()` allocated and initializes a structure.
- `_free()` deallocates structures.
- `_init()` initializes a structure but does not allocate memory for the structure itself.
- `_destroy()`  uninitializes a sturcture, but does not free memory for the structure itself.
- `_copy()` allocates memory and performs a deep copy.

Particle types
--------------

- Massive particles
- Small particles (old testparticle type 1). Small particles only see massive particle, massive particles see small particles. But small particles see no other small particles.
- Test particles (old testparticle type 0). Test particles only see massive particles. Test particles do not see small particles.


