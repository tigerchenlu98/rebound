Naming Convention for REBOUND Version 4
=======================================

|                                    | C                                  | Python                 |
| ---                                | ---                                | ---                    |
| module prefix                      | `reb_`                             |                        |
| Macros, constants                  | `REB_PI`                           | `REB_PI`               | 
| Structs                            | `RebSimulation`, `RebParticle`     |                        | 
| Classes                            |                                    | `TestParticle`         | 
| Functions operating on structs     | `reb_simulation_free()`            |                        | 
| Instance methods                   |                                    | `Simulation.free()`    | 
| Private Functions                  | `_reb_simulation_cache()` (static if possible, then prefix is not needed)  | `Simulation._cache()`  |


Other naming rules
------------------

- If a function operates on a struct, the struct's name appears first: e.g. `reb_simulation_move_to_com()`, **not** `reb_move_simulation_to_com().


