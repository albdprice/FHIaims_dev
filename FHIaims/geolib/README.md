# Library of geometry.in files for performance analyses

Purpose of this small collection of geometry.in files is to have representatives of different domains. The domains can be roughly structured in:

- Number of atoms per unit cell
- Number of basis functions per unit cell
- Number of k-points
- Heavy/light elements
- Organic/Inorganic crystals

Danger of any code optimization is that the optimization adresses only some domains, but may harm other. 

Please note that this library should be follow the idea: as small as possible as big as necessary.

Right now there are two sub-folders: molecules and solids. Each of the directory contains representative-folders with the geometry.in, a control.in, and the output of a `dry_run` for *light*, *intermediate*(if available), and *tight* default settings, which can be used to look up important numerical parameters from those systems to estimate runtime and memory.

If you feel something is missing, please add it to the library, however please note the reasons:

## Solids
- Si: Well, what would be a testset without Si?
- GaAs: And what would be a testset without GaAs?
- LiF: System with only very small number of basis functions per unit cell.
- Urea: Smallest organic crystal with all important "organic elemenets" (HCNO)
- Pentacene: Organic crystal, medium size, with two medium sized molecules per unit cell.
- Cu2BaSnS4: Heavy elements only; large cutoff radii need to see convergence;, medium number of atoms per unit cell
- MAPbI3: hybrid inorganic/organic crystal, but still relatively small. Both light and heavy elements.

## Molecules
- H2: Smallest two-atom system.
- H2O: Actually no specific reason, but just put it here (can be removed).
- Pentacene: Medium sized H-C-only system.

**TODO**:
- The library for molecules is still quite small.
- We have no metalls so far.

