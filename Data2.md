# Update 0

Due to the large break between now and the previous simulations, in addition to the presence of some gaps in the previous simulations, we will fix more experimental parameters and take a few steps back to investigate the confinement problem.

## Default values

*Model parameters*
- R = `CONST_RLIST(1)` = 2.37633
- r = 0.45 mm
- h1 = 0.30 mm
- mem = 0.99
- theta = 1.20 (implicit *π)

*Simulation parameters*
- topography_type = `'circular_well'`
- damping_type = `'none'`
- corral_type = `'none'`
- droplet_collision_type = `'spring'`
- droplet_collision_k = 0.2

# Week 0

## First
1. Make sure frame times are reasonable for the runs planned

## Goals
1. What is the current state of confinement for the experiment-matching configuration? If you just match the experiment, what happens? Do they all leave?
	76. Set all parameters to defaults using standard numerical resolution values. NOTE: not a batch run.
2. Will increasing spatial/temporal resolution improve confinement?
	77. **TODO FIX**: Default, N = 10, 10x temporal resolution. NOTE: not a batch run.
	78. **TODO**
2. Is the domain size large enough to damp waves outside bath for flat topography + rigid damping configuration? For one? Many?
	- **TODO**: Try simulation at 2x the domain size
	79. Default, N = [1, 10], 2x domain size
3. Compare free walking speed from experiments & simulations using many- and single-droplet speed histograms
	- Experiment histogram has peak around 17 mm/s (see experiment **2025-07-22/114** for more)
	- Current parameters put the value around 22 mm/s but will investigate more this time.
	80. Default, N = 1, sending droplet to the right. Sweeping depth and none vs. circular_well.
4. How effective is wave damping in the shallow region and what does the droplet's wavefield look like as the height of the flat bath (h0 == h1) goes from deep to shallow? What about when just h1 approaches `circular_well` topography?
	- **TODO**: Sweep h1 for a single droplet moving across a flat bath
	77. Hello
	78. Hello
	79. Hello
	80. Hello

# For later
- Sweep droplet collision spring constant ('pressure term') 
- **Investigating the 1.2pi -> 1.3pi phase transition is also a worthy endeavor.**
	- Check out behavior of many low-phase particles with rigid boundary
	- "" but with topography