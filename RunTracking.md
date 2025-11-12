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
- droplet_collision_type = `'none'`

# Week 0

## First
1. Make sure frame times are reasonable for the runs planned

## Goals
1. What is the current state of confinement for the experiment-matching configuration? If you just match the experiment, what happens? Do they all leave?
	- 76. **TODO**: Set all parameters to defaults using standard numerical resolution values. NOTE: not a batch run.
2. Will increasing spatial/temporal resolution improve confinement? Is the domain size large enough to damp waves outside bath for flat topography + rigid damping configuration? For one? Many?
	- 81. Default
		res: escape
	- 82. Default, double domain size
		res: blow up
	- 83. Default, double time resolution
		res: no result for 2 days
	- 84. Default, double time resolution, double domain size
		res: blow up
	- 85. Default, double spatial resolution
		res: no result for 2 days
3. Compare free walking speed from experiments & simulations using many- and single-droplet speed histograms
	- Experiment histogram has peak around 17 mm/s (see experiment **2025-07-22/114** for more)
	- Current parameters put the value around 22 mm/s but will investigate more this time.
	- 80. Default, N = 1, sending droplet to the right. Sweeping damper depth and phase.

# For later
- Sweep droplet collision spring constant ('pressure term') 
- **Investigating the 1.2pi -> 1.3pi phase transition is also a worthy endeavor.**
	- Check out behavior of many low-phase particles with rigid boundary
	- "" but with topography



2025-11-08	
Goal: Add extra spatial damping outside the corral before every impact to prevent escape, by simulated_extraDaampingOutside.m file 
RUN 86
	Damping by a Guassian function with damping wavelength = 1.

RUN 87
	Change the Gaussian to Exponentially Decaying type.

RUN 88
	Change the dampling wavelength to 0.2

2025-11-09
Goal: Single Droplet Minimize the outreaching distance by changing parameters 
RUN 89
	Sweep c4 = 0.01:0.05:0.31
	