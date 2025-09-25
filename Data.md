# 2025-08-06/18
Recalculating specific Faraday threshold (0.004950_0.000100_2.376330) using higher-res (256x256) grid

# 2025-08-07/19
Redoing the same threshold (0.004950_0.000100_2.376330) but with 10 gaussians instead of the above which used a Bessel. We find the threshold is lower, 5.60 vs 5.74.

# 2025-08-07/20
We looked at 2025-08-04/17 when we swept h1 and radius for a fixed phase 1.2 (which TODO we now know is too low). There was only one h1 >= 0.2 mm run that blew up - the cavity mode one (R = 2.63) - so we will pick that one and sweep droplet number.

# 2025-08-20/48, 49, 50, 51
Using well topography, spring corral and setting N = 1. Sweeping phase, corral radius, and memory

# 2025-08-20/52
Using flat topography, rigid corral and N = 10. 99% memory. Sweeping phase and corral radius.

# 2025-08-20/53
Using well topography, spring corral. Sweeping N and memory. and setting N = 1.

# 2025-08-20/54
Using well topography, spring corral and setting N = 10. Sweeping phase and memory.

# 2025-08-20/55
Using well topography, no corral and setting N = 10. Sweeping phase and memory.

# 2025-08-20/56
Using well topography, no corral and setting N = 1. Sweeping phase and memory.

# 2025-08-20/57
Using well topography, no corral and setting N = 10 at r = 0.45 mm. Sweeping phase and memory.

# 2025-08-20/58
Using well topography, no corral and setting N = 1 at r = 0.45 mm. Sweeping phase and memory.

# 2025-08-20/59
Special branch - testing arrested bouncers to check wavefield slope

# 2025-08-20/60
Using well topography, no corral and setting N = 10 at r = 0.36 mm. Sweeping h1 and memory. We intentionally picked unusually high h1 values to see if we can convince the slow 1.3pi droplets to leave the middle.

# 2025-08-24/61
Copy of 56 - sending r = 0.36 mm droplet to top boundary and plotting video. Now trying for doubled time resolution.

# 2025-08-24/62
Copy of 58 - sending r = 0.45 mm droplet to top boundary and plotting video. Now trying for doubled time resolution.

# 2025-08-24/63
Free walking + collisions. r = 0.45 mm, 98% memory, fixed corral size, phase. No sweeping.

# 2025-08-24/64
Free walking + collisions. r = 0.36 mm, 98% memory, fixed corral size, phase. No sweeping.

# 2025-08-26/65
Complement to 43. Sweeping memory and phase for free walker but at r = 0.45 mm.

# 2025-08-25/73
After a chat with another working on the code, decided to revisit the very-low-damper case and consider what happens if we accept a high Faraday threshold. Sent off a run for 0.1 mm, sweeping 90 to 99 memory and radius

# 2025-08-25/74
Same as above but tried 0.05 mm and for a larger range of memories beginning at 80%

# 2025-08-25/76
Different goal - while the above run, revisiting artificially-damped droplet motion and reran the simulations but this time shifted the effective corral size in very slightly to 98%.