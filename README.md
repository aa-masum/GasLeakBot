# GasLeakBot
Autonomous gas-leak localization robot with LiDAR occupancy mapping (MATLAB simulation)

A simulated differential-drive robot explores a room with walls, builds an occupancy-grid map from a 72-beam LiDAR, and uses five gas sensors to follow the concentration gradient to a hidden leak while avoiding obstacles.

**Features**
- 72-beam LiDAR simulated by ray marching on a grid world
- Occupancy-grid mapping with a log-odds update
- Five-sensor gas array (centre, front, back, left, right) with Gaussian sensor noise
- Gas-gradient estimation from sensor differences, smoothed with an exponential filter
- Gas-seeking navigation: gas gradient attracts, LiDAR hits repel (potential field)
- Obstacle-aware gas model: the plume spreads around walls (geodesic distance), not through them
- Live three-panel view: true world, LiDAR occupancy map, gas map
- Built-in metrics: localization error, path length, minimum clearance, area explored, map accuracy, obstacle IoU


--How it works:** 
- sense  ->  map  ->  estimate gradient  ->  decide  ->  move   (repeats every dt = 0.1 s)
