# NTL-Asteroid-Tracker/

The program is executed calling from this directory:

`java -jar AsteroidTrackerVisualizer.jar`

Testcases can be changed using option `-testcase`, while output file can be specified using option `-output`. The default test case is `testcase1`. The default output is `output.txt` in this directory.

Visualization can be switched off using option `-no_visualization`. 

The following constants are hardcoded:
```
MAX_ASTEROIDS = 75;
MAX_ANTENNAS = 20;
TURN_LENGTH = 180;
TOTAL_TIME = 7 days
```

Configuration parameters can be changed in config file. However, the algorithm has been optimized for the default values of the parameters.

Build jar locally using:
```
cd executable/AsteroidTracker
ant
```