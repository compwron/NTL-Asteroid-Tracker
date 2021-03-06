/**
 * Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
 */
package asteroidtracker;

import java.awt.Point;
import java.awt.geom.*;
import java.io.Serializable;
import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;
import java.io.FileNotFoundException;
import java.util.List;

class TestCase implements Serializable {
    static public double SPEED_OF_LIGHT = 299792458.0;
    static public double MAX_TRANSMITTING_POWER = 40000.0;
    static public double BACKGROUND_NOISE = 1.0e-30;
    static public double SHIELDING_MULTIPLIER = 1.0e-27;
    static public double RELOCATE_SPEED = Math.PI / (10 * 60);
    static public double RELOCATE_POWER = 5000.0;
    static public double SIMULATION_TIME = 7 * 24 * 60 * 60;
    static public double T_MIN = 0.1;
    static public double CRITICAL_TRACKING_ANGLE = 1.0 * Math.PI / (180.0 * 60);
    static public double ANTENNA_SAFE_RADIUS = 11.0;

    static public double Q_LOST = (24 * 60 * 60) / Math.log(2.0);
    static public double Q_TRAJECTORY = 6.0e3;
    static public double Q_IMAGE = 10.0e5;
    static public double IMAGE_SCORE_MULTIPLIER = 30.0;
    static public double TRAJECTORY_SCORE_MULTIPLIER = 60.0 / SIMULATION_TIME;

    public Point2D.Double antennaPositions[];
    public Asteroid asteroids[];
    public double antennaMinDistance;
    public double antennaMaxDistance;
    public double peakGain;
    public double minDistanceGain[];
    public double maxDistanceGain[];

public void readConfig() throws FileNotFoundException
{
Scanner sc = new Scanner(new File("config"));
while(sc.hasNext())
{
	String param=sc.next();
	sc.next();
	double value=sc.nextDouble();
	sc.nextLine();


	if(param.equals("SPEED_OF_LIGHT"))
		SPEED_OF_LIGHT=value;
	if(param.equals("MAX_TRANSMITTING_POWER"))
{
System.out.println("RRR");
		MAX_TRANSMITTING_POWER=value;
}
	if(param.equals("BACKGROUND_NOISE"))
		BACKGROUND_NOISE=value;
	if(param.equals("SHIELDING_MULTIPLIER"))
		SHIELDING_MULTIPLIER=value;
	if(param.equals("RELOCATE_SPEED"))
		RELOCATE_SPEED=value;
	if(param.equals("RELOCATE_POWER"))
		RELOCATE_POWER=value;
//	if(param.equals("SIMULATION_TIME"))
//		SIMULATION_TIME=value;
	if(param.equals("T_MIN"))
		T_MIN=value;
	if(param.equals("CRITICAL_TRACKING_ANGLE"))
		CRITICAL_TRACKING_ANGLE=value;
	if(param.equals("ANTENNA_SAFE_RADIUS"))
		ANTENNA_SAFE_RADIUS=value;
	if(param.equals("Q_LOST"))
		Q_LOST=value;
	if(param.equals("Q_TRAJECTORY"))
		Q_TRAJECTORY=value;
	if(param.equals("Q_IMAGE"))
		Q_IMAGE=value;
	if(param.equals("IMAGE_SCORE_MULTIPLIER"))
		IMAGE_SCORE_MULTIPLIER=value;
	if(param.equals("TRAJECTORY_SCORE_MULTIPLIER"))
		TRAJECTORY_SCORE_MULTIPLIER=value;

}
System.out.println("Parameters are:");
System.out.println("SPEED_OF_LIGHT: "+SPEED_OF_LIGHT);
System.out.println("MAX_TRANSMITTING_POWER: "+MAX_TRANSMITTING_POWER);
System.out.println("BACKGROUND_NOISE: "+BACKGROUND_NOISE);
System.out.println("SHIELDING_MULTIPLIER: "+SHIELDING_MULTIPLIER);
System.out.println("RELOCATE_SPEED: "+RELOCATE_SPEED);
System.out.println("RELOCATE_POWER: "+RELOCATE_POWER);
//System.out.println("SIMULATION_TIME: "+SIMULATION_TIME);
System.out.println("T_MIN: "+T_MIN);
System.out.println("CRITICAL_TRACKING_ANGLE: "+CRITICAL_TRACKING_ANGLE);
System.out.println("ANTENNA_SAFE_RADIUS: "+ANTENNA_SAFE_RADIUS);
System.out.println("Q_LOST: "+Q_LOST);
System.out.println("Q_TRAJECTORY: "+Q_TRAJECTORY);
System.out.println("Q_IMAGE: "+Q_IMAGE);
System.out.println("IMAGE_SCORE_MULTIPLIER: "+IMAGE_SCORE_MULTIPLIER);
System.out.println("TRAJECTORY_SCORE_MULTIPLIER: "+TRAJECTORY_SCORE_MULTIPLIER);
System.out.println();
}

    public void generate(String path) throws FileNotFoundException {
	
	readConfig();


        Scanner sc = new Scanner(new File(path));

        final int numberOfAntennas = sc.nextInt();

        this.antennaPositions = new Point2D.Double[numberOfAntennas];

        for (int antennaIndex = 0; antennaIndex < numberOfAntennas; ++antennaIndex) {
            final double x = sc.nextDouble();
            final double y = sc.nextDouble();
            this.antennaPositions[antennaIndex] = new Point2D.Double(x, y);
        }

        this.peakGain = sc.nextDouble();

        final int gainLength = sc.nextInt();
        this.minDistanceGain = new double[gainLength];
        this.maxDistanceGain = new double[gainLength];

        {
            for (int i = 0; i < gainLength; ++i) {
                this.minDistanceGain[i] = sc.nextDouble();
            }
        }

        {
            for (int i = 0; i < gainLength; ++i) {
                this.maxDistanceGain[i] = sc.nextDouble();
            }
        }

        final int numberOfAsteroids = sc.nextInt();
        this.asteroids = new Asteroid[numberOfAsteroids];

        for (int asteroidIndex = 0; asteroidIndex < numberOfAsteroids; ++asteroidIndex) {
            Asteroid asteroid = new Asteroid();
            asteroid.scienceScoreMultiplier = sc.nextDouble();
            asteroid.reflectivityMultiplier = sc.nextDouble();
            asteroid.initialImageInformation = sc.nextDouble();
            asteroid.initialTrajectoryInformation = sc.nextDouble();

            final int len = sc.nextInt();

            asteroid.trajectory = new AsteroidSighting[len];

            for (int i = 0; i < len; ++i) {
                AsteroidSighting s = new AsteroidSighting();
                s.pos = new Vector3d();

                s.time = sc.nextDouble();
                s.pos.x = sc.nextDouble();
                s.pos.y = sc.nextDouble();
                s.pos.z = sc.nextDouble();

                asteroid.trajectory[i] = s;
            }

            this.asteroids[asteroidIndex] = asteroid;
        }

        this.antennaMinDistance = Double.MAX_VALUE;
        this.antennaMaxDistance = 0.0;

        for (int i = 0; i < numberOfAntennas; ++i) {
            for (int j = 0; j < numberOfAntennas; ++j) {
                if (i != j) {
                    double dx = antennaPositions[i].x - antennaPositions[j].x;
                    double dy = antennaPositions[i].y - antennaPositions[j].y;
                    double dist = Math.sqrt(dx * dx + dy * dy);

                    this.antennaMinDistance = Math.min(this.antennaMinDistance, dist);
                    this.antennaMaxDistance = Math.max(this.antennaMaxDistance, dist);
                }
            }
        }
    }
}
