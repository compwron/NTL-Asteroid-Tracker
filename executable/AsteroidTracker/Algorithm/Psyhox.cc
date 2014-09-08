bool stage2 = false;

double gUseVis = 1;
double gLastUsed = 1.25;
double gScienceMultiplier = 0.86;
double gMinimumScience = 0.54;
double gImageThreshold = 0.15;
double gPowerMultiplier = 1.1;
double gSacriface = 1.0;

double timer1 = 0;
double timer2 = 0;
double timer3 = 0;
double timer4 = 0;

bool VERBOSE = false;
bool SHOW_ASTEROIDS = false;

#define TIMERS

#include<stdio.h>
#include<stdlib.h>
//#include <bits/stdc++.h>
// C
#ifndef _GLIBCXX_NO_ASSERT
#include <cassert>
#endif
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <ciso646>
#include <climits>
#include <clocale>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
 
#if __cplusplus >= 201103L
#include <ccomplex>
#include <cfenv>
#include <cinttypes>
#include <cstdalign>
#include <cstdbool>
#include <cstdint>
#include <ctgmath>
#include <cwchar>
#include <cwctype>
#endif
 
// C++
#include <algorithm>
#include <bitset>
#include <complex>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>
 
#if __cplusplus >= 201103L
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <forward_list>
#include <future>
#include <initializer_list>
#include <mutex>
#include <random>
#include <ratio>
#include <regex>
#include <scoped_allocator>
#include <system_error>
#include <thread>
#include <tuple>
#include <typeindex>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#endif
#include <sys/time.h>
#include <iostream>
#include <fstream>
//#include <emmintrin.h>

using namespace std;

#define INLINE   inline __attribute__ ((always_inline))
#define NOINLINE __attribute__ ((noinline))

#define ALIGNED __attribute__ ((aligned(16)))

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

#define SSELOAD(a)     _mm_load_si128((__m128i*)&a)
#define SSESTORE(a, b) _mm_store_si128((__m128i*)&a, b)

#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
#define LL          long long
#define ULL         unsigned long long
#define LD          long double
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define PII         pair <int, int>
#define VI          VC < int >
#define VVI         VC < VI >
#define VVVI        VC < VVI >
#define VPII        VC < PII >
#define VD          VC < double >
#define VVD         VC < VD >
#define VS          VC < string >
#define VVS         VC < VS >
#define DB(a)       cerr << #a << ": " << (a) << endl;


template<class T> void print(VC < T > v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]" << endl;}
template<class T> string i2s(T x) {ostringstream o; o << x; return o.str();}
VS splt(string s, char c = ' ') {VS all; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) all.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.S) all.PB(s.substr(p)); return all;}

const double EPS = 1e-6;
const double PI = 2 * acos(0);

double SPEED_OF_LIGHT = 299792458.0;
double MAX_TRANSMITTING_POWER = 40000.0;
double BACKGROUND_NOISE = 1.0e-30;
double SHIELDING_MULTIPLIER = 1.0e-27;
double RELOCATE_SPEED = PI / (10 * 60);
double RELOCATE_POWER = 5000.0;
const double SIMULATION_TIME = 7 * 24 * 60 * 60;
const int TOTAL_TIME = 60 * 60 * 24 * 7;
double T_MIN = 0.1;
double CRITICAL_TRACKING_ANGLE = 1.0 * PI / (180.0 * 60);
double ANTENNA_SAFE_RADIUS = 11.0 + EPS;

double Q_LOST = (24 * 60 * 60) / log(2.0);
double Q_TRAJECTORY = 6.0e3;
double Q_IMAGE = 10.0e5;
double IMAGE_SCORE_MULTIPLIER = 30.0;
double TRAJECTORY_SCORE_MULTIPLIER = 60.0 / SIMULATION_TIME;

const int MAX_ASTEROIDS = 75;
const int MAX_ANTENNAS = 20;

const int TURN_LENGTH = 180;

const int TOTAL_TURNS = TOTAL_TIME / TURN_LENGTH + 1;

// const int MIN_SIGNAL_DIFFERENCE = 50;
const int MIN_SIGNAL_DIFFERENCE = 10;

double expPCT[2000];
INLINE double myexp(double x) {
	int a = (int)(x * 10);
	return expPCT[a];
}

double tanhPCT[3000];
INLINE double mytanh(double x) {
	int a = (int)(x * 1000);
	return a < 3000 ? tanhPCT[a] : 1;
}

double initexp() {
	REP(i, 2000) expPCT[i] = exp(-(i * 0.1 + 0.05) / Q_LOST);
}

double inittanh() {
	REP(i, 3000) tanhPCT[i] = tanh(i * 0.001 + 0.0005);
}

double getTime() {
	timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}

struct RNG {
    unsigned int MT[624];
    int index;
	
	RNG(int seed = 1) {
		init(seed);
	}
    
    void init(int seed = 1) {
        MT[0] = seed;
        FOR(i, 1, 624) MT[i] = (1812433253UL * (MT[i-1] ^ (MT[i-1] >> 30)) + i);
        index = 0;
    }
    
    void generate() {
        const unsigned int MULT[] = {0, 2567483615UL};
        REP(i, 227) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i+397] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        FOR(i, 227, 623) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i-227] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        unsigned int y = (MT[623] & 0x8000000UL) + (MT[0] & 0x7FFFFFFFUL);
        MT[623] = MT[623-227] ^ (y >> 1);
        MT[623] ^= MULT[y&1];
    }
    
    unsigned int rand() {
        if (index == 0) {
            generate();
        }
        
        unsigned int y = MT[index];
        y ^= y >> 11;
        y ^= y << 7  & 2636928640UL;
        y ^= y << 15 & 4022730752UL;
        y ^= y >> 18;
        index = index == 623 ? 0 : index + 1;
        return y;
    }
    
    INLINE int next(int x) {
        return rand() % x;
    }
    
    INLINE int next(int a, int b) {
        return a + (rand() % (b - a));
    }
    
    INLINE double nextDouble() {
        return (rand() + 0.5) * (1.0 / 4294967296.0);
    }
};

static RNG rng;

double peakGain;
VD minDistanceGain;
VD maxDistanceGain;
double antennaMaxDistance;
double antennaMinDistance;


struct Vector3D {
	double x, y, z;
	
	Vector3D(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	
	Vector3D() {
		x = y = z = 0;
	}
	
	Vector3D cross(Vector3D &v) {
		Vector3D rv;
		rv.x = y * v.z - z * v.y;
		rv.y = v.x * z - v.z * x;
		rv.z = x * v.y - y * v.x;
		return rv;
	}
	
	double length() {
		return sqrt(x * x + y * y + z * z);
	}
	
	double lengthSq() {
		return x * x + y * y + z * z;
	}
	
	Vector3D scale(double scale) {
		return Vector3D(x * scale, y * scale, z * scale);
	}
	
	Vector3D normalize() {
		double len = length();
		return scale(1 / len);
	}
	
	Vector3D sub(Vector3D &v) {
		return Vector3D(x - v.x, y - v.y, z - v.z);
	}
	
	double dot(Vector3D &v) {
		return x * v.x + y * v.y + z * v.z;
	}
	
	double angle(Vector3D &v) {
		double cosValue = dot(v) / length() / v.length();
		if (cosValue < -1) cosValue = -1;
		if (cosValue > +1) cosValue = +1;
		return acos(cosValue);
	}
};

struct Asteroid {
	int id;

	bool known;
	double science;
	double reflectivity;
	double initialImageInformation;
	double initialTrajectoryInformation;
	double totalVisibility;
	double relativeVisibility;
	VD trajectory;
	VD positions;
	
	Vector3D currentPosition;
	double imageScore;
	double trajectoryScore;
	double trajectoryInformation;
	queue<pair<double, double> > changeSignalEvents;
	double currentSignal;
	
	Vector3D getPosition(int time) {
		int turn = time / TURN_LENGTH;
		if (!known || turn < 0 || turn >= TOTAL_TURNS || positions[turn * 3] == 0) return Vector3D();
		
		return Vector3D(positions[turn * 3], positions[turn * 3 + 1], positions[turn * 3 + 2]);	
	}
	
	Vector3D getPositionTurn(int turn) {
		return Vector3D(positions[turn * 3], positions[turn * 3 + 1], positions[turn * 3 + 2]);	
	}
	
	Vector3D getPosition(double time) {
		return getPosition((int)(floor(time) + EPS));
	}
	
	double distance(int time) {
		return getPosition(time).length();
	}
	
	double distance(double time) {
		return getPosition(time).length();
	}
	
	Asteroid() {
		known = false;
imageScore = 0;
trajectoryScore = 0;
trajectoryInformation = 0;
currentSignal = 0;
	}
	
	Asteroid(int id, double science, double reflectivity, double initialImageInformation, double initialTrajectoryInformation, VD &trajectory) {
imageScore = 0;
trajectoryScore = 0;
trajectoryInformation = 0;
currentSignal = 0;
		this->id = id;
		this->science = science;
		this->reflectivity = reflectivity;
		this->initialImageInformation = initialImageInformation;
		this->initialTrajectoryInformation = initialTrajectoryInformation;
		this->trajectory = trajectory;
		known = true;
		
		this->positions = VD(TOTAL_TURNS * 3, -1);
		
		int firstTurn = (int)trajectory[0] / TURN_LENGTH;
		
		REP(i, trajectory.S / 4) {
			int turn = (int)trajectory[i * 4 + 0] / TURN_LENGTH;
			if (turn >= TOTAL_TURNS) continue;
			REP(j, 3) positions[turn * 3 + j] = trajectory[i * 4 + 1 + j];
		}
		FOR(i, 1, TOTAL_TURNS) if (positions[i * 3 + 0] == -1) {
			REP(j, 3) positions[i * 3 + j] = positions[i * 3 + j - 3];
		}
		
		imageScore = initialImageInformation * Q_IMAGE;
		trajectoryInformation = initialTrajectoryInformation * Q_TRAJECTORY;
		currentPosition = Vector3D(trajectory[1], trajectory[2], trajectory[3]);
		
		REP(i, TOTAL_TURNS) relativeVisibility += getPositionTurn(i).z > 1e8;
		totalVisibility = relativeVisibility;
		
		totalVisibility /= TOTAL_TURNS;
		relativeVisibility /= (TOTAL_TURNS - firstTurn);
		
	}
	
	double optimalRatio(int time, int n, double power) {
		double dist = currentPosition.lengthSq();
		return peakGain * reflectivity * n * n * n * power / dist / dist * 1e30;
	}
	
	double expectedRatio(int time, int n, double power) {
		double dist = currentPosition.lengthSq();
		return peakGain * reflectivity * max(mytanh(trajectoryInformation / Q_TRAJECTORY), T_MIN) * n * n * n * power / dist / dist * 1e30;
	}
	
	double expectedRatio(Vector3D &position, int n, double power, double d) {
		double dist = position.lengthSq();
		return peakGain * reflectivity * max(mytanh(trajectoryInformation / Q_TRAJECTORY), T_MIN) * n * n * n * power / dist / dist * 1e30;
	}
	
	double expectedInformation(int time, int n, double power) {
		double ratio = expectedRatio(time, n, power);
		return ratio > 1.0 ? 10 * log10(ratio) : 0;
	}
};

struct Antenna {
	int id;
	double posX;
	double posY;
	
	int power;
	int asteroid;
	
	bool relocating;
	Vector3D position;
	
	Antenna() {
power = 0;
asteroid = -1;
	
relocating = false;
position = Vector3D(0.0, 0.0, 1.0);
}
	
	Antenna(int id, double posX, double posY) {
power = 0;
asteroid = -1;
	
relocating = false;
position = Vector3D(0.0, 0.0, 1.0);
		this->id = id;
		this->posX = posX;
		this->posY = posY;
	}
};


Asteroid asteroids[MAX_ASTEROIDS];
Antenna antennas[MAX_ANTENNAS];
int antennasNo;

double asteroidImageScoreCopy[MAX_ASTEROIDS];
double asteroidTrajectoryInformationCopy[MAX_ASTEROIDS];
Vector3D asteroidCurrentPositionCopy[MAX_ASTEROIDS];

void saveAsteroidData() {
	REP(i, MAX_ASTEROIDS) if (asteroids[i].known) {
		asteroidImageScoreCopy[i] = asteroids[i].imageScore;
		asteroidTrajectoryInformationCopy[i] = asteroids[i].trajectoryInformation;
		asteroidCurrentPositionCopy[i] = asteroids[i].currentPosition;
	}
}

void restoreAsteroidData() {
	REP(i, MAX_ASTEROIDS) if (asteroids[i].known) {
		asteroids[i].imageScore = asteroidImageScoreCopy[i];
		asteroids[i].trajectoryInformation = asteroidTrajectoryInformationCopy[i];
		asteroids[i].currentPosition = asteroidCurrentPositionCopy[i];
	}
}


double interpolateGain(double distance, double angle) {
	const int len = minDistanceGain.S;
	
	double u = (len - 1) * angle / PI;
	int index = (int)floor(u);
	u -= index;
	
	if (index == len - 1) {
		index--;
		u += 1;
	}
	
	const double v = (distance - antennaMinDistance) / (antennaMaxDistance - antennaMinDistance);
	
	const double p0 = (1 - v) * minDistanceGain[index + 0] + v * maxDistanceGain[index + 0];
	const double p1 = (1 - v) * minDistanceGain[index + 1] + v * maxDistanceGain[index + 1];
	
	return (1 - u) * p0 + u * p1;
}

double calcInducedNoise(Antenna &a1, Antenna &a2, Asteroid &t, double power) {
	Vector3D antennaDirection = t.currentPosition.normalize();
	double antennaDX = a2.posX - a1.posX;
	double antennaDY = a2.posY - a1.posY;
	double distanceBetweenAntennas = hypot(antennaDX, antennaDY);
	
    double scalarProduct = antennaDX * antennaDirection.x + antennaDY * antennaDirection.y;
	
	if (VERBOSE) cerr << t.id << ' ' << t.currentPosition.x << ' ' << t.currentPosition.y << ' ' << t.currentPosition.z << ' ' << antennaDX << ' ' << antennaDY << ' ' << scalarProduct << ' ' << distanceBetweenAntennas << ' ' << (scalarProduct / distanceBetweenAntennas) << endl;
    double angle = acos(scalarProduct / distanceBetweenAntennas);
	
    return power * interpolateGain(distanceBetweenAntennas, angle) / (distanceBetweenAntennas * distanceBetweenAntennas);
}

double calcInducedNoise(Antenna &a, Asteroid &t, double power) {
	double total = 0;
	REP(i, antennasNo) if (i != a.id) total += calcInducedNoise(a, antennas[i], t, power);
	DB(total);
	return total;
}



bool checkBeamIntersection(Antenna &a1, Antenna &a2, Asteroid &t, int time) {
	Vector3D receiverPosition;
	receiverPosition.x = a2.posX - a1.posX;
	receiverPosition.y = a2.posY - a1.posY;
	
	Vector3D beamDirection = t.getPosition(time).normalize();
	
	double closestBeamDistance = beamDirection.dot(receiverPosition);
	
	if (closestBeamDistance <= 0) return false;
	
	Vector3D closestBeamPosition = beamDirection.scale(closestBeamDistance);
	Vector3D relativeBeamPosition = receiverPosition.sub(closestBeamPosition);
	double distanceFromBeam = relativeBeamPosition.length();
	
	return distanceFromBeam < ANTENNA_SAFE_RADIUS;
}

bool checkBeamIntersection(Antenna &a, Asteroid &t, int time) {
	REP(i, antennasNo) if (i != a.id && checkBeamIntersection(a, antennas[i], t, time)) return true;
	return false;
}

int xplanSize[TOTAL_TURNS + 1000];
int xplan[TOTAL_TURNS + 1000][100];
int planSize[TOTAL_TURNS + 1000];
int plan[TOTAL_TURNS + 1000][100];

void savePlan(int turn) {
	FOR(t, turn, turn + 100) {
		xplanSize[t] = planSize[t];
		REP(i, xplanSize[t]) xplan[t][i] = plan[t][i];
	}
}

void restorePlan(int turn) {
	FOR(t, turn, turn + 100) {
		planSize[t] = xplanSize[t];
		REP(i, xplanSize[t]) plan[t][i] = xplan[t][i];
	}
}

VC<pair<double,int> > vp;

class AsteroidTracker {
public:

AsteroidTracker() {

time = 0;
lastTime = 0;
giveOrders = true;
newAsteroid = false;
lastCmd = "";
timePassed = 0;
deltaTimeSum = 0;
plansSimulated = 0;
turnsSimulated = 0;
lastChange = -100000;
checkScoreAt = 0;

	initexp();
	inittanh();
}

int time;
int lastTime;
bool giveOrders;
bool newAsteroid;
#define PIS pair<int, string>
priority_queue<PIS, VC<PIS>, greater<PIS> > cmdPQ;
string lastCmd;

double timePassed;

VVI getSubarrays() { 
	VVI rv(MAX_ASTEROIDS);

	REP(i, antennasNo) {
		Antenna &a = antennas[i];
		if (a.asteroid != -1 && !a.relocating && asteroids[a.asteroid].currentPosition.z >= 0)
			rv[a.asteroid].PB(i);
	}

	return rv;
}

double deltaTimeSum;

//simulator stuff
void passTime(double deltaTime) {
	// if (deltaTime < 0) DB(deltaTime);

	// if (deltaTime == 0) return;
	if (deltaTime <= 0) return;
	
	deltaTimeSum += deltaTime;

	//follow
	double angle = deltaTime * RELOCATE_SPEED;
	double cosA = cos(angle);
	double oneMinusCosA = 1.0 - cosA;
	double sinA = sin(angle);
	
	double energy = 0;
	
	//relocate
	REP(antennaIndex, antennasNo) {
		Antenna &a = antennas[antennaIndex];
	
		Vector3D antennaDirection = a.position;
		if (a.asteroid == -1 || !a.relocating) continue;

		Vector3D asteroidPosition = asteroids[a.asteroid].currentPosition;

		Vector3D rotationAxis = antennaDirection.cross(asteroidPosition).normalize();
		double ux = rotationAxis.x;
		double uy = rotationAxis.y;
		double uz = rotationAxis.z;

		double m00 = cosA + ux * ux * oneMinusCosA;
		double m01 = ux * uy * oneMinusCosA - uz * sinA;
		double m02 = ux * uz * oneMinusCosA + uy * sinA;
		double m10 = ux * uy * oneMinusCosA + uz * sinA;
		double m11 = cosA + uy * uy * oneMinusCosA;
		double m12 = uy * uz * oneMinusCosA - ux * sinA;
		double m20 = ux * uz * oneMinusCosA - uy * sinA;
		double m21 = uy * uz * oneMinusCosA + ux * sinA;
		double m22 = cosA + uz * uz * oneMinusCosA;

		double newX = m00 * antennaDirection.x + m01 * antennaDirection.y + m02 * antennaDirection.z;
		double newY = m10 * antennaDirection.x + m11 * antennaDirection.y + m12 * antennaDirection.z;
		double newZ = m20 * antennaDirection.x + m21 * antennaDirection.y + m22 * antennaDirection.z;

		a.position = Vector3D(newX, newY, newZ).normalize();
	}
	
	//update energy
	/*
	REP(antennaIndex, antennaNo) {
		
		if (antennaAssignments[antennaIndex] == -1) {
			continue;
		}

		if (isRelocating[antennaIndex]) {
			this.energySpent += deltaTime * testCase.RELOCATE_POWER;
			this.energySpentTurning += deltaTime * testCase.RELOCATE_POWER;
		}
		else {
			this.energySpent += deltaTime * transmittingPowers[antennaIndex];
			this.energySpentTransmitting += deltaTime * testCase.RELOCATE_POWER;
		}
	}
	*/

	VVI subarrays = getSubarrays();
	double decayMultiplier = exp(-deltaTime / Q_LOST);
	REP(asteroidIndex, MAX_ASTEROIDS) {
		Asteroid &a = asteroids[asteroidIndex];
		if (!a.known) continue;
		
		VI &antennaIndices = subarrays[a.id];

		double informationRate = 0;
		if (antennaIndices.S) {
		
			double dist = a.currentPosition.lengthSq();

			double powerReturn = peakGain
               * a.reflectivity
               * max(tanh(a.trajectoryInformation / Q_TRAJECTORY), T_MIN)
               * a.currentSignal
               * (antennaIndices.S * antennaIndices.S)
               / (dist * dist);

			double inducedNoiseSum = 0.0;
			REP(i, antennasNo) for (int receivingAntennaIndex : antennaIndices) if (i != receivingAntennaIndex)
				inducedNoiseSum += calcInducedNoise(antennas[i], antennas[receivingAntennaIndex], a, antennas[i].power);

			double noise = SHIELDING_MULTIPLIER * inducedNoiseSum + antennaIndices.S * BACKGROUND_NOISE;

			double ratio = powerReturn / noise;
			if (ratio > 1.0) informationRate = 10.0 * log10(ratio);
		}

		a.imageScore += deltaTime * informationRate;

		double previousTrajectoryInformation = a.trajectoryInformation;
		double newTrajectoryInformation = decayMultiplier * previousTrajectoryInformation + (1.0 - decayMultiplier) * Q_LOST * informationRate;
		a.trajectoryInformation = newTrajectoryInformation;
		a.trajectoryScore += deltaTime * (tanh(previousTrajectoryInformation / Q_TRAJECTORY) + tanh(newTrajectoryInformation / Q_TRAJECTORY)) / 2.0;
	}
}

void simulate() {
	// cerr << "simulate() " << ' ' << lastTime << ' ' << time << endl;
	double curTime = lastTime;
	while (true) {
		//check signal update
		VVI subarrays = getSubarrays();
		REP(i, MAX_ASTEROIDS) {
			Asteroid &a = asteroids[i];
			if (!a.known) continue;
			
			double amplitudeSum = 0.0;
			for (int antennaIndex : subarrays[i])
				amplitudeSum += sqrt(antennas[antennaIndex].power);
			double combinedOutputPower = amplitudeSum * amplitudeSum;
			
			double lastOutputPower = a.changeSignalEvents.empty() ? a.currentSignal : a.changeSignalEvents.back().Y;
			
			if (lastOutputPower != combinedOutputPower) {
				double travelTime = 2.0 * a.currentPosition.length() / SPEED_OF_LIGHT;
				a.changeSignalEvents.push(MP(lastTime + travelTime, combinedOutputPower));
			}
		}
		
	
		double eventTime = 1e9;
		int eventType = -1; //0 - Antenna Relocation : 1 - Signal Change
		int eventValue = -1; //Antenna / Asteroid id
		
		REP(antennaIndex, antennasNo) {
			Antenna &a = antennas[antennaIndex];
			if (!a.relocating) continue;
			
			Vector3D targetPosition = asteroids[a.asteroid].currentPosition.normalize();
			double timeLeft = a.position.angle(targetPosition) / RELOCATE_SPEED;
			
			if (curTime + timeLeft < eventTime) {
				eventTime = curTime + timeLeft;
				eventType = 0;
				eventValue = antennaIndex;
			}
		}
		
		REP(i, MAX_ASTEROIDS) {
			Asteroid &a = asteroids[i];
			if (!a.known || a.changeSignalEvents.empty()) continue;
			
			if (a.changeSignalEvents.front().X < eventTime) {
				eventTime = a.changeSignalEvents.front().X;
				eventType = 1;
				eventValue = i;
			}
		}
		
		if (eventTime > time) {
			// cerr << time << ' ' << curTime << ' ' << (time - curTime) << endl;
			passTime(time - curTime);
			break;
		}
		
		// cerr << eventType << ' ' << eventValue << ' ' << eventTime << ' ' << curTime << ' ' << (eventTime - curTime) << endl;
		passTime(eventTime - curTime);
		curTime = eventTime;
		
		if (eventType == 0) {
			// cerr << "Done Relocating: " << eventValue << " to " << antennas[eventValue].asteroid << " at " << eventTime << endl;
			antennas[eventValue].relocating = false;
		} else {
			asteroids[eventValue].currentSignal = asteroids[eventValue].changeSignalEvents.front().Y;
			asteroids[eventValue].changeSignalEvents.pop();
		}
		
	}
	
	
}



//adding Commands
void addCommandRelocate(int time, int antenna, int asteroid) {
	char cmd[1000];
	sprintf(cmd, "%d %d R %d", time, antenna, asteroid);
	cmdPQ.push(MP(time, cmd));
}

void addCommandPower(int time, int antenna, int power) {
	char cmd[1000];
	sprintf(cmd, "%d %d P %d", time + (power > 0 ? 1 : 0), antenna, power);
	cmdPQ.push(MP(time, cmd));
}

void addCommandRelocate(int antenna, int asteroid) {
	if (antennas[antenna].asteroid == asteroid) return;
	addCommandRelocate(time, antenna, asteroid);
}

void addCommandPower(int antenna, int power) {
	if (abs(antennas[antenna].power - power) < MIN_SIGNAL_DIFFERENCE) return;
	addCommandPower(time, antenna, power);
}

VI generatePlan(int xtime, set<int> &lastUsed, bool useCurTime = false) {
	int ztime = useCurTime ? time : xtime;
	vp.clear();
	
	REP(i, MAX_ASTEROIDS) if (asteroids[i].known && asteroids[i].getPosition(xtime).z > 1e8 && (asteroids[i].imageScore < Q_IMAGE * 1.8 || asteroids[i].trajectoryInformation < Q_TRAJECTORY * 2) && asteroids[i].expectedRatio(xtime, antennasNo, 60000 / antennasNo) > 0.25) {
		double value = 1;
		if (lastUsed.count(i)) {
			value *= gLastUsed;
			if (ztime > TOTAL_TIME - 3600) value *= 1.5;
		} else {
			double vtime = ztime + 2 * asteroids[i].currentPosition.length() / SPEED_OF_LIGHT;
			if (asteroids[i].getPosition(vtime + 5 * TURN_LENGTH).z < 1e8) value *= 0.85;
			if (asteroids[i].getPosition(vtime + 10 * TURN_LENGTH).z < 1e8) value *= 0.85;
			if (asteroids[i].getPosition(vtime + 15 * TURN_LENGTH).z < 1e8) value *= 0.85;
		}
		
		if (gUseVis) value /= pow(asteroids[i].totalVisibility, 0.1);
		// value /= pow(
		value *= pow(asteroids[i].science, gScienceMultiplier);
		value /= pow(asteroids[i].distance(xtime), 0.6);
		value *= asteroids[i].imageScore > Q_IMAGE * gImageThreshold - 1e-9 ? 6 * (mytanh(sqrt(asteroids[i].imageScore / Q_IMAGE) + .1) - mytanh(sqrt(asteroids[i].imageScore / Q_IMAGE))) : 1;
		value *= (xtime < 4.5 * TOTAL_TIME / 7 && asteroids[i].trajectoryInformation < Q_TRAJECTORY * 1.75) ? 1.3 : 1;
		value *= (xtime < 1 * TOTAL_TIME / 7 && asteroids[i].trajectoryInformation < Q_TRAJECTORY * 1.75) ? 1.2 : 1;
		vp.PB(MP(-value, i));
	}
	
	// if (vp.S > 10) {
		// nth_element(vp.begin(), vp.begin() + 10, vp.end());
		// sort(vp.begin(), vp.begin() + 10);
	// } else {
		sort(ALL(vp));
	// }
	
	int size = min(1, (int)vp.S);
	
	int tries = 0;
	while (size < vp.S && asteroids[vp[size].Y].science > gMinimumScience && antennasNo / (size + 1) >= 2) {
		int oldAmount = antennasNo / size;
		int newAmount = antennasNo / (size + 1);
		bool ok = true;
		size++;
		REP(i, size) {
			double ratioOld = asteroids[vp[i].Y].expectedRatio(xtime, oldAmount, (55000 + size * 10000) / antennasNo);
			double ratioNew = asteroids[vp[i].Y].expectedRatio(xtime, newAmount, (55000 + size * 10000) / antennasNo);
			double optimalRatioNew = asteroids[vp[i].Y].optimalRatio(xtime, newAmount, (55000 + size * 10000) / antennasNo);
			if (i < size - 1 && ratioNew < 3 / gSacriface) ok = false;
			if (ratioNew < 0.5 / gSacriface || optimalRatioNew < 1.25) ok = false;
		}
		if (!ok) {
			size--;
			if (tries <= 3 && size + tries + 1 < vp.S) {
				tries++;
				swap(vp[size], vp[size+tries]);
			} else {
				break;
			}
		} else if (tries) {
			break;
		}		
	}
	
	VI rv(min(10, (int)vp.S) + 1);
	rv[0] = size;
	REP(i, min((int)vp.S, 10)) rv[i+1] = vp[i].Y;
	
	return rv;
}

//TODO: speedup
//TODO: simulate only specified asteroids //speed x15-30
//TODO: don't simulate inactive asteroids
//TODO: add custom mytanh
int plansSimulated;
int turnsSimulated;

double simulatePlan(int p[TOTAL_TURNS + 1000][100], int pno[TOTAL_TURNS + 1000], int turns, int generate = 10000) {
	plansSimulated++;
	
	int lastChangePlan = -10000;
	if (generate < 10000) saveAsteroidData();

	if (turns < 0) turns = 1 << 20;

	//assumes current starting conditions
	queue < pair < double, double > > signalsChange[MAX_ASTEROIDS];
	queue < pair < double, int > > antennaChange[MAX_ASTEROIDS]; //TODO: queue is wrong
	int antennaFollowing[MAX_ASTEROIDS] = {0};
	int antennaFollowingTotal[MAX_ASTEROIDS] = {0};
	int targets[MAX_ANTENNAS];
	double signals[MAX_ASTEROIDS];
	double imageScore[MAX_ASTEROIDS];
	double trajectoryScore[MAX_ASTEROIDS];
	double trajectoryInformation[MAX_ASTEROIDS];
	double power[MAX_ASTEROIDS];
	double newPower[MAX_ASTEROIDS];
	int newTargets[MAX_ANTENNAS];
	int targetsNo[MAX_ASTEROIDS];
	int needed[MAX_ASTEROIDS];
	double relocationTime[MAX_ASTEROIDS];
	int neededCopy[MAX_ASTEROIDS];
	Vector3D antennaDirection[MAX_ASTEROIDS];
	
	int totalAsteroids = 0;
	REP(i, MAX_ASTEROIDS) if (asteroids[i].known) totalAsteroids = i + 1;
	
	REP(i, antennasNo) {
		targets[i] = antennas[i].asteroid;
		if (antennas[i].asteroid == -1) {
			antennaDirection[i] = antennas[i].position;
			continue;
		}
		
		antennaFollowingTotal[targets[i]]++;
		if (antennas[i].relocating) {
			antennaChange[targets[i]].push(MP(time + asteroids[targets[i]].currentPosition.normalize().angle(antennas[i].position) / RELOCATE_SPEED, 1));			
		} else {
			antennaFollowing[targets[i]]++;
		}
	}
	
	REP(i, totalAsteroids) if (asteroids[i].known) {
		signals[i] = asteroids[i].currentSignal;
		if (antennaFollowingTotal[i]) signalsChange[i] = asteroids[i].changeSignalEvents;
		imageScore[i] = asteroids[i].imageScore;
		trajectoryScore[i] = asteroids[i].trajectoryScore;
		trajectoryInformation[i] = asteroids[i].trajectoryInformation;
	}
	
	
	int curTime = time;
	int turn = time / TURN_LENGTH + 1;
	
	double powerUsed = 0;
	
	REP(steps, turns) {
		if (curTime >= TOTAL_TIME) break;
		turnsSimulated++;
		
		int newTime = curTime + TURN_LENGTH;
		
		double timer1Start = getTime();
		if (steps >= generate) {
			// bool regenerate = steps == 0;
			bool regenerate = false;
			
			REP(i, totalAsteroids) if (asteroids[i].known) {
				asteroids[i].imageScore = imageScore[i];
				asteroids[i].trajectoryInformation = trajectoryInformation[i];
				asteroids[i].currentPosition = asteroids[i].getPositionTurn(turn - 1);
			}
			
			double maxDist = 0;
			REP(i, pno[turn-1]) maxDist = max(maxDist, asteroids[p[turn-1][i]].distance(curTime));
			regenerate |= curTime > lastChangePlan + 1 * maxDist / SPEED_OF_LIGHT + 10 * TURN_LENGTH;
			REP(i, pno[turn-1]) regenerate |= asteroids[p[turn-1][i]].getPosition(curTime).z < 1e8;
			REP(i, totalAsteroids) if (asteroids[i].known && asteroids[i].getPosition(curTime).z > 1e8 && asteroids[i].known && asteroids[i].getPosition(curTime - TURN_LENGTH).z < 1e8 && (imageScore[i] < Q_IMAGE * 1.8 || trajectoryInformation[i] < Q_TRAJECTORY * 2) && asteroids[i].expectedRatio(curTime, antennasNo, 60000 / antennasNo) > 0.25)
				regenerate = true;
				
			double timer4Start = getTime();
			if (regenerate) {
				set<int> lastUsed;
				REP(i, pno[turn-1]) lastUsed.insert(p[turn-1][i]);
				VI np = generatePlan(curTime, lastUsed);
				pno[turn] = np[0];
				FOR(i, 1, np.S) p[turn][i-1] = np[i];
				lastChangePlan = curTime;
			} else {
				pno[turn] = pno[turn-1];
				REP(i, totalAsteroids) p[turn][i] = p[turn-1][i];
			}
			timer4 += getTime() - timer4Start;
		}
		timer1 += getTime() - timer1Start;
		
		bool change = false;
		change |= pno[turn] != pno[turn-1];
		REP(i, pno[turn]) change |= p[turn][i] != p[turn-1][i];
		
		double timer2Start = getTime();
		//relocate and stuff
		if (change) {
			if (pno[turn] == 0) {
				REP(i, antennasNo) if (targets[i] != -1) {
					antennaFollowing[targets[i]] = 0;
					antennaFollowingTotal[targets[i]] = 0;
					while (!antennaChange[targets[i]].empty()) antennaChange[targets[i]].pop();
					while (!signalsChange[targets[i]].empty()) signalsChange[targets[i]].pop();
					signals[targets[i]] = 0;
					antennaDirection[i] = asteroids[targets[i]].getPosition(curTime - TURN_LENGTH).normalize();
					newTargets[i] = -1;
				}
			} else {
				REP(i, pno[turn]) targetsNo[i] = 0;
				REP(i, antennasNo) {
					targetsNo[i % pno[turn]]++;
					newTargets[i] = -1;
				}
				REP(i, pno[turn]) needed[i] = targetsNo[i];
			
				REP(i, antennasNo) {
					REP(j, pno[turn]) if (needed[j] > 0 && targets[i] == p[turn][j]) {
						needed[j]--;
						newTargets[i] = targets[i];
					}
				}
				
				REP(i, pno[turn]) {
					neededCopy[i] = needed[i];
					relocationTime[i] = 0;
				}
				
				REP(j, pno[turn]) {
					Vector3D apos = asteroids[p[turn][j]].getPosition(curTime);
					REP(i, antennasNo) if (needed[j] > 0 && newTargets[i] == -1) {
						needed[j]--;
						newTargets[i] = p[turn][j];
						antennaFollowingTotal[newTargets[i]]++;
						if (targets[i] != -1) {
							antennaFollowing[targets[i]]--;
							antennaFollowingTotal[targets[i]]--;
							relocationTime[j] = max(relocationTime[j], asteroids[targets[i]].getPosition(curTime).angle(apos) / RELOCATE_SPEED);
						} else {
							relocationTime[j] = max(relocationTime[j], antennaDirection[i].angle(apos) / RELOCATE_SPEED);
						}
					}
				}
				
				REP(i, pno[turn]) if (neededCopy[i])
					antennaChange[p[turn][i]].push(MP(curTime + relocationTime[i], neededCopy[i]));
			}
			
			REP(i, antennasNo) targets[i] = newTargets[i];
		}
		
		//compute power
		double defPower = (55000 + pno[turn] * 10000) / antennasNo;
		REP(i, pno[turn]) {
			// double ratio = asteroids[p[turn][i]].expectedRatio(curTime, targetsNo[i], newPower[a.id]);
			// if (ratio < 2.5) {
				// double adjustedPower = newPower[a.id] * 2.5 / ratio;
				// power[a.id] = (adjustedPower < 80000 && ratio > 0.25) ? min(MAX_TRANSMITTING_POWER, adjustedPower) : 0;
			// }
			double ratio = asteroids[p[turn][i]].expectedRatio(curTime, targetsNo[i], defPower);
			newPower[p[turn][i]] = ratio < 2.5 * gSacriface ? min(MAX_TRANSMITTING_POWER, defPower * 2.5 * gSacriface / ratio) : defPower;
			newPower[p[turn][i]] = ratio < 0.2 ? 0 : min(MAX_TRANSMITTING_POWER, newPower[p[turn][i]] * gPowerMultiplier);
		}
		
		
		//change power signal
		REP(i, pno[turn]) change |= power[p[turn][i]] != newPower[p[turn][i]];
		
		if (change) {
		
			//kill old targets
			REP(i, pno[turn-1]) {
				int id = p[turn-1][i];
				if (antennaFollowingTotal[id]) continue;
				power[id] = 0;
				signals[id] = 0;
				while (!antennaChange[id].empty()) antennaChange[id].pop();
				while (!signalsChange[id].empty()) signalsChange[id].pop();
			}
			
			//adjust new targets
			//TODO: early exit if power == newPower && old targets == new targets
			REP(i, pno[turn]) {
				int id = p[turn][i];
				power[id] = newPower[id];
				if (antennaFollowing[id]) {
					double newSignal = power[id] * antennaFollowing[id] * antennaFollowing[id];
					double oldSignal = signalsChange[id].empty() ? signals[id] : signalsChange[id].back().Y;
					if (abs(newSignal - oldSignal) > MIN_SIGNAL_DIFFERENCE) {
						signalsChange[id].push(MP(curTime + 2 * asteroids[id].getPosition(curTime).length() / SPEED_OF_LIGHT, newSignal));
					}
				}
			}
		}
		timer2 += getTime() - timer2Start;
		
		
		double timer3Start = getTime();
		double basicDecayMultiplier = myexp(TURN_LENGTH);
		REP(i, totalAsteroids) if (asteroids[i].known) {
			if (antennaFollowingTotal[i] == 0) {
				trajectoryScore[i] += TURN_LENGTH * (mytanh(trajectoryInformation[i] / Q_TRAJECTORY) + mytanh(basicDecayMultiplier * trajectoryInformation[i] / Q_TRAJECTORY)) / 2.0;
				trajectoryInformation[i] *= basicDecayMultiplier;
				continue;
			}
		
			double t = curTime;
			while (true) {
				double t0 = signalsChange[i].empty() ? 1e9 : signalsChange[i].front().X;
				double t1 = antennaChange[i].empty() ? 1e9 : antennaChange[i].front().X;
				double nt = min(t0, t1);
				
				if (nt >= newTime)
					nt = newTime;
				
				if (nt > t) {
					double deltaTime = nt - t;
					
					//TODO: add smaller steps?
					if (signals[i] && antennaFollowing[i] > 0) {
						double dist = asteroids[i].getPosition(curTime).lengthSq(); //TODO: precompute distances?
						double powerReturn = peakGain
						   * asteroids[i].reflectivity
						   * max(mytanh(trajectoryInformation[i] / Q_TRAJECTORY), T_MIN)
						   * signals[i]
						   * (antennaFollowing[i] * antennaFollowing[i])
						   / (dist * dist);
						double noise = antennaFollowing[i] * BACKGROUND_NOISE * 1.04; //TODO: adjust noise multiplier
						double ratio = powerReturn / noise;
						double informationRate = ratio > 1 ? 10.0 * log10(ratio) : 0;
						
						imageScore[i] += deltaTime * informationRate;

						double decayMultiplier = myexp(deltaTime);
						double previousTrajectoryInformation = trajectoryInformation[i];
						double newTrajectoryInformation = decayMultiplier * previousTrajectoryInformation + (1.0 - decayMultiplier) * Q_LOST * informationRate;
						trajectoryInformation[i] = newTrajectoryInformation;
						trajectoryScore[i] += deltaTime * (mytanh(previousTrajectoryInformation / Q_TRAJECTORY) + mytanh(newTrajectoryInformation / Q_TRAJECTORY)) / 2.0;
					} else {
						double decayMultiplier = myexp(deltaTime);
						trajectoryScore[i] += deltaTime * (mytanh(trajectoryInformation[i] / Q_TRAJECTORY) + mytanh(decayMultiplier * trajectoryInformation[i] / Q_TRAJECTORY)) / 2.0;
						trajectoryInformation[i] *= decayMultiplier;
					}
					
					if (antennaFollowing[i] > 0) powerUsed += power[i] * antennaFollowing[i] * deltaTime;
					if (antennaFollowingTotal[i] != antennaFollowing[i]) {
						// assert(antennaFollowingTotal[i] > antennaFollowing[i]);
						powerUsed += 5000 * (antennaFollowingTotal[i] - antennaFollowing[i]) * deltaTime;
					}
					
					t = nt;
				}
				
				if (nt == newTime) break;
				
				// assert(t0 >= curTime);
				// assert(t1 >= curTime);
				
				if (t0 < t1) {
					signals[i] = signalsChange[i].front().Y;
					signalsChange[i].pop();
				} else {
					antennaFollowing[i] += antennaChange[i].front().Y;
					antennaChange[i].pop();
					double signal = power[i] * antennaFollowing[i] * antennaFollowing[i];
					signalsChange[i].push(MP(t + 2 * asteroids[i].getPosition(t).length() / SPEED_OF_LIGHT, signal));
				}
				
			}
		}
		timer3 += getTime() - timer3Start;
		
		turn++;
		curTime = newTime;
	}
	
	if (generate < 10000) restoreAsteroidData();
	
	double totalScore = 0;
	REP(i, totalAsteroids) {
		Asteroid &a = asteroids[i];
		if (!a.known) continue;
		totalScore += a.science * (IMAGE_SCORE_MULTIPLIER * mytanh(imageScore[i] / Q_IMAGE) + TRAJECTORY_SCORE_MULTIPLIER * trajectoryScore[i]);;
	}
	
	powerUsed *= 1e-9;
	if (VERBOSE && stage2 && turns > 1) cerr << time << ' ' << powerUsed << ' ' << totalScore << ' ' << (totalScore - powerUsed) << endl;
	return totalScore - powerUsed;
}

/*
void modifyPlan(int p[TOTAL_TURNS + 1000][100], int pno[TOTAL_TURNS + 1000], int turns, int tries) {
	double curScore = simulatePlan(p, pno, turns);
	while (tries--) {
		
	}
}
*/


int initialize(VD antennaPositions, double peakGain, VD minDistanceGain, VD maxDistanceGain) {
	double startTime = getTime();
	
	antennasNo = antennaPositions.S / 2;
	REP(i, antennasNo) antennas[i] = Antenna(i, antennaPositions[i * 2], antennaPositions[i * 2 + 1]);
	
	::minDistanceGain = minDistanceGain;
	::maxDistanceGain = maxDistanceGain;
	::peakGain = peakGain;
	
	antennaMinDistance = 1e9;
	antennaMaxDistance = 0;
	REP(i, antennasNo) REP(j, i) {
		double d = hypot(antennas[i].posX - antennas[j].posX, antennas[i].posY - antennas[j].posY);
		antennaMinDistance = min(antennaMinDistance, d);
		antennaMaxDistance = max(antennaMaxDistance, d);
	}
	
	timePassed += getTime() - startTime;
	return 0;
}

	
int asteroidAppearance(int asteroidIndex, double scienceScoreMultiplier, double reflectivityMultiplier, double initialImageInformation, double initialTrajectoryInformation, VD trajectory) {
	double startTime = getTime();
	
	asteroids[asteroidIndex] = Asteroid(asteroidIndex, scienceScoreMultiplier, reflectivityMultiplier, initialImageInformation, initialTrajectoryInformation, trajectory);
	newAsteroid = true;
	
	if (VERBOSE) cerr << time << ' ' << asteroidIndex << endl;
	
	if (SHOW_ASTEROIDS) {
		double startDistance = asteroids[asteroidIndex].getPosition(trajectory[0]).length();
		double endDistance = asteroids[asteroidIndex].getPosition(min(TOTAL_TIME, (int)trajectory[trajectory.S - 4])).length();
		double aboveHorizon = 0;
		REP(i, TOTAL_TURNS) aboveHorizon += asteroids[asteroidIndex].getPosition(i * TURN_LENGTH).z > 0;
		aboveHorizon /= TOTAL_TURNS;
		cerr << asteroidIndex << ' ' << (1.0 * time / TOTAL_TIME) << ' ' << scienceScoreMultiplier << ' ' << reflectivityMultiplier << ' ' << startDistance << ' ' << endDistance << ' ' << aboveHorizon << ' ' << trajectory.S << endl;
	}
	
	timePassed += getTime() - startTime;
	return 0;
}

int lastChange;
int checkScoreAt;

string nextCommand(double currentTime) {
	double startTime = getTime();
	lastTime = time;
	time = (int)currentTime;
	
	bool newAsteroidCopy = newAsteroid;
	if (newAsteroid) {
		while (!cmdPQ.empty()) cmdPQ.pop();
		giveOrders = true;
		newAsteroid = false;
		checkScoreAt = 0;
	} else if (lastCmd.S) {
		int time;
		int antenna;
		char type;
		int value;
		sscanf(lastCmd.c_str(), "%d %d %c %d", &time, &antenna, &type, &value);
		if (VERBOSE && antenna != antennasNo - 1) cerr << time << ' ' << lastCmd << endl;
		if (type == 'R') {
			antennas[antenna].asteroid = value;
			antennas[antenna].relocating = value == -1 ? false : true;
		} else if (type == 'P') {
			antennas[antenna].power = value;
		} else {
			assert(false);
		}
	}
	
	//simulate stuff
	if (lastTime != time) {
		if (VERBOSE) REP(i, MAX_ASTEROIDS) if (asteroids[i].known)
			cerr << i << ' ' << asteroids[i].trajectoryInformation / Q_TRAJECTORY << ' ' << tanh(asteroids[i].imageScore / Q_IMAGE) << endl;
			
		simulate();
		if (time % TURN_LENGTH == 0) {
			//update asteroid position + follow
			REP(i, MAX_ASTEROIDS) asteroids[i].currentPosition = asteroids[i].getPosition(time);
			REP(i, antennasNo) {
				Antenna &a = antennas[i];
				if (a.asteroid == -1 || a.relocating || asteroids[a.asteroid].currentPosition.z <= 0) continue;
				a.position = asteroids[a.asteroid].currentPosition.normalize();
			}
			
		}
	}
	
	
	//TODO: energy saving
	
	//executePlan
	if (giveOrders) {
		// if (time == checkScoreAt) cerr << simulatePlan(plan, planSize, 0) << endl;
	
		int turn = time / TURN_LENGTH + 1;
		
		//create plan
		bool change = false;
		
		set<int> lastVisible;
		REP(i, MAX_ASTEROIDS) if (asteroids[i].known && asteroids[i].getPosition(time - TURN_LENGTH).z > 1e8) lastVisible.insert(i);
		
		set<int> visible;
		REP(i, MAX_ASTEROIDS) if (asteroids[i].known && asteroids[i].getPosition(time).z > 1e8) visible.insert(i);
		
		set<int> lastUsed;
		REP(i, planSize[turn-1]) lastUsed.insert(plan[turn-1][i]);
		
		double maxDist = 0;
		REP(i, planSize[turn-1]) maxDist = max(maxDist, asteroids[plan[turn-1][i]].distance(time));
		change |= !stage2 && time > lastChange + 1 * maxDist / SPEED_OF_LIGHT + 10 * TURN_LENGTH;
		REP(i, planSize[turn-1]) change |= !stage2 && !visible.count(plan[turn-1][i]);
		REP(i, MAX_ASTEROIDS) if (asteroids[i].known && asteroids[i].getPosition(time).z > 1e8 && (asteroids[i].imageScore < Q_IMAGE * 1.8 || asteroids[i].trajectoryInformation < Q_TRAJECTORY * 2) && asteroids[i].expectedRatio(time, antennasNo, 60000 / antennasNo) > 0.25 && asteroids[i].getPosition(time + 5 * TURN_LENGTH).z > 1e8)
			change |= !stage2 && !lastVisible.count(i);
			// change |= !lastVisible.count(i);
			
		change |= newAsteroidCopy && stage2;
		change |= time > lastChange + 120 * TURN_LENGTH;
		// if (stage2) {
			// change = newAsteroidCopy;
			// change |= time > lastChange + 50 * TURN_LENGTH;
		// }
		
		if (change) {
			if (turn > TOTAL_TURNS - 2500 && !stage2) stage2 = true;
			if (stage2) {
				static int third = 3;
				third = (third + 1) % 4;
				if (timePassed < 27 && third != 3) {
					int turnsLeft = TOTAL_TURNS - turn;
					double bv = 0;
					double bestImageThreshold = gImageThreshold;
					double bestScienceMultiplier = gScienceMultiplier;
					double bestLastUsed = gLastUsed;
					double bestPowerMultiplier = gPowerMultiplier;
					double bestSacriface = gSacriface;
					double bestUseVis = gUseVis;
					// savePlan(turn);
					// for (gSacriface = 1.0; gSacriface < 1.5 + EPS; gSacriface += 0.5)
					for (gUseVis = 0; gUseVis <= 1; gUseVis++)
					for (gPowerMultiplier = 1.1; gPowerMultiplier < 1.7 + EPS; gPowerMultiplier += 0.6) 
					for (gLastUsed = 1.25; gLastUsed < 1.8 + EPS; gLastUsed += 0.55) 
					for (gImageThreshold = 0.0; gImageThreshold < 0.30 + EPS; gImageThreshold += 0.15)
					for (gScienceMultiplier = 0.76; gScienceMultiplier < 0.96 + EPS; gScienceMultiplier += 0.1) {
						if (gUseVis == 0 && turnsLeft > 1500) continue;
						double av = simulatePlan(plan, planSize, -1, 0);
						if (av > bv) {
							bestPowerMultiplier = gPowerMultiplier;
							bestLastUsed = gLastUsed;
							bestImageThreshold = gImageThreshold;
							bestScienceMultiplier = gScienceMultiplier;
							bestSacriface = gSacriface;
							bestUseVis = gUseVis;
							bv = av;
						}
					}
					// restorePlan(turn);
					gPowerMultiplier = bestPowerMultiplier;
					gImageThreshold = bestImageThreshold;
					gScienceMultiplier = bestScienceMultiplier;
					gLastUsed = bestLastUsed;
					gSacriface = bestSacriface;
					gUseVis = bestUseVis;
				}
				// VERBOSE = true;
				// if (VERBOSE) cerr << bestImageThreshold << ' ' << bestScienceMultiplier << ' ' << bestLastUsed << ' ' << bestPowerMultiplier << endl;
				simulatePlan(plan, planSize, 480, 0);
				// VERBOSE = false;
				// lastChange = TOTAL_TIME;
				lastChange = time;
			} else {
				simulatePlan(plan, planSize, 1, 0);
				lastChange = time;
			}
		} else {
			if (!stage2) {
				planSize[turn] = planSize[turn-1];
				REP(i, MAX_ASTEROIDS) plan[turn][i] = plan[turn-1][i];
			}
		}
		
		// cerr << time << ' ' << turn << ' ' << planSize[turn] << ' ' << plan[turn][0] << ' ' << plan[turn][1] << ' ' << endl;
		
		double defPower = min((int)MAX_TRANSMITTING_POWER, (55000 + 10000 * planSize[turn]) / antennasNo);
		
		//set targets
		set<int> asteroidsFollowed;
		REP(i, antennasNo) if (antennas[i].asteroid != -1) asteroidsFollowed.insert(antennas[i].asteroid);
		
		set<int> curAsteroidsFollowed;
		REP(i, planSize[turn]) curAsteroidsFollowed.insert(plan[turn][i]);
		
		VI targets(antennasNo, -1);
		VI targetsNo(planSize[turn]);
		
		if (asteroidsFollowed != curAsteroidsFollowed && planSize[turn]) {
			VI t(planSize[turn]);
			REP(i, planSize[turn]) t[i] = plan[turn][i];
			
			REP(i, antennasNo) targetsNo[i % planSize[turn]]++;

			VI n = targetsNo;
			
			REP(i, antennasNo) {
				REP(j, planSize[turn]) if (antennas[i].asteroid == t[j] && n[j] > 0) {
					n[j]--;
					targets[i] = t[j];
				}
			}
			REP(i, antennasNo) {
				REP(j, planSize[turn]) if (targets[i] == -1 && n[j] > 0) {
					n[j]--;
					targets[i] = t[j];
					break;
				}
			}
		} else {
			REP(i, antennasNo) {
				targets[i] = antennas[i].asteroid;
				if (targets[i] == -1) continue;
				int pos = -1;
				REP(j, planSize[turn]) if (targets[i] == plan[turn][j]) pos = j;
				if (pos != -1) targetsNo[pos]++;
			}
		}
		
		//set power
		VD power(MAX_ASTEROIDS, defPower);
		
		
		// DB(time);
		REP(i, planSize[turn]) {
			Asteroid &a = asteroids[plan[turn][i]];
			
			// cerr << time << ' ' << a.id << ' ' << (a.trajectoryInformation / Q_TRAJECTORY) << ' ' << endl;
			
			assert(a.currentPosition.z > 0);
			double ratio = a.expectedRatio(time, targetsNo[i], power[a.id]);
			if (ratio < 2.5) {
				double newPower = power[a.id] * 2.5 * gSacriface / ratio;
				power[a.id] = (newPower < 80000 && ratio > 0.25) ? min(MAX_TRANSMITTING_POWER, newPower) : 0;
			}
			
			bool found = false;
			if (stage2) {
				int futureTurn = turn + (int)(2 * a.currentPosition.length() / SPEED_OF_LIGHT / TURN_LENGTH);
				REP(j, planSize[futureTurn]) found |= plan[futureTurn][j] == a.id;
			} else {
				VI futurePlan = generatePlan((int)(time + 2 * a.currentPosition.length() / SPEED_OF_LIGHT - TURN_LENGTH), lastUsed, true);
				REP(j, max(planSize[turn], futurePlan[0])) if (j + 1 < futurePlan.S) found |= futurePlan[j + 1] == a.id;
			}
			if (!found) power[a.id] = 0;
		}
		
		REP(i, MAX_ASTEROIDS) power[i] = min(MAX_TRANSMITTING_POWER, power[i] * gPowerMultiplier);
		
		//add commands
		REP(i, antennasNo) {
			bool use = targets[i] != -1 && power[targets[i]] > 0;
		
			if (use) {
				double totalNoise = 0;
				REP(j, antennasNo) if (j != i && !antennas[j].relocating && targets[j] != -1 && (asteroids[targets[j]].currentSignal > 0)) totalNoise += calcInducedNoise(antennas[i], antennas[j], asteroids[targets[i]], power[targets[i]]);
				if (checkBeamIntersection(antennas[i], asteroids[targets[i]], time) || totalNoise > 1.5e-3) use = false;
			}
			
			addCommandRelocate(i, targets[i]);
			addCommandPower(i, use ? power[targets[i]] : 0);
		}
		giveOrders = false;		
		
	}
	
	if (cmdPQ.S) {
		string rv = cmdPQ.top().Y; cmdPQ.pop();
		
		timePassed += getTime() - startTime;
		return lastCmd = rv;
	}
	
	const int STEP = 90;
	
	if (time == TOTAL_TIME - STEP - TURN_LENGTH) {
	//	DB(timePassed);
	//	DB(timer1);
	//	DB(timer2);
	//	DB(timer3);
	//	DB(timer4);
		// DB(plansSimulated);
		// DB(turnsSimulated);
		// double totalScore = 0;
		// REP(i, MAX_ASTEROIDS) {
			// Asteroid &a = asteroids[i];
			// if (!a.known) continue;
			// double scoreA = tanh(a.imageScore / Q_IMAGE);
			// double scoreB = a.trajectoryScore * TRAJECTORY_SCORE_MULTIPLIER / 60;
			// double score = a.science * (IMAGE_SCORE_MULTIPLIER * tanh(a.imageScore / Q_IMAGE) + TRAJECTORY_SCORE_MULTIPLIER * a.trajectoryScore);
			// totalScore += score;
			// cerr << i << ' ' << scoreA << ' ' << scoreB << ' ' << score << endl;
		// }
		// DB(totalScore);
		// DB(deltaTimeSum);
		// fflush(stderr);
		// REP(i, 2500000) timer2 += timer1;
	}
	
	char cmd[1000];
	sprintf(cmd, "%d %d P %d", (time + STEP) - time % STEP, antennasNo - 1, antennas[antennasNo-1].power); //Empty Command
	giveOrders = time % 180 >= 180 - STEP;
	
	timePassed += getTime() - startTime;
	return lastCmd = string(cmd);
}

};


void readConfig()
{
  string line;
  ifstream myfile ("config");
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
	istringstream iss(line);
	string command;
	iss >> command;

if(command=="SPEED_OF_LIGHT" || command=="MAX_TRANSMITTING_POWER" || command=="BACKGROUND_NOISE"|| command=="SHIELDING_MULTIPLIER"
|| command=="RELOCATE_SPEED" || command=="RELOCATE_POWER" || command=="T_MIN" || command=="CRITICAL_TRACKING_ANGLE"
|| command=="ANTENNA_SAFE_RADIUS" || command=="Q_LOST" || command=="Q_TRAJECTORY" || command=="Q_IMAGE"
|| command=="IMAGE_SCORE_MULTIPLIER" || command=="TRAJECTORY_SCORE_MULTIPLIER")
{
	string blank;
	iss >> blank;
	double value;
	iss >> value;
	if(command=="SPEED_OF_LIGHT")
		SPEED_OF_LIGHT=value;
	if(command=="MAX_TRANSMITTING_POWER")
		MAX_TRANSMITTING_POWER=value;
	if(command=="BACKGROUND_NOISE")
		BACKGROUND_NOISE=value;
	if(command=="SHIELDING_MULTIPLIER")
		SHIELDING_MULTIPLIER=value;
	if(command=="RELOCATE_SPEED")
		RELOCATE_SPEED=value;
	if(command=="RELOCATE_POWER")
		RELOCATE_POWER=value;
/*	if(command=="SIMULATION_TIME")
	{
		TOTAL_TIME=SIMULATION_TIME=value;
		TOTAL_TURNS=TOTAL_TIME / TURN_LENGTH + 1;
	}
*/	if(command=="T_MIN")
		T_MIN=value;
	if(command=="CRITICAL_TRACKING_ANGLE")
		CRITICAL_TRACKING_ANGLE=value;
	if(command=="ANTENNA_SAFE_RADIUS")
		ANTENNA_SAFE_RADIUS=value;
	if(command=="Q_LOST")
		Q_LOST=value;
	if(command=="Q_TRAJECTORY")
		Q_TRAJECTORY=value;
	if(command=="Q_IMAGE")
		Q_IMAGE=value;
	if(command=="IMAGE_SCORE_MULTIPLIER")
		IMAGE_SCORE_MULTIPLIER=value;
	if(command=="TRAJECTORY_SCORE_MULTIPLIER")
		TRAJECTORY_SCORE_MULTIPLIER=value;
}
    }
    myfile.close();
  }

  else cerr << "Unable to open file config";

cerr<<"Algorithm params:"<<endl;
cerr<<"SPEED_OF_LIGHT: "<<SPEED_OF_LIGHT<<endl;
cerr<<"MAX_TRANSMITTING_POWER: "<<MAX_TRANSMITTING_POWER<<endl;
cerr<<"BACKGROUND_NOISE: "<<BACKGROUND_NOISE<<endl;
cerr<<"SHIELDING_MULTIPLIER: "<<SHIELDING_MULTIPLIER<<endl;
cerr<<"RELOCATE_SPEED: "<<RELOCATE_SPEED<<endl;
cerr<<"RELOCATE_POWER: "<<RELOCATE_POWER<<endl;
cerr<<"T_MIN: "<<T_MIN<<endl;
cerr<<"CRITICAL_TRACKING_ANGLE: "<<CRITICAL_TRACKING_ANGLE<<endl;
cerr<<"ANTENNA_SAFE_RADIUS: "<<ANTENNA_SAFE_RADIUS<<endl;
cerr<<"Q_LOST: "<<Q_LOST<<endl;
cerr<<"Q_TRAJECTORY: "<<Q_TRAJECTORY<<endl;
cerr<<"Q_IMAGE: "<<Q_IMAGE<<endl;
cerr<<"IMAGE_SCORE_MULTIPLIER: "<<IMAGE_SCORE_MULTIPLIER<<endl;
cerr<<"TRAJECTORY_SCORE_MULTIPLIER: "<<TRAJECTORY_SCORE_MULTIPLIER<<endl;
}

int main()
{

readConfig();



AsteroidTracker at;

int numberOfAntennas;
cin >> numberOfAntennas;
VD positions;
for(int i=0;i<numberOfAntennas;i++)
{
	double a,b;
	cin >> a >> b;
	positions.push_back(a);
	positions.push_back(b);
}
double peakGain;
cin >> peakGain;
int len;
cin >> len;
VD minDG;
for(int i=0;i<len;i++)
{
	double a;
	cin >> a;
	minDG.push_back(a);
}
VD maxDG;
for(int i=0;i<len;i++)
{
	double a;
	cin >> a;
	maxDG.push_back(a);
}

at.initialize(positions, peakGain, minDG, maxDG);

string line;
getline (cin,line);
while(true)
{
getline (cin,line);
istringstream iss(line);

string com;
iss >> com;


if(com=="C")
{
	double ct;
	iss >> ct; 
	cout << at.nextCommand(ct) << endl;
}
else
{
	int index, tl;
	double ssm, rm, iii, iti;
	iss >> index >> ssm >> rm >> iii >> iti >> tl;
	VD tra;
	for(int i=0;i<tl;i++)
	{
		string temp;
		getline (cin,temp);
		istringstream iss2(temp);

		double t,x,y,z;
		iss2 >> t >> x >> y >> z;
		tra.push_back(t);
		tra.push_back(x);
		tra.push_back(y);
		tra.push_back(z);
	}

	at.asteroidAppearance(index, ssm, rm, iii, iti, tra);
}
}

return(0);
}


