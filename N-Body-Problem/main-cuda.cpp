#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>

using namespace std;

struct Particle {
	static mt19937 eng;
	double pos[3];
	double vel[3];
	Particle();
};

mt19937 Particle::eng(static_cast<mt19937::result_type>(12345));

__global__ void integrateKernel(Particle *particles,
		int NumParticles, double GM2, double InverseMass, double Dt);

void SphericToCartesian(double cartesianVector[3], const double sphericVector[3])
{
	cartesianVector[0] = sphericVector[0] * cos(sphericVector[1]) * sin(sphericVector[2]);
	cartesianVector[1] = sphericVector[0] * sin(sphericVector[1]) * sin(sphericVector[2]);
	cartesianVector[2] = sphericVector[0] * cos(sphericVector[2]);
}

Particle::Particle()
{
	const double maxVelocity = 5000.0;

	uniform_real_distribution<double> dist_radius(0.0, 1.0);
	uniform_real_distribution<double> dist_psi(0.0, 2 * M_PI);
	uniform_real_distribution<double> dist_phi(0.0, M_PI);

	double positionSpheric[3] = {dist_radius(eng), dist_psi(eng), dist_phi(eng)};

	SphericToCartesian(pos, positionSpheric);

	vel[0] = pos[0];
       	vel[1] = pos[1];
       	vel[2] = pos[2];
	double mod = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
	for (auto& v : vel)
		v *= maxVelocity / mod;
}

struct Problem {
	const double mG = 6.674e-11;
	const double mMass;
	const double mInverseMass;
	const double mDt;
	const unsigned mNumParticles;
	const double mGM2 = mG * mMass * mMass;
	int mNumSMs;
	Particle *mParticles;
	Particle *d_mParticles; // copy in the GPU
	Problem(double Mass, double dt, unsigned numParticles);
	~Problem();
	void mMoveToDevice();
	void mMoveToHost();

	void integrate();

	void write_csv(unsigned ts);
	void write_console(unsigned ts);

};


Problem::Problem(double Mass, double dt, unsigned numParticles)
: 
	mMass(Mass),
       	mInverseMass(1 / Mass),
	mGM2(mG * mMass * mMass),
       	mDt(dt),
       	mNumParticles(numParticles)
{
	mParticles = new Particle[mNumParticles];

	// create the particles in the device
	//
	cudaError_t code;
	code = cudaMalloc((void **)&d_mParticles, mNumParticles * sizeof(Particle));
	if (code != cudaSuccess) {
		cout << cudaGetErrorString(code) << endl;
	}

	cudaDeviceGetAttribute(&mNumSMs, cudaDevAttrMultiProcessorCount, 0);
	cout << "GPU" << endl;
	cout << "Number of Stream Multiprocessor SMs : " << mNumSMs << endl;
}

Problem::~Problem()
{ 
	// delete from host
	delete [] mParticles;

	// delete from device
	cudaError_t code;
	code = cudaFree(d_mParticles);
	if (code != cudaSuccess) {
		cout << cudaGetErrorString(code) << endl;
	}
}

void Problem::mMoveToDevice()
{
	cudaError_t code;
	code = cudaMemcpy(d_mParticles, mParticles, mNumParticles * sizeof(Particle), cudaMemcpyHostToDevice);
	if (code != cudaSuccess) {
		cout << cudaGetErrorString(code) << endl;
	}
}

void Problem::mMoveToHost()
{
	cudaMemcpy(mParticles, d_mParticles, mNumParticles * sizeof(Particle), cudaMemcpyDeviceToHost);
}

void Problem::write_csv(unsigned ts)
{
	std::stringstream ss;
	ss << "particles_" << ts << ".csv";
	ofstream outFile;
	outFile.open(ss.str());

	outFile << "x,y,z,vx,vy,vz" << endl;
	outFile << scientific ;
	for (int pi = 0; pi < mNumParticles; pi++) {
		outFile << mParticles[pi].pos[0] << ","
			<< mParticles[pi].pos[1] << ","
			<< mParticles[pi].pos[2] << "," 
			<< mParticles[pi].vel[0] << ","
			<< mParticles[pi].vel[1] << ","
			<< mParticles[pi].vel[2] << endl; 
	}
	outFile.close();
}

void Problem::write_console(unsigned ts)
{
	cout << scientific ;
	for (int pi = 0; pi < mNumParticles; pi++) {
		for (auto& p : mParticles[pi].pos)
			cout << setw(15) << p << ",";
		for (auto& v : mParticles[pi].vel)
			cout << setw(15) << v << ",";
		cout << endl;
	}
}

void Problem::integrate()
{
	integrateKernel<<<8 * mNumSMs, 64>>>(d_mParticles, mNumParticles, mGM2, mInverseMass, mDt);
	if (cudaSuccess != cudaGetLastError()) {
		cout << "Error in kernel execution" << endl;
	}
}

__global__
void integrateKernel(Particle *particles, int numParticles, double GM2, double inverseMass, double dt)
{

	for (int pi = blockIdx.x * blockDim.x + threadIdx.x; pi < numParticles; pi += blockDim.x * gridDim.x) {

		double force[3] = {0, 0, 0};

		// Calculate total force
		for (int pj = 0; pj < numParticles; pj++) {

			if (pj != pi) {
				const double dij[3] = {
					particles[pj].pos[0] - particles[pi].pos[0],
					particles[pj].pos[1] - particles[pi].pos[1],
					particles[pj].pos[2] - particles[pi].pos[2]};

				const double dist2 = dij[0] * dij[0] + dij[1] * dij[1] + dij[2] * dij[2];
				const double ConstDist2 = GM2 / dist2;
				const double idist = 1 / sqrt(dist2);
				//const double idist = 1 / dist2;

				// F = C * m * m / ||x2 - x1||^2 * (x2 - x1) / ||x2 - x1||
				force[0] += ConstDist2 * dij[0] * idist;
				force[1] += ConstDist2 * dij[1] * idist;
				force[2] += ConstDist2 * dij[2] * idist;
			}
		}

		// dv / dt = a = F / m
		particles[pi].vel[0] += force[0] * inverseMass * dt;
		particles[pi].vel[1] += force[1] * inverseMass * dt;
		particles[pi].vel[2] += force[2] * inverseMass * dt;
	}
	__syncthreads();

	//
	// Update pos this should be done after all forces/velocities have being computed
	//
	for (int pi = blockIdx.x * blockDim.x + threadIdx.x; pi < numParticles; pi += blockDim.x * gridDim.x) {
		// dx / dt = v
		particles[pi].pos[0] += particles[pi].vel[0] * dt;
		particles[pi].pos[1] += particles[pi].vel[1] * dt;
		particles[pi].pos[2] += particles[pi].vel[2] * dt;
	}
	__syncthreads();
}

int main()
{
	const int nTimeSteps = 5;
	const double Mass = 1e12;
	const double dt = 1e-4;
	const unsigned numParticles = 100000;
	Problem problem(Mass, dt, numParticles);
	int tsp = 0;

	problem.mMoveToDevice();

	for (int ts = 0; ts < nTimeSteps; ts++) {

		cout << ts << endl;

		//if (ts % 1 == 0) {
		//	problem.mMoveToHost();
		//	problem.write_csv(tsp++);
		//	problem.write_console(ts);
		//}

		problem.integrate();
	}
	return 0;
}
