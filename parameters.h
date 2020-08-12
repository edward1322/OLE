#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <complex>

# define PI      3.14159265358979323846  /* pi */
# define ONE_MINUS_COSZERO 1.0e-6	/* using in angle */
# define LIGHTSPEED 3*1e8	/* light speed */
# define CHANCE 0.1    /* used in roulette */
# define FACTOR 1.0/180.0 * PI

enum STATUS { ALIVE, DEAD };
enum bio_or_opt { bio, opt };
enum HIT { NOHIT, HITSEASURF, HITSEABOT };
enum spftyp { H_G, MHG, TTHG, FF, MIE_CAL,S_CUSTOME};
enum icdftyp { I_HG, I_MIE,I_CUSTOME };
enum  Phase_DTheta { Phase_DTheta0_1, Phase_DTheta0_5, Phase_DTheta_1};
enum seabottype { ruser_def, rsand,	rgreen,	rbrown	, rred};

class parameters
{
public:
    parameters(void);
    ~parameters(void);

    double *radArray;           // radArray: radius Array
    double *numDensityArray;    // numDensityArray: Number Density array
    double meanRadius;          // meanRadius: "mean radius" in Poly disperse and "radius" in mono disperse
    double minRadius;           // minRadius: minimum radius
    double maxRadius;           // maxRadius: maximum radius
    double stdDev;              // stdDev: Std. deviation
    int nRadius;                // nRadius: number of radii

    double *scatRefRealArray;   // scatRefRealArray: refractive index of scatterer -Real  data array
    double *scatRefImagArray;   // scatRefImagArray: refractive index of scatterer -Imag  data array
    double scatRefReal;         // scatRefReal: refractive index of scatterer - Real
    double scatRefImag;         // scatRefImg: refractive index of scatterer - Imaginary ((n-ik): negative sign convention as Scott Prahl's Mie calculator)
    double medRef;              // medRef: refractive index of medium
    double volFraction;         // volFraction: Volume fraction of sphere volume
    double sphNumDensity;       // sphNumDensity: Sphere concentration/volume (Number Density)

    double startWavel;          // startWavel: starting wavelength
    double endWavel;            // endWavel: starting wavelength
    double stepWavel;           // stepWavel: starting wavelength
    double *wavelArray;         // wavel: wavelength array
    int nWavel;                 // nWavel: number of wavelength
	int laserWavel;

    double minTheta;            // minTheta: minimum angle
    double maxTheta;            // maxTheta: maximum angle
    double stepTheta;           // stepTheta: maximum angle
    int nTheta;                 // nTheta: number of angles

    double **phaseFunctionAve;  // phaseFunctionAve: Average phase function
    double **phaseFunctionPara; // phaseFunctionPara: parallel phase function
    double **phaseFunctionPerp; // phaseFunctionPerp: perpendicularphase function
    double *curPolarAng;        // curPolarAng: current polar angles for phase function
    double *curPhaseFunc;       // curPhaseFunc: current Phase function
    double *cSca;               // cSca: Scattering cross section
    double *cExt;               // cExt: Extinction cross section
    double *cBack;              // cBack: Backscattering cross section
    double *SizePara;           // SizePara: Size Parameter
    double *mus;                // mus: scattering coefficient
    double *g;                  // g: average cosine of phase function
    std::complex<double> **S1;  // S1: amplitude matrix component S1
    std::complex<double> **S2;  // S2: amplitude matrix component S2
    double *forward;            // Total intensity from forward hemisphere
    double *backward;           // Total intensity from backward hemisphere

    double minPolarPtheta;      // minR: Minimum radius for polar plot
    double maxPolarPtheta;      // rMax: Maximum radius for polar plot

    double fRay;                //fRay: fitting parameter fRay
    double bMie;                //bMie: fitting parameter bMie
    double fittedA;             //FittedA:  A value
    double muspFittingError;    //muspfittingError: fitting error

	//mie parameter
	int Disperse;  //MonoDisperse , PolyDisperse
	int distIndex ;   //Log_Normal, Guassian, Junge(plan distribution), CustomData
	int m_type ;                         //Conc_mm3, VolFrac
	const char* custom_psd ;
	const char* custom_spf;
	const char* custom_cdf;
	Phase_DTheta DTheta;


	const char* custom_in;
	const char* custom_out;
	// seawater parameters
	double* mua0;
	double* mus0;
	double g0;
	double ffn;
	double ffmu;
	double seaindex;
	double airindex;
	double* muat;
	double seadepth;
	double fseabott; //the reflection coefficient of the sea bottom
	int fseabotttype;
	const char* seabotfile;
	spftyp st;
	icdftyp it;
	bio_or_opt bo;
	double interface[4];
	double* i_m_rnd;
	double* i_angle;
	double* imie_m_rnd;
	double* imie_angle;
	int i_len;
	double* s_m_rnd;
	double* s_angle;
	int s_len;

	//lidar system parameters
	double planeheight;
	double viewfield;
	double receivearea;

	int Nphotons;
	int ns_per_bin; //the digitizing interleaves
	int BssArrayLength;
	double inc_ang;
	double light_R;
	int light_typ;
	double windspeed;
	double THRESHOLD;

	//output parameters
	double* bss;
	double* bss1;
	double* bss2;
	double* bss3;
	double* bss3plus;

	//photon parameters
	double x;
	double y;
	double z;
	double ux;
	double uy;
	double uz;
	STATUS status;
	HIT hit;
	double totalpath;
	double sleft;
	double weight;
	int collision;
	double steplength;

};

#endif // PARAMETERS_H
