#ifndef ENZO_ENZO_INITIAL_ISOLATED_GALAXY_HPP
#define ENZO_ENZO_INITIAL_ISOLATED_GALAXY_HPP

//#define USE_NANOFLANN

#ifdef USE_NANOFLANN
  //#include "nanoflann.hpp"

struct Point
{
  double x,y,z;

/*
  inline void pup(PUP::er &p){
    p | x;
    p | y;
    p | z;
  }

  PUPable_decl(Point);
*/
};

struct PointCloud
{

  std::vector<Point>  pts;

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return pts.size(); }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const
  {
          if (dim == 0) return pts[idx].x;
          else if (dim == 1) return pts[idx].y;
          else return pts[idx].z;
  }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

/*
  inline void pup(PUP::er &p){
    p | pts;
  }

  PUPable_decl(PointCloud);
*/
};


inline void operator|(PUP::er &p, Point &pt){
  p | pt.x;
  p | pt.y;
  p | pt.z;
}

inline void operator|(PUP::er &p, PointCloud &ptc){
  p | ptc.pts;
}


#endif


class EnzoInitialIsolatedGalaxy : public Initial {


private:

  int  ntypes_;
  int  nparticles_;
  int  ndim_;
  // Dimensions for pos / vel array :
  //    1        : particle type
  //    2        : dimension (x,y,z)
  //    3        : particle(s)
  //    mass is 2D and has just dimensions 1 and 3

#ifdef USE_NANOFLANN
  std::vector<PointCloud> particlePointCloud;
#else
  enzo_float *** particleIcPosition;
#endif

  enzo_float *** particleIcVelocity;
  enzo_float **  particleIcMass;
  enzo_float **  particleIcLifetime;
  enzo_float **  particleIcCreationTime;
  int * particleIcTypes;
  std::vector<std::string> particleIcFileNames;

  // Utility deallocation routine for particle ICs
  void allocateParticles(void){

    if (ntypes_ == 0){
#ifndef USE_NANOFLANN
      particleIcPosition     = NULL;
#endif
      particleIcVelocity     = NULL;
      particleIcTypes        = NULL;
      particleIcMass         = NULL;
      particleIcLifetime     = NULL;
      particleIcCreationTime = NULL;
      return;
    }

#ifdef USE_NANOFLANN
    particlePointCloud.resize(ntypes_);
#else
    particleIcPosition  = new enzo_float **[ntypes_];
#endif
    particleIcVelocity  = new enzo_float **[ntypes_];
    particleIcMass      = new enzo_float  *[ntypes_];
    particleIcCreationTime = new enzo_float *[ntypes_];
    particleIcLifetime     = new enzo_float *[ntypes_];

    if (!(particleIcTypes)) particleIcTypes = new int[ntypes_];

    for (int k = 0; k < ntypes_; k++){
#ifndef USE_NANOFLANN
      particleIcPosition[k] = new enzo_float*[ndim_];
#endif
      particleIcVelocity[k]  = new enzo_float*[ndim_];

      for(int j = 0; j < ndim_; j++){
#ifndef USE_NANOFLANN
        particleIcPosition[k][j] = new enzo_float[nparticles_];
#endif
        particleIcVelocity[k][j] = new enzo_float[nparticles_];

        for (int i = 0; i < nparticles_; i++){
#ifndef USE_NANOFLANN
          particleIcPosition[k][j][i] = -1.0;
#endif
          particleIcVelocity[k][j][i] = -1.0;
        }
      }

      particleIcMass[k]         = new enzo_float[nparticles_];
      particleIcCreationTime[k] = new enzo_float[nparticles_];
      particleIcLifetime[k]     = new enzo_float[nparticles_];
      for (int i = 0; i < nparticles_; i++){
        particleIcMass[k][i] = -1.0;
        particleIcCreationTime[k][i] = 0.0;
        particleIcLifetime[k][i]     = -999999.0;
      }
    }
    return;
  }

  void freeParticles(void){

    for (int k=0; k < ntypes_; k++){
      for (int j = 0; j < ndim_; j ++){
#ifndef USE_NANOFLANN
        delete [] particleIcPosition[k][j];
#endif
        delete [] particleIcVelocity[k][j];
      }
#ifndef USE_NANOFLANN
      particleIcPosition[k] = NULL;
#endif
      particleIcVelocity[k] = NULL;

      delete [] particleIcMass;
      delete [] particleIcLifetime;
      delete [] particleIcCreationTime;
    }

    delete [] particleIcTypes;

#ifdef USE_NANOFLANN
    particlePointCloud.clear();
    //particlePointCloud.erase();
#else
    particleIcPosition = NULL;
#endif
    particleIcVelocity = NULL;
    particleIcMass     = NULL;
    particleIcTypes    = NULL;
    particleIcLifetime = NULL;
    particleIcCreationTime = NULL;

    return;
  }


public: // interface

  /// CHARM++ constructor
  EnzoInitialIsolatedGalaxy(const EnzoConfig * enzo_config) throw();

  // do I need to init things here????
  PUPable_decl(EnzoInitialIsolatedGalaxy); // do i need this??? AE

  /// Charm++ PUP::able migration constructor (AE unsure of this?)
  EnzoInitialIsolatedGalaxy (CkMigrateMessage *m)
    : Initial (m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize block
  virtual void enforce_block
  ( Block * block,
    const Hierarchy * hierarchy ) throw();

  void InitializeParticles(Block * block,
                           Particle * particle);

  void InitializeExponentialGasDistribution(Block * block);
  void InitializeGasFromParticles(Block * block);

  void ReadParticles(void);

  /// Read in particle data (DM and stars)
  void ReadParticlesFromFile_(const int&nl, const int& ipt);

  void ReadParticlesFromFile(const int& nl,
                      #ifdef USE_NANOFLANN
                             const int& pcindex,
                      #else
                             enzo_float *position[],
                      #endif
                             enzo_float *velocity[],
                             enzo_float *mass, const std::string& filename);

  /// Read in circular velocity table
  void    ReadInVcircData(void);

  /// Intarpolate circular velocity table
  double InterpolateVcircTable(double radius);

  /// Compute cell mass with gaussian interpolation
  double gauss_mass(const double rho_zero,
                    const double xpos,
                    const double ypos,
                    const double zpos,
                    const double dx);

  /// Destructor
  virtual ~EnzoInitialIsolatedGalaxy(void) throw() {
    if (particleIcVelocity) freeParticles();
  };

private: // attributes

  double center_position_[3];
  double bfield_[3];
  double scale_length_;
  double scale_height_;
  double disk_mass_;
  double gas_fraction_;
  double disk_temperature_;
  double disk_metal_fraction_;
  double gas_halo_mass_;
  double gas_halo_temperature_;
  double gas_halo_metal_fraction_;

  double gas_halo_density_;
  double gas_halo_radius_;

  int dual_energy_;
  double uniform_density_;
  double gamma_;
  double mu_;
  double tiny_number_;

  bool analytic_velocity_;

  bool use_gas_particles_;
  bool live_dm_halo_;
  bool stellar_bulge_;
  bool stellar_disk_;

  const int VCIRC_TABLE_LENGTH = 10000;

  double vcirc_radius[10000];
  double vcirc_velocity[10000];

  bool include_recent_SF;
  double recent_SF_start;
  double recent_SF_end;
  double recent_SF_bin_size;
  double recent_SF_SFR;
  int recent_SF_seed;

};


#endif /* ENZO_ENZO_INITIAL_ISOLATED_GALAXY_HPP */
