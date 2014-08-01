#include <libmints/mints.h>
#include <libfock/jk.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include "mex.h"
#include "class_handle.hpp"
#include <boost/shared_array.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libscf_solver/rhf.h>
#include <read_options.cc>

using namespace std;
using namespace psi;
using namespace boost;

namespace psi {
    //FILE* outfile = fopen("/dev/null", "w");
    FILE* outfile = stdout;
    char* psi_file_prefix = "matpsi";
    std::string outfile_name = "";
    extern int read_options(const std::string &name, Options & options, bool suppress_printing = false);
}

class MatPsi {
protected:
    //~ Options options_; // for now, just use an independent Options object; 
                      // in the future may be we should use independent Enviroment object as well 
    Process::Environment process_environment_;
    boost::shared_ptr<worldcomm> worldcomm_;
    
    std::string fakepid_; // serves as a fake pid 
    boost::shared_ptr<PSIO> psio_;
    
    std::string path_;
    std::string molstring_;
    std::string basisname_;
    
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
	boost::shared_ptr<IntegralFactory> intfac_;
    boost::shared_ptr<TwoBodyAOInt> eri_;
	boost::shared_ptr<MatrixFactory> matfac_;
    boost::shared_ptr<DirectJK> directjk_;
    boost::shared_ptr<scf::RHF> rhf_;
    
    // common initializer for constructors 
    void common_init(std::string molstring, std::string basisname, int ncores = 4, unsigned long int memory = 1000000000L);
    
    // create basis object 
    void create_basis();
    
    // create one & two electron integral factories and directjk object 
    void create_integral_factories();
    
    // initialize the directjk object 
    void init_directjk(double cutoff = 1.0E-12);
    
public:
    // constructor; takes in 2 strings and parse them 
    MatPsi(std::string molstring, std::string basisname, std::string path);
	
    // destructor 
	~MatPsi();
    
    // global molecule is shared among all instances, this method is for debugging 
    void testmol() {Process::environment.molecule()->print();}
    
    // enable DirectJK object 
    void UseDirectJK(double cutoff = 1.0E-12);
    
    // the string describing the molecule 
    std::string molecule_string() { return molstring_; }
    
    // basis set name string 
    std::string basis_name() { return basisname_; }
    
    // Molecule operations 
    // fix the molecule
    void fix_mol();
    
    // free the molecule
    void free_mol();
    
    // Molecule properties 
    // number of atoms 
    int natom() { return molecule_->natom(); }
    
    // geometry 
    SharedMatrix geom() { return molecule_->geometry().clone(); }
    
    // set geometry 
    void set_geom(SharedMatrix newGeom);
    
    // nuclear repulsion energy 
    double Enuc() { return molecule_->nuclear_repulsion_energy(); }
    
    // Z list 
    SharedVector Zlist();
    
    // number of electrons 
    int nelec();
    
    // Basis set properties 
    // number of basis functions 
    int nbasis() { return basis_->nbf(); }
    
    // map basis number to the number of atom it is centred on 
    SharedVector func2center();
    
    // map basis number to its angular momentum 
    SharedVector func2am();
    
	// One-electron integrals 
    // compute the overlap matrix S 
	SharedMatrix overlap();
    
    // compute the kinetic energy matrix KE 
    SharedMatrix kinetic();
    
    // compute the total potential energy matrix EN 
    SharedMatrix potential();
    
    // compute the dipole integrals 
    std::vector<SharedMatrix> dipole();
    
    // compute the atom-separated potential energy matrix ENI 
    boost::shared_array<SharedMatrix> potential_sep();
    
    // compute from a given point charge list the environment potential energy matrix ENVI
    SharedMatrix potential_Zxyz(SharedMatrix Zxyz_list);
    
    // Two-electron integrals 
    // compute the 4-indexed two-electron integral H2(i, j, k, l) 
    double tei_ijkl(int i, int j, int k, int l);
    
    // number of unique two-electron integrals 
    int tei_uniqN();
    
    // compute all unique two-electron integrals and put them in a vector; be careful as it costs a huge amount of memory 
    void tei_alluniq(double* matpt);
    
    // compute all, full, nbasis by nbasis by nbasis by nbasis two-electron integrals and put them in a vector; be careful as it costs a super huge amount of memory 
    void tei_allfull(double* matpt);
    
    // compute all unique two-electron integrals and pre-arrange them for the forming of J and K 
    void tei_alluniqJK(double* matptJ, double* matptK);
    
    // SCF related 
    
    SharedMatrix Density2J(SharedMatrix Density);
    
    // for restricted Hartree Fock, compute 2-electron Coulomb interaction J matrix from occupied molecular orbital coefficient matrix, direct algorithm, consider no geometrical symmetry 
    SharedMatrix OccMO2J(SharedMatrix OccMO);
    
    // for restricted Hartree Fock, compute 2-electron exchange interaction K matrix from occupied molecular orbital coefficient matrix, direct algorithm, consider no geometrical symmetry 
    SharedMatrix OccMO2K(SharedMatrix OccMO);
    
    // for restricted Hartree Fock, compute 2-electron G matrix from occupied molecular orbital coefficient matrix, direct algorithm, consider no geometrical symmetry 
    SharedMatrix OccMO2G(SharedMatrix OccMO);
    
    // restricted Hartree-Fock; quick but uses our filesystem and causes a lot of risky issues 
    double RHF();
    
    // restricted Hartree-Fock with environment potential; quick but uses our filesystem and causes a lot of risky issues 
    double RHF(SharedMatrix EnvMat);
    
    // release memory and clean temporary files for Hartree-Fock 
    void RHF_finalize();
    
    // restricted Hartree-Fock energy 
    double RHF_EHF();
    
    // restricted Hartree-Fock molecular orbitals 
    SharedMatrix RHF_C();
    
    // restricted Hartree-Fock molecular orbital energies 
    SharedVector RHF_EMO();
    
    // restricted Hartree-Fock density matrix 
    SharedMatrix RHF_D();
    
    // restricted Hartree-Fock one-electron (core) Hamiltonian matrix 
    SharedMatrix RHF_Ha();
    
    // restricted Hartree-Fock two-electron Coulomb matrix 
    SharedMatrix RHF_J();
    
    // restricted Hartree-Fock two-electron exchange matrix 
    SharedMatrix RHF_K();
    
    // restricted Hartree-Fock Fock matrix 
    SharedMatrix RHF_F();
    
    void switch_worldcomm() {
        WorldComm = worldcomm_;
    }
    
};
