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
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace psi;
using namespace boost;

namespace psi {
#ifdef PSIDEBUG
    FILE* outfile = stdout;
#else
    FILE* outfile = fopen("/dev/null", "w");
#endif
    char* psi_file_prefix = "matpsi";
    std::string outfile_name = "";
    extern int read_options(const std::string &name, Options & options, bool suppress_printing = false);
}

class MatPsi {
protected:

    int ncores_;
    unsigned long int memory_;
    
    std::string molstring_;
    std::string basisname_;
    std::string path_;

    Process::Environment process_environment_;
    boost::shared_ptr<worldcomm> worldcomm_;
    
    std::string matpsi_id;
    std::string matpsi_tempdir_str;
    boost::shared_ptr<PSIO> psio_;
    
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
    boost::shared_ptr<IntegralFactory> intfac_;
    boost::shared_ptr<TwoBodyAOInt> eri_;
    boost::shared_ptr<MatrixFactory> matfac_;
    boost::shared_ptr<JK> jk_;
    boost::shared_ptr<scf::RHF> rhf_;
    
    // common initializer for constructors 
    void common_init();
    
    // create basis object 
    void create_basis();
    
    // create one & two electron integral factories and jk object 
    void create_integral_factories();
    
public:
    // constructor; takes in 2 strings and parse them 
    MatPsi(const std::string& molstring, const std::string& basisname, int ncores, const std::string& memory_str, const std::string& path);
    
    // destructor 
    ~MatPsi();
    
    // the string describing the molecule 
    std::string& molecule_string() { return molstring_; }
    
    // basis set name string 
    std::string& basis_name() { return basisname_; }
    
    void set_basis(const std::string& basisname);
    
    void set_ncores(int ncores);
    
    void set_memory(std::string);
    
    //*** Molecule operations 
    // fix the molecule
    void fix_mol();
    
    // free the molecule 
    void free_mol();
    
    
    //*** Molecule properties 
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
    
    
    //*** Basis set properties 
    // number of basis functions 
    int nbasis() { return basis_->nbf(); }
    
    // map basis number to the number of atom it is centred on 
    SharedVector func2center();
    
    // map basis number to its angular momentum 
    SharedVector func2am();
    
    
    //*** One-electron integrals 
    // compute the overlap matrix S 
    SharedMatrix overlap();
    
    // compute the kinetic energy matrix KE 
    SharedMatrix kinetic();
    
    // compute the total potential energy matrix EN 
    SharedMatrix potential();
    
    // compute the dipole integrals 
    std::vector<SharedMatrix> dipole();
    
    // compute the atom-separated potential energy matrix ENI 
    //~ boost::shared_array<SharedMatrix> potential_sep();
    std::vector<SharedMatrix> potential_sep();
    
    // compute from a given point charge list the environment potential energy matrix ENVI
    SharedMatrix potential_Zxyz(SharedMatrix Zxyz_list);
    
    
    //*** Two-electron integrals 
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
    
    
    //*** JK related
    // enables different types of JK  
    void UseDirectJK();
    void UsePKJK();
    
    // for restricted Hartree Fock, compute 2-electron Coulomb interaction J matrix from density matrix, consider no geometrical symmetry 
    SharedMatrix Density2J(SharedMatrix Density);
    
    // for restricted Hartree Fock, compute 2-electron exchange interaction K matrix from density matrix, consider no geometrical symmetry 
    SharedMatrix Density2K(SharedMatrix Density);
    
    // for restricted Hartree Fock, compute 2-electron total interaction G matrix from density matrix, consider no geometrical symmetry 
    SharedMatrix Density2G(SharedMatrix Density);
    
    // for restricted Hartree Fock, compute 2-electron Coulomb interaction J matrix from occupied molecular orbital coefficient matrix, consider no geometrical symmetry 
    SharedMatrix OccMO2J(SharedMatrix OccMO);
    
    // for restricted Hartree Fock, compute 2-electron exchange interaction K matrix from occupied molecular orbital coefficient matrix, consider no geometrical symmetry 
    SharedMatrix OccMO2K(SharedMatrix OccMO);
    
    // for restricted Hartree Fock, compute 2-electron exchange interaction G matrix from occupied molecular orbital coefficient matrix, consider no geometrical symmetry 
    SharedMatrix OccMO2G(SharedMatrix OccMO);
    
    
    //*** SCF related
    // RHF engine for MSQC 
    double RHF_msqc(SharedMatrix given_H_in, SharedMatrix Jmodifier_in, SharedMatrix Kmodifier_in);
    
    // create/reset RHF object 
    void RHF_reset();
    
    // restricted Hartree-Fock 
    double RHF();
    
    // restricted Hartree-Fock with environment potential 
    double RHFenv(SharedMatrix EnvMat);
    
    // restricted Hartree-Fock starting from a given density matrix 
    double RHF_fromC(SharedMatrix C_in);
    
    // restricted Hartree-Fock with environment potential starting from a given density matrix 
    double RHFenv_fromC(SharedMatrix EnvMat, SharedMatrix C_in);
    
    // enable MOM in restricted Hartree-Fock to solve convergence issue 
    void RHF_EnableMOM(int mom_start);
    
    // enable Damping in restricted Hartree-Fock to solve convergence issue 
    void RHF_EnableDamping(double damping_percentage);
    
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
    
    
    //*** used at the beginning of mex file to let a global pointer pointing to a MatPsi class member property 
    void switch_worldcomm() { WorldComm = worldcomm_; }
    
};
