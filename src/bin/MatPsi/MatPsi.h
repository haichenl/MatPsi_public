#include <libmints/mints.h>
#include <libfock/jk.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <boost/shared_array.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libscf_solver/rhf.h>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>



using namespace std;
using namespace psi;
using namespace boost;

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
    
    // create basis object and one & two electron integral factories 
    void create_basis_and_integral_factories();
    
public:
    // constructor; takes in 2 strings and parse them 
    MatPsi(const std::string& molstring, const std::string& basisname, 
        int ncores, const std::string& memory_str, const std::string& path);
    
    // destructor 
    virtual ~MatPsi();
    
    // Construcing properties 
    std::string& molecule_string() { return molstring_; } // the string describing the molecule 
    std::string& basis_name() { return basisname_; } // basis set name string 
    
    // CPU and memory controll 
    void set_ncores(int ncores);
    void set_memory(std::string);
    
    
    
    //*** Molecule properties 
    int natom() { return molecule_->natom(); } // number of atoms 
    int nelec(); // number of electrons 
    SharedMatrix geom() { return molecule_->geometry().clone(); } // geometry 
    void set_geom(SharedMatrix newGeom); // set a new geometry 
    double Enuc() { return molecule_->nuclear_repulsion_energy(); } // nuclear repulsion energy 
    SharedVector Zlist(); // protonic number list vector 
    
    //*** Molecule operations 
    void fix_mol();
    void free_mol();
    
    
    //*** Basis set properties 
    void set_basis(const std::string& basisname); // set a new basis set 
    int nbf() { return basis_->nbf(); } // number of basis functions 
    int nshell() { return basis_->nshell(); }
    SharedVector shellTypes();
    SharedVector shellNprims();
    SharedVector func2center(); // map basis function number to the number of atom it is centred on 
    SharedVector func2am(); // map basis function number to its angular momentum 
    SharedVector primExps();
    SharedVector primCoefs();
    
    
    //*** One-Electron Integrals 
    SharedMatrix overlap(); // overlap matrix S <i|j>
    SharedMatrix kinetic(); // kinetic energy matrix KE 
    SharedMatrix potential(); // total potential energy matrix EN <i|sum(1/R)|j>
    std::vector<SharedMatrix> dipole(); // dipole matrices <i|x|j>, <i|y|j>, <i|z|j>
    std::vector<SharedMatrix> potential_sep(); // atom-separated EN 
    SharedMatrix potential_Zxyz(SharedMatrix Zxyz_list); // compute from a given point charge list 
                                                         // the environment potential energy matrix ENVI 
    
    
    //*** Two-Electron Integrals (TEIs) 
    int tei_uniqN(); // number of unique TEIs 
    double tei_ijkl(int i, int j, int k, int l); // (ij|kl), chemist's notation 
    // ## HIGH MEMORY COST METHODS ## 
    void tei_alluniq(double* matpt); // all unique TEIs in a vector 
    void tei_allfull(double* matpt); // all (repetitive) TEIs in a 4D-array 
    void tei_alluniqJK(double* matptJ, double* matptK); // pre-arrange TEI vectors for forming J/K 
    // ## HIGH MEMORY COST METHODS ## 
    
    
    //*** JK related
    // use different types of JK 
    void InitJK(std::string jktype);
    void UseDirectJK();
    void UsePKJK();
    void UseICJK();
    const std::string& JKtype();
    // ### EXPERT ### 
    void UseMatlabJK();
    void SetMatlabJK(boost::shared_array<double*> Jcell_ptr, boost::shared_array<double*> Kcell_ptr);
    void DisableMatlabJK();
    // ### EXPERT ### 
    
    // methods computing J/K/G 
    SharedMatrix Density2J(SharedMatrix Density);
    SharedMatrix Density2K(SharedMatrix Density); 
    SharedMatrix Density2G(SharedMatrix Density);
    SharedMatrix OccMO2J(SharedMatrix OccMO);
    SharedMatrix OccMO2K(SharedMatrix OccMO);
    SharedMatrix OccMO2G(SharedMatrix OccMO);
    
    
    //*** SCF related
    // create/reset RHF object 
    void RHF_reset();
    
    // RHF engine for MSQC 
    double RHF_msqc(SharedMatrix given_H_in, SharedMatrix Jmodifier_in, 
        SharedMatrix Kmodifier_in);
    double RHF_msqc_fromC(SharedMatrix given_H_in, SharedMatrix Jmodifier_in, 
        SharedMatrix Kmodifier_in, SharedMatrix C_in);
    
    // method of doing RHF calculations 
    double RHF(); // "regular" restricted Hartree-Fock 
    double RHFenv(SharedMatrix EnvMat);  // with an environment potential 
    double RHF_fromC(SharedMatrix C_in); // starting from a given molecular orbital matrix 
    double RHFenv_fromC(SharedMatrix EnvMat, SharedMatrix C_in);
    
    // methods controlling RHF algorithm 
    void RHF_EnableMOM(int mom_start);
    void RHF_EnableDamping(double damping_percentage);
    void RHF_EnableDIIS();
    void RHF_DisableDIIS();
    void RHF_GuessSAD();
    void RHF_GuessCore();
    
    // methods extracting restricted Hartree-Fock results
    double RHF_EHF();       // restricted Hartree-Fock energy 
    SharedMatrix RHF_C();   // molecular orbital coefficients
    SharedVector RHF_EMO(); // molecular orbital engenvalues 
    SharedMatrix RHF_D();   // density matrix 
    SharedMatrix RHF_Ha();  // core Hamiltonian 
    SharedMatrix RHF_J();   // Coulomb interaction matrix J 
    SharedMatrix RHF_K();   // exchange interaction matrix K
    SharedMatrix RHF_F();   // entire Fock matrix 
    
};
