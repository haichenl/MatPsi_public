
#include "MatPsi.h"
#include "matlabjk.hpp"
#include <read_options.cc>

namespace psi {
#ifdef PSIDEBUG
    FILE* outfile = stdout;
#else
    FILE* outfile = fopen("/dev/null", "w");
#endif
    char* psi_file_prefix = "matpsi";
    std::string outfile_name = "";
    extern int read_options(const std::string &name, Options & options, 
            bool suppress_printing = false);
}


unsigned long int parse_memory_str(const std::string& memory_str) {
    std::string memory_str_ = memory_str;
    boost::algorithm::to_lower(memory_str_);
    boost::cmatch cm;
    boost::regex_search(memory_str_.c_str(), cm, boost::regex("([+]?[0-9]*.?[0-9]+)"));
    double memory = boost::lexical_cast<double>(std::string(cm[1].first, cm[1].second));
    boost::regex_search(memory_str_.c_str(), cm, boost::regex("([a-z])"));
    std::string unit_str(cm[1].first, cm[1].second);
    unsigned long int unit;
    if(boost::iequals(unit_str, "b"))
        unit = 1L;
    else if(boost::iequals(unit_str, "k"))
        unit = 1000L;
    else if(boost::iequals(unit_str, "m"))
        unit = 1000000L;
    else
        unit = 1000000000L;
    if((unsigned long int)(memory*unit) < 100000000L) // less than 100 mb 
        return 100000000L; // if less than 100mb then return 100 mb 
    else
        return (unsigned long int)(memory*unit);
}

// Constructor
MatPsi::MatPsi(const std::string& molstring, const std::string& basisname, int ncores, const std::string& memory_str, const std::string& path)
    : molstring_(molstring), basisname_(basisname), path_(path)
{
    ncores_ = ncores;
    memory_ = parse_memory_str(memory_str);
    common_init();
}

void MatPsi::common_init() {
    // some necessary initializations
    process_environment_.initialize();
    
    // set cores and memory 
    process_environment_.set_n_threads(ncores_);
    process_environment_.set_memory(memory_);
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    switch_worldcomm();
    
    // read in options 
    process_environment_.options.set_read_globals(true);
    read_options("", process_environment_.options, true);
    process_environment_.options.set_read_globals(false);
    process_environment_.set("PSIDATADIR", path_);
    process_environment_.options.set_global_int("MAXITER", 100);
    
    Wavefunction::initialize_singletons();
    
    // initialize psio 
    boost::filesystem::path uniqname = boost::filesystem::unique_path();
    matpsi_id = uniqname.string();
    boost::filesystem::path tempdir = boost::filesystem::temp_directory_path();
    matpsi_tempdir_str = tempdir.string();
    matpsi_tempdir_str += "/matpsi.temp.";
    matpsi_tempdir_str += matpsi_id;
    //~ boost::filesystem::create_directories(matpsi_tempdir_str);
    psio_ = boost::shared_ptr<PSIO>(new PSIO);
    psio_->set_pid(matpsi_id);
    for (int i=1; i<=PSIO_MAXVOL; ++i) {
        char kwd[20];
        sprintf(kwd, "VOLUME%u", i);
        psio_->filecfg_kwd("DEFAULT", kwd, PSIF_CHKPT, matpsi_tempdir_str.c_str());
        psio_->filecfg_kwd("DEFAULT", kwd, -1, matpsi_tempdir_str.c_str());
    }
    psio_->_psio_manager_->set_default_path(matpsi_tempdir_str);
    
    // create molecule object and set its basis set name 
    molecule_ = psi::Molecule::create_molecule_from_string(process_environment_, molstring_);
    molecule_->set_basis_all_atoms(basisname_);
    process_environment_.set_molecule(molecule_);
    
    // create basis object and one & two electron integral factories & rhf 
    create_basis_and_integral_factories();
    RHF_reset();
    rhf_->extern_finalize(); // after this finalize() rhf_ seems to be stable and does not crash after RHF() or set_basis() 
    rhf_.reset();
    
    // create matrix factory object 
    int nbf[] = { basis_->nbf() };
    matfac_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory);
    matfac_->init_with(1, nbf, nbf);
}

void MatPsi::set_basis(const std::string& basisname) {
    if(jk_ != NULL)
        jk_->finalize();
    psio_->_psio_manager_->psiclean();
    jk_.reset();
    
    basisname_ = basisname;
    molecule_->set_basis_all_atoms(basisname_);
    
    // create basis object and one & two electron integral factories & rhf 
    create_basis_and_integral_factories();
    RHF_reset();
    rhf_->extern_finalize();
    rhf_.reset();
    
    // create matrix factory object 
    int nbf[] = { basis_->nbf() };
    matfac_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory);
    matfac_->init_with(1, nbf, nbf);
}

void MatPsi::create_basis() {
    // create basis object 
    boost::shared_ptr<PointGroup> c1group(new PointGroup("C1"));
    molecule_->set_point_group(c1group); 
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    basis_ = BasisSet::construct(process_environment_, parser, molecule_, "BASIS");  
    
    molecule_->set_point_group(c1group); // creating basis set object change molecule's point group, for some reasons 
}

void MatPsi::create_basis_and_integral_factories() {
    
    create_basis();
    
    // create integral factory object 
    intfac_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basis_, basis_, basis_, basis_));
    
    // create two electron integral generator
    eri_ = boost::shared_ptr<TwoBodyAOInt>(intfac_->eri());
}

// destructor 
MatPsi::~MatPsi() {
    if(rhf_ != NULL)
        rhf_->extern_finalize();
    if(jk_ != NULL)
        jk_->finalize();
    psio_->_psio_manager_->psiclean();
    boost::filesystem::remove_all(matpsi_tempdir_str);
}

void MatPsi::set_ncores(int ncores) {
    // set cores and update worldcomm_ 
    ncores_ = ncores;
    process_environment_.set_n_threads(ncores_);
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    switch_worldcomm();
}

void MatPsi::set_memory(std::string memory_str) {
    // set memory and update worldcomm_ 
    memory_ = parse_memory_str(memory_str);
    process_environment_.set_memory(memory_);
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    switch_worldcomm();
}

void MatPsi::fix_mol() {
    molecule_->set_orientation_fixed();
    molecule_->set_com_fixed();
    molecule_->set_reinterpret_coordentry(false);
}

void MatPsi::free_mol() {
    // done by re-generating a new molecule object (and basis object etc.) while retaining the geometry 
    SharedMatrix oldgeom = geom();
    molecule_ = psi::Molecule::create_molecule_from_string(process_environment_, molstring_);
    molecule_->set_basis_all_atoms(basisname_);
    process_environment_.set_molecule(molecule_);
    create_basis(); // the molecule object isn't complete before we create the basis object, according to psi4 documentation 
    molecule_->set_geometry(*(oldgeom.get()));
    create_basis_and_integral_factories();
    RHF_reset();
    rhf_->extern_finalize();
    rhf_.reset();
}

void MatPsi::set_geom(SharedMatrix newGeom) {
    
    // store the old geometry
    Matrix oldgeom = molecule_->geometry();
    molecule_->set_geometry(*(newGeom.get()));
    
    // determine whether the new geometry will cause a problem (typically 2 atoms are at the same point) 
    Matrix distmat = molecule_->distance_matrix();
    bool nonbreak = true;
    for(int i = 0; i < natom() && nonbreak; i++) {
        for(int j = 0; j < i - 1; j++) {
            if(distmat.get(i, j) == 0) {
                molecule_->set_geometry(oldgeom);
                cout<<"set_geom: Geometry is too crazy; keeping the old one."<<endl;
                nonbreak = false;
                break;
            }
        }
    }
    
    // update other objects 
    if(nonbreak) {
        if(jk_ != NULL)
            jk_->finalize();
        psio_->_psio_manager_->psiclean();
        jk_.reset();
        create_basis_and_integral_factories();
        RHF_reset();
        rhf_->extern_finalize();
        rhf_.reset();
    }
}

SharedVector MatPsi::Zlist() {
    SharedVector zlistvec(new Vector(molecule_->natom()));
    for(int i = 0; i < molecule_->natom(); i++) {
        zlistvec->set(i, (double)molecule_->Z(i));
    }
    return zlistvec;
}

int MatPsi::nelec() {
    int charge = molecule_->molecular_charge();
    int nelectron  = 0;
    for(int i = 0; i < molecule_->natom(); i++)
        nelectron += (int)molecule_->Z(i);
    nelectron -= charge;
    return nelectron;
}

SharedVector MatPsi::func2center() {
    SharedVector func2centerVec(new Vector(basis_->nbf()));
    for(int i = 0; i < basis_->nbf(); i++) {
        func2centerVec->set(i, (double)basis_->function_to_center(i));
    }
    return func2centerVec;
}

SharedVector MatPsi::func2am() {
    SharedVector func2amVec(new Vector(basis_->nbf()));
    for(int i = 0; i < basis_->nbf(); i++) {
        func2amVec->set(i, (double)basis_->shell(basis_->function_to_shell(i)).am());
    }
    return func2amVec;
}

SharedMatrix MatPsi::overlap() {
    SharedMatrix sMat(matfac_->create_matrix("Overlap"));
    boost::shared_ptr<OneBodyAOInt> sOBI(intfac_->ao_overlap());
    sOBI->compute(sMat);
    sMat->hermitivitize();
    return sMat;
}

SharedMatrix MatPsi::kinetic() {
    SharedMatrix tMat(matfac_->create_matrix("Kinetic"));
    boost::shared_ptr<OneBodyAOInt> tOBI(intfac_->ao_kinetic());
    tOBI->compute(tMat);
    tMat->hermitivitize();
    return tMat;
}

SharedMatrix MatPsi::potential() {
    SharedMatrix vMat(matfac_->create_matrix("Potential"));
    boost::shared_ptr<OneBodyAOInt> vOBI(intfac_->ao_potential());
    vOBI->compute(vMat);
    vMat->hermitivitize();
    return vMat;
}

std::vector<SharedMatrix> MatPsi::dipole() {
    std::vector<SharedMatrix> ao_dipole;
    SharedMatrix dipole_x(matfac_->create_matrix("Dipole x"));
    SharedMatrix dipole_y(matfac_->create_matrix("Dipole y"));
    SharedMatrix dipole_z(matfac_->create_matrix("Dipole z"));
    ao_dipole.push_back(dipole_x);
    ao_dipole.push_back(dipole_y);
    ao_dipole.push_back(dipole_z);
    boost::shared_ptr<OneBodyAOInt> dipoleOBI(intfac_->ao_dipole());
    dipoleOBI->compute(ao_dipole);
    ao_dipole[0]->hermitivitize();
    ao_dipole[1]->hermitivitize();
    ao_dipole[2]->hermitivitize();
    return ao_dipole;
}

std::vector<SharedMatrix> MatPsi::potential_sep() {
    int natom_ = molecule_->natom();
    std::vector<SharedMatrix> viMatVec;
    boost::shared_ptr<OneBodyAOInt> viOBI(intfac_->ao_potential());
    boost::shared_ptr<PotentialInt> viPtI = boost::static_pointer_cast<PotentialInt>(viOBI);
    SharedMatrix Zxyz = viPtI->charge_field();
    SharedMatrix Zxyz_rowi(new Matrix(1, 4));
    for( int i = 0; i < natom_; i++) {
        SharedVector Zxyz_rowi_vec = Zxyz->get_row(0, i);
        Zxyz_rowi->set_row(0, 0, Zxyz_rowi_vec);
        viPtI->set_charge_field(Zxyz_rowi);
        viMatVec.push_back(matfac_->create_shared_matrix("potential_sep"));
        viOBI->compute(viMatVec[i]);
        viMatVec[i]->hermitivitize();
    }
    return viMatVec;
}

SharedMatrix MatPsi::potential_Zxyz(SharedMatrix Zxyz_list) {
    boost::shared_ptr<OneBodyAOInt> viOBI(intfac_->ao_potential());
    boost::shared_ptr<PotentialInt> viPtI = boost::static_pointer_cast<PotentialInt>(viOBI);
    viPtI->set_charge_field(Zxyz_list);
    SharedMatrix vZxyzListMat(matfac_->create_matrix("Potential_ZxyzList"));
    viOBI->compute(vZxyzListMat);
    vZxyzListMat->hermitivitize();
    return vZxyzListMat;
}

double MatPsi::tei_ijkl(int i, int j, int k, int l) {
    int ish = basis_->function_to_shell(i);
    int jsh = basis_->function_to_shell(j);
    int ksh = basis_->function_to_shell(k);
    int lsh = basis_->function_to_shell(l);
    int ii = i - basis_->shell_to_basis_function(ish);
    int jj = j - basis_->shell_to_basis_function(jsh);
    int kk = k - basis_->shell_to_basis_function(ksh);
    int ll = l - basis_->shell_to_basis_function(lsh);
    int ni = basis_->shell(ish).nfunction();
    int nj = basis_->shell(jsh).nfunction();
    int nk = basis_->shell(ksh).nfunction();
    int nl = basis_->shell(lsh).nfunction();
    eri_->compute_shell(ish, jsh, ksh, lsh);
    const double *buffer = eri_->buffer();
    return buffer[ll+nl*(kk+nk*(jj+nj*ii))];
}

inline int ij2I(int i, int j) {
    if(i < j) {
        int tmp = i;
        i = j;
        j = tmp;
    }
    return i * ( i + 1 ) / 2 + j;
}

int MatPsi::tei_uniqN() {
    return ( basis_->nbf() * ( basis_->nbf() + 1 ) * ( basis_->nbf() * basis_->nbf() + basis_->nbf() + 2 ) ) / 8;
}

void MatPsi::tei_alluniq(double* matpt) {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nuniq = tei_uniqN();
    const double *buffer = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            matpt[ ij2I( ij2I(intIter.i(), intIter.j()), ij2I(intIter.k(), intIter.l()) ) ] = buffer[intIter.index()];
        }
    }
}

void MatPsi::tei_allfull(double* matpt) {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nbf_ = basis_->nbf();
    const double *buffer = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            int i = intIter.i();
            int j = intIter.j();
            int k = intIter.k();
            int l = intIter.l();
            matpt[ l+nbf_*(k+nbf_*(j+nbf_*i)) ] = buffer[intIter.index()];
            matpt[ l+nbf_*(k+nbf_*(i+nbf_*j)) ] = buffer[intIter.index()];
            matpt[ k+nbf_*(l+nbf_*(j+nbf_*i)) ] = buffer[intIter.index()];
            matpt[ k+nbf_*(l+nbf_*(i+nbf_*j)) ] = buffer[intIter.index()];
            matpt[ j+nbf_*(i+nbf_*(l+nbf_*k)) ] = buffer[intIter.index()];
            matpt[ j+nbf_*(i+nbf_*(k+nbf_*l)) ] = buffer[intIter.index()];
            matpt[ i+nbf_*(j+nbf_*(l+nbf_*k)) ] = buffer[intIter.index()];
            matpt[ i+nbf_*(j+nbf_*(k+nbf_*l)) ] = buffer[intIter.index()];
        }
    }
}

void MatPsi::tei_alluniqJK(double* matptJ, double* matptK) {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    const double *buffer = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            int i = intIter.i();
            int j = intIter.j();
            int k = intIter.k();
            int l = intIter.l();
            matptJ[ij2I( ij2I(i, j), ij2I(k, l) )] = buffer[intIter.index()];
        }
    }
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            int i = intIter.i();
            int j = intIter.j();
            int k = intIter.k();
            int l = intIter.l();
            matptK[ij2I( ij2I(i, j), ij2I(k, l) )] = matptJ[ij2I( ij2I(i, l), ij2I(k, j) )] + matptJ[ij2I( ij2I(i, k), ij2I(j, l) )];
        }
    }
}

void MatPsi::UseDirectJK() {
    // create DirectJK object
    DirectJK* jk = new DirectJK(process_environment_, basis_);
    if (process_environment_.options["INTS_TOLERANCE"].has_changed())
        jk->set_cutoff(process_environment_.options.get_double("INTS_TOLERANCE"));
    if (process_environment_.options["PRINT"].has_changed())
        jk->set_print(process_environment_.options.get_int("PRINT"));
    if (process_environment_.options["DEBUG"].has_changed())
        jk->set_debug(process_environment_.options.get_int("DEBUG"));
    if (process_environment_.options["BENCH"].has_changed())
        jk->set_bench(process_environment_.options.get_int("BENCH"));
    if (process_environment_.options["DF_INTS_NUM_THREADS"].has_changed())
        jk->set_df_ints_num_threads(process_environment_.options.get_int("DF_INTS_NUM_THREADS"));
    jk_ = boost::shared_ptr<JK>(jk);
    jk_->set_memory((ULI)(process_environment_.options.get_double("SCF_MEM_SAFETY_FACTOR")*(process_environment_.get_memory() / 8L)));
    jk_->initialize();
}

void MatPsi::UsePKJK() {
    if(jk_ != NULL)
        jk_->finalize();
    if ( !boost::filesystem::exists(matpsi_tempdir_str) ) {
        boost::filesystem::create_directories(matpsi_tempdir_str);
    }
    // create PKJK object
    PKJK* jk = new PKJK(process_environment_, basis_, psio_);

    if (process_environment_.options["INTS_TOLERANCE"].has_changed())
        jk->set_cutoff(process_environment_.options.get_double("INTS_TOLERANCE"));
    if (process_environment_.options["PRINT"].has_changed())
        jk->set_print(process_environment_.options.get_int("PRINT"));
    if (process_environment_.options["DEBUG"].has_changed())
        jk->set_debug(process_environment_.options.get_int("DEBUG"));
    jk_ = boost::shared_ptr<JK>(jk);
    jk_->set_memory((ULI)(process_environment_.options.get_double("SCF_MEM_SAFETY_FACTOR")*(process_environment_.get_memory() / 8L)));
    jk_->initialize();
}

void MatPsi::UseICJK() {
    if(jk_ != NULL)
        jk_->finalize();
    // create ICJK object
    ICJK* jk = new ICJK(process_environment_, basis_);

    if (process_environment_.options["INTS_TOLERANCE"].has_changed())
        jk->set_cutoff(process_environment_.options.get_double("INTS_TOLERANCE"));
    if (process_environment_.options["PRINT"].has_changed())
        jk->set_print(process_environment_.options.get_int("PRINT"));
    if (process_environment_.options["DEBUG"].has_changed())
        jk->set_debug(process_environment_.options.get_int("DEBUG"));
    jk_ = boost::shared_ptr<JK>(jk);
    jk_->set_memory(process_environment_.get_memory());
    jk_->initialize();
}

void MatPsi::UseMatlabJK() {
    if(jk_ != NULL)
        jk_->finalize();
    // create MatlabJK object
    MatlabJK* jk = new MatlabJK(process_environment_, basis_);

    if (process_environment_.options["INTS_TOLERANCE"].has_changed())
        jk->set_cutoff(process_environment_.options.get_double("INTS_TOLERANCE"));
    if (process_environment_.options["PRINT"].has_changed())
        jk->set_print(process_environment_.options.get_int("PRINT"));
    if (process_environment_.options["DEBUG"].has_changed())
        jk->set_debug(process_environment_.options.get_int("DEBUG"));
    jk_ = boost::shared_ptr<JK>(jk);
    jk_->set_memory(process_environment_.get_memory());
    jk_->initialize();
}

const std::string& MatPsi::JKtype() {
    if(jk_ == NULL)
        throw PSIEXCEPTION("JKtype: JK object has not been initialized.");
    return jk_->JKtype();
}

void MatPsi::SetMatlabJK(boost::shared_array<double*> Jcell_ptr_in, boost::shared_array<double*> Kcell_ptr_in) {
    boost::static_pointer_cast<MatlabJK>(jk_)->set_JKcell_ptrs(Jcell_ptr_in, Kcell_ptr_in);
}

void MatPsi::DisableMatlabJK() {
    boost::static_pointer_cast<MatlabJK>(jk_)->disable();
}

SharedMatrix MatPsi::Density2J(SharedMatrix Density) {
    if(jk_ == NULL) {
        UsePKJK();
    }
    jk_->set_do_K(false);
    jk_->C_left().clear();
    jk_->D().clear();
    jk_->D().push_back(Density);
    
    jk_->compute_from_D();
    SharedMatrix Jnew = jk_->J()[0];
    Jnew->hermitivitize();
    jk_->set_do_K(true);
    return Jnew;
}

SharedMatrix MatPsi::Density2K(SharedMatrix Density) {
    if(jk_ == NULL) {
        UsePKJK();
    }
    jk_->set_do_J(false);
    jk_->C_left().clear();
    jk_->D().clear();
    jk_->D().push_back(Density);
    
    jk_->compute_from_D();
    SharedMatrix Knew = jk_->K()[0];
    Knew->hermitivitize();
    jk_->set_do_J(true);
    return Knew;
}

SharedMatrix MatPsi::Density2G(SharedMatrix Density) {
    if(jk_ == NULL) {
        UsePKJK();
    }
    
    jk_->C_left().clear();
    jk_->D().clear();
    jk_->D().push_back(Density);
    
    jk_->compute_from_D();
    SharedMatrix Gnew = jk_->J()[0];
    Gnew->scale(2.0);
    Gnew->subtract(jk_->K()[0]);
    Gnew->hermitivitize();
    return Gnew;
}

SharedMatrix MatPsi::OccMO2J(SharedMatrix OccMO) {
    if(jk_ == NULL) {
        UsePKJK();
    }
    jk_->set_do_K(false);
    jk_->C_left().clear();
    jk_->C_left().push_back(OccMO);
    jk_->compute();
    SharedMatrix Jnew = jk_->J()[0];
    Jnew->hermitivitize();
    jk_->set_do_K(true);
    return Jnew;
}

SharedMatrix MatPsi::OccMO2K(SharedMatrix OccMO) {
    if(jk_ == NULL) {
        UsePKJK();
    }
    jk_->set_do_J(false);
    jk_->C_left().clear();
    jk_->C_left().push_back(OccMO);
    jk_->compute();
    SharedMatrix Knew = jk_->K()[0];
    Knew->hermitivitize();
    jk_->set_do_J(true);
    return Knew;
}

SharedMatrix MatPsi::OccMO2G(SharedMatrix OccMO) {
    if(jk_ == NULL) {
        UsePKJK();
    }
    jk_->C_left().clear();
    jk_->C_left().push_back(OccMO);
    jk_->compute();
    SharedMatrix Gnew = jk_->J()[0];
    Gnew->scale(2.0);
    Gnew->subtract(jk_->K()[0]);
    Gnew->hermitivitize();
    return Gnew;
}

void MatPsi::RHF_reset() {
    //~ if(rhf_ != NULL)
        //~ rhf_->extern_finalize(); // really tired of trying when we should do finalize, just be it, seems to work 
    rhf_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, process_environment_.options, jk_, psio_));
    process_environment_.set_wavefunction(rhf_);
}

void MatPsi::RHF_EnableMOM(int mom_start) {
    process_environment_.options.set_global_int("MOM_START", mom_start);
}

void MatPsi::RHF_EnableDamping(double damping_percentage) {
    process_environment_.options.set_global_double("DAMPING_PERCENTAGE", damping_percentage);
}

void MatPsi::RHF_DisableDIIS() {
    process_environment_.options.set_global_bool("DIIS", false);
    process_environment_.options.set_global_int("MAXITER", 500);
}

void MatPsi::RHF_EnableDIIS() {
    process_environment_.options.set_global_int("DIIS", true);
    process_environment_.options.set_global_int("MAXITER", 100);
}

void MatPsi::RHF_GuessSAD() {
    process_environment_.options.set_global_str("GUESS", "SAD");
}

void MatPsi::RHF_GuessCore() {
    process_environment_.options.set_global_str("GUESS", "Core");
}

double MatPsi::RHF() {
    if(jk_ == NULL)
        UsePKJK();
    rhf_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, process_environment_.options, jk_, psio_));
    process_environment_.set_wavefunction(rhf_);
    return rhf_->compute_energy();
}

double MatPsi::RHFenv(SharedMatrix EnvMat) {
    if(jk_ == NULL)
        UsePKJK();
    rhf_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, process_environment_.options, jk_, psio_));
    process_environment_.set_wavefunction(rhf_);
    rhf_->set_EnvMat(EnvMat);
    return rhf_->compute_energy();
}

double MatPsi::RHF_fromC(SharedMatrix C_in) {
    if(jk_ == NULL)
        UsePKJK();
    rhf_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, process_environment_.options, jk_, psio_));
    process_environment_.set_wavefunction(rhf_);
    rhf_->set_StartingC(C_in);
    return rhf_->compute_energy();
}

double MatPsi::RHFenv_fromC(SharedMatrix EnvMat, SharedMatrix C_in) {
    if(jk_ == NULL)
        UsePKJK();
    rhf_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, process_environment_.options, jk_, psio_));
    process_environment_.set_wavefunction(rhf_);
    rhf_->set_EnvMat(EnvMat);
    rhf_->set_StartingC(C_in);
    return rhf_->compute_energy();
}

double MatPsi::RHF_msqc(SharedMatrix given_H_in, SharedMatrix Jmodifier_in, SharedMatrix Kmodifier_in) {
    if(jk_ == NULL)
        throw PSIEXCEPTION("RHF_msqc: Please explictly initialize JK first.");
    if(rhf_ == NULL || rhf_->jk() != jk_) // always use our jk_ 
        RHF_reset();
    // now we don't need to cancel given_H_ or JKmodifiers_ since we re-generate an rhf_ each time we call 
    rhf_->set_print(0);
    rhf_->set_given_H(given_H_in);
    rhf_->set_JKmodifiers(Jmodifier_in, Kmodifier_in);
    return rhf_->compute_energy_minIO();
}

double MatPsi::RHF_msqc_fromC(SharedMatrix given_H_in, SharedMatrix Jmodifier_in, SharedMatrix Kmodifier_in, SharedMatrix C_in) {
    if(jk_ == NULL)
        UseICJK();
    if(rhf_ == NULL)
        RHF_reset();
    rhf_->set_print(0);
    rhf_->set_given_H(given_H_in);
    rhf_->set_JKmodifiers(Jmodifier_in, Kmodifier_in);
    rhf_->set_StartingC(C_in);
    return rhf_->compute_energy_minIO();
}

double MatPsi::RHF_EHF() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_EHF: Hartree-Fock calculation has not been done.");
    }
    return rhf_->EHF(); 
}

SharedMatrix MatPsi::RHF_C() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_C: Hartree-Fock calculation has not been done.");
    }
    return rhf_->Ca(); 
}

SharedVector MatPsi::RHF_EMO() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_EMO: Hartree-Fock calculation has not been done.");
    }
    return rhf_->epsilon_a(); 
}

SharedMatrix MatPsi::RHF_D() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_D: Hartree-Fock calculation has not been done.");
    }
    return rhf_->Da(); 
}

SharedMatrix MatPsi::RHF_Ha() { // RHF_H leads to naming issues due to "ifdef RHF_H" or something 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_H: Hartree-Fock calculation has not been done.");
    }
    return rhf_->H(); 
}

SharedMatrix MatPsi::RHF_J() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_J: Hartree-Fock calculation has not been done.");
    }
    return rhf_->J(); 
}

SharedMatrix MatPsi::RHF_K() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_K: Hartree-Fock calculation has not been done.");
    }
    return rhf_->K(); 
}

SharedMatrix MatPsi::RHF_F() { 
    if(rhf_ == NULL) {
        throw PSIEXCEPTION("RHF_F: Hartree-Fock calculation has not been done.");
    }
    return rhf_->Fa(); 
}

