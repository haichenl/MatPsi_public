// Please directly include this file in MatPsi_mex.cpp to make compilation easier... 

std::string tempname() {
    // an ugly trick to get a temp file name.. subject to change 
    char* trickbuffer = tempnam("/tmp/", "tmp");
    std::string tempname_ = std::string(trickbuffer + 5); 
    delete trickbuffer;
    return tempname_;
}

// Constructor
MatPsi::MatPsi(std::string molstring, std::string basisname, std::string path) {
    fakepid_ = tempname();
    path_ = path + "/share/";
    common_init(molstring, basisname);
}

void MatPsi::common_init(std::string molstring, std::string basisname, int ncores, unsigned long int memory) {
    // some necessary initializations
    process_environment_.initialize();
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    switch_worldcomm();
    process_environment_.options.set_read_globals(true);
    
    read_options("", process_environment_.options, true);
    process_environment_.options.set_read_globals(false);
    process_environment_.set("PSIDATADIR", path_);
    //~ process_environment_.options.set_current_module("MatPsi");
    process_environment_.options.set_global_int("MAXITER", 100);
    
    Wavefunction::initialize_singletons();
    
    // initialize psio 
    psio_ = boost::shared_ptr<PSIO>(new PSIO);
    psio_->set_pid(fakepid_);
    
    // create molecule object and set its basis set name 
    molstring_ = molstring;
    basisname_ = basisname;
    molecule_ = psi::Molecule::create_molecule_from_string(process_environment_, molstring);
    molecule_->set_basis_all_atoms(basisname);
    
    process_environment_.set_molecule(molecule_);
    
    // set cores and memory 
    process_environment_.set_n_threads(ncores); // these values are shared among all instances, so better be constants 
    process_environment_.set_memory(memory);
    
    // create basis object and one & two electron integral factories 
	create_basis();
    create_integral_factories();
    
    // create matrix factory object 
    int nbf[] = { basis_->nbf() };
    matfac_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory);
    matfac_->init_with(1, nbf, nbf);
      
}

void MatPsi::create_basis() {
    // create basis object 
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    basis_ = BasisSet::construct(process_environment_, parser, molecule_, "BASIS");  
    boost::shared_ptr<PointGroup> c1group(new PointGroup("C1"));
    molecule_->set_point_group(c1group); // creating basis set object change molecule's point group, for some reasons 
}

void MatPsi::create_integral_factories() {
    // create integral factory object 
    intfac_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basis_, basis_, basis_, basis_));
    
    // create two electron integral generator
    eri_ = boost::shared_ptr<TwoBodyAOInt>(intfac_->eri());
}

void MatPsi::UseDirectJK(double cutoff) {
    // create directJK object
    directjk_ = boost::shared_ptr<DirectJK>(new DirectJK(process_environment_, basis_));
    directjk_->set_cutoff(cutoff);
    directjk_->initialize();
    directjk_->remove_symmetry();
}

// destructor 
MatPsi::~MatPsi() {
    if(rhf_ != NULL)
        RHF_finalize();
    if(directjk_ != NULL)
        directjk_->finalize();
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
    create_basis();
    create_integral_factories();
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
        create_basis();
        create_integral_factories();
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

boost::shared_array<SharedMatrix> MatPsi::potential_sep() {
    int natom_ = molecule_->natom();
    boost::shared_array<SharedMatrix> viMatArray(new SharedMatrix [natom_]);
    boost::shared_ptr<OneBodyAOInt> viOBI(intfac_->ao_potential());
    boost::shared_ptr<PotentialInt> viPtI = boost::static_pointer_cast<PotentialInt>(viOBI);
    SharedMatrix Zxyz = viPtI->charge_field();
    SharedMatrix Zxyz_rowi(new Matrix(1, 4));
    for( int i = 0; i < natom_; i++) {
        SharedVector Zxyz_rowi_vec = Zxyz->get_row(0, i);
        Zxyz_rowi->set_row(0, 0, Zxyz_rowi_vec);
        viPtI->set_charge_field(Zxyz_rowi);
        viMatArray[i] = matfac_->create_shared_matrix("potential_sep");
        viOBI->compute(viMatArray[i]);
        viMatArray[i]->hermitivitize();
    }
    return viMatArray;
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

void MatPsi::init_directjk(double cutoff) {
    directjk_->set_cutoff(cutoff);
    directjk_->initialize();
    directjk_->remove_symmetry();
}

SharedMatrix MatPsi::Density2J(SharedMatrix Density) {
    if(directjk_ == NULL) {
        throw PSIEXCEPTION("Density2J: DirectJK object has not been enabled.");
    }
    
    std::vector<SharedMatrix>& D = directjk_->D();
    std::vector<SharedMatrix>& C_left = directjk_->C_left();
    C_left.clear();
    D.clear();
    D.push_back(Density);
    
    directjk_->compute_from_D();
    SharedMatrix Jnew = directjk_->J()[0];
    Jnew->hermitivitize();
    //~ directjk_->finalize();
    return Jnew;
}

SharedMatrix MatPsi::OccMO2J(SharedMatrix OccMO) {
    if(directjk_ == NULL) {
        throw PSIEXCEPTION("OccMO2J: DirectJK object has not been enabled.");
    }
    std::vector<SharedMatrix>& C_left = directjk_->C_left();
    C_left.clear();
    C_left.push_back(OccMO);
    directjk_->compute();
    SharedMatrix Jnew = directjk_->J()[0];
    Jnew->hermitivitize();
    //~ directjk_->finalize();
    return Jnew;
}

SharedMatrix MatPsi::OccMO2K(SharedMatrix OccMO) {
    if(directjk_ == NULL) {
        throw PSIEXCEPTION("OccMO2K: DirectJK object has not been enabled.");
    }
    std::vector<SharedMatrix>& C_left = directjk_->C_left();
    C_left.clear();
    C_left.push_back(OccMO);
    directjk_->compute();
    SharedMatrix Knew = directjk_->K()[0];
    Knew->hermitivitize();
    //~ directjk_->finalize();
    return Knew;
}

SharedMatrix MatPsi::OccMO2G(SharedMatrix OccMO) {
    if(directjk_ == NULL) {
        throw PSIEXCEPTION("OccMO2G: DirectJK object has not been enabled.");
    }
    std::vector<SharedMatrix>& C_left = directjk_->C_left();
    C_left.clear();
    C_left.push_back(OccMO);
    directjk_->compute();
    SharedMatrix Gnew = directjk_->J()[0];
    Gnew->scale(2);
    Gnew->subtract(directjk_->K()[0]); // 2 J - K 
    Gnew->hermitivitize();
    //~ directjk_->finalize();
    return Gnew;
}

double MatPsi::RHF() {
    boost::shared_ptr<PointGroup> c1group(new PointGroup("C1"));
    molecule_->set_point_group(c1group); // for safety 
    
    //~ Process::environment.options.set_global_int("MOM_START", 20);
    //~ Process::environment.options.set_global_double("DAMPING_PERCENTAGE", 20.0);
    
    rhf_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, process_environment_.options, psio_));
    process_environment_.set_wavefunction(rhf_);
    try {
        double Ehf = rhf_->compute_energy();
        rhf_->J()->scale(0.5);
        RHF_finalize();
        return Ehf;
    }
    catch (...) {
        RHF_finalize();
        //~ rhf_.reset();
        throw PSIEXCEPTION("RHF: Hartree-Fock possibly not converged.");
    }
}

double MatPsi::RHF(SharedMatrix EnvMat) {
    boost::shared_ptr<PointGroup> c1group(new PointGroup("C1"));
    molecule_->set_point_group(c1group); // for safety 
    rhf_ = boost::shared_ptr<scf::RHF>(new scf::RHF(process_environment_, process_environment_.options, psio_));
    process_environment_.set_wavefunction(rhf_);
    try {
        double Ehf = rhf_->compute_energy(EnvMat);
        rhf_->J()->scale(0.5);
        RHF_finalize();
        return Ehf;
    }
    catch (...) {
        RHF_finalize();
        //~ rhf_.reset();
        throw PSIEXCEPTION("RHF(env): Hartree-Fock possibly not converged.");
    }
    
}

void MatPsi::RHF_finalize() {
    rhf_->extern_finalize();
    psio_->_psio_manager_->psiclean();
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

