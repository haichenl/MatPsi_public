#include "mex.h"
#include "class_handle.hpp"
#include "MatPsi.h"

using namespace std;
using namespace psi;
using namespace boost;

SharedMatrix InputMatrix(const mxArray*& Mat_m) {
    int nrow = mxGetM(Mat_m);
    int ncol = mxGetN(Mat_m);
    SharedMatrix Mat_c(new Matrix(nrow, ncol));
    double* Mat_m_pt = mxGetPr(Mat_m);
    double** Mat_c_pt = Mat_c->pointer();
    for(int i = 0; i < ncol; i++)
        for(int j = 0; j < nrow; j++)
            Mat_c_pt[j][i] = *Mat_m_pt++; // Matlab loops over a column first, but C++ loops over a row first 
    return Mat_c;
}

double InputScalar(const mxArray*& Mat_m) {
    double* Mat_m_pt = mxGetPr(Mat_m);
    return *Mat_m_pt;
}

void OutputMatrix(mxArray*& Mat_m, SharedMatrix Mat_c) {
    int nrow = Mat_c->nrow();
    int ncol = Mat_c->ncol();
    Mat_m = mxCreateDoubleMatrix( nrow, ncol, mxREAL);
    double* Mat_m_pt = mxGetPr(Mat_m);
    double** Mat_c_pt = Mat_c->pointer();
    for(int i = 0; i < ncol; i++)
        for(int j = 0; j < nrow; j++)
            *Mat_m_pt++ = Mat_c_pt[j][i];
}

void OutputVector(mxArray*& Mat_m, SharedVector Vec_c) {
    int dim = Vec_c->dim();
    Mat_m = mxCreateDoubleMatrix( 1, dim, mxREAL);
    double* Mat_m_pt = mxGetPr(Mat_m);
    double* Vec_c_pt = Vec_c->pointer();
    for(int i = 0; i < dim; i++) {
        *Mat_m_pt++ = *Vec_c_pt++;
    }
}

void OutputScalar(mxArray*& Mat_m, double scalar) {
    Mat_m = mxCreateDoubleMatrix( 1, 1, mxREAL);
    double* Mat_m_pt = mxGetPr(Mat_m);
    *Mat_m_pt = scalar;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    
    // Constructor
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("MatPsi Constructor: One output expected.");
        if ( (nrhs!=4 && nrhs!=6) || !mxIsChar(prhs[1]) || !mxIsChar(prhs[2]))
            mexErrMsgTxt("MatPsi Constructor: MatPsi(mol_string, basis_name) input expected.");
        if (nrhs == 4){
            plhs[0] = convertPtr2Mat<MatPsi>(new MatPsi((std::string)mxArrayToString(prhs[1]) , (std::string)mxArrayToString(prhs[2]), 1, "1000mb", (std::string)mxArrayToString(prhs[3]) + "/share/"));
            return;
        }
        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<MatPsi>(new MatPsi((std::string)mxArrayToString(prhs[1]) , (std::string)mxArrayToString(prhs[2]), (int)InputScalar(prhs[3]), (std::string)mxArrayToString(prhs[4]), (std::string)mxArrayToString(prhs[5]) + "/share/" ));
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
        mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<MatPsi>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    MatPsi* MatPsi_obj = convertMat2Ptr<MatPsi>(prhs[1]);
    
    //*** Call the various class methods 
    
    // molecule_string 
    if (!strcmp("molecule_string", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->molecule_string()).c_str());
        return;
    }
    
    // basis_name 
    if (!strcmp("basis_name", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->basis_name()).c_str());
        return;
    }
    
    // set_basis
    if (!strcmp("set_basis", cmd)) {
        if ( nrhs!=3 || !mxIsChar(prhs[2]))
            mexErrMsgTxt("set_basis(\"basis\"): String input expected.");
        MatPsi_obj->set_basis((std::string)mxArrayToString(prhs[2]));
        return;
    }
    
    // set_ncores 
    if (!strcmp("set_ncores", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->set_ncores(1);
            return;
        }
        if (nrhs!=3 ||  mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("set_ncores(ncores): Integer input expected.");
        MatPsi_obj->set_ncores((int)InputScalar(prhs[2]));
        return;
    }
    
    // set_memory 
    if (!strcmp("set_memory", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->set_memory("1000mb");
            return;
        }
        if (nrhs!=3 || !mxIsChar(prhs[2]))
            mexErrMsgTxt("set_memory(\"memory\"): String input expected.");
        MatPsi_obj->set_memory((std::string)mxArrayToString(prhs[2]));
        return;
    }
    
    
    //*** Molecule operations 
    // fix_mol 
    if (!strcmp("fix_mol", cmd)) {
        MatPsi_obj->fix_mol();
        return;
    }
    
    // free_mol 
    if (!strcmp("free_mol", cmd)) {
        MatPsi_obj->free_mol();
        return;
    }
    
    
    //*** Molecule properties 
    // natom 
    if (!strcmp("natom", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->natom());
        return;
    }
    
    // nelec 
    if (!strcmp("nelec", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->nelec());
        return;
    }
    
    // geom 
    if (!strcmp("geom", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->geom());
        return;
    }
    
    // set_geom 
    if (!strcmp("set_geom", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != MatPsi_obj->natom() || mxGetN(prhs[2]) != 3)
            mexErrMsgTxt("set_geom(newGeom): natom by 3 matrix input expected.");
        // Call the method
        MatPsi_obj->set_geom(InputMatrix(prhs[2]));
        return;
    }
    
    // Zlist 
    if (!strcmp("Zlist", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->Zlist());
        return;
    }
    
    // Enuc
    if (!strcmp("Enuc", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->Enuc());
        return;
    }
    
    
    //*** Basis set properties 
    // nbf
    if (!strcmp("nbasis", cmd) || !strcmp("nbf", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->nbf());
        return;
    }
    
    // nshell
    if (!strcmp("nshell", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->nshell());
        return;
    }
    
    // shellTypes 
    if (!strcmp("shellTypes", cmd)) {
        SharedVector shellTypesVec = MatPsi_obj->shellTypes();
        OutputVector(plhs[0], shellTypesVec);
        return;
    }
    
    // shellNprims 
    if (!strcmp("shellNprims", cmd)) {
        SharedVector shellNprimsVec = MatPsi_obj->shellNprims();
        OutputVector(plhs[0], shellNprimsVec);
        return;
    }
    
    // func2center 
    if (!strcmp("func2center", cmd)) {
        SharedVector func2centerVec = MatPsi_obj->func2center();
        for(int i = 0; i < func2centerVec->dim(); i++)
            func2centerVec->add(i, 1.0); // + 1 convert C++ convention to Matlab convention 
        OutputVector(plhs[0], func2centerVec);
        return;
    }
    
    // func2am 
    if (!strcmp("func2am", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->func2am());
        return;
    }
    
    // primExps 
    if (!strcmp("primExps", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->primExps());
        return;
    }
    
    // primCoefs 
    if (!strcmp("primCoefs", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->primCoefs());
        return;
    }
    
    
    //*** One-electron integrals 
    // overlap(nbasis, nbasis)
    if (!strcmp("overlap", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->overlap());
        return;
    }
    
    // kinetic 
    if (!strcmp("kinetic", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->kinetic());
        return;
    }
    
    // potential 
    if (!strcmp("potential", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->potential());
        return;
    }
    
    // potential_sep 
    if (!strcmp("potential_sep", cmd)) {
        //~ boost::shared_array<SharedMatrix> viMatArray = MatPsi_obj->potential_sep();
        std::vector<SharedMatrix> viMatArray = MatPsi_obj->potential_sep();
        int ncol = viMatArray[0]->ncol();
        int nrow = viMatArray[0]->nrow();
        int natom = MatPsi_obj->natom();
        mwSize dims[3] = {ncol, nrow, natom};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        double* matlab_pt = mxGetPr(plhs[0]);
        for(int iatom = 0; iatom < natom; iatom++) {
            double* tmp_pt = viMatArray[iatom]->get_pointer();
            for(int i = 0; i < ncol * nrow; i++) {
                matlab_pt[iatom*ncol*nrow + i] = tmp_pt[i];
            }
        }
        return;
    }
    
    // potential_Zxyz 
    if (!strcmp("potential_Zxyz", cmd)) {
        // Check parameters
        if (nrhs!=3)
            mexErrMsgTxt("potential_Zxyz(Zxyz_mat): (number of point charges) by 4 matrix input expected.");
        if (mxGetN(prhs[2]) != 4)
            mexErrMsgTxt("potential_Zxyz: Zxyz list matrix dimension does not agree.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->potential_Zxyz(InputMatrix(prhs[2])));
        return;
    }
    
    // dipole 
    if (!strcmp("dipole", cmd)) {
        std::vector<SharedMatrix> dipole = MatPsi_obj->dipole();
        if(nlhs == 3) {
            OutputMatrix(plhs[0], dipole[0]);
            OutputMatrix(plhs[1], dipole[1]);
            OutputMatrix(plhs[2], dipole[2]);
            return;
        }
        int ncol = dipole[0]->ncol();
        int nrow = dipole[0]->nrow();
        mwSize dims[3] = {ncol, nrow, 3};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        double* matlab_pt = mxGetPr(plhs[0]);
        for(int idim = 0; idim < 3; idim++) {
            double* tmp_pt = dipole[idim]->get_pointer();
            for(int i = 0; i < ncol * nrow; i++) {
                matlab_pt[idim*ncol*nrow + i] = tmp_pt[i];
            }
        }
        return;
    }
    
    
    //*** Two-electron integrals 
    // tei_ijkl
    if (!strcmp("tei_ijkl", cmd)) {
        // Check parameters
        if (nrhs!=6 || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]) || !mxIsDouble(prhs[5]))
            mexErrMsgTxt("tei_ijkl(i, j, k, l): 4 integers input expected.");
        // Call the method
        int ind[4];
        for(int i = 0; i < 4; i++) {
            ind[i] = (int)mxGetScalar(prhs[2+i]) - 1; // -1 convert Matlab convention to C++ convention
        if(ind[i] < 0 || ind[i] >= MatPsi_obj->nbf())
                mexErrMsgTxt("tei_ijkl: Required index not within scale.");
        }
        OutputScalar(plhs[0], MatPsi_obj->tei_ijkl(ind[0], ind[1], ind[2], ind[3]));
        return;
    }
    
    // tei_uniqN 
    if (!strcmp("tei_uniqN", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->tei_uniqN());
        return;
    }
    
    // tei_alluniq 
    if (!strcmp("tei_alluniq", cmd)) {
        plhs[0] = mxCreateDoubleMatrix( 1, MatPsi_obj->tei_uniqN(), mxREAL);
        double* matpt = mxGetPr(plhs[0]);
        MatPsi_obj->tei_alluniq(matpt);
        return;
    }
    
    // tei_allfull 
    if (!strcmp("tei_allfull", cmd)) {
        mwSize dims[4] = {MatPsi_obj->nbf(), MatPsi_obj->nbf(), MatPsi_obj->nbf(), MatPsi_obj->nbf()};
        plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
        double* matpt = mxGetPr(plhs[0]);
        MatPsi_obj->tei_allfull(matpt);
        return;
    }
    
    // tei_alluniqJK 
    if (!strcmp("tei_alluniqJK", cmd)) {
        plhs[0] = mxCreateDoubleMatrix( 1, MatPsi_obj->tei_uniqN(), mxREAL);
        double* matptJ = mxGetPr(plhs[0]);
        plhs[1] = mxCreateDoubleMatrix( 1, MatPsi_obj->tei_uniqN(), mxREAL);
        double* matptK = mxGetPr(plhs[1]);
        MatPsi_obj->tei_alluniqJK(matptJ, matptK);
        return;
    }
    
    
    //*** JK related 
    // InitJK 
    if (!strcmp("InitJK", cmd)) {
        if ( nrhs!=3 || !mxIsChar(prhs[2]))
            mexErrMsgTxt("InitJK(\"jktype\"): String input expected.");
        MatPsi_obj->InitJK((std::string)mxArrayToString(prhs[2]));
        return;
    }
    
    // UseDirectJK 
    if (!strcmp("UseDirectJK", cmd)) {
        MatPsi_obj->UseDirectJK();
        return;
    }
    
    // UsePKJK 
    if (!strcmp("UsePKJK", cmd)) {
        MatPsi_obj->UsePKJK();
        return;
    }
    
    // UseICJK 
    if (!strcmp("UseICJK", cmd)) {
        MatPsi_obj->UseICJK();
        return;
    }
    
    // UseMatlabJK 
    if (!strcmp("UseMatlabJK", cmd)) {
        MatPsi_obj->UseMatlabJK();
        return;
    }
    
    // JKtype 
    if (!strcmp("JKtype", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->JKtype()).c_str());
        return;
    }
    
    // SetMatlabJK 
    if (!strcmp("SetMatlabJK", cmd)) {
        int nbf = MatPsi_obj->nbf();
        if (nrhs!=4 || mxGetM(prhs[2]) != nbf || mxGetN(prhs[2]) != nbf || mxGetM(prhs[3]) != nbf || mxGetN(prhs[3]) != nbf)
            mexErrMsgTxt("SetMatlabJK(Jcell, Kcell): 2 nbasis by nbasis cell arrays input expected.");
        boost::shared_array<double*> Jcell_ptr_in(new double* [nbf * nbf]);
        boost::shared_array<double*> Kcell_ptr_in(new double* [nbf * nbf]);
        mxArray* Jcell = (mxArray*)prhs[2];
        mxArray* Kcell = (mxArray*)prhs[3];
        for( int i = 0; i < nbf * nbf; i++) {
            mxArray* Jelement = mxGetCell(Jcell, i);
            mxArray* Kelement = mxGetCell(Kcell, i);
            if (mxGetM(Jelement) != nbf || mxGetN(Jelement) != nbf || mxGetM(Kelement) != nbf || mxGetN(Kelement) != nbf)
                mexErrMsgTxt("SetMatlabJK(Jcell, Kcell): Element not having nbasis by nbasis dimension detected.");
            Jcell_ptr_in[i] = mxGetPr(Jelement);
            Kcell_ptr_in[i] = mxGetPr(Kelement);
        }
        MatPsi_obj->SetMatlabJK(Jcell_ptr_in, Kcell_ptr_in);
        return;
    }
    
    // DisableMatlabJK 
    if (!strcmp("DisableMatlabJK", cmd)) {
        MatPsi_obj->DisableMatlabJK();
        return;
    }
    
    // Density2J 
    if (!strcmp("Density2J", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != MatPsi_obj->nbf() || mxGetN(prhs[2]) != MatPsi_obj->nbf())
            mexErrMsgTxt("Density2J(MOmat): nbasis by nbasis matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->Density2J(InputMatrix(prhs[2])));
        return;
    }
    
    // Density2K 
    if (!strcmp("Density2K", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != MatPsi_obj->nbf() || mxGetN(prhs[2]) != MatPsi_obj->nbf())
            mexErrMsgTxt("Density2K(MOmat): nbasis by nbasis matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->Density2K(InputMatrix(prhs[2])));
        return;
    }
    
    // Density2G 
    if (!strcmp("Density2G", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != MatPsi_obj->nbf() || mxGetN(prhs[2]) != MatPsi_obj->nbf())
            mexErrMsgTxt("Density2G(MOmat): nbasis by nbasis matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->Density2G(InputMatrix(prhs[2])));
        return;
    }
    
    // OccMO2J 
    if (!strcmp("OccMO2J", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != MatPsi_obj->nbf())
            mexErrMsgTxt("OccMO2J(MOmat): nbasis by noccupy matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->OccMO2J(InputMatrix(prhs[2])));
        return;
    }
    
    // OccMO2K 
    if (!strcmp("OccMO2K", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != MatPsi_obj->nbf())
            mexErrMsgTxt("OccMO2K(MOmat): nbasis by noccupy matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->OccMO2K(InputMatrix(prhs[2])));
        return;
    }
    
    // OccMO2G 
    if (!strcmp("OccMO2G", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != MatPsi_obj->nbf())
            mexErrMsgTxt("OccMO2G(MOmat): nbasis by noccupy matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->OccMO2G(InputMatrix(prhs[2])));
        return;
    }
    
    
    //*** SCF related 
    // RHF_msqc 
    if (!strcmp("RHF_msqc", cmd)) {
        // Check parameters
        if (nrhs!=5)
            mexErrMsgTxt("RHF_msqc(H, Jmod, Kmod): nbasis by nbasis matrices input expected.");
        SharedMatrix given_H_in = InputMatrix(prhs[2]);
        SharedMatrix Jmodifier_in = InputMatrix(prhs[3]);
        SharedMatrix Kmodifier_in = InputMatrix(prhs[4]);
        // matrix dimension checks will be done inside 
        // Call the method
        OutputScalar(plhs[0], MatPsi_obj->RHF_msqc(given_H_in, Jmodifier_in, Kmodifier_in));
        return;
    }
    
    // RHF_msqc_fromC 
    if (!strcmp("RHF_msqc_fromC", cmd)) {
        // Check parameters
        if (nrhs!=6)
            mexErrMsgTxt("RHF_msqc_fromC(H, Jmod, Kmod, Cin): nbasis by nbasis matrices input expected.");
        SharedMatrix given_H_in = InputMatrix(prhs[2]);
        SharedMatrix Jmodifier_in = InputMatrix(prhs[3]);
        SharedMatrix Kmodifier_in = InputMatrix(prhs[4]);
        SharedMatrix C_in = InputMatrix(prhs[5]);
        // matrix dimension checks will be done inside 
        // Call the method
        OutputScalar(plhs[0], MatPsi_obj->RHF_msqc_fromC(given_H_in, Jmodifier_in, Kmodifier_in, C_in));
        return;
    }
    
    // RHF
    if (!strcmp("RHF", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->RHF());
        return;
    }
    
    // RHFenv 
    if (!strcmp("RHFenv", cmd)) {
        // Check parameters
        if (nrhs!=3)
            mexErrMsgTxt("RHFenv(EnvMat): nbasis by nbasis matrix input expected.");
        SharedMatrix EnvMat = InputMatrix(prhs[2]);
        if(EnvMat->nirrep() != 1 || EnvMat->nrow() != MatPsi_obj->nbf() || EnvMat->ncol() != MatPsi_obj->nbf())
            mexErrMsgTxt("RHFenv(EnvMat): Environment potential energy matrix dimensions do not agree.");
        // Call the method
        OutputScalar(plhs[0], MatPsi_obj->RHFenv(EnvMat));
        return;
    }
    
    // RHF_fromC 
    if (!strcmp("RHF_fromC", cmd)) {
        // Check parameters
        if (nrhs!=3)
            mexErrMsgTxt("RHF_fromC(MO): nbasis by nbasis matrix input expected.");
        SharedMatrix C_in = InputMatrix(prhs[2]);
        if(C_in->nirrep() != 1 || C_in->nrow() != MatPsi_obj->nbf() || C_in->ncol() != MatPsi_obj->nbf())
            mexErrMsgTxt("RHF_fromC(MO): MO matrix dimensions do not agree.");
        // Call the method
        OutputScalar(plhs[0], MatPsi_obj->RHF_fromC(C_in));
        return;
    }
    
    // RHFenv_fromC 
    if (!strcmp("RHFenv_fromC", cmd)) {
        // Check parameters
        if (nrhs!=4)
            mexErrMsgTxt("RHFenv_fromC(EnvMat, MO): nbasis by nbasis matrix input expected.");
        SharedMatrix EnvMat = InputMatrix(prhs[2]);
        SharedMatrix C_in = InputMatrix(prhs[3]);
        if(EnvMat->nirrep() != 1 || EnvMat->nrow() != MatPsi_obj->nbf() || EnvMat->ncol() != MatPsi_obj->nbf()
            || C_in->nirrep() != 1 || C_in->nrow() != MatPsi_obj->nbf() || C_in->ncol() != MatPsi_obj->nbf())
            mexErrMsgTxt("RHFenv_fromC(EnvMat, MO): MO matrix dimensions do not agree.");
        // Call the method
        OutputScalar(plhs[0], MatPsi_obj->RHFenv_fromC(EnvMat, C_in));
        return;
    }
    
    // RHF_reset 
    if (!strcmp("RHF_reset", cmd)) {
        MatPsi_obj->RHF_reset();
        return;
    }
    
    // RHF_EnableMOM 
    if (!strcmp("RHF_EnableMOM", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->RHF_EnableMOM(20);
            return;
        }
        if (nrhs!=3 || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("RHF_EnableMOM(mom_start): Integer input expected.");
        MatPsi_obj->RHF_EnableMOM((int)InputScalar(prhs[2]));
        return;
    }
    
    // RHF_DisableMOM 
    if (!strcmp("RHF_DisableMOM", cmd)) {
        MatPsi_obj->RHF_EnableMOM(0);
        return;
    }
    
    // RHF_EnableDamping 
    if (!strcmp("RHF_EnableDamping", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->RHF_EnableDamping(20.0);
            return;
        }
        if (nrhs!=3 || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("RHF_EnableDamping(damping_percentage): 1 double input expected.");
        MatPsi_obj->RHF_EnableDamping(InputScalar(prhs[2]));
        return;
    }
    
    // RHF_DisableDamping 
    if (!strcmp("RHF_DisableDamping", cmd)) {
        MatPsi_obj->RHF_EnableDamping(0.0);
        return;
    }
    
    // RHF_EnableDIIS 
    if (!strcmp("RHF_EnableDIIS", cmd)) {
        MatPsi_obj->RHF_EnableDIIS();
        return;
    }
    
    // RHF_DisableDIIS 
    if (!strcmp("RHF_DisableDIIS", cmd)) {
        MatPsi_obj->RHF_DisableDIIS();
        return;
    }
    
    // RHF_GuessSAD 
    if (!strcmp("RHF_GuessSAD", cmd)) {
        MatPsi_obj->RHF_GuessSAD();
        return;
    }
    
    // RHF_GuessCore 
    if (!strcmp("RHF_GuessCore", cmd)) {
        MatPsi_obj->RHF_GuessCore();
        return;
    }
    
    // RHF_EHF 
    if (!strcmp("RHF_EHF", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->RHF_EHF());
        return;
    }
    
    // RHF_C 
    if (!strcmp("RHF_C", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->RHF_C());
        return;
    }
    
    // RHF_EMO 
    if (!strcmp("RHF_EMO", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->RHF_EMO());
        return;
    }
    
    // RHF_D 
    if (!strcmp("RHF_D", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->RHF_D());
        return;
    }
    
    // RHF_H 
    if (!strcmp("RHF_H", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->RHF_Ha());
        return;
    }
    
    // RHF_J 
    if (!strcmp("RHF_J", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->RHF_J());
        return;
    }
    
    // RHF_K 
    if (!strcmp("RHF_K", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->RHF_K());
        return;
    }
    
    // RHF_F 
    if (!strcmp("RHF_F", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->RHF_F());
        return;
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}

