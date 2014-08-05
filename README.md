MatPsi: Quantum Chemistry in Matlab based on PSI4
======

##Compilation

configure, make, find useable stuffs in $objdir/MatPsi

Refer to PSI4 compilation when in trouble

##Usage 

Usage is (nearly) all of Matlab convention.

    >> [output1, output2, ...] = matpsi.WhatEverFunction(input1, input2, ...);

This work is based on _Example MATLAB class wrapper for a C++ class_ by Oliver Woodford. 
See http://www.mathworks.com/matlabcentral/fileexchange/38964-example-matlab-class-wrapper-for-a-c++-class for more details. 

####Constructor 

Usually, we construct a MatPsi object using 2 strings, one describes the molecule's geometry (as well as charge and spin), 
and one sets the name of the basis set we use. 

```
mol_string = 

O
H 1 R
H 1 R 2 A

R = .9
A = 104.5

basis_name = 

6-31g*

>> matpsi = MatPsi(mol_string, basis_name);
```

####Copy Constructor

Construct from an existing MatPsi object. 

    >> matpsi2 = matpsi.MatPsiCopy();

####Constructor related properties
1. molecule_string: Output the string used in constructing the MatPsi object. 

    ```
    >> matpsi.molecule_string(); 
    ```

2. basis_name: Output the basis set name used in constructing the MatPsi object. 

    ```
    >> matpsi.basis_name(); 
    ```

####Molecule operations 
1. fix_mol: Fix the molecule at its current frame, i.e. disable re-orientation, moving to center-of-mass, etc. 

    ```
    >> matpsi.fix_mol(); 
    ```

2. free_mol: Free the molecule, i.e. re-enable re-orientation, moving to center-of-mass, etc. 

    ```
    >> matpsi.free_mol(); 
    ```

####Molecule properties 
1. natom: Number of atoms. 

    ```
    >> matpsi.natom(); 
    ```

2. nelec: Number of electrons. 

    ```
    >> matpsi.nelec(); 
    ```

3. geom: Geometry or Cartesian coordinates of all atoms. __Unit is a.u. (i.e. Bohr radius) rather than Angstrom!__  

    ```
    >> matpsi.geom(); 
    ```

4. set_geom: Set geometry. Requires an natom by 3 matrix as input. __Unit is a.u. (i.e. Bohr radius) rather than Angstrom!__  

    ```
    >> matpsi.set_geom(newGeom); 
    ```

5. Enuc:  Nuclear repulsion energy. 

    ```
    >> matpsi.Enuc(); 
    ```

####Basis set properties 

1. nbasis: Number of basis functions. 

    ```
    >> matpsi.nbasis(); 
    ```

2. func2center: A vector mapping basis function number to the number of atom it is centred on, a.k.a "basisAtom" in MSQC. 

    ```
    >> matpsi.func2center(); 
    ```

3. func2am(): A vector mapping basis function number to the corresponding angular momentum number, a.k.a "basisType" in MSQC. 

    ```
    >> matpsi.func2am(); 
    ```

####One-electron integrals 

1. overlap: Atomic orbital overlap matrix (S). 

    ```
    >> matpsi.overlap(); 
    ```

2. kinetic: Kinetic energy matrix (KE). 

    ```
    >> matpsi.kinetic(); 
    ```

3. potential: One-electron potential energy matrix (EN). Summed over all atoms. 

    ```
    >> matpsi.potential(); 
    ```

4. potential_sep: Atom-separated One-electron potential energy matrices (ENI). Return a 3-D, (nbasis by nbasis by natom) array. 

    ```
    >> matpsi.potential_sep(); 
    ```

5. potential_Zxyz: Environment potential energy matrix (ENVI) for a list of point charges in the format of ZxyzMatrix = [Z1 x1 y1 z1; Z2 x2 y2 z2; ...]; Zi stands for 
charge magnitude; xi, yi, zi are Cartesian coordinates of point charges. Use for the summed environment potential matrix. 

    ```
    >> matpsi.potential_Zxyz(ZxyzMatrix); 
    ```

6. dipole: Dipole integrals. __Unit is a.u. (i.e. Bohr radius) rather than Angstrom!__ 

    ```
    >> [x y z] = matpsi.dipole(); 
    ```

####Two-electron integrals 

1. tei_ijkl: 4-indexed two-electron interaction integral. Need four indices as input arguments. Return only one integral value. 

    ```
    >> matpsi.tei_ijkl(i, j, k, l); 
    ```

2. tei_uniqN: Total number of unique two-electron interaction integrals. No geometrical symmetry is considered. 

    ```
    >> matpsi.tei_uniqN(); 
    ```

3. tei_alluniq: All unique two-electron interaction integrals. No geometrical symmetry is considered. Be careful as it consumes 
a huge amount of memory. 

    ```
    >> matpsi.tei_alluniq(); 
    ```

4. tei_alluniqJK: Still all unique two-electron interaction integrals but ordered and summed in quick forming of 
the exchange energy matrix K. See Developer's Note for detailed discussions. 

    ```
    >> [Jvec, Kvec] = matpsi.tei_alluniqJK(); 
    ```

5. tei_allfull: Full, 4-D two-electron interaction integrals. No geometrical symmetry is considered. Be careful as it consumes 
a really huge amount of memory. 

    ```
    >> matpsi.tei_allfull(); 
    ```

####SCF related 

1. OccMO2J: For restricted Hartree Fock theory, form the 2-electron Coulomb interaction J matrix from occupied molecular orbital coefficients matrix. 
Direct algorithm. No geometrical symmetry is considered. 

    ```
    >> matpsi.OccMO2J(OccMO); 
    ```

2. OccMO2K: For restricted Hartree Fock theory, form the 2-electron exchange interaction K matrix from occupied molecular orbital coefficients matrix. 
Direct algorithm. No geometrical symmetry is considered. 

    ```
    >> matpsi.OccMO2K(OccMO); 
    ```

3. OccMO2G: For restricted Hartree Fock theory, form the 2-electron G = 2 * J - K matrix from occupied molecular orbital coefficients matrix. 
Direct algorithm. No geometrical symmetry is considered. 

    ```
    >> matpsi.OccMO2G(OccMO); 
    ```

4. RHF: Solve restricted Hartree-Fock functions and returns the final Hartree-Fock energy 

    ```
    >> matpsi.RHF(); 
    ```

5. RHFenv: Solve restricted Hartree-Fock functions with a given environment potential energy matrix EnvMat and returns the final Hartree-Fock energy. 
EnvMat must be an nbasis by nbasis Hermitian(symmetric) matrix. 

    ```
    >> matpsi.RHFenv(EnvMat); 
    ```

6. RHF_finalize: Free some memory and disk space required by ```RHF()```. Right now executed automatically thus no need to call. 

    ```
    >> matpsi.RHF_finalize(); 
    ```

7. RHF_EHF: Get the final restricted Hartree-Fock energy. Executable after RHF. 

    ```
    >> matpsi.RHF_EHF(); 
    ```

8. RHF_C: Restricted Hartree-Fock molecular orbital coefficients. Executable after RHF. 

    ```
    >> matpsi.RHF_C(); 
    ```

9. RHF_EMO: Restricted Hartree-Fock molecular orbital energies (eigenvalues). Executable after RHF. 

    ```
    >> matpsi.RHF_EMO(); 
    ```

10. RHF_D: Restricted Hartree-Fock density matrix. Executable after RHF. 

    ```
    >> matpsi.RHF_D(); 
    ```

11. RHF_H: One-electron (core) Hamiltonian matrix. Executable after RHF. 

    ```
    >> matpsi.RHF_H(); 
    ```

12. RHF_J: 2-electron Coulomb interaction J matrix. Executable after RHF. 

    ```
    >> matpsi.RHF_J(); 
    ```

13. RHF_K: 2-electron exchange interaction K matrix. Executable after RHF. 

    ```
    >> matpsi.RHF_K(); 
    ```

14. RHF_F: Fock matrix. Executable after RHF. 

    ```
    >> matpsi.RHF_F(); 
    ```

##TODO 

remove python dependence

incorejk

##Developer's Note 

Mar. 29: To make two-electron integrals work, the original PSI4 source code has been slightly changed. Below is where and why. 

From "$psi4dir/src/lib/libmins/twobody.h", we know that the virtual class `TwoBodyAOInt` has some properties probably related with 
Psi4's python interface: 

    class TwoBodyAOInt
    {
    protected:
    ...
        /// The PyBuffer object used for sharing the target_ buffer without copying data
        PyBuffer<double> target_pybuffer_;
        /// Whether or not to use the PyBuffer
        bool enable_pybuffer_;
    ...
    }; // class end

and in "$psi4dir/src/lib/libmins/twobody.cc", `TwoBodyAOInt` class virtual constructor: 

    TwoBodyAOInt::TwoBodyAOInt(const IntegralFactory* intsfactory, int deriv) :
        integral_(intsfactory),
        original_bs1_(integral_->basis1()),
        original_bs2_(integral_->basis2()),
        original_bs3_(integral_->basis3()),
        original_bs4_(integral_->basis4()),
        deriv_(deriv),
        target_(0),
        target_pybuffer_(&target_, true)

the last line shown, `target_pybuffer_(&target_, true)`, is known causing some strange segmentation fault errors. 
Eliminating this line fixes it, but as a result, we probably need to set `bool enable_pybuffer_` as `false` forever. 

__(Jul. 27 added) Deleted all things related with `pybuffer` from "eribase.cc", "twobody.cc", "twobody.h", and "export_mints.cc".__

Apr. 02: A new method, void remove_symmetry(), has been added to JK class. It allows us to get rid of the geometrical symmetry automatically 
imposed (but not used in real computation for some weird reasons) by JK constructor. 

    class JK {
        ...
    public:
        ...
        void remove_symmetry() {
            AO2USO_ = SharedMatrix(new Matrix("AO->SO matrix", primary_->nbf(), primary_->nbf()));
            AO2USO_->identity();
        }
    };

Apr. 05: About forming the exchange energy matrix K 

In order to (efficiently) form K, we loop over all unique two-electron integrals and store them as: 

```
I = i*(i-1)/2 + j;
J = k*(k-1)/2 + l;
Jvec( I*(I-1)/2 + J ) = H2(i, j, k, l);
Kvec( I*(I-1)/2 + J ) = H2(i, l, k, j) + H2(i, k, j, l);
```

where i >= j, k >= l, and I >= J. Then we can form K from teiVecForK and the density matrix in the same manner as forming J, that is, 
restore this vector into an (upper case) N by N (N = n*(n+1)/2, n is the number of basis functions) Hermitian matrix, contract 
the density matrix into an N-length column vector and scale it by C(J) = 2 or C(K) = -1/2 and scale again all of the diagonal elements by 1/2, 
then time the N by N matrix by this density vector. In the real algorithm we can avoid the "restore", or Hermitianize, step and 
save some memory as well as sacrifice some speed, though. 

Apr. 11: To correctly input/output matrices, we have to be aware of the following. 

If we are going to use ```double* pointer; pointer++;``` to loop over all elements of a matrix, the pointer get from the Matlab mex function 
```double* mxGetPr(mxArray*)``` will first go through a column, i.e. mat(1, 1)->mat(2, 1)->mat(3, 1)... etc.; but the C++ double* pointer will 
first go through a row, i.e. mat[0][0]->mat[0][1]->mat[0][2]... etc. See ```SharedMatrix InputMatrix(mxArray*& Mat_m)``` in MatPsi_mex.cpp 
for an appropriate way to deal with that. 

Apr. 13: ```<libscf_solver/hf.h>``` and ```<libscf_solver/rhf.h>``` are slightly changed, with some matrix wiping 
turned off and ```finalize()``` method wrapped in ```public extern_finalize()```, as well as canceling a python external field interface. 


