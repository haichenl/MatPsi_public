% MatPsi: An interface between Matlab and Psi4 
classdef MatPsi < handle
    properties (SetAccess = private)
        path;
    end
    properties (SetAccess = private, Hidden = true, Transient = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance  
        function this = MatPsi(varargin)
            if(exist('./share/basis', 'file'))
                this.path = pwd;
            elseif(exist('./@MatPsi/share/basis', 'file'))
                this.path = [pwd, '/@MatPsi'] ;
            else
                currpath = path();
                paths_num = length(regexp(currpath, ':', 'match')) + 1;
                for i = 1:paths_num
                    toppath = regexp(currpath, '[^:]*', 'match', 'once');
                    if(exist([toppath, '/share/basis'], 'file'))
                        this.path = toppath;
                        break;
                    elseif(exist([toppath, '/@MatPsi/share/basis'], 'file'))
                        this.path = [toppath, '/@MatPsi'];
                        break;
                    else
                        currpath = currpath(length(toppath)+2:end);
                    end
                end
                if(i>=paths_num)
                    throw(MException('MatPsi:MatPsi','MatPsi cannot find basis set files.'));
                end
            end
            if(nargin > 3 && isfloat(varargin{4}))
                varargin{4} = num2str(varargin{4});
            end
            this.objectHandle = MatPsi_mex('new', varargin{:}, this.path);
        end
        
        %% Copy Constructor
        function this2 = MatPsiCopy(this, varargin)
            this2 = MatPsi(this.molecule_string(), this.basis_name());
        end
        
        %% Destructor - Destroy the C++ class instance 
        function delete(this)
            if(~isempty(this.objectHandle))
                MatPsi_mex('delete', this.objectHandle);
            end
        end
        
        %% Constructor related properties
        % string = matpsi.molecule_string(); 
        % The molecule description string used in constructor 
        function varargout = molecule_string(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('molecule_string', this.objectHandle, varargin{:});
        end
        
        % string = matpsi.basis_name(); 
        % The basis set name used in constructor  
        function varargout = basis_name(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('basis_name', this.objectHandle, varargin{:});
        end
        
        % matpsi.set_basis(integer); 
        % Set basis set 
        function varargout = set_basis(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('set_basis', this.objectHandle, varargin{:});
        end
        
        % matpsi.set_ncores(integer); 
        % Set the number of CPU cores this MatPsi instance can use 
        function varargout = set_ncores(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('set_ncores', this.objectHandle, varargin{:});
        end
        
        % matpsi.set_memory(string); 
        % Set the amount of memory this MatPsi instance can use; default unit GB 
        function varargout = set_memory(this, varargin)
            if(isfloat(varargin{1}))
                varargin{1} = num2str(varargin{1});
            end
            [varargout{1:nargout}] = MatPsi_mex('set_memory', this.objectHandle, varargin{:});
        end
        
        %% Molecule operations 
        % matpsi.fix_mol(); 
        % Fix the molecule at its current frame 
        function varargout = fix_mol(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('fix_mol', this.objectHandle, varargin{:});
        end
        
        % matpsi.free_mol(); 
        % Free the molecule so that the program can move/reorient the molecule 
        function varargout = free_mol(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('free_mol', this.objectHandle, varargin{:});
        end
        
        %% Molecule properties 
        % integer = matpsi.natom(); 
        % Number of atoms in the molecule 
        function varargout = natom(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('natom', this.objectHandle, varargin{:});
        end
        
        % integer = matpsi.nelec(); 
        % Number of electrons in the molecule 
        function varargout = nelec(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('nelec', this.objectHandle, varargin{:});
        end
        
        % matrix(natom, 3) = matpsi.geom(); 
        % The molecule's geometry in Cartesian coordinates 
        function varargout = geom(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('geom', this.objectHandle, varargin{:});
        end
        
        % matpsi.set_geom( matrix(natom, 3) ); 
        % Set the molecule's geometry to the given Cartesian coordinates; 
        % if matpsi.fix_mol() has not been executed the geometry will only be "relative" 
        function varargout = set_geom(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('set_geom', this.objectHandle, varargin{:});
        end
        
        % vector(natom) = matpsi.Zlist(); 
        % Return a vector containing every atom's proton number 
        function varargout = Zlist(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('Zlist', this.objectHandle, varargin{:});
        end
        
        % double = matpsi.Enuc(); 
        % Return the molecule's nuclear Coulomb repulsion energy in Hartree 
        function varargout = Enuc(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('Enuc', this.objectHandle, varargin{:});
        end
        
        
        %% Basis set properties 
        % integer = matpsi.nbasis(); 
        % Number of basis functions 
        function varargout = nbasis(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('nbasis', this.objectHandle, varargin{:});
        end
        
        % vector(nbasis) = matpsi.func2center(); 
        % A mapping from every basis function to the number of atom it is centered on 
        function varargout = func2center(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('func2center', this.objectHandle, varargin{:});
        end
        
        % vector(nbasis) = matpsi.func2am; 
        % A mapping from every basis function to its angular momentum quantum number 
        function varargout = func2am(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('func2am', this.objectHandle, varargin{:});
        end
        
        
        %% One-electron integrals 
        % matrix(nbasis, nbasis) = matpsi.overlap(); 
        % Overlap matrix 
        function varargout = overlap(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('overlap', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.kinetic(); 
        % Kinetic energy matrix; one part of the core Hamiltonian H1 
        function varargout = kinetic(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('kinetic', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.potential(); 
        % Potential energy matrix; the other part of H1 
        function varargout = potential(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential', this.objectHandle, varargin{:});
        end
        
        % 3-D_array(nbasis, nbasis, natom) = matpsi.potential_sep(); 
        % Single-atom potential energy matrices; use mat(:, :, iatom) to extract one
        function varargout = potential_sep(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential_sep', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.potential_Zxyz(); 
        % Coulomb potential for a list of point charges 
        function varargout = potential_Zxyz(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential_Zxyz', this.objectHandle, varargin{:});
        end
        
        % dipole, 3 (nbasis, nbasis) matrix 
        function varargout = dipole(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('dipole', this.objectHandle, varargin{:});
        end
        
        %% Two-electron integrals 
        % double = matpsi.tei_ijkl(i ,j ,k ,l); 
        % Two-electron integral (ij|kl) 
        function varargout = tei_ijkl(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_ijkl', this.objectHandle, varargin{:});
        end
        
        % integer = matpsi.tei_uniqN(): 
        % The number of unique (without spatial symmetry) two-electron integrals 
        function varargout = tei_uniqN(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_uniqN', this.objectHandle, varargin{:});
        end
        
        % vector(nuniq) = matpsi.tei_alluniq(); 
        % A vector of all unique two-electron integrals 
        function varargout = tei_alluniq(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_alluniq', this.objectHandle, varargin{:});
        end
        
        % 4-D_array(nbasis, nbasis, nbasis, nbasis) = matpsi.tei_allfull(); 
        % All (with redundancy) two-electron integrals in a 4-D array 
        function varargout = tei_allfull(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_allfull', this.objectHandle, varargin{:});
        end
        
        % [vector(nuniq), vector(nuniq)] = matpsi.tei_alluniqJK(); 
        % Two vectors containing all unique two-electron integrals; 
        % 1st is the "regular" one for computing Coulomb matrix J 
        % 2nd is specifically arranged for computing exchange matrix K 
        function varargout = tei_alluniqJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_alluniqJK', this.objectHandle, varargin{:});
        end
        
        %% JK related 
        % matpsi.UseDirectJK(); 
        % Enable DirectJK object that computes J and K with direct algoritm 
        function varargout = UseDirectJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('UseDirectJK', this.objectHandle, varargin{:});
        end
        
        % matpsi.UsePKJK(); 
        % Enable PKJK object that computes J and K with PK algoritm 
        % Psi4's default JK method, uses disk space to store integrals 
        function varargout = UsePKJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('UsePKJK', this.objectHandle, varargin{:});
        end
        
        % matpsi.UseICJK(); 
        % Enable ICJK object that computes J and K in core 
        % In-Core JK 
        function varargout = UseICJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('UseICJK', this.objectHandle, varargin{:});
        end
        
        % matpsi.UseMatlabJK(); 
        % Enable MatlabJK object that computes J and K in core 
        % ### EXPERT ###  
        function varargout = UseMatlabJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('UseMatlabJK', this.objectHandle, varargin{:});
        end
        
        function varargout = SetMatlabJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('SetMatlabJK', this.objectHandle, varargin{:});
        end
        
        function varargout = DisableMatlabJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('DisableMatlabJK', this.objectHandle, varargin{:});
        end
        % ### EXPERT ###  
        
        % matpsi.JKtype(); 
        % See type of JK currently using 
        function varargout = JKtype(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('JKtype', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.Density2J( matrix(nbasis, nbasis) ); 
        % Compute J from density matrix 
        function varargout = Density2J(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('Density2J', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.Density2K( matrix(nbasis, nbasis) ); 
        % Compute K from density matrix 
        function varargout = Density2K(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('Density2K', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.Density2G( matrix(nbasis, nbasis) ); 
        % Compute G from density matrix 
        function varargout = Density2G(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('Density2G', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.OccMO2J( matrix(nbasis, nocc) ); 
        % Compute J from occupied molecular orbital matrix 
        function varargout = OccMO2J(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('OccMO2J', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.OccMO2K( matrix(nbasis, nocc) ); 
        % Compute K from occupied molecular orbital matrix 
        function varargout = OccMO2K(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('OccMO2K', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.OccMO2G( matrix(nbasis, nocc) ); 
        % Compute G from occupied molecular orbital matrix 
        function varargout = OccMO2G(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('OccMO2G', this.objectHandle, varargin{:});
        end
        
        %% SCF related 
        % double = matpsi.RHF_msqc( matrix(nbasis, nbasis), matrix(nbasis, nbasis), matrix(nbasis, nbasis) ); 
        % Hartree-Fock engine for MSQC 
        function varargout = RHF_msqc(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_msqc', this.objectHandle, varargin{:});
        end
        
        % double = matpsi.RHF_msqc_fromC( matrix(nbasis, nbasis), matrix(nbasis, nbasis), matrix(nbasis, nbasis), matrix(nbasis, nbasis) ); 
        % Hartree-Fock engine for MSQC 
        function varargout = RHF_msqc_fromC(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_msqc_fromC', this.objectHandle, varargin{:});
        end
        
        % double = matpsi.RHF(); 
        % Run a Hartree-Fock calculation and returns the final total energy 
        function varargout = RHF(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF', this.objectHandle, varargin{:});
        end
        
        % double = matpsi.RHFenv( matrix(nbasis, nbasis) ); 
        % Run a Hartree-Fock with a given external potential field 
        function varargout = RHFenv(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHFenv', this.objectHandle, varargin{:});
        end
        
        % double = matpsi.RHF_fromC( matrix(nbasis, nbasis) ); 
        % Run a Hartree-Fock starting from a given molecular orbital matrix 
        function varargout = RHF_fromC(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_fromC', this.objectHandle, varargin{:});
        end
        
        % double = matpsi.RHFenv_fromC( matrix(nbasis, nbasis) ); 
        % Run a Hartree-Fock with a given external potential field starting from a given molecular orbital matrix 
        function varargout = RHFenv_fromC(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHFenv_fromC', this.objectHandle, varargin{:});
        end
        
        % matpsi.RHF_reset(); 
        % Reset the RHF object 
        function varargout = RHF_reset(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_reset', this.objectHandle, varargin{:});
        end
        
        % matpsi.RHF_EnableMOM(integer); 
        % Start MOM from the given number of iteration, default 20 
        % may help converge 
        function varargout = RHF_EnableMOM(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_EnableMOM', this.objectHandle, varargin{:});
        end
        
        % matpsi.RHF_DisableMOM(); 
        % Disable MOM 
        function varargout = RHF_DisableMOM(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_DisableMOM', this.objectHandle, varargin{:});
        end
        
        % matpsi.RHF_EnableDamping(double); 
        % Damp each Hartree-Fock iteration with the given rate, default 20.0 
        % may help converge; minimum 0.0, maximum 100.0 
        function varargout = RHF_EnableDamping(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_EnableDamping', this.objectHandle, varargin{:});
        end
        
        % matpsi.RHF_DisableDamping(); 
        % Disable damping 
        function varargout = RHF_DisableDamping(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_DisableDamping', this.objectHandle, varargin{:});
        end
        
        % matpsi.RHF_EnableDIIS(); 
        % Enable DIIS 
        % may help converge 
        function varargout = RHF_EnableDIIS(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_EnableDIIS', this.objectHandle, varargin{:});
        end
        
        % matpsi.RHF_DisableDIIS(); 
        % Enable DIIS 
        % may help converge 
        function varargout = RHF_DisableDIIS(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_DisableDIIS', this.objectHandle, varargin{:});
        end
        
        % matpsi.RHF_GuessSAD(); 
        % Start from SAD guess density 
        % may help converge 
        function varargout = RHF_GuessSAD(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_GuessSAD', this.objectHandle, varargin{:});
        end
        
        % matpsi.RHF_GuessCore(); 
        % Start from core guess density 
        % default 
        function varargout = RHF_GuessCore(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_GuessCore', this.objectHandle, varargin{:});
        end
        
        %~ The below RHF_X functions must be executed after matpsi.RHF(); 
        % double = matpsi.RHF_EHF(); 
        % Hartree-Fock energy 
        function varargout = RHF_EHF(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_EHF', this.objectHandle, varargin{:});
        end
        
        % matrix(nbasis, nbasis) = matpsi.RHF_C(); 
        % Molecular orbital coefficients, ascending sorted by their eigenvalues 
        function varargout = RHF_C(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_C', this.objectHandle, varargin{:});
        end
        
        % vector(nbasis) = matpsi.RHF_EMO(); 
        % Molecular orbital energies, ascending sorted by their eigenvalues 
        function varargout = RHF_EMO(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_EMO', this.objectHandle, varargin{:});
        end
        
        % vector(nbasis) = matpsi.RHF_D(); 
        % Final density matrix 
        function varargout = RHF_D(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_D', this.objectHandle, varargin{:});
        end
        
        % vector(nbasis) = matpsi.RHF_H(); 
        % Core 1-electron Hamiltonian H1 
        function varargout = RHF_H(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_H', this.objectHandle, varargin{:});
        end
        
        % vector(nbasis) = matpsi.RHF_J();
        % Coulomb matrix J 
        function varargout = RHF_J(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_J', this.objectHandle, varargin{:});
        end
        
        % vector(nbasis) = matpsi.RHF_K(); 
        % Exchange matrix K 
        function varargout = RHF_K(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_K', this.objectHandle, varargin{:});
        end
        
        % vector(nbasis) = matpsi.RHF_F(); 
        % Total Fock matrix F 
        function varargout = RHF_F(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_F', this.objectHandle, varargin{:});
        end

    end
end
