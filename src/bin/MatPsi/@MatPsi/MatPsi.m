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
                toppath = regexp(path(), '[^:]*', 'match', 'once');
                if(exist([toppath, '/share/basis'], 'file'))
                    this.path = toppath;
                elseif(exist([toppath, '/@MatPsi/share/basis'], 'file'))
                    this.path = [toppath, '/@MatPsi'];
                else
                    throw(MException('MatPsi:MatPsi','MatPsi cannot find basis set files.'));
                end
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
        % molecule_string, 1 string 
        function varargout = molecule_string(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('molecule_string', this.objectHandle, varargin{:});
        end
        
        % basis_name, 1 string 
        function varargout = basis_name(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('basis_name', this.objectHandle, varargin{:});
        end

        % ignore this 
        function varargout = testmol(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('testmol', this.objectHandle, varargin{:});
        end
        
        %% Molecule operations 
        % fix_mol 
        function varargout = fix_mol(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('fix_mol', this.objectHandle, varargin{:});
        end
        
        % free_mol 
        function varargout = free_mol(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('free_mol', this.objectHandle, varargin{:});
        end
        
        %% Molecule properties 
        % natom, 1 double 
        function varargout = natom(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('natom', this.objectHandle, varargin{:});
        end
        
        % nelec, 1 double 
        function varargout = nelec(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('nelec', this.objectHandle, varargin{:});
        end
        
        % geom, (natom, 3) matrix 
        function varargout = geom(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('geom', this.objectHandle, varargin{:});
        end
        
        % set_geom, (natom, 3) matrix 
        function varargout = set_geom(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('set_geom', this.objectHandle, varargin{:});
        end
        
        % Zlist, (natom, 1) matrix 
        function varargout = Zlist(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('Zlist', this.objectHandle, varargin{:});
        end
        
        % Enuc, 1 double 
        function varargout = Enuc(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('Enuc', this.objectHandle, varargin{:});
        end
        
        %% Basis set properties 
        % nbasis, 1 double 
        function varargout = nbasis(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('nbasis', this.objectHandle, varargin{:});
        end
        
        % func2center, (nbasis, 1) matrix 
        function varargout = func2center(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('func2center', this.objectHandle, varargin{:});
        end
        
        % func2am, (nbasis, 1) matrix 
        function varargout = func2am(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('func2am', this.objectHandle, varargin{:});
        end
        
        %% One-electron integrals 
        % overlap, (nbasis, nbasis) matrix 
        function varargout = overlap(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('overlap', this.objectHandle, varargin{:});
        end
        
        % kinetic, (nbasis, nbasis) matrix 
        function varargout = kinetic(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('kinetic', this.objectHandle, varargin{:});
        end
        
        % potential, (nbasis, nbasis) matrix 
        function varargout = potential(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential', this.objectHandle, varargin{:});
        end
        
        % atom-separated potential, (nbasis, nbasis, natom) 3-d matrix 
        function varargout = potential_sep(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential_sep', this.objectHandle, varargin{:});
        end
        
        % environment potential for a list containing a lot of point charges, (nbasis, nbasis) matrix 
        function varargout = potential_Zxyz(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('potential_Zxyz', this.objectHandle, varargin{:});
        end
        
        % dipole, 3 (nbasis, nbasis) matrix 
        function varargout = dipole(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('dipole', this.objectHandle, varargin{:});
        end
        
        %% Two-electron integrals 
        % tei_ijkl, 1 double 
        function varargout = tei_ijkl(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_ijkl', this.objectHandle, varargin{:});
        end
        
        % tei_uniqN, 1 double  
        function varargout = tei_uniqN(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_uniqN', this.objectHandle, varargin{:});
        end
        
        % tei_alluniq, (nuniq, 1) vector 
        function varargout = tei_alluniq(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_alluniq', this.objectHandle, varargin{:});
        end
        
        % tei_allfull, (nbasis, nbasis, nbasis, nbasis) array 
        function varargout = tei_allfull(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_allfull', this.objectHandle, varargin{:});
        end
        
        % tei_alluniqJK, (nuniq, 1) vector 
        function varargout = tei_alluniqJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('tei_alluniqJK', this.objectHandle, varargin{:});
        end
        
        %% SCF related 
        % UseDirectJK 
        function varargout = UseDirectJK(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('UseDirectJK', this.objectHandle, varargin{:});
        end
        
        % OccMO2J, (nbasis, nbasis) matrix 
        function varargout = OccMO2J(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('OccMO2J', this.objectHandle, varargin{:});
        end
        
        % OccMO2K, (nbasis, nbasis) matrix 
        function varargout = OccMO2K(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('OccMO2K', this.objectHandle, varargin{:});
        end
        
        % OccMO2G, (nbasis, nbasis) matrix 
        function varargout = OccMO2G(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('OccMO2G', this.objectHandle, varargin{:});
        end
        
        % RHF, 1 double 
        function varargout = RHF(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF', this.objectHandle, varargin{:});
        end
        
        % RHFenv, 1 double 
        function varargout = RHFenv(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHFenv', this.objectHandle, varargin{:});
        end
        
        % RHF_finalize, nothing  
        function varargout = RHF_finalize(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_finalize', this.objectHandle, varargin{:});
        end 
        
        % RHF_EHF, 1 double 
        function varargout = RHF_EHF(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_EHF', this.objectHandle, varargin{:});
        end
        
        % RHF_C, (nbasis, nbasis) matrix  
        function varargout = RHF_C(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_C', this.objectHandle, varargin{:});
        end
        
        % RHF_EMO, (nbasis, 1) vector  
        function varargout = RHF_EMO(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_EMO', this.objectHandle, varargin{:});
        end
        
        % RHF_D, (nbasis, nbasis) matrix  
        function varargout = RHF_D(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_D', this.objectHandle, varargin{:});
        end
        
        % RHF_H, (nbasis, nbasis) matrix  
        function varargout = RHF_H(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_H', this.objectHandle, varargin{:});
        end
        
        % RHF_J, (nbasis, nbasis) matrix  
        function varargout = RHF_J(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_J', this.objectHandle, varargin{:});
        end
        
        % RHF_K, (nbasis, nbasis) matrix  
        function varargout = RHF_K(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_K', this.objectHandle, varargin{:});
        end
        
        % RHF_F, (nbasis, nbasis) matrix  
        function varargout = RHF_F(this, varargin)
            [varargout{1:nargout}] = MatPsi_mex('RHF_F', this.objectHandle, varargin{:});
        end

    end
end
