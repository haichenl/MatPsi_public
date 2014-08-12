
/* 2D Column-Major-Order array lookup pattern. */
#define CMO2(i, j, I) ((j) * (I) + (i))

namespace psi {
/**
 * Class MatlabJK
 *
 * Receives Matlab cell array pointers then compute J and K 
 */
class MatlabJK : public JK {
protected:
    boost::shared_array<double*> Jcell_ptr_;
    boost::shared_array<double*> Kcell_ptr_;
    
    virtual bool C1() const { return false; }
    virtual void preiterations();
    virtual void compute_JK();
    virtual void postiterations();
    
public:
    MatlabJK(Process::Environment& process_environment_in, boost::shared_ptr<BasisSet> primary);
    virtual ~MatlabJK();
    
    virtual void print_header() const;
    
    void set_JKcell_ptrs(boost::shared_array<double*> Jcell_ptr_in, boost::shared_array<double*> Kcell_ptr_in);
    
    void disable();

};

MatlabJK::MatlabJK(Process::Environment& process_environment_in, boost::shared_ptr<BasisSet> primary) : JK(process_environment_in, primary) {
    JKtype_ = "MatlabJK";
}

MatlabJK::~MatlabJK() {
}

void MatlabJK::postiterations() {
}

void MatlabJK::print_header() const {
}

void MatlabJK::preiterations() {
}

void MatlabJK::set_JKcell_ptrs(boost::shared_array<double*> Jcell_ptr_in, boost::shared_array<double*> Kcell_ptr_in) {
    Jcell_ptr_ = Jcell_ptr_in;
    Kcell_ptr_ = Kcell_ptr_in;
}

inline void loop_compute(double* const target_ptr, double* const density_ptr, boost::shared_array<double*> cell_ptr, int& nbf) {
    for (int j = 0; j < nbf; j++) {
        for (int i = j; i < nbf; i++) {
            double element = 0.0;
            
            double* c_mat_ptr = cell_ptr[CMO2(i, j, nbf)];
            
            /* Build the elements of J or K. */
            for (int k = 0; k < nbf * nbf; k++) {
                element += c_mat_ptr[k] * density_ptr[k];
            }
            
            target_ptr[CMO2(i, j, nbf)] = element;
            
            /* Taking advantage of symmetry. */
            if (i != j) {
                target_ptr[CMO2(j, i, nbf)] = element;
            }
        }
    }
}

void MatlabJK::compute_JK() {
    if(Jcell_ptr_ == NULL || Kcell_ptr_ == NULL)
        throw PSIEXCEPTION("MatlabJK::compute_JK(): J and K cell arrays not initialized yet.");
    int nbf = primary_->nbf();
    if(do_J_ && J_.size()) {
        for (int N = 0; N < D_.size(); ++N) { 
            loop_compute(J_[N]->get_pointer(), D_[N]->get_pointer(), Jcell_ptr_, nbf);
        }
    }
    
    if(do_K_ && K_.size()) {
        for (int N = 0; N < D_.size(); ++N) { 
            loop_compute(K_[N]->get_pointer(), D_[N]->get_pointer(), Kcell_ptr_, nbf);
        }
    }
    
    if(do_wK_ || wK_.size()) {
        throw PSIEXCEPTION("MatlabJK::compute_JK(): wK not supported for MatlabJK now.");
    }
}

void MatlabJK::disable() {
    Jcell_ptr_.reset();
    Kcell_ptr_.reset();
}

} // end of namespace psi 

