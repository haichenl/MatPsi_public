/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*
 *  rhf.h
 *  matrix
 *
 *  Created by Justin Turney on 4/10/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef RHF_H
#define RHF_H

#include <libpsio/psio.hpp>
#include "hf.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class TwoBodySOInt;
class PSIO;
class Chkpt;
class Matrix;
class Vector;

namespace scf {

class RHF : public HF {
protected:
    SharedMatrix D_;
    SharedMatrix Dold_;
    SharedMatrix G_;
    SharedMatrix J_;
    SharedMatrix K_;
    
    // control whether we manually specify a starting molecular orbital matrix 
    bool StartingC_enabled_;
    SharedMatrix StartingC_;
    
    // control whether we modify J_ and K_ 
    bool JKmodifiers_enabled_;
    SharedMatrix Jmodifier_;
    SharedMatrix Kmodifier_;

    void form_C();
    void form_D();
    virtual void damp_update();
    double compute_initial_E();
    virtual double compute_E();
    virtual void stability_analysis();

    //Some stuff for Ed Hohenstein's SAPT code
    // TODO: This must be removed for a conforming SCF module
    // The SAPT driver should save the three references and extract info from
    // That point
    void save_sapt_info();

    virtual void form_F();
    virtual void form_G();
    virtual void compute_orbital_gradient(bool save_fock);

    bool diis();

    bool test_convergency();
    void save_information();

    void common_init();

    // Finalize memory/files
    virtual void finalize();

    void save_density_and_energy();
    
    // Automatically judge whether we start from a manually specified starting molecular orbital 
    void whether_to_use_StartingC();

public:
    RHF(Process::Environment& process_environment_in, Options& options, boost::shared_ptr<JK> jk_in, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    RHF(Process::Environment& process_environment_in, Options& options, boost::shared_ptr<JK> jk_in, boost::shared_ptr<PSIO> psio);
    virtual ~RHF();
    
    // Finalize memory/files
    void extern_finalize();

    virtual SharedMatrix Da() const;
    virtual SharedMatrix J() const { return J_; }
    virtual SharedMatrix K() const { return K_; }

    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
    
    // Control whether we manually specify a starting molecular orbital matrix 
    void set_StartingC(SharedMatrix StartingC_in);
    void disable_StartingC();
    
    // Control whether we modify J_ and K_ 
    void set_JKmodifiers(SharedMatrix Jmodifier_in, SharedMatrix Kmodifier_in);
    void disable_JKmodifiers();
};

}}

#endif
