/******************************************
 *           ENERGY FACTORY               *
 * Interface for creating an object       *
 * that represents an energy function.    *
 ******************************************/

#ifndef COCOS_ENERGY_FACTORY__
#define COCOS_ENERGY_FACTORY__

#include "globals.h"
#include "energy.h"
#include "potential_energy.h"
#include "rmsd_energy.h"
#include "contact_energy.h"

class EnergyFactory {
public:
  enum EnergyType {
    PotentialEnergy,
    RmsdEnergy,
    ContactEnergy,
    ContactDecayEnergy
  };
  
  static Energy*
  getEnergyFunction( EnergyType energyType, agent_type agentType=structure ) {
    switch ( energyType ) {
      case PotentialEnergy:
      {
        real w_a = gh_params.str_weights[ 0 ];
        real w_b = gh_params.str_weights[ 1 ];
        real w_c = gh_params.str_weights[ 2 ];
        if ( agentType == coordinator ) {
          w_a = gh_params.crd_weights[ 0 ];
          w_b = gh_params.crd_weights[ 1 ];
          w_c = gh_params.crd_weights[ 2 ];
        }
        return ((new PotentialEnergy::PotentialEnergy ())->set_parameters ( gd_params.secondary_s_info,
                                                                            gd_params.h_distances,
                                                                            gd_params.h_angles,
                                                                            gd_params.contact_params,
                                                                            gd_params.aa_seq,
                                                                            gd_params.tors, gd_params.tors_corr,
                                                                            w_a, w_b, w_c ));
      }
      case RmsdEnergy:
        return ((new RmsdEnergy::RmsdEnergy())->set_parameters ( gd_params.known_prot ));
      case ContactEnergy:
        return ((new ContactEnergy::ContactEnergy())->set_parameters ( gd_params.aa_seq ));
      case ContactDecayEnergy:
        return ((new ContactEnergy::ContactEnergy())->set_parameters ( gd_params.aa_seq ));
    }
    throw "Invalid Energy function.";
  }//getEnergyFunction
};


#endif
