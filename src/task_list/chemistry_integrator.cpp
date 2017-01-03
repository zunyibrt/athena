//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file chemistry_integrator.cpp
//  \brief derived class for chemistry integrator task list.
//======================================================================================

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#ifdef INCLUDE_CHEMISTRY
#include "../chemistry/species.hpp" 
#include "../chemistry/network/network.hpp" 
#endif

// this class header
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
//  ChemistryIntegratorTaskList constructor
ChemistryIntegratorTaskList::ChemistryIntegratorTaskList(ParameterInput *pin, Mesh *pm)
  : TaskList(pm)
{
  nsub_steps = 1;
  // Now assemble list of tasks for each step of chemistry integrator
  {using namespace ChemistryIntegratorTaskNames;
    AddChemistryIntegratorTask(START_SPEC_RECV, NONE);
    AddChemistryIntegratorTask(INT_CHEM_SRC,NONE);
    //MPI boundary
    AddChemistryIntegratorTask(SEND_SPEC, INT_CHEM_SRC);
    AddChemistryIntegratorTask(RECV_SPEC, START_SPEC_RECV);
    AddChemistryIntegratorTask(CLEAR_SPEC_RECV, RECV_SPEC);

    //add advection term here
  } // end of using namespace block
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  
void ChemistryIntegratorTaskList::AddChemistryIntegratorTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;
  using namespace ChemistryIntegratorTaskNames;
  switch((id)) {
    case (INT_CHEM_SRC):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemistryIntegratorTaskList::IntegrateSourceTerm);
      break;
    case (START_SPEC_RECV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemistryIntegratorTaskList::StartSpeciesReceive);
      break;
    case (CLEAR_SPEC_RECV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemistryIntegratorTaskList::ClearSpeciesReceive);
      break;
    case (SEND_SPEC):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemistryIntegratorTaskList::SpeciesSend);
      break;
    case (RECV_SPEC):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemistryIntegratorTaskList::SpeciesReceive);
      break;
    //add advection here
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in Add Chemistry Task" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

//--------------------------------------------------------------------------------------
// Functions to integrate chemistry
enum TaskStatus ChemistryIntegratorTaskList::IntegrateSourceTerm(MeshBlock *pmb,
                                                                 int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->pspec->podew->Integrate();
#endif
  return TASK_SUCCESS;
}

//MPI boundary
enum TaskStatus ChemistryIntegratorTaskList::StartSpeciesReceive(MeshBlock *pmb,
                                                                 int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->pbval->StartReceivingSpecies();
#endif
  return TASK_SUCCESS;
}

enum TaskStatus ChemistryIntegratorTaskList::ClearSpeciesReceive(MeshBlock *pmb,
                                                                 int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->pbval->ClearBoundarySpecies();
#endif
  return TASK_SUCCESS;
}

enum TaskStatus ChemistryIntegratorTaskList::SpeciesSend(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->pbval->SendSpeciesBoundaryBuffers(pmb->pspec->s);
#endif
  return TASK_SUCCESS;
}

enum TaskStatus ChemistryIntegratorTaskList::SpeciesReceive(MeshBlock *pmb, int step)
{
  bool ret=true;
#ifdef INCLUDE_CHEMISTRY
  ret = pmb->pbval->ReceiveSpeciesBoundaryBuffers(pmb->pspec->s);
#endif
  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

//add advection here
