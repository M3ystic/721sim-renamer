#include "renamer.h"
//#define DEBUG

#ifdef DEBUG
#define print(x) std::cout << x
#else
#define print(x)
#endif



	/////////////////////////////////////////////////////////////////////
	// This is the constructor function.
	// When a renamer object is instantiated, the caller indicates:
	// 1. The number of logical registers (e.g., 32).
	// 2. The number of physical registers (e.g., 128).
	// 3. The maximum number of unresolved branches.
	//    Requirement: 1 <= n_branches <= 64.
	// 4. The maximum number of active instructions (Active List size).
	//
	// Tips:
	//
	// Assert the number of physical registers > number logical registers.
	// Assert 1 <= n_branches <= 64.
	// Assert n_active > 0.
	// Then, allocate space for the primary data structures.
	// Then, initialize the data structures based on the knowledge
	// that the pipeline is intially empty (no in-flight instructions yet).
	/////////////////////////////////////////////////////////////////////
	renamer::renamer(uint64_t n_log_regs,
		uint64_t n_phys_regs,
		uint64_t n_branches,
		uint64_t n_active){
		
		this->n_phys_regs = n_phys_regs;
		this->n_log_regs = n_log_regs;
		this->n_branches = n_branches;
		this->n_active = n_active;

		 free_list_size = this->n_phys_regs - this->n_log_regs;

		// assert the input parameters are valid
		assert(n_phys_regs > n_log_regs);
		assert(n_active > 0);
		assert (n_branches >=1 && n_branches <= 64);
		// reserve space for RMT
		rename_map_table.resize(n_log_regs);
		// initialize RMT
		for (uint64_t i = 0; i < n_log_regs; i++){
			rename_map_table[i] = i;
		}
		// reserve space for AMT
		arch_map_table.resize(n_log_regs);
		// initialize AMT
		arch_map_table = rename_map_table; // initially, RMT is same as AMT
		// reserve space for free list
		free_list.phys_reg.reserve(n_phys_regs - n_log_regs);
		// initialiaze free list
		for (uint64_t i = n_log_regs; i < n_phys_regs; i++){
			free_list.phys_reg.emplace_back(i);
		}

		PRF.resize(n_phys_regs, 0);
		//make sure ready bits of commited state (those originally set in amt/rmt above are 1)
		PRF_ready.resize(n_phys_regs, 0);
		 for (uint64_t i = 0; i < n_log_regs; i++){
		 	PRF_ready.at(i) = 1;
		 }

		active_list.resize(n_active);

		gbm_mask.resize(n_branches);

			for (uint64_t i = 0; i < n_branches; i++){
				gbm_mask.at(i).shadow_table.resize(n_log_regs);
				
			}
		}

	/////////////////////////////////////////////////////////////////////
	// This is the destructor, used to clean up memory space and
	// other things when simulation is done.
	// I typically don't use a destructor; you have the option to keep
	// this function empty.
	/////////////////////////////////////////////////////////////////////
	renamer::~renamer(){
	}

	
	void renamer::write_gbm(uint64_t branch_ID) { GBM = GBM | (1ULL << branch_ID); }
	void renamer::clear_gbm(uint64_t branch_ID) { GBM = GBM & ~(1ULL << branch_ID);}

	uint64_t renamer::count_gbm_zeroes(){
		print("enter count gbm zeroes function \n");
		uint64_t zereo_count = 0;
		for (uint64_t i = 0; i < this->n_branches; i++){
			if ((GBM & (1ULL << i)) == 0){
				zereo_count++;
			}
		}
		print("count gbm zeroes function :: free checkpoints = " << zereo_count << "\n");
		print("exit count gbm zeroes function \n");
		return zereo_count;
	}

	uint64_t renamer::get_branch_id(){
		print("enter get branch id function \n");
		for (uint64_t i = 0; i < this->n_branches; i++){
			if (!(GBM & ( 1ULL << i))){
				print(" get branch id function :: branch id = " << i << "\n");
				print("exit get branch id function \n");
				return i;
			}
		}
		print("exit get branch id function :: no free branch id found\n");
		assert(false); // should not happen
	}
	


	//////////////////////////////////////////
	// Functions related to Rename Stage.   //
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// The Rename Stage must stall if there aren't enough free physical
	// registers available for renaming all logical destination registers
	// in the current rename bundle.
	//
	// Inputs:
	// 1. bundle_dst: number of logical destination registers in
	//    current rename bundle
	//
	// Return value:
	// Return "true" (stall) if there aren't enough free physical
	// registers to allocate to all of the logical destination registers
	// in the current rename bundle.
	/////////////////////////////////////////////////////////////////////
	bool renamer::stall_reg(uint64_t bundle_dst){
		print("enter stall reg function \n");
		bool space_left;
		uint64_t regs;
			if (free_list.head_phase == free_list.tail_phase){
				regs = (free_list.tail - free_list.head + free_list_size) % free_list_size;
			}
			else{
				// this means the free list is full (no registers are in used)
				regs = free_list_size - ((free_list.head - free_list.tail + free_list_size) % free_list_size);
			}
		
			space_left = regs >= bundle_dst;

		print("stall reg function :: free list head pointer = " << free_list.head << " free list tail pointer = " << free_list.tail << "\n");
		print("stall reg function :: free phy regs = " << regs << "\n");
		if (!space_left){
			print("**stall reg function :: not enough registers so stall**\n");
			print("exit stall reg function \n");
				return true;
		}
		print("exit stall reg function \n");
		return false;
	}	

	/////////////////////////////////////////////////////////////////////
	// The Rename Stage must stall if there aren't enough free
	// checkpoints for all branches in the current rename bundle.
	//
	// Inputs:
	// 1. bundle_branch: number of branches in current rename bundle
	//
	// Return value:
	// Return "true" (stall) if there aren't enough free checkpoints
	// for all branches in the current rename bundle.
	/////////////////////////////////////////////////////////////////////
	bool renamer::stall_branch(uint64_t bundle_branch){
		print("enter stall branch function \n");
		bool space_left = count_gbm_zeroes() >= bundle_branch;
		print("stall branch function :: free checkpoints = " << count_gbm_zeroes() << "\n");
		if (!space_left){
			print("**stall branch function :: not enough branch checkpoints so stall**\n");
			print("exit stall branch function \n");
			return true;
		}
		print("exit stall branch function \n");
		return false;
	}

	/////////////////////////////////////////////////////////////////////
	// This function is used to get the branch mask for an instruction.
	/////////////////////////////////////////////////////////////////////
	uint64_t renamer::get_branch_mask(){
		print("enter get branch mask function \n");
		uint64_t branch_mask = GBM;
		print("get branch mask function :: branch mask = " << branch_mask << "\n");
		print("exit get branch mask function \n");
		return branch_mask;
		
	}

	/////////////////////////////////////////////////////////////////////
	// This function is used to rename a single source register.
	//
	// Inputs:
	// 1. log_reg: the logical register to rename
	//
	// Return value: physical register name
	/////////////////////////////////////////////////////////////////////
	uint64_t renamer::rename_rsrc(uint64_t log_reg){
		print("enter rename src function \n");
		print("rename src function :: log reg : " << log_reg << " phys reg : " << rename_map_table.at(log_reg) << "\n");
		print("exit rename src function \n");
		return rename_map_table.at(log_reg);
		
	}

	/////////////////////////////////////////////////////////////////////
	// This function is used to rename a single destination register.
	//
	// Inputs:
	// 1. log_reg: the logical register to rename
	//
	// Return value: physical register name
	/////////////////////////////////////////////////////////////////////
	uint64_t renamer::rename_rdst(uint64_t log_reg){
		/// dont know if need to check if phase bits are the same and if head pointer has wrapped around to tail pointer
		print("enter rename dst function \n");
		if ((free_list.head || free_list.tail) > free_list.phys_reg.size()){
			print("assertion failed in rename dst function :: free list head pointer = " << free_list.head << " free list tail pointer = " << free_list.tail << "\n");
		}

		assert(free_list.head < free_list.phys_reg.size());
		assert(free_list.tail < free_list.phys_reg.size());

		assert(free_list.head != free_list.tail || free_list.head_phase != free_list.tail_phase); // assert that there is a free physical register to allocate
		rename_map_table.at(log_reg) = free_list.phys_reg.at(free_list.head);
		print("renamed log reg:" << log_reg << " to phys reg: " <<  free_list.phys_reg.at(free_list.head) << "\n");
		// update the head pointer and phase bits if necessary
		free_list.head = (free_list.head + 1) % (this->n_phys_regs - n_log_regs);
		print("rename dst function :: free list head pointer after update = " << free_list.head << "\n");
		// check if head pointer wrapped around the free list and update phase bits if necessary
		if (free_list.head == 0){
			print("rename dst function :: free list head pointer wrapped around, so update phase bit\n");
			free_list.head_phase = !free_list.head_phase;
		}
		print("exit rename dst function \n");
		return rename_map_table.at(log_reg);
	}

	/////////////////////////////////////////////////////////////////////
	// This function creates a new branch checkpoint.
	//
	// Inputs: none.
	//
	// Output:
	// 1. The function returns the branch's ID. When the branch resolves,
	//    its ID is passed back to the renamer via "resolve()" below.
	//
	// Tips:
	//
	// Allocating resources for the branch (a GBM bit and a checkpoint):
	// * Find a free bit -- i.e., a '0' bit -- in the GBM. Assert that
	//   a free bit exists: it is the user's responsibility to avoid
	//   a structural hazard by calling stall_branch() in advance.
	// * Set the bit to '1' since it is now in use by the new branch.
	// * The position of this bit in the GBM is the branch's ID.
	// * Use the branch checkpoint that corresponds to this bit.
	// 
	// The branch checkpoint should contain the following:
	// 1. Shadow Map Table (checkpointed Rename Map Table)
	// 2. checkpointed Free List head pointer and its phase bit
	// 3. checkpointed GBM
	/////////////////////////////////////////////////////////////////////
	uint64_t renamer::checkpoint(){
		  print("check_point function :: checkpointing a new branch\n");
		// uint64_t branch_id = get_branch_id();
		// print("check_point function :: branch id = " << branch_id << "\n");
		// gbm_mask.at(branch_id).shadow_table = rename_map_table;
		// gbm_mask.at(branch_id).fl_head = free_list.head;
		// gbm_mask.at(branch_id).fl_head_phase = free_list.head_phase;
		// gbm_mask.at(branch_id).gbm_snapshot = GBM;
		// write_gbm(branch_id);
		print("exit check_point function \n");
		return 0;
		// return branch_id;
	}

	//////////////////////////////////////////
	// Functions related to Dispatch Stage. //
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// The Dispatch Stage must stall if there are not enough free
	// entries in the Active List for all instructions in the current
	// dispatch bundle.
	//
	// Inputs:
	// 1. bundle_inst: number of instructions in current dispatch bundle
	//
	// Return value:
	// Return "true" (stall) if the Active List does not have enough
	// space for all instructions in the dispatch bundle.
	/////////////////////////////////////////////////////////////////////
	bool renamer::stall_dispatch(uint64_t bundle_inst){
		print("enter stall dispatch function \n");
		bool space_left;
		if (al_head == al_tail){
			if (al_head_phase == al_tail_phase){
				//active list is empty
				space_left = true;
			}
			else{
				// this means the active list is full
				space_left = false;
			}	
		} else 
		{
			space_left = (this->n_active - ((al_tail - al_head + this->n_active) % this->n_active) >= bundle_inst);
		}

		print("stall dispatch function :: active list head pointer = " << al_head << " active list tail pointer = " << al_tail << "\n");
		print("stall dispatch function :: free entries in active list = " << ((this->n_active - ((al_tail - al_head + this->n_active) % this->n_active))) << "\n");
		if(!space_left){
			print("**stall dispatch function :: not enough space in active list so stall**\n");
			print("exit stall dispatch function \n");
			return true;
		}
		print("exit stall dispatch function \n");
		return false;
	}

	/////////////////////////////////////////////////////////////////////
	// This function dispatches a single instruction into the Active
	// List.
	//
	// Inputs:
	// 1. dest_valid: If 'true', the instr. has a destination register,
	//    otherwise it does not. If it does not, then the log_reg and
	//    phys_reg inputs should be ignored.
	// 2. log_reg: Logical register number of the instruction's
	//    destination.
	// 3. phys_reg: Physical register number of the instruction's
	//    destination.
	// 4. load: If 'true', the instr. is a load, otherwise it isn't.
	// 5. store: If 'true', the instr. is a store, otherwise it isn't.
	// 6. branch: If 'true', the instr. is a branch, otherwise it isn't.
	// 7. amo: If 'true', this is an atomic memory operation.
	// 8. csr: If 'true', this is a system instruction.
	// 9. PC: Program counter of the instruction.
	//
	// Return value:
	// Return the instruction's index in the Active List.
	//
	// Tips:
	//
	// Before dispatching the instruction into the Active List, assert
	// that the Active List isn't full: it is the user's responsibility
	// to avoid a structural hazard by calling stall_dispatch()
	// in advance.
	/////////////////////////////////////////////////////////////////////
	uint64_t renamer::dispatch_inst(bool dest_valid,
	                       uint64_t log_reg,
	                       uint64_t phys_reg,
	                       bool load,
	                       bool store,
	                       bool branch,
	                       bool amo,
	                       bool csr,
	                       uint64_t PC){
	
	print("enter dispatch inst function\n");
	assert(al_head != al_tail || al_head_phase == al_tail_phase); // assert active list isn't full

	active_list.at(al_tail).has_dest = dest_valid;
	active_list.at(al_tail).log_reg = log_reg;
	active_list.at(al_tail).phys_reg = phys_reg;
	active_list.at(al_tail).is_load = load;
	active_list.at(al_tail).is_store = store;
	active_list.at(al_tail).is_branch = branch;
	active_list.at(al_tail).is_amo = amo;
	active_list.at(al_tail).is_csr = csr;
	active_list.at(al_tail).PC = PC;

	active_list.at(al_tail).completed = false;
	active_list.at(al_tail).threw_exception = false;
	active_list.at(al_tail).load_violation = false;
	active_list.at(al_tail).b_mispred = false;
	active_list.at(al_tail).v_mispred = false;

	
	uint64_t current_tail = al_tail;
	al_tail = (al_tail + 1) % this->n_active;
	print("dispatch inst function :: active list tail pointer after update = " << al_tail << "\n");
	if (al_tail == 0){
		al_tail_phase = !al_tail_phase;
		print("dispatch inst function :: active list tail pointer wrapped around, so update phase bit\n");
	}
	print("exit dispatch inst function\n");
	return current_tail;
	}

	/////////////////////////////////////////////////////////////////////
	// Test the ready bit of the indicated physical register.
	// Returns 'true' if ready.
	/////////////////////////////////////////////////////////////////////
	bool renamer::is_ready(uint64_t phys_reg){
		print("enter is_ready function \n");
		print("PRF index: " << phys_reg << " ready test: " << PRF_ready.at(phys_reg) << "\n");
		bool result = PRF_ready.at(phys_reg);
		print("exit is_ready function \n");
		return result;
	}

	/////////////////////////////////////////////////////////////////////
	// Clear the ready bit of the indicated physical register.
	/////////////////////////////////////////////////////////////////////
	void renamer::clear_ready(uint64_t phys_reg){
		print("enter clear_ready function \n");
		print("PRF index: " << phys_reg << " cleared ready bit\n");
		PRF_ready.at(phys_reg) = 0;
		print("exit clear_ready function \n");
	}


	//////////////////////////////////////////
	// Functions related to the Reg. Read   //
	// and Execute Stages.                  //
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// Return the contents (value) of the indicated physical register.
	/////////////////////////////////////////////////////////////////////
	uint64_t renamer::read(uint64_t phys_reg){
		print("enter read function \n");
		print("PRF index: " << phys_reg << " read with value: " << PRF.at(phys_reg) << "\n");
		uint64_t result = PRF.at(phys_reg);
		print("exit read function \n");
		return result;
	}

	/////////////////////////////////////////////////////////////////////
	// Set the ready bit of the indicated physical register.
	/////////////////////////////////////////////////////////////////////
	void renamer::set_ready(uint64_t phys_reg){
		print("enter set_ready function \n");
		print("PRF index: " << phys_reg << " set ready bit\n");
	    PRF_ready.at(phys_reg) = true;
		print("exit set_ready function \n");
	}


	//////////////////////////////////////////
	// Functions related to Writeback Stage.//
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// Write a value into the indicated physical register.
	/////////////////////////////////////////////////////////////////////
	void renamer::write(uint64_t phys_reg, uint64_t value){
		print("enter write function \n");
		print("PRF index: " << phys_reg << " write with value: " << value << "\n");
		PRF.at(phys_reg) = value;
		print("exit write function \n");
	}

	/////////////////////////////////////////////////////////////////////
	// Set the completed bit of the indicated entry in the Active List.
	/////////////////////////////////////////////////////////////////////
	void renamer::set_complete(uint64_t AL_index){
		print("enter set complete function \n");
		print("active list index: " << AL_index << " set completed bit\n");
		active_list.at(AL_index).completed = true;
		print("exit set complete function \n");
	}

	/////////////////////////////////////////////////////////////////////
	// This function is for handling branch resolution.
	//
	// Inputs:
	// 1. AL_index: Index of the branch in the Active List.
	// 2. branch_ID: This uniquely identifies the branch and the
	//    checkpoint in question.  It was originally provided
	//    by the checkpoint function.
	// 3. correct: 'true' indicates the branch was correctly
	//    predicted, 'false' indicates it was mispredicted
	//    and recovery is required.
	//
	// Outputs: none.
	//
	// Tips:
	//
	// While recovery is not needed in the case of a correct branch,
	// some actions are still required with respect to the GBM and
	// all checkpointed GBMs:
	// * Remember to clear the branch's bit in the GBM.
	// * Remember to clear the branch's bit in all checkpointed GBMs.
	//
	// In the case of a misprediction:
	// * Restore the GBM from the branch's checkpoint. Also make sure the
	//   mispredicted branch's bit is cleared in the restored GBM,
	//   since it is now resolved and its bit and checkpoint are freed.
	// * You don't have to worry about explicitly freeing the GBM bits
	//   and checkpoints of branches that are after the mispredicted
	//   branch in program order. The mere act of restoring the GBM
	//   from the checkpoint achieves this feat.
	// * Restore the RMT using the branch's checkpoint.
	// * Restore the Free List head pointer and its phase bit,
	//   using the branch's checkpoint.
	// * Restore the Active List tail pointer and its phase bit
	//   corresponding to the entry after the branch's entry.
	//   Hints:
	//   You can infer the restored tail pointer from the branch's
	//   AL_index. You can infer the restored phase bit, using
	//   the phase bit of the Active List head pointer, where
	//   the restored Active List tail pointer is with respect to
	//   the Active List head pointer, and the knowledge that the
	//   Active List can't be empty at this moment (because the
	//   mispredicted branch is still in the Active List).
	// * Do NOT set the branch misprediction bit in the Active List.
	//   (Doing so would cause a second, full squash when the branch
	//   reaches the head of the Active List. We don’t want or need
	//   that because we immediately recover within this function.)
	/////////////////////////////////////////////////////////////////////
	void renamer::resolve(uint64_t AL_index,
		     uint64_t branch_ID,
		     bool correct){
					 print("enter resolve function \n");
				// if (correct){
				// 	print("resolve function :: branch was predicted correct" << "\n");
				// 	clear_gbm(branch_ID);

				// 	for (uint64_t i=0; i < this->n_branches; i++)
				// 	{
				// 		gbm_mask.at(i).gbm_snapshot = (gbm_mask.at(i).gbm_snapshot & ~(1ULL << branch_ID));
				// 	}
					

				// }
				// else{
				// 	print("resolve function :: branch was predicted *incorrect*\n");
				// 	GBM = (gbm_mask.at(branch_ID).gbm_snapshot & ~(1ULL << branch_ID));
				// 	assert(!(GBM & (1ULL << branch_ID)));

				// 	rename_map_table = gbm_mask.at(branch_ID).shadow_table;
				// 	free_list.head = gbm_mask.at(branch_ID).fl_head;
				// 	free_list.head_phase = gbm_mask.at(branch_ID).fl_head_phase;

				// 	al_tail= (AL_index + 1) % n_active;
					
				// 	if (al_head == al_tail  || al_tail < al_head)
				// 	{
				// 		print("resolve function :: active list tail phase updated to be different from active list head phase\n");
				// 		al_tail_phase = !al_head_phase;
				// 	} 
				// 	else {
				// 		assert(!(al_head >= al_tail));
				// 		print("resolve function :: active list tail phase updated to be same as active list head phase\n");
				// 		al_tail_phase = al_head_phase;
				// 	}
				// }
				//end of outside if else statement
				print("exit resolve function \n");
			 }

	//////////////////////////////////////////
	// Functions related to Retire Stage.   //
	//////////////////////////////////////////

	///////////////////////////////////////////////////////////////////
	// This function allows the caller to examine the instruction at the head
	// of the Active List.
	//
	// Input arguments: none.
	//
	// Return value:
	// * Return "true" if the Active List is NOT empty, i.e., there
	//   is an instruction at the head of the Active List.
	// * Return "false" if the Active List is empty, i.e., there is
	//   no instruction at the head of the Active List.
	//
	// Output arguments:
	// Simply return the following contents of the head entry of
	// the Active List.  These are don't-cares if the Active List
	// is empty (you may either return the contents of the head
	// entry anyway, or not set these at all).
	// * completed bit
	// * exception bit
	// * load violation bit
	// * branch misprediction bit
	// * value misprediction bit
	// * load flag (indicates whether or not the instr. is a load)
	// * store flag (indicates whether or not the instr. is a store)
	// * branch flag (indicates whether or not the instr. is a branch)
	// * amo flag (whether or not instr. is an atomic memory operation)
	// * csr flag (whether or not instr. is a system instruction)
	// * program counter of the instruction
	/////////////////////////////////////////////////////////////////////
	bool renamer::precommit(bool &completed,
                       bool &exception, bool &load_viol, bool &br_misp, bool &val_misp,
	               bool &load, bool &store, bool &branch, bool &amo, bool &csr,
		       uint64_t &PC){

			 print("enter precommit function \n");
			 print("al_head = " << al_head << "al_tail = " << al_tail << " al_head_phase = " << al_head_phase << " al_tail_phase = " << al_tail_phase << "\n");
			 print("completed=" << completed  << " exception=" << exception  << " load_viol=" << load_viol << " br_misp=" << br_misp << "\n");
	
		
			if (al_head == al_tail){
				if (al_head_phase == al_tail_phase){
					//active list is empty
					// dont cares
					    completed = false;
						exception = false;
						load_viol = false;
						br_misp = false;
						val_misp = false;
						load = false;
						store = false;
						branch = false;
						amo = false;
						csr = false;
						PC = 0;
					print("precommit function :: active list is empty so return false\n");
					print("exit precommit function \n");
					return false;
				}
				else{
					// this means the active list is full
					completed = active_list.at(al_head).completed;
					exception = active_list.at(al_head).threw_exception;
					load_viol = active_list.at(al_head).load_violation;
					br_misp = active_list.at(al_head).b_mispred;
					val_misp = active_list.at(al_head).v_mispred;
					load = active_list.at(al_head).is_load;
					store = active_list.at(al_head).is_store;
					branch = active_list.at(al_head).is_branch;
					amo = active_list.at(al_head).is_amo;
					csr = active_list.at(al_head).is_csr;
					PC = active_list.at(al_head).PC;
					 print("precommit function :: active list is full so return true\n");
					 print("exit precommit function \n");
					return true;
				}
			}
			else {
					// active list is not empty
					completed = active_list.at(al_head).completed;
					exception = active_list.at(al_head).threw_exception;
					load_viol = active_list.at(al_head).load_violation;
					br_misp = active_list.at(al_head).b_mispred;
					val_misp = active_list.at(al_head).v_mispred;
					load = active_list.at(al_head).is_load;
					store = active_list.at(al_head).is_store;
					branch = active_list.at(al_head).is_branch;
					amo = active_list.at(al_head).is_amo;
					csr = active_list.at(al_head).is_csr;
					PC = active_list.at(al_head).PC;
					//  std::cout << "precommit function :: active list is not empty so return true\n";
					//  std::cout << "exit precommit function \n";
				    return true;

			}
			

	}

	/////////////////////////////////////////////////////////////////////
	// This function commits the instruction at the head of the Active List.
	//
	// Tip (optional but helps catch bugs):
	// Before committing the head instruction, assert that it is valid to
	// do so (use assert() from standard library). Specifically, assert
	// that all of the following are true:
	// - there is a head instruction (the active list isn't empty)
	// - the head instruction is completed
	// - the head instruction is not marked as an exception
	// - the head instruction is not marked as a load violation
	// It is the caller's (pipeline's) duty to ensure that it is valid
	// to commit the head instruction BEFORE calling this function
	// (by examining the flags returned by "precommit()" above).
	// This is why you should assert() that it is valid to commit the
	// head instruction and otherwise cause the simulator to exit.
	/////////////////////////////////////////////////////////////////////
	void renamer::commit(){
		print("enter commit function \n");
		assert(al_head != al_tail || al_head_phase != al_tail_phase); // assert active list isn't empty
		assert(active_list.at(al_head).completed);
	    assert(!active_list.at(al_head).threw_exception);
		assert(!active_list.at(al_head).load_violation);

		if (active_list.at(al_head).has_dest){
			free_list.phys_reg.at(free_list.tail) = arch_map_table.at(active_list.at(al_head).log_reg);
			free_list.tail = (free_list.tail + 1) % free_list_size;
			
			if (free_list.tail == 0){
			print(" commit function :: free list tail pointer wrapped around, so update phase bit\n");
			free_list.tail_phase = !free_list.tail_phase;
			}

			arch_map_table.at(active_list.at(al_head).log_reg) = active_list.at(al_head).phys_reg;
	  }

		al_head = (al_head + 1) % this->n_active;
		if (al_head == 0){
		print(" commit function :: active list head pointer wrapped around, so update phase bit\n");
		al_head_phase = !al_head_phase;
		}
		print("exit commit function \n");
	}

	//////////////////////////////////////////////////////////////////////
	// Squash the renamer class.
	//
	// Squash all instructions in the Active List and think about which
	// sructures in your renamer class need to be restored, and how.
	//
	// After this function is called, the renamer should be rolled-back
	// to the committed state of the machine and all renamer state
	// should be consistent with an empty pipeline.
	/////////////////////////////////////////////////////////////////////
	void renamer::squash(){
		print("enter squash function \n");
		rename_map_table = arch_map_table; 

		free_list.head = free_list.tail;
		free_list.head_phase = !free_list.tail_phase;


		al_tail = al_head;
		al_tail_phase = al_head_phase;
		for (uint64_t i = 0; i < this->n_active; i++){
			active_list.at(i).has_dest = false;
			active_list.at(i).log_reg = 0;
			active_list.at(i).phys_reg = 0;
			active_list.at(i).completed = false;
			active_list.at(i).threw_exception = false;
			active_list.at(i).load_violation = false;
			active_list.at(i).b_mispred = false;
			active_list.at(i).v_mispred = false;
			active_list.at(i).is_load = false;
			active_list.at(i).is_store = false;
			active_list.at(i).is_branch = false;
			active_list.at(i).is_amo = false;
			active_list.at(i).is_csr = false;
			active_list.at(i).PC = 0;
		}

		// for (uint64_t i = 0; i < this->n_phys_regs; i++){
		// 	PRF_ready.at(i) = 0;
		// }

		//assert that the commited state are ready 
		 for (uint64_t i = 0; i < n_log_regs; i++){
		 	PRF_ready.at(arch_map_table.at(i)) = 1;
		 }

		GBM = 0ULL;
		for (uint64_t i = 0; i < this->n_branches; i++){
			gbm_mask.at(i).gbm_snapshot = 0ULL;
			gbm_mask.at(i).fl_head = 0ULL;
			gbm_mask.at(i).fl_head_phase = 0ULL;
			gbm_mask.at(i).shadow_table.clear();
		}

		print("exit squash function \n");
	}

	//////////////////////////////////////////
	// Functions not tied to specific stage.//
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// Functions for individually setting the exception bit,
	// load violation bit, branch misprediction bit, and
	// value misprediction bit, of the indicated entry in the Active List.
	/////////////////////////////////////////////////////////////////////
	void renamer::set_exception(uint64_t AL_index){active_list.at(AL_index).threw_exception = true; print("set exception function\n"); }
	void renamer::set_load_violation(uint64_t AL_index){active_list.at(AL_index).load_violation = true; print("set load violation function\n"); }
	void renamer::set_branch_misprediction(uint64_t AL_index){active_list.at(AL_index).b_mispred = true; print("set branch misprediction function\n"); }
	void renamer::set_value_misprediction(uint64_t AL_index){active_list.at(AL_index).v_mispred = true; print("set value misprediction function\n"); }

	/////////////////////////////////////////////////////////////////////
	// Query the exception bit of the indicated entry in the Active List.
	/////////////////////////////////////////////////////////////////////
	bool renamer::get_exception(uint64_t AL_index){
		print("enter get exception function \n");
		return active_list.at(AL_index).threw_exception;
	}