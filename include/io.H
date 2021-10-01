#include <string>
#include <vector>
#include <list>
#include "dspaces.h"

class RunData;
class GridData;
class PhysicsData;
class RTS;

// extern var for 2d slices use
extern dspaces_client_t ds_client;

// struct for IO logs
struct log_entry {
  char name[128];
	std::vector<int> rank_root;
  std::vector<int> iter;
  std::vector<double> api_time;
  std::vector<double> wait_time;
  std::vector<double> time;
  double total_api_time;
  double total_wait_time;
  double total_time;
  double avg_api_time;
  double avg_wait_time;
  double avg_time;
};

struct log {
	char name[128];
  struct log_entry* eos;
  struct log_entry* diag;
  struct log_entry* tau;
  struct log_entry* yz;
  struct log_entry* xy;
  struct log_entry* xz;
  struct log_entry* corona;
  struct log_entry* analyze_vp;
};

struct root_list {
  std::list<int> xy_root;
  std::list<int> xz_root;
  std::list<int> yz_root;
  std::list<int> xcol_root;
  std::list<int> ycol_root;
  std::list<int> zcol_root;
};

extern int  HasBackupFile(const char*);
extern void OutputSolution(const RunData&, const GridData&, const PhysicsData&);
extern void BackupSolution(const RunData&, const GridData&, const PhysicsData&);
extern void RestoreSolution(RunData&, GridData&, const PhysicsData&);
extern void OutputDerivatives(const RunData&, const GridData&, const PhysicsData&);
extern void RestoreDerivatives(RunData&, GridData&, const PhysicsData&);
extern void IO_Init(const GridData&, const RunData&);
extern void IO_Finalize();

extern void Iout(const RunData&, const GridData&, const PhysicsData&, RTS*);
extern void tau_slice(const RunData&, const GridData&, const PhysicsData&, RTS*);
extern void xy_slice(const RunData&, const GridData&, const PhysicsData&);
extern void yz_slice(const RunData&, const GridData&, const PhysicsData&, RTS*);
extern void xz_slice(const RunData&, const GridData&, const PhysicsData&);

extern void diag_output(const RunData&,const GridData&,const PhysicsData&, RTS*);
extern void eos_output(const RunData&,const GridData&,const PhysicsData&, RTS*);
extern void eos_output_gpu(const RunData&,const GridData&,const PhysicsData&, RTS*);
extern void corona_emission_dem_xyz(const RunData&  Run, const GridData& Grid, 
		               const PhysicsData& Physics);

extern void log_entry_init(struct log_entry* le, std::string name);

