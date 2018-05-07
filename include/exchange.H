class GridData;
class PhysicsData;

extern void exchange_grid(GridData&, const PhysicsData&, const int);
extern void exchange_B(GridData&);
extern void exchange_single(const GridData&,double*);


