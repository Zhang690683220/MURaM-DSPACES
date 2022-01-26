#ifndef __MURAM_DSERV_H
#define __MURAM_DSERV_H

#include<dspaces.h>

#define STEP_DONE -1
#define STEP_ERR -2

enum muram_var {
    EOS,
    DIAG,
    SOLUTION,
    TAU_SLICE,
    YZ_SLICE,
    XY_SLICE,
    XZ_SLICE,
    CORONA,
    ANALYZE_VP
};

struct ds_meta {
    int globiter;
    int nslvar; // use nslvar as DIAG_flag only for DIAG
    
    enum muram_var var;
};

inline void publish_meta(dspaces_client_t dsp, int globiter, enum muram_var var, int nslvar)
{
    static int ver = 0;
    struct ds_meta mdata;

    mdata.globiter = globiter;
    mdata.var = var;
    mdata.nslvar = nslvar;

    dspaces_put_meta(dsp, "muram_meta", ver, &mdata, sizeof(mdata));

    ver++;
}

#endif /* __MURAM_DSERV_H */
