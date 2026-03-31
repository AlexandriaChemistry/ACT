
#include "import.h"

namespace alexandria
{

void importFile(MsgHandler           *msg_handler,
                const ForceField     *pd,
                const std::string    &filenm,
                std::vector<MolProp> *mp,
                const char           *conf,
                JobType               jobtype,
                bool                  userqtot,
                double               *qtot,
                bool                  oneH)
{
    (void)msg_handler;
    (void)pd;
    (void)filenm;
    (void)mp;
    (void)conf;
    (void)jobtype;
    (void)userqtot;
    (void)qtot;
    (void)oneH;
}

} // namespace alexandria
