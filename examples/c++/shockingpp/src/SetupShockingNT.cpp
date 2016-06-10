#include "SetupShockingNT.h"
#include "DataShockingNT.h"
#include "RankineHugoniotNT.h"
#include "ShockingTTv.h"

Data* SetupShockingNT::getDataPreShock(Mutation::Mixture& l_mix, const std::string& l_free_stream_conditions){ return new DataShockingNT(l_mix, l_free_stream_conditions); }

Data* SetupShockingNT::getDataPostShock(Mutation::Mixture& l_mix){ return new DataShockingNT(l_mix); }

ShockRelations* SetupShockingNT::getShockRelations(Mutation::Mixture& l_mix){ return new RankineHugoniotNT(l_mix); }

Problem* SetupShockingNT::getProblem(Mutation::Mixture& l_mix, Data& l_data){ return new ShockingTTv(l_mix, l_data); }

