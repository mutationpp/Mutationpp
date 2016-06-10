#ifndef SHOCKRELATIONS_H
#define SHOCKRELATIONS_H

#include "Data.h"

class ShockRelations {
public:
    ShockRelations(){}
    virtual ~ShockRelations(){}

    virtual void applyShockRelations(const Data& l_data_before, Data& l_data_after ){ l_data_after = l_data_before; };
};

#endif /* SHOCKRELATIONS_H */
