#ifndef GSIREACTION_H
#define GSIREACTION_H

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReaction{

public:    
    GSIReaction(){ }
    virtual ~GSIReaction(){ }

    GSIReaction( const GSIReaction& l_gsi_reaction ){ }

protected:
    std::string m_formula;
    std::vector<int> m_reactants;
    std::vector<int> m_products;

};

//=========================================================

class CatalysisReaction : public GSIReaction {

public:
    CatalysisReaction(){ }
    virtual ~CatalysisReaction(){ }

    CatalysisReaction( const CatalysisReaction& l_catalysis_reaction ){ }

};

//=========================================================

class CatalysisReactionGamma : public CatalysisReaction {

public:
    CatalysisReactionGamma(){ }
    virtual ~CatalysisReactionGamma(){ }

    CatalysisReactionGamma( const CatalysisReactionGamma& l_catalysis_reaction ){ }

};

//=========================================================

class CatalysisReactionFiniteRateChemistry : public CatalysisReaction {

public:
    CatalysisReactionFiniteRateChemistry(){ }
    virtual ~CatalysisReactionFiniteRateChemistry(){ }

    CatalysisReactionFiniteRateChemistry( const CatalysisReactionFiniteRateChemistry& l_catalysis_reaction ){ }

};

//=========================================================

class AblationReaction : public GSIReaction {

public:
    AblationReaction(){ }
    virtual ~AblationReaction(){ }

    AblationReaction( const AblationReaction& l_ablation_reaction ){ }

};

//=========================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation





#endif // GSIREACTION_H
