
namespace Mutation {
    namespace Thermodynamics {

class StatetoState : public ThermoDB
{
public:







}; // class RrhoDB

// Register the RRHO model with the other thermodynamic databases
Utilities::Config::ObjectProvider<StatetoState, ThermoDB> stsDB("StatetoState");

    } // namespace Thermodynamics
} // namespace Mutation
