cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    JacobianManager.cpp
    Kinetics.cpp
    RateLaws.cpp
    RateManager.cpp
    Reaction.cpp
    ReactionType.cpp
    StoichiometryManager.cpp
)

install(FILES JacobianManager.h DESTINATION include/mutation++)
install(FILES Kinetics.h DESTINATION include/mutation++)
install(FILES RateLaws.h DESTINATION include/mutation++)
install(FILES RateLawGroup.h DESTINATION include/mutation++)
install(FILES RateManager.h DESTINATION include/mutation++)
install(FILES Reaction.h DESTINATION include/mutation++)
install(FILES ReactionType.h DESTINATION include/mutation++)
install(FILES StoichiometryManager.h DESTINATION include/mutation++)
install(FILES ThirdBodyManager.h DESTINATION include/mutation++)

