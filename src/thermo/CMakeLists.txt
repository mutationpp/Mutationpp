cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    ChemNonEqStateModel.cpp
    #ChemNonEqTTvStateModel.cpp
    EquilStateModel.cpp
    MultiPhaseEquilSolver.cpp
    Nasa7Polynomial.cpp
    Nasa9Polynomial.cpp
    Nasa7DB.cpp
    Nasa9DB.cpp
    ParticleRRHO.cpp
    RrhoDB.cpp
    Species.cpp
    SpeciesListDescriptor.cpp
    Thermodynamics.cpp
    ThermoDB.cpp
)

install(FILES GfcEquilSolver.h DESTINATION include/mutation++)
install(FILES MultiPhaseEquilSolver.h DESTINATION include/mutation++)
install(FILES Species.h DESTINATION include/mutation++)
install(FILES SpeciesListDescriptor.h DESTINATION include/mutation++)
install(FILES StateModel.h DESTINATION include/mutation++)
install(FILES ThermoDB.h DESTINATION include/mutation++)
install(FILES Thermodynamics.h DESTINATION include/mutation++)
