cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    GfcEquilSolver.cpp
    Nasa7Polynomial.cpp
    Nasa9Polynomial.cpp
    NasaDB.cpp
    ParticleRRHO.cpp
    RrhoDB.cpp
    Species.cpp
    StateModel.cpp
    Thermodynamics.cpp
)

install(FILES GfcEquilSolver.h DESTINATION include/mutation++)
install(FILES Species.h DESTINATION include/mutation++)
install(FILES StateModel.h DESTINATION include/mutation++)
install(FILES ThermoDB.h DESTINATION include/mutation++)
install(FILES Thermodynamics.h DESTINATION include/mutation++)
