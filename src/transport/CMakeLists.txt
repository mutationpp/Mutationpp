cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    CollisionDB.cpp
    ThermalConductivityCG.cpp
    ThermalConductivityWilke.cpp
    ViscosityCG.cpp
    ViscosityGuptaYos.cpp
    ViscosityWilke.cpp
)

install(FILES CollisionDB.h DESTINATION include/mutation++)
install(FILES ThermalConductivityAlgorithm.h DESTINATION include/mutation++)
install(FILES ViscosityAlgorithm.h DESTINATION include/mutation++)
install(FILES Ramshaw.h DESTINATION include/mutation++)
install(FILES Transport.h DESTINATION include/mutation++)
