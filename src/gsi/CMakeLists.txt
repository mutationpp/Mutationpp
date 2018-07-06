cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    DiffusionVelocityCalculator.cpp
    GasSurfaceInteraction.cpp
    GSIReactionAblation.cpp
    GSIReactionCatalysis.cpp
    GSIRateLawGammaConst.cpp
    GSIRateLawGammaT.cpp
    GSIRateLawSublimation.cpp
    GSIRateManagerGamma.cpp
    GSIStoichiometryManager.cpp
    MassBlowingRateAblation.cpp
    MassBlowingRateNull.cpp
    SurfaceBalanceSolverMass.cpp
    SurfacePropertiesNull.cpp
    WallState.cpp
    WallProductionTermsEmpty.cpp
    WallProductionSurfaceChemistry.cpp
    WallProductionPyrolysis.cpp
)

install(FILES DiffusionVelocityCalculator.h DESTINATION include/mutation++)
install(FILES GasSurfaceInteraction.h DESTINATION include/mutation++)
install(FILES GSIRateLaw.h DESTINATION include/mutation++)
install(FILES GSIRateManager.h DESTINATION include/mutation++)
install(FILES GSIReaction.h DESTINATION include/mutation++)
install(FILES GSIStoichiometryManager.h DESTINATION include/mutation++)
install(FILES MassBlowingRate.h DESTINATION include/mutation++)
install(FILES SurfaceBalanceSolver.h DESTINATION include/mutation++)
install(FILES SurfaceProperties.h DESTINATION include/mutation++)
install(FILES WallProductionTerms.h DESTINATION include/mutation++)
install(FILES WallState.h DESTINATION include/mutation++)
