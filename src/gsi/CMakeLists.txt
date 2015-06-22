cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    BuilderGasSurfaceInteraction.cpp
    CatalysisRateLaws.cpp
    CatalysisRateManager.cpp
#    JacobianManagerGSI.cpp
    SurfaceProperties.cpp
    GasSurfaceInteraction.cpp
    GSIReaction.cpp
    GSIStoichiometryManager.cpp
    WallSolver.cpp
    WallSteadyStateCoverageSolver.cpp
)

install(FILES BuilderGasSurfaceInteraction.h DESTINATION include/mutation++)
install(FILES CatalysisRateLaws.h DESTINATION include/mutation++)
install(FILES CatalysisRateManager.h DESTINATION include/mutation++)
#install(FILES JacobianManagerGSI.h DESTINATION include/mutation++)
install(FILES SurfaceProperties.h DESTINATION include/mutation++)
install(FILES GasSurfaceInteraction.h DESTINATION include/mutation++)
install(FILES GSIReaction.h DESTINATION include/mutation++)
install(FILES GSIStoichiometryManager.h DESTINATION include/mutation++)
install(FILES WallSolver.h DESTINATION include/mutation++)
install(FILES WallSteadyStateCoverageSolver.h DESTINATION include/mutation++)
