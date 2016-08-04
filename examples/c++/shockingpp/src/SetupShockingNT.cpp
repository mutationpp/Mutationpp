#include "SetupShockingNT.h"
#include "DataShockingNT.h"
#include "RankineHugoniotNT.h"
#include "ShockingTTv.h"

SetupShockingNT::SetupShockingNT(){
          std::cout << "             _______                                                     " << std::endl;
          std::cout << "            (______()     ___ _  The       _   _               \\W_W/     " << std::endl;
          std::cout << "     \\W/___(______()()   / __| |_  ___  __| |_(_)_ _  __ _     /o o\\     " << std::endl;
          std::cout << "     /M\\   (__TNT_()()   \\__ \\ ' \\/ _ \\/ _| / / | ' \\/ _` |  _( <   )_   " << std::endl;
          std::cout << "            (______()    |___/_||_\\___/\\__|_\\_\\_|_||_\\__, |  \\\\\\ O ///   " << std::endl;
          std::cout << "                                               Code  |___/   _))\\_/((_   " << std::endl;
          std::cout << "                                                            / /     \\ \\  " << std::endl;
          std::cout << "                                                           /_/|     |\\_\\ " << std::endl;
}

Data* SetupShockingNT::getDataPreShock(Mutation::Mixture& l_mix, const std::string& l_free_stream_conditions){ return new DataShockingNT(l_mix, l_free_stream_conditions); }

Data* SetupShockingNT::getDataPostShock(Mutation::Mixture& l_mix){ return new DataShockingNT(l_mix); }

ShockRelations* SetupShockingNT::getShockRelations(Mutation::Mixture& l_mix){ return new RankineHugoniotNT(l_mix); }

Problem* SetupShockingNT::getProblem(Mutation::Mixture& l_mix, Data& l_data){ return new ShockingTTv(l_mix, l_data); }

