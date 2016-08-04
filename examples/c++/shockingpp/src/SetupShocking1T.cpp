#include "SetupShocking1T.h"
#include "DataShocking1T.h"
#include "RankineHugoniot1T.h"
#include "Shocking1T.h"

SetupShocking1T::SetupShocking1T(){
          std::cout << "             _______                                                     " << std::endl;
          std::cout << "            (______()     ___ _  The       _   _               \\W_W/     " << std::endl;
          std::cout << "     \\W/___(______()()   / __| |_  ___  __| |_(_)_ _  __ _     /o o\\     " << std::endl;
          std::cout << "     /M\\   (__TNT_()()   \\__ \\ ' \\/ _ \\/ _| / / | ' \\/ _` |  _( <   )_   " << std::endl;
          std::cout << "            (______()    |___/_||_\\___/\\__|_\\_\\_|_||_\\__, |  \\\\\\ O ///   " << std::endl;
          std::cout << "                                               Code  |___/   _))\\_/((_   " << std::endl;
          std::cout << "                                                            / /     \\ \\  " << std::endl;
          std::cout << "                                                           /_/|     |\\_\\ " << std::endl;
}
Data* SetupShocking1T::getDataPreShock(Mutation::Mixture& l_mix, const std::string& l_free_stream_conditions){ return new DataShocking1T(l_mix, l_free_stream_conditions); }

Data* SetupShocking1T::getDataPostShock(Mutation::Mixture& l_mix){ return new DataShocking1T(l_mix); }

ShockRelations* SetupShocking1T::getShockRelations(Mutation::Mixture& l_mix){ return new RankineHugoniot1T(l_mix); }

Problem* SetupShocking1T::getProblem(Mutation::Mixture& l_mix, Data& l_data){ return new Shocking1T(l_mix, l_data); }

