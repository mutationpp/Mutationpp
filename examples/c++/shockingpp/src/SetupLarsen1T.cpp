#include "SetupLarsen1T.h"
#include "DataLarsen1T.h"
#include "Larsen1T.h"

SetupLarsen1T::SetupLarsen1T(){
      std::cout << "                                                                    " << std::endl;
      std::cout << "              _ _ __ __ ________   __         __ \\  |  /           " << std::endl;
      std::cout << "     _ _ __ __  _______________   /  \\_/\\ /\\_/  \\       -       " << std::endl;
      std::cout << "              _ _ __ __ ________  \\__/   V   \\__/     __          " << std::endl;
      std::cout << "                                                 /   /  \\          " << std::endl;
      std::cout << "   _      _   ___  ___ ___ _  _                      \\__/          " << std::endl;
      std::cout << "  | |    /_\\ | _ \\/ __| __| \\| |  LAgrangian               \\    " << std::endl;
      std::cout << "  | |__ / _ \\|   /\\__ \\ _|| .` |  Reactor for StrEams  \\  \\ \\ " << std::endl;
      std::cout << "  |____/_/ \\_\\_|_\\|___/___|_|\\_|  in Nonequilibrium     \\  \\  " << std::endl;
      std::cout << "                                                            \\      " << std::endl;
      std::cout << "                                                                    " << std::endl;
}

Data* SetupLarsen1T::getData(Mutation::Mixture& l_mix){ return new DataLarsen1T(l_mix); }
Problem* SetupLarsen1T::getProblem(Mutation::Mixture& l_mix, Data& l_data){ return new Larsen1T(l_mix, l_data); }

