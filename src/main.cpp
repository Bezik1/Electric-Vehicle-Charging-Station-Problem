#include <print>
#include "gurobi_c++.h"

int main() 
{
    GRBEnv env = GRBEnv(true);
    env.start();
    GRBModel model = GRBModel(env);
    std::printf("Initiaization of EVCSP project!");

    return 0;
}