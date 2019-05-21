#include "math.hpp"
#include "config.hpp"
#include "tj.hpp"
#include "subdomain.hpp"
#include <iostream>

int main(int argc, char** argv){
    if(argc != 2){
        std::cout<<"Usage: TJ2GBE.out PATH/TO/CONFIG/FILE"<<std::endl;
        return 1;
    }
    TJ myTJ(argv[1]);
    myTJ.cal_idxcell();
    myTJ.write_idxcell(); // debug purpose
    myTJ.find_neighbor();
    myTJ.subdomain.write_cellInfo(); // GB ids in each cell
    myTJ.write_neighborInfo(); // to varify the threshold value is proper
    myTJ.make_A();
    myTJ.write_A(); // input for python script for energy reconstruction
    return 0;
}
