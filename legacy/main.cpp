#include "binning.hpp"
#include <iostream>
#include <string>

int main(){
    TJ myTJ("Cubic.config_bak");
    myTJ.cal_idxcell();
    myTJ.write_idxcell("idxcell_old_cub.txt");
    myTJ.find_neighbor();
    myTJ.write_neighborInfo("Tranbs_old_cub.txt");
    myTJ.makeA();
    myTJ.writeA();
    return 0;
}
