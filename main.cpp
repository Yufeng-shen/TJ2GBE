#include "binning.hpp"
#include <iostream>
#include <string>

int main(){
    TJ myTJ("config_Hex.txt");
    myTJ.cal_idxcell();
    myTJ.find_neighbor();
    myTJ.makeA();
    myTJ.writeA();
//    myTJ.write_neighborInfo("Tranbs_cpp.txt");
//    myTJ.write_idxcell("idxcell_cpp.txt");
    return 0;
}
