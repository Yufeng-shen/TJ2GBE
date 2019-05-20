#include "math.hpp"
#include "config.hpp"
#include "tj.hpp"
#include "subdomain.hpp"

int main(){
    TJ myTJ("../CfgFile/Cubic.config");
    myTJ.cal_idxcell();
    myTJ.write_idxcell(); // debug purpose
    myTJ.find_neighbor();
    myTJ.write_neighborInfo(); // to varify the threshold value is proper
    myTJ.make_A();
    myTJ.write_A(); // input for python script
    return 0;
}
