#include "unap.hpp"
#include "unapMPI.hpp"
using namespace UNAP;


int main(){
    // 切分通信域
    COMM::init(NULL,NULL);

    Communicator &globalCommunicator = COMM::getGlobalComm();

    // unapMPI::initMPI();

    COUT <<"全局通信域的进程的个数: "<< globalCommunicator.getMySize() <<ENDL;


}