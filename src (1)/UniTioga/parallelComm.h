#pragma once
#include"codetypes.h"
#include"mpi.h"
class parallelComm
{
private:
    int nsend;
    int nrecv;
    int* sndMap;
    int* rcvMap;

public:
    int myid;
    int numprocs;
    MPI_Comm scomm;

    parallelComm() { sndMap = NULL; rcvMap = NULL; }

    ~parallelComm() {
        if (sndMap) free(sndMap);
        if (rcvMap) free(rcvMap);
    }

    void sendRecvPacketsAll(PACKET* sndPack, PACKET* rcvPack);

    void sendRecvPackets(PACKET* sndPack, PACKET* rcvPack);

    void sendRecvPacketsCheck(PACKET* sndPack, PACKET* rcvPack);

    void setMap(int ns, int nr, int* snd, int* rcv);

    void getMap(int* ns, int* nr, int** snd, int** rcv);

    void initPackets(PACKET* sndPack, PACKET* rcvPack);

    void clearPackets(PACKET* sndPack, PACKET* rcvPack);
    void clearPackets2(PACKET* sndPack, PACKET* rcvPack);

};
