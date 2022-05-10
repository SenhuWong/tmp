#include"commCell.h"

//Increasing ranking: e1<e2
bool rankWithSource(const CommunicationCell &e1,const CommunicationCell &e2)
{
    //Proc of e1 < e2 or Proc equals but id of e1 < e2.
    bool WithinBool = (e1.d_localId<e2.d_localId) and (e1.d_remoteProc == e2.d_remoteProc);
    bool InbetwBool = e1.d_remoteProc < e2.d_remoteProc;
    return (WithinBool or InbetwBool);
}
bool rankWithDestin(const CommunicationCell &e1,const CommunicationCell &e2)
{
    bool WithinBool = (e1.d_remoteId<e2.d_remoteId) and (e1.d_remoteProc == e2.d_remoteProc);
    bool InbetwBool = e1.d_remoteProc < e2.d_remoteProc;
    return (WithinBool or InbetwBool);
}