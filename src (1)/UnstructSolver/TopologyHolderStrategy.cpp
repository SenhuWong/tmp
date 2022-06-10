#include "TopologyHolderStrategy.h"
#include "UnstructIntegrator.h"
#include "mpi.h"
#include <iostream>
#include <cmath>
TopologyHolderStrategy::TopologyHolderStrategy(UnstructTopologyHolder *hder,int nequ)
    : d_hder(hder)
{
    d_dim = d_hder->d_dim;
    d_nmesh = d_hder->d_nmesh;
    d_NEQU = nequ;
    commSendBuffer = new double *[d_hder->relatedProcs.size()];
    commRecvBuffer = new double *[d_hder->relatedProcs.size()];
    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        commSendBuffer[i] = new double[d_hder->sendCellPtr(i, d_nmesh) * d_NEQU];
        commRecvBuffer[i] = new double[d_hder->recvCellPtr(i, d_nmesh) * d_NEQU];
    }
    std::cout<<"Constructor of TopologyHoldersTRATEGy finished\n";
}

TopologyHolderStrategy::~TopologyHolderStrategy()
{
}

void TopologyHolderStrategy::initializeConsData(double** W)
{
    const double Gamma = 1.4;
    double **U = getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0; k< d_hder->nCells(i); k++)
        {
            double VMSquare = 0;
            for (int l = 0; l < d_dim; l++)
            {
                VMSquare += U[i][l + 2+d_NEQU*k] * U[i][l + 2+d_NEQU*k];
                W[i][l + 2+d_NEQU*k] = U[i][0+d_NEQU*k] * U[i][l + 2+d_NEQU*k];
            }
            W[i][0+d_NEQU*k] = U[i][0+d_NEQU*k];
            W[i][1+d_NEQU*k] = U[i][1+d_NEQU*k] / (Gamma - 1) + 0.5 * U[i][0+d_NEQU*k] * VMSquare;
        }
    }
}

void TopologyHolderStrategy::UpdatePrimitiveVariables(double** W)
{
    const double Gamma = 1.4;
    double** U = getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            U[i][0+d_NEQU*k] = W[i][0+d_NEQU*k];
            double VMSquare = 0;
            for (int l = 0; l < d_dim; l++)
            {
                U[i][l + 2+d_NEQU*k] = W[i][l + 2+d_NEQU*k] / U[i][0+d_NEQU*k];
                VMSquare += U[i][l + 2+d_NEQU*k] * U[i][l + 2+d_NEQU*k];            
            }
            U[i][1+d_NEQU*k] = (Gamma - 1) * (W[i][1+d_NEQU*k] - 0.5 * U[i][0+d_NEQU*k] * VMSquare);
        }
    }


}

void TopologyHolderStrategy::UpdateConservativeForRecvCells(double** W)
{
    double Gamma = 1.4;
    double** U = getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(auto iter = d_hder->recvBegin(i);iter!=d_hder->recvEnd(i);iter++)
        {
            int cellId = iter->d_localId;
            W[i][0+d_NEQU*cellId] = U[i][0+d_NEQU*cellId];
            W[i][2+d_NEQU*cellId] = U[i][0+d_NEQU*cellId]*U[i][2+d_NEQU*cellId];
            W[i][3+d_NEQU*cellId] = U[i][0+d_NEQU*cellId]*U[i][3+d_NEQU*cellId];
            W[i][1+d_NEQU*cellId] = U[i][1+d_NEQU*cellId]/(Gamma - 1) 
            + 0.5*U[i][0+d_NEQU*cellId]*(powf64(U[i][2+d_NEQU*cellId],2)+powf64(U[i][3+d_NEQU*cellId],2));
        }
    }
}

void TopologyHolderStrategy::ConservativeParam2ConvectiveFlux(double* W, double* Fc, GeomElements::vector3d<2,double>& norm_vec)
{
    double Gamma = 1.4;
    GeomElements::vector3d<2,double> velocity(W[2]/W[0],W[3]/W[0]);
    double Vn = velocity.dot_product(norm_vec);
    double DatCell = W[0];
    double PatCell = (Gamma-1)*(W[1] - 0.5*DatCell*velocity.L2Square());
    double rhoEatCell = W[1];
    Fc[0] = W[0]*Vn;
    Fc[1] = (W[1] + PatCell)*Vn;
    Fc[2] = W[2]*Vn + PatCell*norm_vec[0];
    Fc[3] = W[3]*Vn + PatCell*norm_vec[1];
}

void TopologyHolderStrategy::SolveViscousFlux(int curMesh,int curCell,int curEdge,int anotherCell,double* detlaW,double* Fv)
{
    const double Gamma = 1.4;
    const double OMEGAN = 1.5;
    auto& curBlk = d_hder->blk2D[curMesh];
    double rA;
    double** U =getU();
    auto& curE = curBlk.d_localEdges[curEdge];
    auto& curC = curBlk.d_localCells[curCell];
    double edgeAreaSquare = curE.area()*curE.area();
    GeomElements::vector3d<2,double> centerLine;
    //It is impossible to have anotherCell at first negative to be raised amid 0 and nCells 
    if(anotherCell>= d_hder->nCells(curMesh))
    {
        throw std::runtime_error("Oh Look what you have done\n");
    }
    // std::cout<<"Viscous Flux computed\n";
    if(anotherCell >=0)
    {
        auto& anotherC = curBlk.d_localCells[anotherCell];
        centerLine = anotherC.center() - curC.center();
    }
    else
    {
        centerLine = curE.center() - curC.center();
    }
    if(anotherCell<0)
    {
        throw std::runtime_error("Another Cell is below 0\n");
    }
    GeomElements::vector3d<2,double> curVelocity(&U[curMesh][2+d_NEQU*anotherCell]);
    GeomElements::vector3d<2,double> norm = curE.normal_vector();
    double Vn = curVelocity.dot_product(norm);
    double c = sqrtf64(Gamma * U[curMesh][1+d_NEQU*anotherCell] / U[curMesh][0+d_NEQU*anotherCell]);
    
    rA  = OMEGAN*(fabsf64(Vn)+c);
    if(!isInvicid())
    {
        rA += std::max<double>(4/(3*U[curMesh][0+d_NEQU*anotherCell]),Gamma/U[curMesh][0+d_NEQU*anotherCell])*
        (getMuOverPrCell(curMesh,curCell))/sqrtf64(centerLine.L2Square());
    }
    for(int j = 0;j<d_NEQU;j++)
    {
        Fv[j] = rA*detlaW[j];
    }
}

void TopologyHolderStrategy::SolveDiag(double** de_diag)
{
    double** spectrum_convec = getSpectrumConvective();
    double OMEGAN = 1.5;
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            auto& curCell = curBlk.d_localCells[k];
            double diag = curCell.volume()/d_stableDt[i][k] + 0.5*OMEGAN*spectrum_convec[i][k];
            if(!isInvicid())
            {
                double ** spectrum_viscous = getSpectrumViscous();
                // std::cout<<"diag before "<<diag<<'\t';
                diag += spectrum_viscous[i][k];
                // std::cout<<"diag after "<<diag<<'\n';
            }
            for(int j = 0;j<d_NEQU;j++)
            {
                de_diag[i][j+d_NEQU*k] = diag;
            }
        }
    }

}

void TopologyHolderStrategy::SolveDF(int curMesh,int curCell,int curEdge,
            double** W_scratch ,double** deltaW0,double* DF,
            int ForB,int cellOffset)
{
    auto& edg = d_hder->blk2D[curMesh].d_localEdges[curEdge];
    double Wp[d_NEQU];
    double W[d_NEQU];
    double dW[d_NEQU];
    double W0[d_NEQU];
    double Fcp[d_NEQU];
    double Fc[d_NEQU];

    double rA;
    double Fv[d_NEQU];
    double DFc[d_NEQU];
    int lC = edg.lCInd();
    int rC = edg.rCInd();
    if(curCell + cellOffset==lC)//Current cell at left side
    {
        bool dothis = false;
        if(rC>=0 and ((int(rC<lC))==(ForB)))
        {
            //If ForB is 1, then calculate when curCell(lC) >the other cell(rC)
            //If ForB is 0, then calculate when curCell(lC) < the other cell(rC)
            if(rC >= d_hder->nCells(curMesh))
            {
                rC = rC - d_hder->nCells(curMesh);
            }
            for(int j = 0;j<d_NEQU;j++)
            {
                W0[j] = W_scratch[curMesh][j+d_NEQU*rC];
                dW[j] = deltaW0[curMesh][j+d_NEQU*rC];
                Wp[j] =W0[j] + dW[j];
                W[j] = W0[j];
            }
            //Get convectiveFlux
            GeomElements::vector3d<2,double> norm = edg.normal_vector();
            ConservativeParam2ConvectiveFlux(W,Fc,norm);
            ConservativeParam2ConvectiveFlux(Wp,Fcp,norm);
            SolveViscousFlux(curMesh,curCell,curEdge,rC,dW,Fv);
            

            //Newly added:compute Viscous FLux
            //i is lc and j is rc
            //Get left cell's Vn
            

            for(int j = 0;j<d_NEQU;j++)
            {
                DFc[j] = Fcp[j] - Fc[j] ;
            }
            for(int j = 0;j<d_NEQU;j++)
            {
                DF[j] = 0.5*(DFc[j] - Fv[j])*edg.area();
            }
        }
        else
        {
            for(int j = 0;j<d_NEQU;j++)
            {
                DF[j] = 0;
            }
        }
    }
    else if(curCell + cellOffset==rC)
    {
        if((int(lC<rC)) == (ForB))
        {
            if(lC>=d_hder->nCells(curMesh))
            {
                lC = lC - d_hder->nCells(curMesh);
            }
            for(int j = 0;j<d_NEQU;j++)
            {
                W0[j] = W_scratch[curMesh][j+d_NEQU*lC];
                dW[j] = deltaW0[curMesh][j+d_NEQU*lC];
                Wp[j] = W0[j] + dW[j];
                W[j] = W0[j];
            }
            //Get convectiveFlux
            GeomElements::vector3d<2,double> norm = edg.normal_vector();
            ConservativeParam2ConvectiveFlux(W,Fc,norm);
            ConservativeParam2ConvectiveFlux(Wp,Fcp,norm);
            SolveViscousFlux(curMesh,curCell,curEdge,lC,dW,Fv);
            for(int j = 0;j<d_NEQU;j++)
            {
                DFc[j] =  Fc[j] - Fcp[j];
            }
            for(int j = 0;j<d_NEQU;j++)
            {
                DF[j] = 0.5*(DFc[j] - Fv[j])*edg.area();
            }
        }
        else
        {
            for(int j = 0;j<d_NEQU;j++)
            {
                DF[j] = 0;
            }
        }
    }
    else
    {
        std::cout<<curCell<<" with offset "<< cellOffset <<" not found in "<<lC<<" and "<<rC<<'\n';
        throw std::runtime_error("Cur cell is not on either side,impossible\n");
    }
}

void TopologyHolderStrategy::set_fs_variable2(double Ma, double AOA, double density, double pressure, double Re)
{
    fs_mach = Ma;
    std::cout<<"fs_mach "<<fs_mach<<'\n';
    fs_AOA = PI * AOA/180;
    std::cout<<"fs_AOA "<<fs_AOA<<'\n';
    fs_density = density;
    std::cout<<"fs_density "<<fs_density<<'\n';
    fs_pressure = pressure;
    std::cout<<"fs_pressure "<<fs_pressure<<'\n';
    fs_Re = Re;
    std::cout<<"fs_Re "<<fs_Re<<'\n';
    fs_soundSpeed = sqrtf64(Gamma*fs_pressure/fs_density);
    std::cout<<"fs_soundSpeed "<<fs_soundSpeed<<'\n';
    fs_Temperature = fs_pressure/(R*fs_density);
    std::cout<<"fs_Temperature "<<fs_Temperature<<'\n';
    fs_velocity_magnitude = fs_soundSpeed*fs_mach;
    std::cout<<"fs_velocity_magnitude "<<fs_velocity_magnitude<<'\n';
    fs_mu = 1.45e-6*pow(fs_Temperature,1.5)/(fs_Temperature+110);
    std::cout<<"fs_mu "<<fs_mu<<'\n';
    fs_velocity_components[0] = fs_velocity_magnitude * cosf64(fs_AOA);
    fs_velocity_components[1] = fs_velocity_magnitude * sinf64(fs_AOA);
    std::cout<<fs_velocity_components[0]<<" "<<fs_velocity_components[1]<<'\n';
    fs_eigen_Length = fs_Re*fs_mu/(fs_density*fs_velocity_magnitude);
    std::cout<<"fs_eigen_Length "<<fs_eigen_Length<<'\n';
}

void TopologyHolderStrategy::set_fs_variable(double Ma, double AOA, double density, double pressure, double eigenLen)
{
    // Ma is 0.8
    fs_mach = Ma;
    // AOA is 1.25
    fs_AOA = PI * AOA / 180;

    // gamma is 1.4
    // density is 1.225
    // pressure is 101325.0

    fs_density = density;
    fs_pressure = pressure;
    fs_eigen_Length = eigenLen;
    
    fs_Temperature = fs_pressure/(fs_density*R);

    fs_mu = 1.45e-6*pow(fs_Temperature,1.5)/(fs_Temperature+110);
    
    fs_soundSpeed = sqrtf64(Gamma * fs_pressure / fs_density);
    fs_velocity_magnitude = fs_mach * fs_soundSpeed;

    fs_Re = fs_density*fs_velocity_magnitude*fs_eigen_Length/fs_mu;

    fs_velocity_components[0] = fs_velocity_magnitude * cosf64(fs_AOA);
    fs_velocity_components[1] = fs_velocity_magnitude * sinf64(fs_AOA);
}

void TopologyHolderStrategy::update_nondim_variable()
{
    // u v
    fsnd_velocity_components[0] = fs_velocity_components[0] / (fs_soundSpeed / sqrtf64(Gamma));
    fsnd_velocity_components[1] = fs_velocity_components[1] / (fs_soundSpeed / sqrtf64(Gamma));
    // p d e
    fsnd_pressure = fs_pressure / fs_pressure;
    fsnd_density = fs_density / fs_density;
    fsnd_soundSpeed = sqrtf64(Gamma);

    fs_primVar[0] = fsnd_density;
    fs_primVar[1] = fsnd_pressure;
    fs_primVar[2] = fsnd_velocity_components[0];
    fs_primVar[3] = fsnd_velocity_components[1];
}

void TopologyHolderStrategy::AllCellCommunication(double** value)
{
    // if (cur_proc == 0)
    //     std::cout << "Within block communication";
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {

        auto iter = d_hder->sendBegin(i); // cur_ub.d_to_send.begin();
        for (int j = 0; j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->sendCellPtr(j, i); k < d_hder->sendCellPtr(j, i + 1); k++)
            {

                for (int l = 0; l < d_NEQU; l++)
                {

                    commSendBuffer[j][d_NEQU * k + l] = value[i][l+d_NEQU*iter->d_localId];
                }
                iter++;
            }
        }
    }

    std::vector<MPI_Request> recvRequests;
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Status> recvStatus;
    recvRequests.resize(d_hder->relatedProcs.size());
    sendRequests.resize(d_hder->relatedProcs.size());
    recvStatus.resize(d_hder->relatedProcs.size());
    // recorder << "-----------------------------------------\n";

    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {

        MPI_Irecv(commRecvBuffer[i], d_NEQU * d_hder->recvCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(recvRequests[i]));
        MPI_Isend(commSendBuffer[i], d_NEQU * d_hder->sendCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(sendRequests[i]));
    }
    MPI_Waitall(d_hder->relatedProcs.size(), recvRequests.data(), recvStatus.data());

    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        auto iter = d_hder->recvBegin(i);
        for (int j = 0; j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->recvCellPtr(j, i); k < d_hder->recvCellPtr(j, i + 1); k++)
            {
                // recorder << iter->d_localId << '\n';
                for (int l = 0; l < d_NEQU; l++)
                {
                    value[i][l+d_NEQU*iter->d_localId] = commRecvBuffer[j][k * d_NEQU + l];
                    // recorder << commRecvBuffer[j][k * d_NEQU + l] << " ";
                }
                // recorder << '\n';
                iter++;
            }
        }
    }
}

void TopologyHolderStrategy::NearCellCommunication(double** value)
{
    int temp_int = 0;
    for (int i = 0;i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j< d_hder->relatedProcs.size();j++)
        {
            for(int k = d_hder->sendNearCellPtr(j,i);k<d_hder->sendNearCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->nearSend(i,temp_int++);
                for(int l = 0;l<d_NEQU;l++)
                {
                    commSendBuffer[j][d_NEQU *k +l] = value[i][l+d_NEQU*iter->d_localId];
                }
            }
        }
    }
    std::vector<MPI_Request> recvRequests;
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Status> recvStatus;
    recvRequests.resize(d_hder->relatedProcs.size());
    sendRequests.resize(d_hder->relatedProcs.size());
    recvStatus.resize(d_hder->relatedProcs.size());

    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        MPI_Irecv(commRecvBuffer[i],d_NEQU* d_hder->recvNearCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(recvRequests[i]));
        MPI_Isend(commSendBuffer[i],d_NEQU* d_hder->sendNearCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(sendRequests[i]));
    }
    MPI_Waitall(d_hder->relatedProcs.size(), recvRequests.data(), recvStatus.data());
    temp_int = 0;
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->recvNearCellPtr(j,i);k<d_hder->recvNearCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->nearRecv(i,temp_int++);
                for(int l = 0;l<d_NEQU;l++)
                {
                    value[i][l+d_NEQU*iter->d_localId] = commRecvBuffer[j][d_NEQU*k+l];
                }
                
            }
        }
    }

}

void TopologyHolderStrategy::RemoteCellCommunication(double** value)
{
    int temp_int = 0;
    for (int i = 0;i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j< d_hder->relatedProcs.size();j++)
        {
            for(int k = d_hder->sendRemoteCellPtr(j,i);k<d_hder->sendRemoteCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->remoteSend(i,temp_int++);
                for(int l = 0;l<d_NEQU;l++)
                {
                    commSendBuffer[j][d_NEQU *k +l] = value[i][l+d_NEQU*iter->d_localId];
                }
            }
        }
    }
    std::vector<MPI_Request> recvRequests;
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Status> recvStatus;
    recvRequests.resize(d_hder->relatedProcs.size());
    sendRequests.resize(d_hder->relatedProcs.size());
    recvStatus.resize(d_hder->relatedProcs.size());

    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        MPI_Irecv(commRecvBuffer[i],d_NEQU* d_hder->recvRemoteCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(recvRequests[i]));
        MPI_Isend(commSendBuffer[i],d_NEQU* d_hder->sendRemoteCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(sendRequests[i]));
    }
    MPI_Waitall(d_hder->relatedProcs.size(), recvRequests.data(), recvStatus.data());
    temp_int = 0;
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->recvRemoteCellPtr(j,i);k<d_hder->recvRemoteCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->remoteRecv(i,temp_int++);
                for(int l = 0;l<d_NEQU;l++)
                {
                    value[i][l+d_NEQU*iter->d_localId] = commRecvBuffer[j][d_NEQU*k+l];
                }                
            }
        }
    }
}


