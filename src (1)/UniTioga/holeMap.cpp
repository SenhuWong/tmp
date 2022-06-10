#include "tioga.h"
#include <array>
#include <fstream>
#include <iostream>
void fillHoleMap(int dim, int *holeMap, int ix[3], int isym);
using namespace TIOGA;
void tioga::getHoleMap()
{
    std::vector<int> existWall(nblocks);

    double ds[3], dsmax, dsbox;
    //
    // get the local bounding box
    //
    std::vector<std::array<double, 6>> wbox(nblocks);
    int meshtag = -BIGINT; //std::numeric_limits<int>::lowest();
    int maxtag, mtagtmp;
    for (int i = 0; i < nblocks; i++)
    {
        auto &mb = mblocks[i];
        mb->getWallBounds(&mtagtmp, &existWall[i], wbox[i].data()); //If it is 2d:{xmin,ymin,xmax,ymax};
        if (mtagtmp > meshtag)
            meshtag = mtagtmp;
    }
    MPI_Allreduce(&meshtag, &maxtag, 1, MPI_INT, MPI_MAX, scomm);
    //
    if (holeMap)
    {
        for (int i = 0; i < nmesh; i++)
            if (holeMap[i].existWall)
                TIOGA_FREE(holeMap[i].sam);
        delete[] holeMap;
    }
    holeMap = new HOLEMAP[maxtag];
    //
    int *existHoleLocal = new int[maxtag];
    int *existHole = new int[maxtag];

    for (int i = 0; i < maxtag; i++)
        existHole[i] = existHoleLocal[i] = 0;
    //
    for (int i = 0; i < nblocks; i++)
    {
        existHoleLocal[mtags[i] - 1] = existWall[i];
    }
    //
    MPI_Allreduce(existHoleLocal, existHole, maxtag, MPI_INT, MPI_MAX, scomm);
    //
    for (int i = 0; i < maxtag; i++)
        holeMap[i].existWall = existHole[i];
    //
    double *bboxLocal = new double[2 * d_dim * maxtag];
    double *bboxGlobal = new double[2 * d_dim * maxtag];
    //
    for (int i = 0; i < d_dim * maxtag; i++)
        bboxLocal[i] = BIGVALUE;
    for (int i = 0; i < d_dim * maxtag; i++)
        bboxLocal[i + d_dim * maxtag] = -BIGVALUE;
    for (int i = 0; i < d_dim * maxtag; i++)
        bboxGlobal[i] = BIGVALUE;
    for (int i = 0; i < d_dim * maxtag; i++)
        bboxGlobal[i + d_dim * maxtag] = -BIGVALUE;

    //
    for (int n = 0; n < nblocks; n++)
    {
        meshtag = mtags[n];
        for (int i = 0; i < d_dim; i++)
        {
            bboxLocal[d_dim * (meshtag - 1) + i] = wbox[n][i];
            bboxLocal[d_dim * (meshtag - 1) + i + d_dim * maxtag] = wbox[n][i + d_dim];
        }
    }
    //
    // get the global bounding box info across all the
    // partitions for all meshes
    //
    MPI_Allreduce(bboxLocal, bboxGlobal, d_dim * maxtag, MPI_DOUBLE, MPI_MIN, scomm);
    MPI_Allreduce(
        &(bboxLocal[d_dim * maxtag]), &(bboxGlobal[d_dim * maxtag]), d_dim * maxtag, MPI_DOUBLE,
        MPI_MAX, scomm);
    //
    // find the bounding box for each mesh
    // from the globally reduced data
    //
    int bufferSize;
    for (int i = 0; i < maxtag; i++)
    {
        if (holeMap[i].existWall)
        {
            for (int j = 0; j < d_dim; j++)
            {
                holeMap[i].extents[j] = bboxGlobal[d_dim * i + j]; //
                holeMap[i].extents[j + d_dim] = bboxGlobal[d_dim * i + j + d_dim * maxtag];
                ds[j] = holeMap[i].extents[j + d_dim] - holeMap[i].extents[j];
            }
            dsmax = TIOGA_MAX(ds[0], ds[1]);
            dsmax = TIOGA_MAX(dsmax, ds[2]);
            dsbox = dsmax / HOLEMAPSIZE;

            for (int j = 0; j < d_dim; j++)
            {
                holeMap[i].extents[j] -= (2 * dsbox);
                holeMap[i].extents[j + d_dim] += (2 * dsbox);
                holeMap[i].nx[j] = floor(
                    TIOGA_MAX((holeMap[i].extents[j + d_dim] - holeMap[i].extents[j]) / dsbox, 1));
            }
            if (d_dim == 2)
            {
                bufferSize = holeMap[i].nx[0] * holeMap[i].nx[1];
            }
            else if (d_dim == 3)
            {
                bufferSize = holeMap[i].nx[0] * holeMap[i].nx[1] * holeMap[i].nx[2];
            }
            holeMap[i].sam = new int[bufferSize];
            holeMap[i].samLocal = new int[bufferSize];
            for (int j = 0; j < bufferSize; j++)
                holeMap[i].sam[j] = holeMap[i].samLocal[j] = 0;
        }
    }
    //
    // mark the wall boundary cells in the holeMap
    //
    for (int ib = 0; ib < nblocks; ib++)
    {
        auto &mb = mblocks[ib];
        meshtag = mb->getMeshTag();
        if (holeMap[meshtag - 1].existWall)
        {
            mb->markWallBoundary(
                holeMap[meshtag - 1].samLocal, holeMap[meshtag - 1].nx,
                holeMap[meshtag - 1].extents);
        }
    }
    //
    // allreduce the holeMap of each mesh
    //
    for (int i = 0; i < maxtag; i++)
    {
        if (holeMap[i].existWall)
        {
            if (d_dim == 2)
            {
                bufferSize = holeMap[i].nx[0] * holeMap[i].nx[1];
            }
            else if (d_dim == 3)
            {
                bufferSize = holeMap[i].nx[0] * holeMap[i].nx[1] * holeMap[i].nx[2];
            }
            MPI_Allreduce(holeMap[i].samLocal, holeMap[i].sam, bufferSize, MPI_INT, MPI_MAX, scomm);
        }
    }
    //
    for (int i = 0; i < maxtag; i++)
        if (holeMap[i].existWall)
            TIOGA_FREE(holeMap[i].samLocal);
    //
    // set the global number of meshes to maxtag
    //
    nmesh = maxtag;
    //
    // now fill the holeMap
    //
    for (int i = 0; i < maxtag; i++)
        if (holeMap[i].existWall)
            fillHoleMap(d_dim, holeMap[i].sam, holeMap[i].nx, isym);
    //
    // output the hole map
    //
    //this->outputHoleMap();
    //
    // free local memory
    //
    TIOGA_FREE(existHoleLocal);
    TIOGA_FREE(existHole);
    TIOGA_FREE(bboxLocal);
    TIOGA_FREE(bboxGlobal);
}
//2D ONLY NOW
void tioga::outputHoleMap(const std::string &filename)
{
    std::string localFilename;
    std::ofstream fout;
    double ds[3];
    if (d_dim == 2)
    {
        for (int i = 0; i < nmesh; i++)
        {
            if (holeMap[i].existWall)
            {
                localFilename = filename + std::to_string(i) + std::to_string(myid) + ".dat";
                fout.open(localFilename);
                if (!fout.is_open())
                {
                    std::cout << "Error at outputHoleMap\n";
                    return;
                }
                fout << "TITLE =\"Tioga output\"\n";
                fout << "VARIABLES=\"X\",\"Y\",\"IBLANK\"\n";
                int nnodes = (holeMap[i].nx[0] + 1) * (holeMap[i].nx[1] + 1);
                int ncells = (holeMap[i].nx[0]) * (holeMap[i].nx[1]);
                fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells << " ET=QUADRILATERAL, F=FEBLOCK\n";
                fout << "VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED)\n";
                for (int k = 0; k < d_dim; k++)
                    ds[k] = (holeMap[i].extents[k + d_dim] - holeMap[i].extents[k]) / (holeMap[i].nx[k]);
                for (int jj = 0; jj < holeMap[i].nx[1] + 1; jj++)
                    for (int ii = 0; ii < holeMap[i].nx[0] + 1; ii++)
                        fout << ii * ds[0] << '\n';

                for (int jj = 0; jj < holeMap[i].nx[1] + 1; jj++)
                    for (int ii = 0; ii < holeMap[i].nx[0] + 1; ii++)
                        fout << jj * ds[1] << '\n';

                int m = 0;

                for (int jj = 0; jj < holeMap[i].nx[1]; jj++)
                    for (int ii = 0; ii < holeMap[i].nx[0]; ii++)
                    {
                        fout << holeMap[i].sam[m] << '\n';
                        m++;
                    }

                m = 0;
                int ns1 = holeMap[i].nx[0] + 1;

                for (int jj = 0; jj < holeMap[i].nx[1]; jj++)
                    for (int ii = 0; ii < holeMap[i].nx[0]; ii++)
                    {
                        m = jj * ns1 + ii + 1;
                        fout << m << " " << m + 1 << " " << m + 1 + ns1 << " " << m + ns1 << '\n';
                    }
                fout.close();
            }
        }
    }
    else if (d_dim == 3)
    {
        for (int i = 0; i < nmesh; i++)
        {
            if (holeMap[i].existWall)
            {
                localFilename = filename + std::to_string(i) + std::to_string(myid) + ".dat";
                fout.open(localFilename);
                if (!fout.is_open())
                {
                    std::cout << "Error at outputHoleMap\n";
                    return;
                }
                fout << "TITLE =\"Tioga output\"\n";
                fout << "VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\"\n";
                int nnodes = (holeMap[i].nx[0] + 1) * (holeMap[i].nx[1] + 1) * (holeMap[i].nx[2] + 1);
                int ncells = (holeMap[i].nx[0]) * (holeMap[i].nx[1]) * (holeMap[i].nx[2]);
                fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells << " ET=BRICK, F=FEBLOCK\n";
                fout << "VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED)\n";
                for (int k = 0; k < d_dim; k++)
                    ds[k] = (holeMap[i].extents[k + d_dim] - holeMap[i].extents[k]) / (holeMap[i].nx[k]);
                for (int kk = 0; kk < holeMap[i].nx[2] + 1; kk++)
                    for (int jj = 0; jj < holeMap[i].nx[1] + 1; jj++)
                        for (int ii = 0; ii < holeMap[i].nx[0] + 1; ii++)
                            fout << ii * ds[0] << '\n';
                for (int kk = 0; kk < holeMap[i].nx[2] + 1; kk++)
                    for (int jj = 0; jj < holeMap[i].nx[1] + 1; jj++)
                        for (int ii = 0; ii < holeMap[i].nx[0] + 1; ii++)
                            fout << jj * ds[1] << '\n';
                for (int kk = 0; kk < holeMap[i].nx[2] + 1; kk++)
                    for (int jj = 0; jj < holeMap[i].nx[1] + 1; jj++)
                        for (int ii = 0; ii < holeMap[i].nx[0] + 1; ii++)
                            fout << kk * ds[2] << '\n';

                int m = 0;
                for (int kk = 0; kk < holeMap[i].nx[2]; kk++)
                    for (int jj = 0; jj < holeMap[i].nx[1]; jj++)
                        for (int ii = 0; ii < holeMap[i].nx[0]; ii++)
                        {
                            fout << holeMap[i].sam[m] << '\n';
                            m++;
                        }
                m = 0;
                int ns1 = holeMap[i].nx[0] + 1;
                int ns2 = ns1 * (holeMap[i].nx[1] + 1);
                for (int kk = 0; kk < holeMap[i].nx[2] + 1; kk++)
                    for (int jj = 0; jj < holeMap[i].nx[1]; jj++)
                        for (int ii = 0; ii < holeMap[i].nx[0]; ii++)
                        {
                            m = kk * ns2 + jj * ns1 + ii + 1;
                            fout << m << " " << m + 1 << " " << m + 1 + ns1 << " " << m + ns1 << " ";
                            fout << m + ns2 << " " << m + 1 + ns2 << " " << m + 1 + ns1 + ns2 << " " << m + ns1 + ns2 << '\n';
                        }
                fout.close();
            }
        }
    }
}

void tioga::outPutQuery(const std::string &filename)
{
    std::ofstream fout;

    for (int i = 0; i < nblocks; i++)
    {

        std::string localFilename = filename + std::to_string(mtags[i]) + std::to_string(myid) + ".dat";
        fout.open(localFilename);
        if (!fout.is_open())
        {
            std::cout << "OutputQuery Failure\n";
            return;
        }
        auto &mb = mblocks[i];
        std::cout << mb->nsearch << " nodes to search on proc " << myid << " block " << i << "\n";
        for (int j = 0; j < mb->nsearch; j++)
        {
            if (mb->donorId[j] > 0)
            {
                for (int k = 0; k < d_dim; k++)
                {
                    fout << mb->xsearch[d_dim * j + k] << " ";
                }
                fout << '\n';
            }
        }
        fout.close();
    }
}

void tioga::outPutCartQuery(const std::string &filename)
{
    std::ofstream fout;
    for (int i = 0; i < nblocks; i++)
    {

        std::string localFilename = filename + std::to_string(mtags[i]) + std::to_string(myid) + ".dat";
        fout.open(localFilename);
        if (!fout.is_open())
        {
            std::cout << "OutputQuery Failure\n";
            return;
        }
        auto &mb = mblocks[i];
        std::cout << mb->ntotalPointsCart << " nodes to give to cart on proc " << myid << " block " << i << "\n";
        for (int j = 0; j < mb->ntotalPointsCart; j++)
        {
            for (int k = 0; k < d_dim; k++)
            {
                fout << mb->rxyzCart[d_dim * j + k] << " ";
            }
            fout << '\n';
        }
        fout.close();
    }
}

// void tioga::outPutGoodQuery(const std::string& filename)
// {
//     std::ofstream fout;

//     for (int i = 0; i < nblocks; i++)
//     {

//         std::string localFilename = filename + std::to_string(mtags[i]) + std::to_string(myid) + ".dat";
//         fout.open(localFilename);
//         if (!fout.is_open())
//         {
//             std::cout << "OutputGoodQuery Failure\n";
//             return;
//         }
//         auto& mb = mblocks[i];
//         std::cout << mb->nsearch << " nodes to search\n";
//         for (int j = 0; j < mb->nsearch; j++)
//         {
//             if (mb->good_or_bad[j] == -1)
//             {
//                 for (int k = 0; k < d_dim; k++)
//                 {
//                     fout << mb->xsearch[d_dim * j + k] << " ";
//                 }
//                 fout << '\n';
//             }
//         }
//         fout.close();
//     }
// }