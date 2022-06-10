#include"MeshBlock.h"
#include<iostream>
#include<fstream>

void MeshBlock::writeGridFile2(const std::string& filename,int* nodeInfo)//This is 2d only for test
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("WriteGridFile failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if (d_dim == 3)
    {
        fout << "\"Z\",";
    }
    fout << "\"IBLANK\"\n";
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    if (d_dim == 2)
    {
        fout << " ET = QUADRILATERAL, F = FEPOINT\n";
    }
    else if (d_dim == 3)
    {
        fout << " ET = BRICK, F = FEPOINT\n";
    }
    
    //fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells << " ET = QUADRILATERAL, F = FEPOINT\n";
    for (int i = 0; i < nnodes; i++)
    {
        for (int j = 0; j < d_dim; j++)
        {
            fout << x[d_dim * i + j] << " ";
        }
        fout << nodeInfo[i] << '\n';
    }
    int ba = 1 - BASE;
    if (d_dim == 2)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 3:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << "\n";
                }
                break;
            case 4:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    else if (d_dim == 3)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 4://Tetrahedron [1233,4444]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            case 5://[1234 5555]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << "\n";
                }
                break;
            case 6://[1233 4566]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 5] + ba << "\n";
                }
                break;
            case 8:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 6] + ba << " " << vconn[n][nvert * i + 7] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    fout.close();
    return;

}

void MeshBlock::writeGridFile2(const std::string& filename, double* nodeInfo)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("WriteGridFile failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if (d_dim == 3)
    {
        fout << "\"Z\",";
    }
    fout<<"\"IBLANK\"\n";
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    if (d_dim == 2)
    {
        fout << " ET = QUADRILATERAL, F = FEPOINT\n";
    }
    else if (d_dim == 3)
    {
        fout << " ET = BRICK, F = FEPOINT\n";
    }
    for (int i = 0; i < nnodes; i++)
    {
        for (int j = 0; j < d_dim; j++)
        {
            fout << x[d_dim * i + j] << " ";
        }
        fout << nodeInfo[i] << '\n';
    }
    int ba = 1 - BASE;
    if (d_dim == 2)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 3:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << "\n";
                }
                break;
            case 4:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    else if (d_dim == 3)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 4://Tetrahedron [1233,4444]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            case 5://[1234 5555]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << "\n";
                }
                break;
            case 6://[1233 4566]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 5] + ba << "\n";
                }
                break;
            case 8:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 6] + ba << " " << vconn[n][nvert * i + 7] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    fout.close();
    return;

}

void MeshBlock::writeCellFile2(const std::string& filename,int * cellInfo)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("WriteCellFile failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if (d_dim == 3)
    {
        fout << "\"Z\",";
    }
    fout << "\"IBLANK_CELL\",\"meshTag\",\"procId\"\n";
    if (d_dim == 2)
    {
        fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E =" << ncells << " ET = QUADRILATERAL, F = FEBLOCK\n";
    }
    else if (d_dim == 3)
    {
        fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E =" << ncells << " ET = BRICK, F = FEBLOCK\n";
    }
    fout << "VARLOCATION =  (";
    for (int i = 0; i < d_dim; i++)
    {
        fout << i + 1 << "=NODAL,";
    }
    fout << d_dim + 1 << "=CELLCENTERED,"<< d_dim + 2 << "=CELLCENTERED,"<< d_dim + 3 << "=CELLCENTERED)\n";
    for (int j = 0; j < d_dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
        {
            fout << x[d_dim * i + j] << '\n';
        }
    }
    for (int i = 0; i < ncells; i++)
    {
        fout << cellInfo[i] << '\n';
    }
    for (int i = 0; i < ncells; i++)
    {
        fout << meshtag<<'\n';
    }
    for (int i = 0; i < ncells; i++)
    {
        fout << myid<<'\n';
    }

    int ba = 1 - BASE;
    if (d_dim == 2)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 3:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << "\n";
                }
                break;
            case 4:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    else if (d_dim == 3)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 4://Tetrahedron [1233,4444]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            case 5://[1234 5555]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << "\n";
                }
                break;
            case 6://[1233 4566]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 5] + ba << "\n";
                }
                break;
            case 8:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 6] + ba << " " << vconn[n][nvert * i + 7] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    fout.close();
    return;
}

void MeshBlock::writeCellFile2(const std::string& filename, double* cellInfo)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("WriteCellFile failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if (d_dim == 3)
    {
        fout << "\"Z\",";
    }
    fout << "\"IBLANK_CELL\",\"ProcId\"\n";
    if (d_dim == 2)
    {
        fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E =" << ncells << " ET = QUADRILATERAL, F = FEBLOCK\n";
    }
    else if (d_dim == 3)
    {
        fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E =" << ncells << " ET = BRICK, F = FEBLOCK\n";
    }
    fout << "VARLOCATION =  (";
    for (int i = 0; i < d_dim; i++)
    {
        fout << i + 1 << "=NODAL,";
    }
    fout << d_dim + 1 << "=CELLCENTERED,\n";
    fout << d_dim + 2 << "=CELLCENTERED)\n";
    for (int j = 0; j < d_dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
        {
            fout << x[d_dim * i + j] << '\n';
        }
    }
    for (int i = 0; i < ncells; i++)
    {
        fout << cellInfo[i] << '\n';
    }
    for(int i = 0;i<ncells;i++)
    {
        fout <<myid+1<<'\n';
    }
    int ba = 1 - BASE;
    if (d_dim == 2)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 3:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << "\n";
                }
                break;
            case 4:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    else if (d_dim == 3)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 4://Tetrahedron [1233,4444]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            case 5://[1234 5555]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << "\n";
                }
                break;
            case 6://[1233 4566]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 5] + ba << "\n";
                }
                break;
            case 8:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 6] + ba << " " << vconn[n][nvert * i + 7] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    fout.close();
    return;
}
