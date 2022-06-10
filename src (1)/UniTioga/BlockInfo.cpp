#include"BlockInfo.h"
#include<iostream>
#include<fstream>
#define BASE 1

void BlockInfo::make_iblank()
{
	ibl = new int[nnodes];
	for (int i = 0; i < nnodes; i++)
	{
		ibl[i] = 1;
	}
	for (int j = 0; j < nwbc; j++)
	{
		ibl[wbc[j] - 1] = -1;
	}
	for (int j = 0; j < nobc; j++)
	{
		ibl[obc[j] - 1] = 2;
	}
}
void BlockInfo::write_grd(const std::string& filename, int meshtag, int cur_proc)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(cur_proc) + ".dat";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        std::cout << "Open file failure\n";
        return;
    }
    
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Tioga output\"\n";
	fout << "VARIABLES=\"X\",\"Y\",";
	if (d_dim == 3)
	{
		fout << "\"Z\"\n";
	}
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
			fout << rxyz[d_dim * i + j] << " ";
		}
		fout << '\n';
	}
    int ba = 1 - BASE;
	if (d_dim == 2)
	{
		for (int n = 0; n < ntypes; n++)
		{
			int nvert = eachNodeCount[n];
			switch (nvert)
			{
			case 3:
				for (int i = 0; i < eachCellCount[n]; i++)
				{
					fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << "\n";				
				}
				break;
			case 4:
				for (int i = 0; i < eachCellCount[n]; i++)
				{
					fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
				}
				break;
			default:
				break;
			}
		}
	}
	if (d_dim == 3)
	{
		for (int n = 0; n < ntypes; n++)
		{
			int nvert = eachNodeCount[n];
			switch (nvert)
			{
			case 4:
				break;
			case 5:
				break;
			case 6:
				break;
			case 8:
			{
				for(int i=0;i<eachCellCount[n];i++)
				{
					for (int j = 0; j < nvert; j++)
					{
						fout << vconn[n][nvert * i + j] + ba << " ";
					}
					fout << '\n';
					

				}
			}
				break;
			}
		}
	}
    fout.close();
    return;

}

void BlockInfo::write_scatter_wbc(const std::string& filename, int meshtag, int cur_proc)
{
	std::string localFilename = filename + std::to_string(meshtag) + std::to_string(cur_proc) + ".dat";
	std::ofstream fout;

	fout.open(localFilename);
	if (!fout.is_open())
	{
		std::cout << "Open file failure\n";
		return;
	}
   // std::cout << "nwbc " << nwbc << '\n';
	for (int i = 0; i < nwbc; i++)
	{
		int inode = wbc[i] - 1;
		//std::cout << inode << " " << nnodes << '\n';
		for (int j = 0; j < d_dim; j++)
		{
			fout << rxyz[d_dim * inode + j] << '\t';
		}
		fout << '\n';
	}
	fout.close();
}
void BlockInfo::write_scatter_obc(const std::string& filename, int meshtag, int cur_proc)
{
	std::string localFilename = filename + std::to_string(meshtag) + std::to_string(cur_proc) + ".dat";
	std::ofstream fout;

	fout.open(localFilename);
	if (!fout.is_open())
	{
		std::cout << "Open file failure\n";
		return;
	}
    //std::cout << "nobc " << nobc << '\n';
	for (int i = 0; i < nobc; i++)
	{
		int inode = obc[i] - 1;
		//std::cout << inode << " " << nnodes << '\n';
		for (int j = 0; j < d_dim; j++)
		{
			fout << rxyz[d_dim * inode + j] << '\t';
		}
		fout << '\n';
	}
	fout.close();
}

void BlockInfo::write_scatter_gwbc(const std::string& filename, int meshtag, int cur_proc)
{
	std::string localFilename = filename + std::to_string(meshtag) + std::to_string(cur_proc) + ".dat";
	std::ofstream fout;

	fout.open(localFilename);
	if (!fout.is_open())
	{
		std::cout << "Open file failure\n";
		return;
	}
	//std::cout << "nwbc " << nwbc << '\n';
	for (int i = 0; i < gNwbc; i++)
	{
		for (int j = 0; j < d_dim; j++)
		{
			fout << gWxyz[d_dim * i + j] << '\t';
		}
		fout << '\n';
	}
	fout.close();
}
void BlockInfo::write_scatter_gobc(const std::string& filename, int meshtag, int cur_proc)
{
	std::string localFilename = filename + std::to_string(meshtag) + std::to_string(cur_proc) + ".dat";
	std::ofstream fout;

	fout.open(localFilename);
	if (!fout.is_open())
	{
		std::cout << "Open file failure\n";
		return;
	}
	//std::cout << "nobc " << nobc << '\n';
	for (int i = 0; i < gNobc; i++)
	{
		for (int j = 0; j < d_dim; j++)
		{
			fout << gOxyz[d_dim * i + j] << '\t';
		}
		fout << '\n';
	}
	fout.close();
}

void BlockInfo::write_info(const std::string& filename, int meshtag, int cur_proc)
{
	std::string localFilename = filename + std::to_string(meshtag) + std::to_string(cur_proc)+"info" + ".dat";
	std::ofstream fout;

	fout.open(localFilename);
	if (!fout.is_open())
	{
		std::cout << "Open file failure\n";
		return;
	}
	fout << ncellDistributed<<'\n';
	for (int i = 0; i < ncellDistributed; i++)
	{
		fout << cellDistributed[i] << '\n';
	}
	fout.close();
}