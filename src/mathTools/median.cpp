#include<math.h>
#include<iostream>
#define BASE 1

void find_median(int* ix, double* x, int n, int& ilow, int& ihigh)
{
    double xtemp;
    int itemp;
    int nleft = ihigh - ilow;
    //Now the nleft+1 >= 3 is the number of elements;
    int low = ilow;
    int high = ihigh;
    //when n is odd, that gives the middle;
    //when n is even that gives the left of the middle;
    int mid = (low + high) / 2;
    int nby2 = (n - 1) / 2;//When n is odd ,that gives
    bool is_odd = n % 2 == 1;
    int nby2p1 = (n - 1) / 2 + 1;
    //make sure low < med < high
    double xlow = x[low];
    double xhigh = x[high];
    double xmed = x[mid];

    if (xlow > xhigh)
    {
        xtemp = xlow;
        xlow = xhigh;
        xhigh = xtemp;
    }
    if (xmed > xhigh)
    {
        xmed = xhigh;
    }
    else if (xmed < xlow)
    {
        xmed = xlow;
    }
    //move all values less than xmed to the left

    int i = low;
    int j = high;
    bool signal = true;
    while (signal)
    {
        signal = false;
        while (true)
        {
            if (x[i] >= xmed) break;
            ++i;
        }
        while (true)
        {
            if (x[j] <= xmed) break;
            --j;
        }
        if (i < j)//there is a pair of misfits
        {
            xtemp = x[i];
            x[i] = x[j];
            x[j] = xtemp;
            itemp = ix[i];
            ix[i] = ix[j];
            ix[j] = itemp;
            //move to the next
            i = i + 1;
            j = j - 1;
            if (i <= j)
            {
                signal = true;
            }
        }
    }

    //now we got all values less than xmed at left side ,others at right side.It's time for discussion
    //Two situations in general: 
    // i==j happens when i and j exchange around their center
    // and i,j neibouring each other at other situations
    if (!is_odd)
    {
        if (j == nby2 and i == nby2p1)
        {
            //special case where it happens to be equally divided into two
            //return;
            int idx = j;
            for (int k = 0; k < nby2; k++)
            {
                if (x[k] > x[idx])
                {
                    idx = k;
                }
            }
            if (idx != j)
            {
                xtemp = x[idx];
                x[idx] = x[j];
                x[j] = xtemp;
                itemp = ix[idx];
                ix[idx] = ix[j];
                ix[j] = itemp;
            }
            return;
        }
        if (j < nby2)
        {
            low = i;
        }
        if (i > nby2p1)
        {
            high = j;
        }
        if (i != j)
        {
            if (low < high - 1)
            {
                find_median(ix, x, n, low, high);
            }
            /*else
            {
                if (x[n / 2 - 1] > x[n / 2])
                {
                    double xtemp = x[n / 2 - 1];
                    x[n / 2 - 1] = x[n / 2];
                    x[n / 2] = xtemp;
                    int itemp = ix[n / 2 - 1];
                    ix[n / 2 - 1] = ix[n / 2];
                    ix[n / 2] = itemp;
                }

            }*/
            //ilow = low;
            //ihigh = high;
            //Call myself again
            return;
        }
        if (i == nby2)
        {
            low = nby2;
        }
        if (j == nby2p1)
        {
            high = nby2p1;
        }
    }
    else
    {
        if (j < nby2) low = i;
        if (i > nby2) high = j;
        if (i != j)
        {
            if (low < high - 1)
            {
                find_median(ix, x, n, low, high);
            }
            //else if (high - low == 1)
            //{
            //    if (x[low] > x[high])
            //    {
            //        xtemp = x[low];
            //        x[low] = x[high];
            //        x[high] = xtemp;
            //        itemp = ix[low];
            //        ix[low] = ix[high];
            //        ix[high] = itemp;
            //    }
            //}
            ilow = low;
            ihigh = high;
            return;

        }
        if (i == nby2)
        {
            ilow = low;
            ihigh = high;
            return;
        }
    }
    if (low < high - 1)
    {
        ilow = low;
        ihigh = high;
        find_median(ix, x, n, low, high);
    }
    return;

}

void find_median_wrap(int* ix, double* x, int n)
{
    if (n <= 2)
    {
        if (n == 1)//Ilow and Ihigh are at the same point
        {
        }
        else if (n == 2)//Ilow and Ihigh are neighbours
        {
            if (x[0] > x[1])
            {
                double xtemp = x[0];
                x[0] = x[1];
                x[1] = xtemp;
                int itemp = ix[0];
                ix[0] = ix[1];
                ix[1] = itemp;
            }
        }
        return;
    }
    int ilow = 0;
    int ihigh = n - 1;
    int nby2 = (n - 1) / 2;
    int nby2p1 = nby2 + 1;

    find_median(ix, x, n, ilow, ihigh);
    if (n % 2 == 0)
    {
        if (x[nby2] > x[nby2p1])
        {
            double xtemp = x[nby2];
            x[nby2] = x[nby2p1];
            x[nby2p1] = xtemp;
            int itemp = ix[nby2];
            ix[nby2] = ix[nby2p1];
            ix[nby2p1] = itemp;
        }
    }
    else
    {
        if (x[ilow] > x[ihigh])
        {
            double xtemp = x[ilow];
            x[ilow] = x[ihigh];
            x[ihigh] = xtemp;
            int itemp = ix[ilow];
            ix[ilow] = ix[ihigh];
            ix[ihigh] = itemp;
        }
    }

    //if (n % 2 == 0)
    //{
    //    if (x[nby2] > x[nby2p1])
    //    {
    //        //exchange
    //    }
    //}
    //else
    //{

    //}
    return;
}