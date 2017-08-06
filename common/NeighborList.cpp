///
/// Neighbor List Implementation
/// 
/// Molecular Dynamics Simulation on GPU
///
/// Written by Vadim Kuras. 2009-2010.
///

#include "NeighborList.h"

///
/// Neighbor list generation function
///
int * buildNeighborList(real4 * posArray, listSettings * listsettings, real boxSize, int buckets, int NumberOfParticles)
{
	int * counter;
	int3 * listbuckets;
	int * list;

	//calculate amount of cells, half the box size & the cellsize
	int cells = (int) ((boxSize / (listsettings->rs)) + 0.99f);

	if (cells < 1) // if the results was < 1 then set cell count to be 1
		cells = 1;

    real npos = boxSize / 2; //half the boxSize
    real cellSize = listsettings->rs; //cellsize is RS

	//allocate the memory for the buckets
	TC(counter = new int[buckets]);
	TC(listbuckets = new int3[buckets]);
    int largest=0;

    int cx;
    int cy;
    int cz;
    int tmp;
   
    int * nlist;

	//reset the allocated memory
	memset(listbuckets,-1,sizeof(int3) * buckets);
	memset(counter,0,sizeof(int) * buckets);

	//build the hash and counter the amount of particles for each cell
	for(int i=0; i<NumberOfParticles; i++)
    {
		//get cell x,y,z
		cx = (int) ((posArray[i].x + npos) / cellSize);
		cy = (int) ((posArray[i].y + npos) / cellSize);
		cz = (int) ((posArray[i].z + npos) / cellSize);

		//get a bucket for the cell
		tmp = computeBucket(cx, cy, cz, listbuckets, buckets, false);

		//count particle in the cell
		counter[tmp]++;

		//get the maximum particle per cell
		if (counter[tmp]>largest)
           largest = counter[tmp];
    }

	//allocate the particle list per bucket
    TC(list = new int[largest*buckets]);

	//reset the allocated memory
	memset(list,-1,sizeof(int) * largest*buckets);
	memset(counter,0,sizeof(int) * buckets);

	//collect the particles to the buckets
    for(int i=0; i<NumberOfParticles; i++)
    {
		//get cell x,y,z
		cx = (int) ((posArray[i].x + npos) / cellSize);
		cy = (int) ((posArray[i].y + npos) / cellSize);
		cz = (int) ((posArray[i].z + npos) / cellSize);

		//get cell bucket
		tmp = computeBucket(cx, cy, cz, listbuckets, buckets, true);

		//add the particle to the corresponding cell bucket
		list[tmp*largest+counter[tmp]] = i;
		counter[tmp]++;
    }

	//calculate neighbor list per particle size
	listsettings->nlistsize = largest*27;
	if (listsettings->nlistsize>NumberOfParticles-1) 
		listsettings->nlistsize = NumberOfParticles-1;

	//allocate neighborlist memory
    TC(nlist = new int[listsettings->nlistsize*NumberOfParticles]);

    int count;
	int largestlist=0;

	//collect particle neighbours for each particle
    for(int i=0; i<NumberOfParticles; i++)
    {
		//set count to 0
		count = 0;

		//get cell x,y,z
		cx = (int) ((posArray[i].x + npos) / cellSize);
		cy = (int) ((posArray[i].y + npos) / cellSize);
		cz = (int) ((posArray[i].z + npos) / cellSize);

		//all the cells surrounding the particle's cell
		for(int icx=-1; icx<2; icx++)
			for(int icy=-1; icy<2; icy++)
				for(int icz=-1; icz<2; icz++)
				{
					//if the cell is outside the box
					if ((cx+icx)>=cells || (cy+icy)>=cells || (cz+icz)>=cells || (cx+icx)<0 || (cy+icy)<0 || (cz+icz)<0)
						continue;

					//find the bucket
					tmp = computeBucket(cx+icx,cy+icy,cz+icz,listbuckets, buckets, true);

					//no bucket - means cell is empty
					if (tmp==-1)
						continue;

					//for each particle in the cell - add it to the particle's neighbour list
					for(int k=0;k<largest && list[tmp*largest+k]!=-1; k++)
					{
						//if it is the same particle
						if (i==list[tmp*largest+k])
							continue;

						//add to the list and count
						nlist[listsettings->nlistsize*i+count] = list[tmp*largest+k];
						count++;
					}
				}

		//largest list
		if (count>largestlist)
			largestlist = count;

		//if the particle's list isn't empty - finish it with a -1
		while (count<listsettings->nlistsize)
			nlist[listsettings->nlistsize*i+count++] = -1;
    }

	listsettings->nlargestsize = largestlist;
	
	//free allocated memory
	delete [] listbuckets;
	delete [] list;
	delete [] counter;

	//return the neighbour list
    return nlist;
}

///
/// Computes the bucket
///
int computeBucket(int cx, int cy, int cz, int3 * bucketlist, int buckets, bool pop)
{
	//first function
	int n = H1*cx + H2*cy + H3*cz;
	n = n%buckets;
	if (n<0) n+=buckets;

	if (pop && bucketlist[n].x==-1)
		return -1;

	if (bucketlist[n].x == -1 || (bucketlist[n].x == cx && bucketlist[n].y == cy && bucketlist[n].z == cz))
	{
		if (!pop)
		{
			bucketlist[n].x = cx;
			bucketlist[n].y = cy;
			bucketlist[n].z = cz;
		}

		return n;	
	}

	//if there is a collision - try second function
	n = H1*cx ^ H2*cy ^ H3*cz;
	n = n%buckets;
	if (n<0) n+=buckets;

	if (pop && bucketlist[n].x==-1)
		return -1;
	
	//if there is a collision again - stupidly jump to the next cell
	while (bucketlist[n].x != -1 && (bucketlist[n].x != cx || bucketlist[n].y != cy || bucketlist[n].z != cz))
	{
		n+=1;
		n = n%buckets;
		if (n<0) n+=buckets;
	}
		
	if (!pop)
	{
		bucketlist[n].x = cx;
		bucketlist[n].y = cy;
		bucketlist[n].z = cz;
	}
	
	return n;		
}