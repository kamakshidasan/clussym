#include "Sampler.hpp"
#include "BD.hpp"
#include <iostream>
Sampler::Sampler(BD* pbd, std::vector<unsigned int> & sidx, std::vector<Vertex> & vlist)
	: bd(pbd), sadidx(sidx), m_vlist(vlist)
{}
void Sampler::PickValues(std::vector<float> & isovals)
{
	std::sort(sadidx.begin(), sadidx.end(), FnCmp(&m_vlist));
	unsigned int sz = sadidx.size();
	float minf = m_vlist[sadidx[0]].w;
	float maxf = m_vlist[sadidx[sz-1]].w;
	float curf = minf, nextf, delta = (maxf - minf)/50;
	for(unsigned int i = 1; i < sz; i++)
	{
		nextf = m_vlist[sadidx[i]].w;
		if(nextf - curf > delta)
		{
			float f = nextf - delta/2;
			isovals.push_back(f);
			std::cout<<"Isovalues: "<<f<<std::endl;
		}
		curf = nextf;
	}
}
void Sampler::Sample(std::vector<float> & isovals)
{
	PickValues(isovals);	
}
void Sampler::RestrictSamples(std::vector<float> & isovals)
{
	for(unsigned int i = 1; i < bd->bridsarr.size(); i++)
	{
		SymBranch* b = bd->bridsarr[i];
		if(i > 1 && (m_vlist[b->ext].w > m_vlist[b->sad].w))
			continue;
		std::vector<unsigned int> canidx;
		unsigned int n = isovals.size();
		for(int i = n-1; i <= 0; i++)
		{
			if(b->sad > isovals[i] && b->ext < isovals[i])
			{
				canidx.push_back(i);
			}
		}
		unsigned int csz = canidx.size();
		if(csz)
		{
			unsigned int selcnt = 1;
			std::vector<unsigned int> selarr(isovals.size(), 0);
			unsigned int cidx0 = canidx[0];
			selarr[cidx0] = 1;
			b->topcomp = cidx0;
			b->comps[cidx0] = -1;


			std::list<SymBranch*>::iterator bit = b->ch.begin();
			for(; bit != b->ch.end() && selcnt < csz; bit++)
			{
				for(unsigned int i = 1; i < csz; i++)
				{
					unsigned int cidx = canidx[i];
					if((*bit)->sad > isovals[cidx])
					{
						if(!selarr[canidx[i]])
						{
							selcnt++;
							selarr[cidx] = 1;
							b->comps[cidx] = -1;
						}
						break;	
					}
				}
			}
		}
	}
}
