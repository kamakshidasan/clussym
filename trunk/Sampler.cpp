#include "Sampler.hpp"
#include "BD.hpp"
#include <iostream>
extern unsigned int fsz;

Sampler::Sampler(BD* pbd, std::vector<unsigned int> & sidx, std::vector<Vertex> & vlist)
	: bd(pbd), sadidx(sidx), m_vlist(vlist)
{}
void Sampler::PickValues(std::vector<float> & isovals, float orgalpha)
{
	std::sort(sadidx.begin(), sadidx.end(), FnCmp(&m_vlist));
	unsigned int sz = sadidx.size();
	float maxf = m_vlist[bd->bridsarr[1]->ext].w;
	float minf = m_vlist[bd->bridsarr[1]->sad].w;
	float curf = minf, nextf;
	float alpha = orgalpha*(maxf - minf);
	std::cout<<"orgalpha range alpha "<<orgalpha<<" "<<maxf-minf<<" "<<alpha<<std::endl;
	std::cout<<"saddle list sz fsz"<<sadidx.size()<<" "<<fsz<<std::endl;
	for(unsigned int i = 0; i < sz; i++)
	{
//		std::cout<<"Saddle feature sz value: "<<sadidx[i]<<" "<<m_vlist[sadidx[i]].feature<<" "<<m_vlist[sadidx[i]].w<<std::endl;
		nextf = m_vlist[sadidx[i]].w;
//		std::cout<<"nextf curf"<<nextf<<" "<<curf<<std::endl;
		if(nextf - curf > alpha)
		{
			float f = nextf - alpha;
			isovals.push_back(f);
			std::cout<<"Iso "<<f<<std::endl;
		}
		curf = nextf;
	}
}
void Sampler::Sample(std::vector<float> & isovals, float alpha)
{
	isovals.push_back(-64);
/*	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);



	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
	isovals.push_back(-8.38);
*/
//	PickValues(isovals, alpha);	
//	RestrictSamples(isovals);
}
void Sampler::RestrictSamples(std::vector<float> & isovals)
{
	float max = m_vlist[bd->bridsarr[1]->ext].w;
	float min = m_vlist[bd->bridsarr[1]->sad].w;
	for(unsigned int j = 1; j <= bd->numbr; j++)
	{
		SymBranch* b = bd->bridsarr[j];
		std::vector<unsigned int> canidx;
		unsigned int n = isovals.size();
		if(j == 1)
		{
			for(int i = n-1; i >= 0; i--)
			{
				canidx.push_back(i);
			}
		}
		else if(m_vlist[b->ext].w > m_vlist[b->sad].w)
			continue;
		else
		{ 
		  	for(int i = n-1; i >= 0; i--)
		  	{
		  		if(m_vlist[b->sad].w > isovals[i] && m_vlist[b->ext].w < isovals[i])
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
					unsigned int pcidx = canidx[i-1];
					if((m_vlist[(*bit)->sad].w >= isovals[cidx] && m_vlist[(*bit)->sad].w < isovals[pcidx]))
					{
						if(selarr[canidx[i]] < fsz && selarr[canidx[i]]+(*bit)->csz >= fsz)
						{
							selcnt++;
							b->comps[cidx] = -1;
						}
						selarr[canidx[i]] += (*bit)->csz;
						break;	
					}
				}
			}
		}
	}
}
void Sampler::PropogateValues(std::vector<float> & isovals, float alpha, float min, float max)
{
	PickValues(isovals, alpha);	
	float curmax = m_vlist[bd->bridsarr[1]->ext].w;
	float curmin = m_vlist[bd->bridsarr[1]->sad].w;
//	isovals.push_back(-95.2415365);
/*	for(unsigned int i = 0; i < fvals.size(); i++)
	{
		float f = (fvals[i] - min)*(curmax-curmin)/(max - min) + curmin;
		isovals.push_back(f);
	}*/
//	RestrictSamples(isovals);
}
