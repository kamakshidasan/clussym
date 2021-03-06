#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include <vector>
#include "Elements.hpp"

class Sampler
{
	public:
		Sampler(class BD* pbd, std::vector<unsigned int> & sidx, std::vector<Vertex> & vlist);
		void Sample(std::vector<float> & isovals, float alpha);
		void PropogateValues(std::vector<float> & fvals, float alpha, float min, float max);
	private:
		void PickValues(std::vector<float> & isovals, float alpha);
		void RestrictSamples(std::vector<float> & isovals);
		BD* bd;
		std::vector<unsigned int> & sadidx;
		std::vector<Vertex> & m_vlist;

};
#endif
