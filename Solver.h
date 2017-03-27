#ifndef SOLVER_H_
#define SOLVER_H_


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <memory.h>


#include "Quadrature.h"


class Solver {

private:

	/* Geometry */
	int _mesh;
	double _width;
	double _delta;

	/* Material */	
	double _sigma_t;
	double _sigma_s;
	double _nu_sigma_f;
	
	/* Quadrature */
	int _n; /* Number of directions */
	Quadrature* _quadrature;

	/* Scalar flux */
	double* _scalar_flux;
	double* _old_scalar_flux;
	
	/* Angular flux */
	double* _angular_flux;
	
	/* Boundary flux */
	double* _boundary_flux;

	/* Source */
	double* _tot_source;
	double* _fiss_source;

	double _keff;
	
	/* Iteration */
	int _max_iter;
	double _converge_thresh;
	double _res_keff;
	double _res_flux;
	
public:

	Solver(int mesh=100, int n=2);
	virtual ~Solver();

	void setQuadrature(Quadrature* quad);

	void startIteration();

	void InitializeFlux();
	void InitializeSource();

	void computeTotalSource();
	void computeFissionSource();

	void sweepMuMinus(int n);
	void sweepMuPlus(int n);

	void storeFlux();
	void zeroFlux();

	void normalize();

	void computeKeff();

	void computeResidual();

	double getResidual() { return _res_flux; }

};

Solver::Solver(int mesh, int n) {

	_mesh = mesh;
	_n = n;

	_width = 66.0053;
	_delta = _width/2.0/_mesh;
	_sigma_t = 0.050;
	_sigma_s = 0.030;
	_nu_sigma_f = 0.0225;
	_keff =1.0;
	_max_iter = 1000;

	_converge_thresh = 0.00001;
	_res_keff = 1.0;
	_res_flux = 1.0;
	_quadrature = NULL;

	/* Scalar flux */
	_scalar_flux = NULL;
	_old_scalar_flux = NULL;

	/* Angular flux */
	_angular_flux = NULL;

	/* Boundary flux */
	_boundary_flux = NULL;

	/* Source */
	_tot_source = NULL;
	_fiss_source = NULL;
	
}

Solver::~Solver() {
	if (_quadrature != NULL)
		delete _quadrature;

	if (_scalar_flux != NULL)
		delete [] _scalar_flux;

	if (_old_scalar_flux != NULL)
		delete [] _old_scalar_flux;

	if (_angular_flux != NULL)
		delete [] _angular_flux;

	if (_boundary_flux != NULL)
		delete [] _boundary_flux;

	if (_tot_source != NULL)
		delete [] _tot_source;

	if (_fiss_source != NULL)
		delete [] _fiss_source;
}

void Solver::setQuadrature(Quadrature* quad) {
	_quadrature = quad;
}

void Solver::startIteration() {
	int num = 0;
	int i = 0;
	_keff = 1.0;
	InitializeFlux();
	InitializeSource();
//	computeFissionSource();
//	normalize();
	storeFlux();
	0();
	computeTotalSource();
// 	while (_res_keff > _converge_thresh) {
// 		zeroFlux();
// //		for (i=0; i<_max_iter; i++) {
// 		while (i<_max_iter && _res_flux > _converge_thresh) {

// 			for (int j=0; j<_n/2; j++)
// 				sweepMuMinus(j);
// 			for (int j=_n/2-1; j>=0; j--)
// 				sweepMuPlus(j);
// //			normalize();
// 			computeResidual();
// 			computeTotalSource();

// //			if (_res_flux < _converge_thresh)
// //				break;
// 			if (_res_flux > _converge_thresh)
// 				storeFlux();
			
// 			i++;
// 		}
// //		computeKeff();
// //		normalize();
// //		computeTotalSource();
// //		storeFlux();
// //		num++;
// //		std::cout << "# " << num << '\t' << "Keff = " << _keff << '\t' << "Inner = " << i << std::endl;
// //		i = 0;
// //	}
	for (int i=0; i<_max_iter; i++) {

		normalize();
		computeTotalSource();
		for (int j=0; j<_n/2; j++)
			sweepMuMinus(j);
		for (int j=_n/2-1; j>=0; j--)
			sweepMuPlus(j);	
		computeKeff();
		computeResidual();
		storeFlux();
		num++;
		if (i>1 && _res_flux < _converge_thresh)
			break;
		std::cout << "# " << num << '\t' << "Keff = " << _keff << std::endl;
	}

	std::cout << "Scalar Flux: " << std::endl;
	for (int i=0; i<_mesh; i++) {
		std::cout << "Mesh #" << i+1 << " = " << _scalar_flux[i] << '\t';
		if ((i+1)%4 ==0)
			std::cout << '\n';
	}
}

void Solver::InitializeFlux() {
	std::cout << "Flux array created..." << std::endl;
	_scalar_flux = new double [_mesh];
	_old_scalar_flux = new double [_mesh];
	_angular_flux = new double [_mesh+1];
	_boundary_flux = new double [_n/2];
	std::cout << "Flux array initialized..." << std::endl;
	memset(_scalar_flux, 1.0, sizeof(double)*_mesh);
	memset(_old_scalar_flux, 0.0, sizeof(double)*_mesh);
	memset(_angular_flux, 0.0, sizeof(double)*(_mesh+1));
	memset(_boundary_flux, 0.0, sizeof(double)*(_n/2));
}

void Solver::InitializeSource() {
	std::cout << "Source array created..." << std::endl;
	_tot_source = new double [_mesh];
	_fiss_source = new double [_mesh];
}

void Solver::computeTotalSource() {

//	std::cout << "computing total source..." << std::endl;
	for (int i=0; i<_mesh; i++) {
		_tot_source[i] = 0.5*(_sigma_s + _nu_sigma_f/_keff)*_scalar_flux[i];
	}
}

void Solver::computeFissionSource() {

//	std::cout << "computing fission source..." << std::endl;
	for (int i=0; i<_mesh; i++) {
		_fiss_source[i] = _nu_sigma_f*_scalar_flux[i];
	}
}

void Solver::sweepMuMinus(int n) {
	double L = _sigma_t*_delta/_quadrature->getMu(n)*(-1.0);
	_angular_flux[_mesh] = 0.0;
	for (int i=_mesh-1; i>=0; i--) {
		_angular_flux[i] = (2.0 + L)/(2.0 - L)*_angular_flux[i+1] - 2.0*L/(2.0 - L)*_tot_source[i]/_sigma_t;
		if (_angular_flux[i] < 0)
			_angular_flux[i] = 0.0;
	}
	for (int i=0; i<_mesh; i++) {
		_scalar_flux[i] += (_angular_flux[i] + _angular_flux[i+1])/2*_quadrature->getWeight(n);
	}
	_boundary_flux[n] = _angular_flux[0];
}
	
void Solver::sweepMuPlus(int n) {
	double L = _sigma_t*_delta/_quadrature->getMu(n);
	_angular_flux[0] = _boundary_flux[n];
	for (int i=0; i<_mesh; i++) {
		_angular_flux[i+1] = (2.0 - L)/(2.0 + L)*_angular_flux[i] + 2.0*L/(2.0 + L)*_tot_source[i]/_sigma_t;
		if (_angular_flux[i+1] < 0)
			_angular_flux[i+1] = 0.0;
	}
	for (int i=0; i<_mesh; i++) {
		_scalar_flux[i] += (_angular_flux[i] + _angular_flux[i+1])/2.0*_quadrature->getWeight(n);
	}
}

void Solver::storeFlux() {
	for (int i=0; i<_mesh; i++)
		_old_scalar_flux[i] = _scalar_flux[i];
}

void Solver::zeroFlux() {

//	std::cout << "reset angular and scalar flux..." << std::endl;
	for (int i=0; i<_mesh; i++)
		_scalar_flux[i] = 0.0;
	for (int i=0; i<=_mesh; i++)
		_angular_flux[i] = 0.0;
}

void Solver::normalize() {
//	computeFissionSource();
//	std::cout << "normalizing fission source..." << std::endl;
	double tot_fiss_source =0;
	for (int i=0; i<_mesh; i++)
		tot_fiss_source += _fiss_source[i];
	for (int i=0; i<_mesh; i++) {
		_scalar_flux[i] /= tot_fiss_source;
		_old_scalar_flux[i] /= tot_fiss_source;
	}
	for (int i=0; i<_n/2; i++) {
		_boundary_flux[i] /= tot_fiss_source;
	}
}

void Solver::computeKeff() {
//	std::cout << "computing keff..." << std::endl;
	double old_keff = _keff;
	double old_fiss_source = 0.0;
	double new_fiss_source = 0.0;
	for (int i=0; i<_mesh; i++)
		old_fiss_source += _fiss_source[i];
	computeFissionSource();
	for (int i=0; i<_mesh; i++)
		new_fiss_source += _fiss_source[i];
	_keff = _keff*new_fiss_source/old_fiss_source;
	_res_keff = fabs((_keff-old_keff)/_keff);
}
	
void Solver::computeResidual() {

//	std::cout << "computing residual..." << std::endl;
//	double old_flux = 0.0;
//	double new_flux = 0.0;
	double max_scalarflux = 0.0;
	double mm;

	for (int i=0; i<_mesh; i++) {

		mm = fabs((_scalar_flux[i]-_old_scalar_flux[i])/_scalar_flux[i]);
		if (mm > max_scalarflux)
			max_scalarflux = mm;
//		old_flux += _old_scalar_flux[i]*_old_scalar_flux[i];
//		new_flux += _scalar_flux[i]*_scalar_flux[i];
	}
	_res_flux = max_scalarflux;//fabs(sqrt(new_flux)-sqrt(old_flux))/sqrt(new_flux);
}

#endif
