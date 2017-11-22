#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <Python.h>

#define CALL (0)
#define PUT  (1)
#define MAX_IV_LOOP (30)

double _bs_d1(double S, double K, double r, double q, double t, double sigma) {
	//(log(_S/_X) + (_r - _b + (_v*_v)/2.0) * _t)/(_v * sqrt(_t));
	return (log(S / K) + (r - q + (sigma * sigma) / 2) * t) / (sigma * sqrt(t));
}

double _bs_d2(double d1, double t, double sigma) {
	//d1 - _v * (sqrt(_t));     
	return d1 - sigma * sqrt(t);
}


double norm_dist(double d) {
	return exp(-(d * d) / 2) / sqrt(6.2831853071795862);
}

double norm_sdist(double d) {
	double d1 = norm_dist(d);
	double d2 = 1.0 / (1.0 + 0.23164190000000001 * fabs(d));
	double d3 = d1
			* (0.31938153000000002 * d2 + -0.356563782 * pow(d2, 2)
					+ 1.781477937 * pow(d2, 3)
					+ -1.8212559779999999 * pow(d2, 4) + 1.3302744289999999 * pow(d2, 5));
	if (d >= 0.0)
		d3 = 1.0 - d3;
	return d3;
}



double _bs_prem_put(double S, double K, double r, double q, double t, double sigma) {
	//exp(-_r * _t) * _X * pnorm(-d2,0,1,1,0) - exp(-_b * _t) * _S * pnorm(-d1,0,1,1,0); 
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return exp(-r * t) * K * norm_sdist(-d2) - exp(-q * t) * S * norm_sdist(-d1);
}

double _bs_prem_call(double S, double K, double r, double q, double t, double sigma) {
	//exp(-_b * _t) * _S * pnorm(d1,0,1,1,0) - exp(-_r * _t) * _K * pnorm(d2,0,1,1,0); 
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return exp(-q * t) * S * norm_sdist(d1) - exp(-r * t) * K * norm_sdist(d2);
	return 0;
}


double _bs_delta_put(double S, double K, double r, double q, double t, double sigma) {
	//-exp(-_b * _t) * pnorm(-d1,0,1,1,0);
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	return -exp(-r * t) * norm_sdist(-d1);
}

double _bs_delta_call(double S, double K, double r, double q, double t, double sigma) {
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	return exp(-r * t) * norm_sdist(d1);
}


double _bs_vega(double S, double K, double r, double q, double t, double sigma) {
	//_S * exp(-_b * _t) * dnorm(d1,0,1,0) * sqrt(_t);
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	return S * exp(-q * t) * norm_dist(d1) * sqrt(t);
}


double _bs_theta_put(double S, double K, double r, double q, double t, double sigma) {
	//-exp(-_b * _t) * (_S * dnorm(d1,0,1,0) * _v)/ (2.0 * sqrt(_t)) + 
	// (_r * _X * exp(-_r * _t) * pnorm(-d2,0,1,1,0));
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return -exp(-q * t) * (S * norm_dist(d1) * sigma)/ (2.0 * sqrt(t)) +
			(r * K * exp(-r * t) * norm_sdist(-d2))
			- q * S * exp(-q * t) * norm_sdist(-d1);
}

double _bs_theta_call(double S, double K, double r, double q, double t, double sigma) {
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return -exp(-q * t) * (S * norm_dist(d1) * sigma)/ (2.0 * sqrt(t)) -
			(r * K * exp(-r * t) * norm_sdist(d2))
			+ q * S * exp(-q * t) * norm_sdist(d1);
}


double _bs_rho_put(double S, double K, double r, double q, double t, double sigma) {
	//-_X * _t * exp(-_r * _t) * pnorm(-d2,0,1,1,0);    
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return -K * t * exp(-r * t) * norm_sdist(-d2);
}

double _bs_rho_call(double S, double K, double r, double q, double t, double sigma) {
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return K * t * exp(-r * t) * norm_sdist(d2);
}


double _bs_gamma(double S, double K, double r, double q, double t, double sigma) {
	//exp(-_b * _t) * ( dnorm(d1,0,1,0)/( _S * _v * sqrt(_t)) );
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	return exp(-q * t) * (norm_dist(d1)/(S * sigma * sqrt(t)));
}


double _bs_vanna(double S, double K, double r, double q, double t, double sigma) {
	// -exp(-_b * _t) * dnorm(d1,0,1,0) * (d2/_v);
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return -exp(-q * t) * norm_dist(d1) * (d2/sigma);
}


double _bs_charm_put(double S, double K, double r, double q, double t, double sigma) {
	// _b * exp(-_b * _t) * pnorm(-d1,0,1,1,0) + exp(-_b * _t) * dnorm(d1,0,1,0) * 
	//	( (2.0 * (_r - _b) * _t - d2 * _v * sqrt(_t))/(2.0 * _t * _v * sqrt(_t)) );
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return -q * exp(-q * t) * norm_sdist(-d1) - exp(-q * t) * norm_dist(d1) *
		( (2.0 * (r - q) * t - d2 * sigma * sqrt(t))/(2.0 * t * sigma * sqrt(t)) );
}

double _bs_charm_call(double S, double K, double r, double q, double t, double sigma) {
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return q * exp(-q * t) * norm_sdist(d1) - exp(-q * t) * norm_dist(d1) *
		( (2.0 * (r - q) * t - d2 * sigma * sqrt(t))/(2.0 * t * sigma * sqrt(t)) );
}


double _bs_speed(double S, double K, double r, double q, double t, double sigma) {
	// -(gamma/_S) * (d1/(_v * sqrt(_t)) + 1.0);
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double g = _bs_gamma(S, K, r, q, t, sigma);
	return -(g / S) * (d1 / (sigma * sqrt(t)) + 1.0);
}


double _bs_zomma(double S, double K, double r, double q, double t, double sigma) {
	// gamma * ( (d1*d2 - 1.0)/_v );    
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	double g = _bs_gamma(S, K, r, q, t, sigma);
	return g * ( (d1*d2 - 1.0) / sigma );
}


double _bs_color(double S, double K, double r, double q, double t, double sigma) {
	// -exp(-_b * _t) * (dnorm(d1,0,1,0)/(2.0*_S*_t*_v*sqrt(_t))) *
	// (2.0 * _b * _t + 1.0 + ((2.0 * (_r-_b) * _t - d2 * _v * sqrt(_t))/ (_v * sqrt(_t))) * d1);
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return -exp(-q * t) * (norm_dist(d1)/(2.0*S*t*sigma*sqrt(t))) *
		(2.0 * q * t + 1.0 + ((2.0 * (r-q) * t - d2 * sigma * sqrt(t))/ (sigma * sqrt(t))) * d1);
}


double _bs_DvegaDtime(double S, double K, double r, double q, double t, double sigma) {
	//_S * exp(-_b * _t) * dnorm(d1,0,1,0) * sqrt(_t) * (_b + ((_r-_b)*d1)/(_v*sqrt(_t))-((1.0+d1*d2)/(2.0*_t))); 
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return S * exp(-q * t) * norm_dist(d1) * sqrt(t) * (q + ((r-q)*d1)/(sigma*sqrt(t))-((1.0+d1*d2)/(2.0*t))); 
}


double _bs_vomma(double S, double K, double r, double q, double t, double sigma) {
	//vega * ((d1*d2) / _v);
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	double v = _bs_vega(S, K, r, q, t, sigma);
	return v * ((d1*d2) / sigma);
}


double _bs_ultima(double S, double K, double r, double q, double t, double sigma) {
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	double v = _bs_vega(S, K, r, q, t, sigma);
	return (-v / (sigma*sigma)) * ((d1*d2)*(1.0-d1*d2) + d1*d1 + d2*d2);
}


double _bs_dualdelta_put(double S, double K, double r, double q, double t, double sigma) {
	//exp(-_r * _t) * pnorm(-d2,0,1,1,0);    
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return exp(-r * t) * norm_sdist(-d2);
}

double _bs_dualdelta_call(double S, double K, double r, double q, double t, double sigma) {
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return -exp(-r * t) * norm_sdist(d2);
}


double _bs_dualgamma(double S, double K, double r, double q, double t, double sigma) {
	//exp(-_r * _t) * ( dnorm(d2,0,1,0)/(_X*_v*sqrt(_t)) );
	double d1 = _bs_d1(S, K, r, q, t, sigma);
	double d2 = _bs_d2(d1, t, sigma);
	return exp(-r * t) * ( norm_dist(d2)/(K*sigma*sqrt(t)) );
}

double _bs_ivolat(double S, double K, double r, double q, double t, double p, int pc) {
	const double e = 0.0001;
	double sigma = sqrt(fabs(log(S/K)+r*t)*2/t);
	double bsprem = 0;
	double v;
	int cnt = 0;

	if(sigma < e)
		sigma = e;

	while (fabs(bsprem - p)> e && cnt < MAX_IV_LOOP){
		if(pc == CALL){
			bsprem = _bs_prem_call(S, K, r, q, t, sigma);
		}else{
			bsprem = _bs_prem_put(S, K, r, q, t, sigma);
		}
		v = _bs_vega(S, K, r, q, t, sigma);
		sigma = sigma - (bsprem-p) / v;
		cnt++;
	}

	if(cnt < MAX_IV_LOOP)
		return sigma;
	return 0;

}

double _bs_ivolat_call(double S, double K, double r, double q, double t, double p) {
	return _bs_ivolat(S, K, r, q, t, p, CALL);
}

double _bs_ivolat_put(double S, double K, double r, double q, double t, double p) {
	return _bs_ivolat(S, K, r, q, t, p, PUT);
}



static PyObject *
bs_prem_put(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_prem_put(S, K, r, q, t, sigma));
}

static PyObject *
bs_prem_call(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_prem_call(S, K, r, q, t, sigma));
}

static PyObject *
bs_delta_put(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_delta_put(S, K, r, q, t, sigma));
}

static PyObject *
bs_delta_call(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_delta_call(S, K, r, q, t, sigma));
}

static PyObject *
bs_vega(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_vega(S, K, r, q, t, sigma));
}

static PyObject *
bs_theta_put(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_theta_put(S, K, r, q, t, sigma));
}

static PyObject *
bs_theta_call(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_theta_call(S, K, r, q, t, sigma));
}

static PyObject *
bs_rho_put(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_rho_put(S, K, r, q, t, sigma));
}

static PyObject *
bs_rho_call(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_rho_call(S, K, r, q, t, sigma));
}

static PyObject *
bs_gamma(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_gamma(S, K, r, q, t, sigma));
}

static PyObject *
bs_vanna(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_vanna(S, K, r, q, t, sigma));
}

static PyObject *
bs_charm_put(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_charm_put(S, K, r, q, t, sigma));
}

static PyObject *
bs_charm_call(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_charm_call(S, K, r, q, t, sigma));
}

static PyObject *
bs_speed(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_speed(S, K, r, q, t, sigma));
}

static PyObject *
bs_zomma(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_zomma(S, K, r, q, t, sigma));
}

static PyObject *
bs_color(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_color(S, K, r, q, t, sigma));
}

static PyObject *
bs_DvegaDtime(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_DvegaDtime(S, K, r, q, t, sigma));
}

static PyObject *
bs_vomma(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_vomma(S, K, r, q, t, sigma));
}

static PyObject *
bs_ultima(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_ultima(S, K, r, q, t, sigma));
}

static PyObject *
bs_dualdelta_put(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_dualdelta_put(S, K, r, q, t, sigma));
}

static PyObject *
bs_dualdelta_call(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_dualdelta_call(S, K, r, q, t, sigma));
}

static PyObject *
bs_dualgamma(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &sigma))
        return NULL;
    return Py_BuildValue("d", _bs_dualgamma(S, K, r, q, t, sigma));
}

static PyObject *
bs_ivolat(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, sigma;
    int pc;
    if (!PyArg_ParseTuple(args, "ddddddi", &S, &K, &r, &q, &t, &sigma, &pc))
        return NULL;
    return Py_BuildValue("d", _bs_ivolat(S, K, r, q, t, sigma, pc));
}

static PyObject *
bs_ivolat_call(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, p;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &p))
        return NULL;
    return Py_BuildValue("d", _bs_ivolat_call(S, K, r, q, t, p));
}

static PyObject *
bs_ivolat_put(PyObject *self, PyObject *args)
{
    double S, K, r, q, t, p;
    if (!PyArg_ParseTuple(args, "dddddd", &S, &K, &r, &q, &t, &p))
        return NULL;
    return Py_BuildValue("d", _bs_ivolat_put(S, K, r, q, t, p));
}


PyDoc_STRVAR(bs_doc, "Python3 extending.\n");

static PyMethodDef methods[] = {
    {"bs_prem_put", bs_prem_put, METH_VARARGS, ".\n"},
    {"bs_prem_call", bs_prem_call, METH_VARARGS, ".\n"},
    {"bs_delta_put", bs_delta_put, METH_VARARGS, ".\n"},
    {"bs_delta_call", bs_delta_call, METH_VARARGS, ".\n"},
    {"bs_vega", bs_vega, METH_VARARGS, ".\n"},
    {"bs_theta_put", bs_theta_put, METH_VARARGS, ".\n"},
    {"bs_theta_call", bs_theta_call, METH_VARARGS, ".\n"},
    {"bs_rho_put", bs_rho_put, METH_VARARGS, ".\n"},
    {"bs_rho_call", bs_rho_call, METH_VARARGS, ".\n"},
    {"bs_gamma", bs_gamma, METH_VARARGS, ".\n"},
    {"bs_vanna", bs_vanna, METH_VARARGS, ".\n"},
    {"bs_charm_put", bs_charm_put, METH_VARARGS, ".\n"},
    {"bs_charm_call", bs_charm_call, METH_VARARGS, ".\n"},
    {"bs_speed", bs_speed, METH_VARARGS, ".\n"},
    {"bs_zomma", bs_zomma, METH_VARARGS, ".\n"},
    {"bs_color", bs_color, METH_VARARGS, ".\n"},
    {"bs_DvegaDtime", bs_DvegaDtime, METH_VARARGS, ".\n"},
    {"bs_vomma", bs_vomma, METH_VARARGS, ".\n"},
    {"bs_ultima", bs_ultima, METH_VARARGS, ".\n"},
    {"bs_dualdelta_put", bs_dualdelta_put, METH_VARARGS, ".\n"},
    {"bs_dualdelta_call", bs_dualdelta_call, METH_VARARGS, ".\n"},
    {"bs_dualgamma", bs_dualgamma, METH_VARARGS, ".\n"},
    {"bs_ivolat", bs_ivolat, METH_VARARGS, ".\n"},
    {"bs_ivolat_call", bs_ivolat_call, METH_VARARGS, ".\n"},
    {"bs_ivolat_put", bs_ivolat_put, METH_VARARGS, ".\n"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef bsmodule = {
	PyModuleDef_HEAD_INIT,
	"bs",   /* name of module */
	bs_doc, /* module documentation, may be NULL */
	-1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
	methods
};

PyInit_bs(void)
{
	return PyModule_Create(&bsmodule);
}

#if 0
int main(){
	double iv = _bs_ivolat(8654.30, 8500, (double)0.5/100, 0, (double)24/365, (double)270.1, CALL);
	printf("iv = %f\n", iv);

}
#endif
