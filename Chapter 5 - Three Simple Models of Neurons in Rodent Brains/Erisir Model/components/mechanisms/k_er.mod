COMMENT
Erisir Model of an Inhibitory 
Neuron in Rat Somatosensory Cortex

Potassium Channel

Reference: Borgers - An Introduction to Modeling Neuronal Dynamics Chapter 5
.mod by Matthew Stroud
ENDCOMMENT

NEURON {
	SUFFIX k_er
	USEION k READ ek WRITE ik
	RANGE gbar, g
	RANGE inf, tau
	RANGE ik
    RANGE ninfvhalf,ninfk
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gbar (siemens/cm2)
    ninfvhalf =1.57
    ninfk = -8.38
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)
	g (siemens/cm2)
	inf
	tau (ms)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*n*n*n*n
	ik = g*(v-ek)
}

INITIAL {
	rate(v)
	n = inf
}

DERIVATIVE states {
	rate(v)
	n' = (inf-n)/tau
}

COMMENT

alpha_n=(95-v)/(exp((95-v)/11.8)-1)
beta_n=0.025/exp(v/22.222)
n_e_inf=alpha_n./(alpha_n+beta_n)


Regression fit INF
ninf = 1.0/(1.0+(exp((v+1.57)/(-8.38))))

Calculated TAU
ntau = (exp(0.0847457627118644*v) - 3136.4518553778)*exp(0.0450004500045*v)/((v - 95)*exp(0.129746212716364*v) + 0.025*exp(0.0847457627118644*v) - 78.4112963844449)

ENDCOMMENT

PROCEDURE rate(v (mV)) {
	UNITSOFF
	inf = 1.0/(1.0+(exp((v+ninfvhalf)/(ninfk))))   
	:inf = (v - 95)*exp(0.129746212716364*v)/((v - 95)*exp(0.129746212716364*v) + 0.025*exp(0.0847457627118644*v) - 78.4112963844449)    
	tau = (exp(5/59*v) - 3136.4518)*exp(500/11111*v)/((v - 95)*exp(0.129746*v) + 0.025*exp(5/59*v) - 78.411296)
	UNITSON
}