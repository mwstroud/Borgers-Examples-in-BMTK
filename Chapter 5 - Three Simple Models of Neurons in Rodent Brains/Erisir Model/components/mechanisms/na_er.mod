COMMENT
Erisir Model of an Inhibitory 
Neuron in Rat Somatosensory Cortex

Sodium Channel

Reference: Borgers - An Introduction to Modeling Neuronal Dynamics Chapter 5
.mod by Matthew Stroud
ENDCOMMENT

NEURON {
	SUFFIX na_er
	USEION na READ ena WRITE ina
	RANGE gbar, g
	RANGE minf, hinf, mtau, htau
	RANGE ina
    RANGE minfvhalf,minfk,hinfvhalf,hinfk
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gbar (siemens/cm2)
    minfvhalf = 25.58
    minfk = -11.4
    hinfvhalf = 58.25
    hinfk = 6.46
}

ASSIGNED {
	v (mV)
	ena (mV)
	ina (mA/cm2)
	minf
	hinf
	mtau (ms)
	htau (ms)
	g (siemens/cm2)
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	m = minf : See Borgers Page 32 Figure 5.1 for explaination
	g = gbar*m*m*m*h
	ina = g*(v-ena)
}

INITIAL {
	rate(v)
	m = minf
	h = hinf
}

DERIVATIVE states {
	rate(v)
	:m' = (minf-m)/mtau
	h' = (hinf-h)/htau
}

COMMENT

alpha_m=40*(75.5-v)./(exp((75.5-v)/13.5)-1)
beta_m=1.2262./exp(v/42.248)

alpha_h = 0.0035/exp(v/24.186)
beta_h = -0.017*(v+51.25)/(exp(-(v+51.25)/5.2)-1)

Regression fit INF
minf = 1.0/(1.0+(exp((v+25.58)/(-11.4))))
hinf = 1.0/(1.0+(exp((v+58.25)/(6.46))))

Calculated TAU
mtau = (exp(0.0740740740740741*v) - 268.430649673547)*exp(0.0236697595152433*v)/((40*v - 3020.0)*exp(0.0977438335893174*v) + 1.2262*exp(0.0740740740740741*v) - 329.149662629703)
htau = (exp(0.192307692307692*v) - 5.24437584203063e-5)*exp(0.0413462333581411*v)/((0.017*v + 0.87125)*exp(0.233653925665833*v) + 0.0035*exp(0.192307692307692*v) - 1.83553154471072e-7)

ENDCOMMENT

PROCEDURE rate(v (mV)) {
	UNITSOFF
	minf = 1.0/(1.0+(exp((v+minfvhalf)/(minfk))))
	:minf = (40*v - 3020.0)*exp(0.0977438335893174*v)/((40*v - 3020.0)*exp(0.0977438335893174*v) + 1.2262*exp(0.0740740740740741*v) - 329.149662629703)
	mtau = (exp(2/27*v) - 268.430)*exp(125/5281*v)/((40*v - 3020.0)*exp(0.0977*v) + 1.2262*exp(2/27*v) - 329.14966)  
	hinf = 1.0/(1.0+(exp((v+hinfvhalf)/(hinfk))))  
	:hinf = (0.0035*exp(0.192307692307692*v) - 1.83553154471072e-7)/((0.017*v + 0.87125)*exp(0.233653925665833*v) + 0.0035*exp(0.192307692307692*v) - 1.83553154471072e-7)
	htau = (exp(5/26*v) - 5.2443e-5)*exp(500/12093*v)/((0.017*v + 697/800)*exp(0.23365*v) + 0.0035*exp(5/26*v) - 1.8355e-7)
	UNITSON
}