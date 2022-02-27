#include "m_pd.h"
#include <stdlib.h>
#include <math.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif
#define WAVETABLESIZE 1024
#define PI 3.141592653589793
#define SAMPLERATE 44000

static t_class *wavefold_allpass_class;

typedef struct _wavefold_allpass
{
    t_object x_obj; 	    // obligatory header
    t_float x_f;    	    // frequency
	t_float *wavetable;     // holds the values for the lookup table oscillator
	t_float phase;          // phase
	t_float samplerate;     // rate of samples

    t_float gain;           // gain value for wavefold
    t_float offset;         // offset value for wavefold
    t_float toggle_wavefold;// wavefold will be on/off if toggle is 0/1

    t_float cutoff;         // cutoff value for allpass filter
    t_float bandwidth;      // bandwidth value for allpass filter
    t_float toggle_allpass; // allpass filter will be on/off if toggle is 0/1

    t_float toggle_wave;    // wave will be sine/square if toggle is 0/1
} t_wavefold_allpass;


// quad interpolation routine
static inline float quad_interpolate(t_wavefold_allpass *x)
{
	int truncphase = (int) (x->phase * WAVETABLESIZE);
	float fr = (x->phase * WAVETABLESIZE) - ((float) truncphase);
	float inm1 = x->wavetable[(truncphase - 1) % WAVETABLESIZE];
	float in   = x->wavetable[(truncphase + 0) % WAVETABLESIZE];
	float inp1 = x->wavetable[(truncphase + 1) % WAVETABLESIZE];
	float inp2 = x->wavetable[(truncphase + 2) % WAVETABLESIZE];

	return in + 0.5 * fr * (inp1 - inm1 +
	 fr * (4.0 * inp1 + 2.0 * inm1 - 5.0 * in - inp2 +
	 fr * (3.0 * (in - inp1) - inm1 + inp2)));
}

// function for wavefolding: If toggle is 0, the previous input is passed through to 
// the output. If toggle is 1, wavefolding is applied to the input. In this case, the 
// input is multiplied by the gain (inlet) which is then added to the offset (inlet)
// to make the ingain. Then, the ingain is plugged into a triangle wave (4 harmonic
// approximation) and outputted as a folded wave.
static inline float wave_folding(float input, float gain, float offset, float toggle){
    
    //handles the case where toggle is not a boolean value
    if(toggle != 0){
        toggle = 1;
    }

    //if toggle is off (= 0), don't add wavefolding
    if(!toggle){
        return input;
    }

    //Note: there aren't any illegal values for gain or offset as far as I can see
    //HOWEVER, high values for gain produce "undesirable" output (although still legal)

    //add gain and offset to input value
    float ingain = (gain * input) + offset;

    //plug input with gain and offset into triangle wave (4 harmonic approximation)
    return cos(0.5 * PI * ingain)
        - 1.0/9.0 * cos(1.5 * PI * ingain)
        + 1.0/25.0 * cos(2.5 * PI * ingain)
        - 1.0/49.0 * cos(3.5 * PI * ingain);
}

//global variables to store the value of allpass filters between sample iterations
double x_0 = 0.0; // input
double x_1 = 0.0; // delayed input
double x_2 = 0.0; // delayed input
double y_0 = 0.0; // allpass output
double y_1 = 0.0; // delayed output
double y_2 = 0.0; // delayed output

// function for allpass filter: If toggle is 0, the previous input is passed through to 
// the output. If toggle is 1, a second order allpass filter is applied to the input. 
// The apf is calculated by storing the values of the previous output and applying it 
// to an equation that puts 2 iterations of this together. Bandwidth and cutoff are
// inlets that shape the allpass filter
static inline float allpass(t_wavefold_allpass *x, float input, float toggle){

    //handles the case where toggle is not a boolean value
    if(toggle != 0){
        toggle = 1;
    }

    //if toggle is off (= 0), don't add apf
    if(!toggle){
        return input;
    }

    //rest of function handles calculations for apf

    float cutoff = x->cutoff/x->samplerate;
    float bandwidth = x->bandwidth/x->samplerate;

    // keep cutoff nonzero and between -0.5 and 0.5 to satisfy calculations later
    if(cutoff == 0){
        cutoff = 0.001;
    }else if(cutoff > 0.5){
        cutoff = 0.5;
    }else if(cutoff < -0.5){
        cutoff = -0.5;
    }

    // keep bandwidth nonzero and between -0.5 and 0.5 to satisfy calculations later
    if(bandwidth == 0){
        bandwidth = 0.001;
    }else if(bandwidth >= 0.5){
        bandwidth = 0.499;
    }else if(bandwidth <= -0.5){
        bandwidth = -0.499;
    }

    // calculate valus for d, tf, and c
    double d = -cos(2.0 * PI * (cutoff));
    double tf = tan(PI * (bandwidth)); 
    double c = (tf - 1.0)/(tf + 1.0); 

    // store values for next function call
    // also, calculate y_0 aka the output
    x_0 = input;   
    y_0 = -c * x_0 + (d - d * c) * x_1 + x_2 - (d - d * c) * y_1 + c * y_2;
    x_2 = x_1;
    x_1 = x_0;
    y_2 = y_1;
    y_1 = y_0;

    return y_0;
}

// function to fill the wave table. If x->toggle_wave is 0, sine wave is put in.
// If it is 1, square wave is put in. The function iterates through the wavetable 
// using its size as an input and fills in each value one by one.
static void fill_wave_table(t_wavefold_allpass *x, float size){
    int i;

    for(i = 0; i < WAVETABLESIZE; i++){

        if(!x->toggle_wave){
            //Sine Wave
            *(x->wavetable+i) = sinf(2 * PI * (float)i/size);
        }else{
            //Square Wave
            if(i <= WAVETABLESIZE / 2) *(x->wavetable+i) = -0.9;
            else *(x->wavetable+i) = 0.9;
        }
    }
}

/* this is the actual performance routine which acts on the samples.
    It's called with a single pointer "w" which is our location in the
    DSP call list.  We return a new "w" which will point to the next item
    after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
    (no use to us), w[1] and w[2] are the input and output vector locations,
    and w[3] is the number of points to calculate. */
int toggle = 0;

static t_int *wavefold_allpass_perform(t_int *w)
{
	t_wavefold_allpass *x = (t_wavefold_allpass *)(w[1]);
    t_float *freq = (t_float *)(w[2]);
    t_float *out = (t_float *)(w[3]);
    int n = (int)(w[4]);

	// count from 0
	int blocksize = n;
	int i, sample = 0;
	float phaseincrement;
    float wave_fold_in;
    float allpass_in;

    double in, out1 = 0.0;

    // iterates by sample to create the wave
    while (n--)
    {
        if(toggle != x->toggle_wave){
            toggle = !toggle;
            fill_wave_table(x, (float)WAVETABLESIZE);
        }

		// calculate the phase increment from the frequency
		// and sample rate - this is the number of cycles per sample
		// freq = cyc/sec, sr = samp/sec, phaseinc = cyc/samp = freq/sr
		phaseincrement = *(freq+sample)/x->samplerate;
		
		// increment the phase and make sure it doesn't go over 1.0
		x->phase += phaseincrement;
		while(x->phase >= 1.0f)
			x->phase -= 1.0f;
		while(x->phase < 0.0f)
			x->phase += 1.0f;

        // get value from interpolation, pass that through the wavefolding function which
        // returns another value. Then pass the new value through the allpass filter function
        // and output it.
        wave_fold_in = quad_interpolate(x);
        allpass_in = wave_folding(wave_fold_in, x->gain, x->offset, x->toggle_wavefold);
        *(out+sample) = allpass(x, allpass_in, x->toggle_allpass);

        // increment the sample
		sample++;
    }
    return (w+5);
}

// function to display oscillator in the external
static void wavefold_allpass_dsp(t_wavefold_allpass *x, t_signal **sp)
{
	// we'll initialize samplerate when starting up
	x->samplerate = sp[0]->s_sr;
    dsp_add(wavefold_allpass_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

static void *wavefold_allpass_new(void)
{
	float size;
    t_wavefold_allpass *x = (t_wavefold_allpass *)pd_new(wavefold_allpass_class);
    outlet_new(&x->x_obj, gensym("signal"));
	// initialize variables
    x->x_f = 0.0f;
	x->phase = 0.0f;
	size = (float)WAVETABLESIZE;

    // create inlets
    // these 3 are for wavefolding
    floatinlet_new(&x->x_obj, &x->gain);
    floatinlet_new(&x->x_obj, &x->offset);
    floatinlet_new(&x->x_obj, &x->toggle_wavefold);

    // these three are for filter
    floatinlet_new(&x->x_obj, &x->cutoff);
    floatinlet_new(&x->x_obj, &x->bandwidth);
    floatinlet_new(&x->x_obj, &x->toggle_allpass);

    // this one is to change the wave
    floatinlet_new(&x->x_obj, &x->toggle_wave);
	
	// space for WAVETABLESIZE samples
	x->wavetable = (t_float *)malloc(WAVETABLESIZE * sizeof(t_float));
	
    // fills the wave table with values for either sine or square wave
    fill_wave_table(x, size);
	
    return (x);
}

// delete allocated memory
static void wavefold_allpass_free(t_wavefold_allpass *x)
{
	free(x->wavetable);
}

// setup external
void wavefold_allpass_tilde_setup(void)
{
    wavefold_allpass_class = class_new(gensym("wavefold_allpass~"), (t_newmethod)wavefold_allpass_new, (t_method)wavefold_allpass_free,
    	sizeof(t_wavefold_allpass), 0, A_DEFFLOAT, 0);
    CLASS_MAINSIGNALIN(wavefold_allpass_class, t_wavefold_allpass, x_f);
    class_addmethod(wavefold_allpass_class, (t_method)wavefold_allpass_dsp, gensym("dsp"), 0);
}
