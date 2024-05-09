//
//  main.c
//  CReverb
//
//  Created by Hans on 7/2/16.
//
//  This file may be used, distributed and modified freely by anyone,
//  for any purpose, without restrictions.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/sysctl.h>
#include "BMVelocityFilter.h"
#include "BMGetOSVersion.h"
#include "BMUpsampler.h"
#include "BMDownsampler.h"
#include "BMFIRFilter.h"
#include "BMMultiLevelBiquad.h"
#include "BMVelvetNoiseDecorrelator.h"
#include "BMRoundRobin.h"
#include "BMBiquadArray.h"
#include "BMMultiTapDelay.h"
#include "BMReverb.h"
#include "TestReverb.h"
#include "BMLagrangeInterpolation.h"
#include "BMMultiLevelSVF.h"
#include "fastlog.h"
#include "fastpow.h"
#include "BMCompressor.h"
#include <string.h>
#include "BMStereoLagTime.h"
#include "BMBinauralSynthesis.h"
#include "BMWetDryMixer.h"
#include "BMWaveshaping.h"
#include "BMIIRUpsampler2x.h"
#include "BMIIRDownsampler2x.h"
#include "BMVAStateVariableFilter.h"
#include "BMNoiseGate.h"
#include "BMSimpleFDN.h"
#include "BMSFM.h"
#include "BMBetaPDF.h"
#include "BMCrossover.h"
#include "BMMonoToStereo.h"
#include "BMAsymptoticLimiter.h"
#include "BMAudioStreamConverter.h"
#include "BMCloudReverb.h"
#include "BMExportWavFile.h"
#include "BMLongLoopFDN.h"
#include "BMSincUpsampler.h"
#include "BMSpectrogram.h"
#include "BMSincDownsampler.h"
#include "BMStaticDelay.h"
#include "BMDPWOscillator.h"
#include "BMCDBlepOscillator.h"
#include "BMBlipOscillator.h"
#include "BMFFT.h"

#define TESTBUFFERLENGTH 128
#define FFTSIZE 4096





void arrayToFileWithName(float* array, char* name, size_t length){
	// open a file for writing
	FILE* audioFile;
	audioFile = fopen(name, "w+");
	
	// print out the entire frame in .csv format
	//fprintf(audioFile,"{");
	for (size_t i=0; i<length; i++) {
		fprintf(audioFile, "%f,%f\n",(float)i/(float)length,array[i]);
		
	}
	
	fclose(audioFile);
	
	system("pwd");
	printf(name);
	printf("\n");
}



void printRGBA(FILE *f, uint8_t *p){
	fprintf(f, "{%i, %i, %i, %i}",p[0], p[1], p[2], p[3]);
}



void printRGBAArray(FILE *f, uint8_t *column, size_t length){
	fprintf(f,"{");
	for(size_t i=0; i<length-4; i+=4){
		printRGBA(f,&column[i]);
		fprintf(f,", ");
	}
	printRGBA(f,&column[length-4]);
	fprintf(f, "}");
}




void imageToFileWithName(uint8_t *image, char *name, size_t width, size_t height){
	// open a file for writing
	FILE* af;
	af = fopen(name, "w+");
	
    size_t bytesPerPixel = 4;
    
	// print out the image in array-as-text format
	fprintf(af, "{");
	for(size_t i=0; i<width-1; i++){
		printRGBAArray(af, &image[i*bytesPerPixel], height*bytesPerPixel);
		fprintf(af,", ");
	}
	printRGBAArray(af, &image[(width-1)*bytesPerPixel], height*bytesPerPixel);
	fprintf(af, "}");

	fclose(af);
	
	system("pwd");
	printf(name);
	printf("\n");
}



double frequencyToPhaseIncrement(double frequency, double sampleRate){
	return (frequency * 2.0 * M_PI) / sampleRate;
}


/*!
 *generateSineSweep
 */
void generateSineSweep(float* output, double startFrequency, double endFrequency, double sampleRate, size_t numSamples){
	/*
	 * logarithmic sweep, direct method
	 */
	double phase = 0.0f;
	//double frequency = startFrequency;
	double startPhaseIncrement = frequencyToPhaseIncrement(startFrequency, sampleRate);
	double endPhaseIncrement = frequencyToPhaseIncrement(endFrequency, sampleRate);
	double phaseIncrement = startPhaseIncrement;
	double phaseIncrementCommonRatio = pow(endPhaseIncrement/startPhaseIncrement,1.0/(double)numSamples);
	for(size_t i=0; i<numSamples; i++){
		output[i] = (float)sin(phase);
		phase += phaseIncrement;
		phaseIncrement *= phaseIncrementCommonRatio;
		if(phase > 2.0 * M_PI) phase -= 2.0 * M_PI;
	}
}





void testDelay(){
	float sampleRate = 44100.;
	size_t delayPerTap = 6;
	size_t numTaps = 50;
	size_t maxNumTaps = 500;
	size_t maxDelayTime = 3000;
	float maxDelayTimeS = maxDelayTime/sampleRate;
	
	size_t* delayTimeL = malloc(sizeof(size_t)*maxNumTaps);
	size_t* delayTimeR = malloc(sizeof(size_t)*maxNumTaps);
	float* gainL = malloc(sizeof(float)*maxNumTaps);
	float* gainR = malloc(sizeof(float)*maxNumTaps);
	
	//Init delayTime
	
	for(int i=0;i<numTaps;i++){
		delayTimeL[i] = i * delayPerTap;
		delayTimeR[i] = i * delayPerTap;
		gainL[i] = BMReverbDelayGainFromRT30(maxDelayTimeS, delayTimeL[i]/sampleRate);
		gainR[i] = BMReverbDelayGainFromRT30(maxDelayTimeS, delayTimeR[i]/sampleRate);
	}
	
	BMMultiTapDelay delay;
	BMMultiTapDelay_Init(&delay, true, delayTimeL, delayTimeR, maxDelayTime, gainL, gainR, numTaps, maxNumTaps);
	
	BMMultiTapDelay_impulseResponse(&delay);
	
	delayPerTap = 16;
	numTaps = 100;
	for(int i=0;i<numTaps;i++){
		delayTimeL[i] = i * delayPerTap;
		delayTimeR[i] = i * delayPerTap;
		gainL[i] = BMReverbDelayGainFromRT30(maxDelayTimeS, delayTimeL[i]/sampleRate);
		gainR[i] = BMReverbDelayGainFromRT30(maxDelayTimeS, delayTimeR[i]/sampleRate);
	}
	
	
	BMMultiTapDelay_setDelayTimeNumTap(&delay, delayTimeL, delayTimeR, numTaps);
	BMMultiTapDelay_setGains(&delay, gainL, gainR);
	
	
}

void testReverb(){
	// open a file for writing
	FILE* audioFile;
	audioFile = fopen("./rvImpulse.csv", "w+");
	
	
	float testBufferInL [TESTBUFFERLENGTH];
	float testBufferInR [TESTBUFFERLENGTH];
	float testBufferOutL [TESTBUFFERLENGTH];
	float testBufferOutR [TESTBUFFERLENGTH];
	
	
	// create the initial impulse followed by zeros
	testBufferInL[0] = 1.0f;
	testBufferInR[0] = 1.0f;
	memset(testBufferInL+1, 0, sizeof(float)*(TESTBUFFERLENGTH-1));
	memset(testBufferInR+1, 0, sizeof(float)*(TESTBUFFERLENGTH-1));
	
	
	// process the first frame twice (the first time to trigger an update)
	//BMCReverbProcessBuffer(&rv, testBufferInL, testBufferInR, testBufferOutL, testBufferOutR, TESTBUFFERLENGTH);
	//BMCReverbProcessBuffer(&rv, testBufferInL, testBufferInR, testBufferOutL, testBufferOutR, TESTBUFFERLENGTH);
	
	// print out the entire frame in .csv format
	for (size_t i=0; i<TESTBUFFERLENGTH; i++) {
		fprintf(audioFile, "%f,%f\n",testBufferOutL[i],testBufferOutR[i]);
	}
	
	// set the input buffers to all zeros (only the first value was non-zero)
	testBufferInL[0] = testBufferInR[0] = 0.0;
	
	
	// start a timer
	clock_t begin, end;
	double time_spent;
	begin = clock();
	
	
	
	
	// print the time taken to process reverb
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time: %f\n", time_spent);
	
	
	fclose(audioFile);
}


void testFIRFilter(){
	float kernel [17] =  {0.0, -0.00130751,
		0.00478209, -0.00725887, 0.0, 0.0255661, -0.0672996,
		0.107966, 0.375, 0.107966, -0.0672996, 0.0255661,
		0.0, -0.00725887,
		0.00478209, -0.00130751, 0.0};
	
	BMFIRFilter f;
	BMFIRFilter_init(&f, kernel, 17);
	
	float* IR = malloc(sizeof(float)*f.length);
	BMFIRFilter_impulseResponse(&f, IR);
	
	printf("\nFIR filter impulse response:\n{");
	for(size_t i=0; i < f.length-1; i++)
		printf("%f,",IR[i]);
	printf("%f}\n",IR[f.length-1]);
}

void testButterworthFilter(){
	float IR [1000];
	
	memset(IR,0,sizeof(float)*1000);
	IR[0]=1.0;
	
	BMMultiLevelBiquad filter;
	// init a filter at 48000 Hz sample rate with four biquad sections
	BMMultiLevelBiquad_init(&filter, 4, 48000, false, true,false);
	
	// init a butterworth lowpass with fc at 7000 hz,
	// starting at the first section of the filter
	// and consuming all four sections
	// this will be a 4*2=8th order filter
	BMMultiLevelBiquad_setHighOrderBWLP(&filter, 7000.0, 0, 4);
	
	// generate the impulse response
	BMMultiLevelBiquad_processBufferMono(&filter, IR, IR, 1000);
	
	// print the impulse response
	printf("\nButterworth filter test output:\n{");
	for(size_t i=0; i < 999; i++)
		printf("%f,",IR[i]);
	printf("%f}\n",IR[999]);
}

void testVND(){
	BMVelvetNoiseDecorrelator vndFilter;
	
	float sampleRate = 48000.0f;
	size_t testLength = sampleRate * 0.5;
	float* inputL = calloc(testLength*4,sizeof(float));
	float* inputR = inputL + testLength;
	float* outputL = inputR + testLength;
	float* outputR = outputL + testLength;
	
	float maxDelay = 1.0 / 3.0;
	size_t numTaps = 16;
	float rt60 = 100.0f;
	bool hasDryTap = false;
	BMVelvetNoiseDecorrelator_initWithEvenTapDensity(&vndFilter, maxDelay, numTaps, rt60, hasDryTap, sampleRate);
	
	inputL[0] = 1.0f;
	inputR[0] = 1.0f;
	
	for(size_t i=0; i<10000; i++)
	BMVelvetNoiseDecorrelator_processBufferStereo(&vndFilter,
												  inputL, inputR,
												  outputL, outputR,
												  testLength);
	
	
	char* filename = "./impulseResponse.csv";
	arrayToFileWithName(outputL, filename, testLength);
	
	free(inputL);
	inputL = NULL;
	BMVelvetNoiseDecorrelator_free(&vndFilter);
}


void testMonoToStereo(){
	float sampleRate = 48000.0f;
	BMMonoToStereo m2s;
	BMMonoToStereo_init(&m2s, sampleRate, true);
	
	
	size_t testLength = sampleRate * 5;
	float* inputL = calloc(testLength*4,sizeof(float));
	float* inputR = inputL + testLength;
	float* outputL = inputR + testLength;
	float* outputR = outputL + testLength;
	
	inputL[(size_t)(sampleRate * 0.5)] = 1.0f;
	inputR[(size_t)(sampleRate * 0.5)] = 1.0f;
	
	BMMonoToStereo_processBuffer(&m2s,
								 inputL, inputR,
								 outputL, outputR,
								 testLength);
	
	char* filename = "./impulseResponse.csv";
	arrayToFileWithName(outputL, filename, testLength);
	
	float* sineSweep = malloc(sizeof(float)*testLength);
	
	generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
	
	BMMonoToStereo_processBuffer(&m2s,
								 sineSweep, sineSweep,
								 outputL, outputR,
								 testLength);
	
	char* filename2 = "./sineSweepResponse.csv";
	arrayToFileWithName(outputL, filename2, testLength);
	
	free(inputL);
	BMMonoToStereo_free(&m2s);
}


void testRRStartNote(BMRoundRobin rrb){
	//Test
	BMRoundRobin_NewNote(&rrb);
	BMRoundRobin_testImpulseResponse(&rrb);
	//Test again
	BMRoundRobin_NewNote(&rrb);
	BMRoundRobin_testImpulseResponse(&rrb);
	//3rd times
	BMRoundRobin_NewNote(&rrb);
	BMRoundRobin_testImpulseResponse(&rrb);
}

void testRoundRobin(){
	BMRoundRobin rrb;
	float freq = 1200;
	for(freq=300;freq<=1800;freq+=300){
		printf("freq %f ==============\n",freq);
		BMRoundRobin_init(&rrb, 48000, freq, 6);
		testRRStartNote(rrb);
	}
}



void testBiquadArray(){
	BMBiquadArray bqa;
	float delayTimes = .2;
	float fc = 5000.0;
	float unfilteredRT60 = 1.0;
	float filteredRT60 = 0.5;
	size_t numChannels = 1;
	float sampleRate = 48000;
	
	BMBiquadArray_init(&bqa, numChannels, sampleRate);
	BMBiquadArray_setLowDecayFDN(&bqa, &delayTimes, fc, unfilteredRT60, filteredRT60, numChannels);
	
	BMBiquadArray_impulseResponse(&bqa);
	
	BMBiquadArray_free(&bqa);
}


void testTestReverb(){
	TestReverb rv;
	TestReverb_init(&rv, 16, 44100.);
	
	TestReverb_impulseResponse(&rv, 4410);
}

void biquadIR(vDSP_biquadm_Setup setup, float* input,float* output, int length){
	const float* inputP [1];
	float* outputP [1];
	inputP[0] = input;
	outputP[0] = output;
	vDSP_biquadm(setup, inputP, 1, outputP, 1, length);
	printf("%f %f %f\n",output[0],output[1],output[2]);
	printf("------------------\n");
}

void testBiquad(){
	vDSP_Length Channels = 1;
	vDSP_Length Sections = 1;
	double *F = malloc(5*Channels*Sections*sizeof *F);
	for (vDSP_Length j = 0; j < Sections; ++j)
	{
		F[j*Channels*5 + 0*5 + 0] = 0.1;    //b0
		F[j*Channels*5 + 0*5 + 1] = 0.1;    //b1
		F[j*Channels*5 + 0*5 + 2] = 0.1;    //b2
		F[j*Channels*5 + 0*5 + 3] = 0;    //a1
		F[j*Channels*5 + 0*5 + 4] = 0;    //a2
	}
	
	int length = 5;
	float* input = malloc(length*sizeof(input));
	for(int i=0;i<length;i++){
		if(i==0)
			input[i] = 1;
		else
			input[i] = 0;
	}
	float* output = malloc(length*sizeof(output));
	
	vDSP_biquadm_Setup Setup = vDSP_biquadm_CreateSetup(F, Sections, Channels);
	
	//New F
	for (vDSP_Length j = 0; j < Sections; ++j)
	{
		F[j*Channels*5 + 0*5 + 0] = 0.04;    //b0
		F[j*Channels*5 + 0*5 + 1] = 0.04;    //b1
		F[j*Channels*5 + 0*5 + 2] = 0.04;    //b2
		F[j*Channels*5 + 0*5 + 3] = 0;    //a1
		F[j*Channels*5 + 0*5 + 4] = 0;    //a2
	}
	
	float rate = 0.995;
	float threshold = 0.05;
	
	vDSP_biquadm_SetTargetsDouble(Setup, F, rate, threshold, 0, 0, Sections, Channels);
	//    vDSP_biquadm_SetCoefficientsSingle(Setup, (const float* _Nonnull)F, 0, 0, Sections, Channels);
	//    vDSP_biquadm_SetCoefficientsDouble(Setup, F, 0, 0, Sections, Channels);
	int count = 0;
	for(int i=0;i<100;i++){
		count++;
		biquadIR(Setup, input, output, length);
		//        if(i==5){
		//            for (vDSP_Length j = 0; j < Sections; ++j)
		//            {
		//                F[j*Channels*5 + 0*5 + 0] = 1.;    //b0
		//                F[j*Channels*5 + 0*5 + 1] = 1.;    //b1
		//                F[j*Channels*5 + 0*5 + 2] = 1.;    //b2
		//                F[j*Channels*5 + 0*5 + 3] = 0;    //a1
		//                F[j*Channels*5 + 0*5 + 4] = 0;    //a2
		//            }
		//            vDSP_biquadm_SetTargetsDouble(Setup, F, rate, threshold, 0, 0, Sections, Channels);
		//        }
	}
	printf("samples %d",count*length);
}

#define LI_inputLength 50
#define LI_order 10
#define LI_Stride_Step 0.1

#define LI_outputLength (LI_inputLength - LI_order) * LI_order
void testLagrangeInterpol(){
	BMLagrangeInterpolation lagInter;
	BMLagrangeInterpolation_init(&lagInter, LI_order);
	
	float* input = malloc(sizeof(float) * LI_inputLength);
	float* stride = malloc(sizeof(float) * LI_outputLength);
	for(int i=0;i<LI_inputLength;i++){
		input[i] = cosf(i/5.);
	}
	//    for(int i=0;i<LI_inputLength;i++){
	//        printf("%f\n",input[i]);
	//    }
	
	for(int i=0;i<LI_outputLength;i++){
		stride[i] = i*LI_Stride_Step + LI_order * 0.5;
		//        printf("%f\n",stride[i]);
	}
	
	float output[LI_outputLength];
	BMLagrangeInterpolation_processUpSample(&lagInter, input, stride, output, 1, LI_inputLength, LI_outputLength);
	for(int i=0;i<LI_outputLength;i++){
		printf("%f\n",output[i]);
	}
	
	BMLagrangeInterpolation_processDownSample(&lagInter, output, 1, input, 1, LI_outputLength, LI_inputLength-LI_order);
	
	for(int i=0;i<LI_inputLength-LI_order;i++){
		printf("%f\n",input[i]);
	}
	
}

void testMultiLevelSVF(){
	BMMultiLevelSVF svf;
	BMMultiLevelSVF_init(&svf, 1, 44100, true);
	
	BMMultiLevelSVF_setBell(&svf, 2000, -40, 2, 0);
	
	BMMultiLevelSVF_impulseResponse(&svf, 1000);
	
	BMMultiLevelSVF_free(&svf);
	
	simd_float4* a = malloc(sizeof(simd_float4));
	a[0].x = 0;
	a[0].y = 1;
	a[0].z = 2;
	a[0].a = 3;
	simd_float2* a2 = (simd_float2*)a;
	printf("%f %f\n",a2[1].x,a2[1].y);
}

void testFastLog(){
	// set up input
	float t [128];
	for(size_t i=0; i<128; i++)
		t[i] = (Float32)(i+1)/128.0;
	
	// run the fast log 2
	vector_fastlog2(t, t, 128);
	
	printf("\nFast log2 test output:\n{");
	for(size_t i=0; i < 128-1; i++)
		printf("%f,",t[i]);
	printf("%f}\n",t[128-1]);
}

void testFasterLog(){
	// set up input
	float t [128];
	for(size_t i=0; i<128; i++)
		t[i] = (Float32)(i+1)/128.0f;
	
	// run the fast log 2
	vector_fasterLog2(t, t, 128);
	
	printf("\nFast log2 test output:\n{");
	for(size_t i=0; i < 128-1; i++)
		printf("%f,",t[i]);
	printf("%f}\n",t[128-1]);
}


void testUserDefinedVectorTypes(){
	vUint32_8 xx;
	vFloat32_8 yy;
	uint32 xi [8] = {0,1,2,3,4,5,6,7};
	Float32 yi [8] = {0.f,1.f,2.f,3.f,4.f,5.f,6.f,7.f};
	xx = *(vUint32_8 *)xi;
	yy = *(vFloat32_8 *)yi;
	
	// print the values in xx using pointer arithmetic
	uint32* p = (uint32*)&xx;
	for(size_t i=0; i<8; i++){
		printf("%u, ",*p++);
	}
	printf("\n");
	
	// print the values in yy using pointer arithmetic
	Float32* pf = (Float32*)&yy;
	for(size_t i=0; i<8; i++){
		printf("%f, ",*pf++);
	}
	printf("\n");
	
	// print the values in yy as ints using pointer arithmetic
	uint32* pfi = (uint32*)&yy;
	for(size_t i=0; i<8; i++){
		printf("%u, ",*pfi++);
	}
	printf("\n");
	
	//    // do a float+uint union on scalars
	//    Float32 testValue = 0.98993;
	//    union { float f; uint32_t i; } vx = { testValue };
	//    union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
	//    printf("vx.f: %f\n",vx.f);
	//    printf("vx.i: %u\n",vx.i);
	//    printf("mx.f: %f\n",mx.f);
	//    printf("mx.i: %u\n",mx.i);
	
	//    // do a float+uint union on vectors
	//    vFloat32_8 vTestValue;
	//    for(size_t i=0; i<8; i++) *((Float32*)&vTestValue+i) = testValue;
	//    union { vFloat32_8 f; vUint32_8 i; } vxv = { vTestValue };
	//    // print the values in vxv using pointer arithmetic
	//    Float32* vxvp = (Float32*)&vxv;
	//    for(size_t i=0; i<8; i++){
	//        printf("%f, ",*vxvp++);
	//    }
	//    printf("\n");
	//    UInt32* vxvpi = (UInt32*)&vxv;
	//    for(size_t i=0; i<8; i++){
	//        printf("%u, ",*vxvpi++);
	//    }
	//    printf("\n");
	
	// init a float vector like this
	Float32 foo [8] = {1.19,1.19,1.19,1.19,1.19,1.19,1.19,1.19};
	vFloat32_8 foov = *(vFloat32_8 *)&foo;
	Float32* foop = (Float32*)&foov;
	for(size_t i=0; i<8; i++){
		printf("%f, ",*foop++);
	}
	printf("\n");
	
	
	/*************************************************************
	 * We test the code from the vector_fastLog2 function itself *
	 *************************************************************/
	
	// set up input and output
	Float32 x = 0.5;
	Float32 log2x;
	Float32 X [8] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
	Float32 log2X [8];
	
	// set up some magic constants
	uint32 c1 [8] = {0x007FFFFF, 0x007FFFFF, 0x007FFFFF, 0x007FFFFF,
		0x007FFFFF, 0x007FFFFF, 0x007FFFFF, 0x007FFFFF};
	vUint32_8 c1v = *(vUint32_8 *)&c1;
	UInt32 c1s = c1[0];
	uint32 c2 [8] = {0x3f000000,0x3f000000,0x3f000000,0x3f000000,
		0x3f000000,0x3f000000,0x3f000000,0x3f000000};
	vUint32_8 c2v = *(vUint32_8 *)&c2;
	UInt32 c2s = c2[0];
	Float32 c3 [8] = {1.1920928955078125e-7f,1.1920928955078125e-7f,
		1.1920928955078125e-7f,1.1920928955078125e-7f,
		1.1920928955078125e-7f,1.1920928955078125e-7f,
		1.1920928955078125e-7f,1.1920928955078125e-7f};
	vFloat32_8 c3v = *(vFloat32_8 *)&c3;
	Float32 c3s = c3[0];
	Float32 c4 [8] = {124.22551499f,124.22551499f,124.22551499f,124.22551499f,
		124.22551499f,124.22551499f,124.22551499f,124.22551499f};
	vFloat32_8 c4v = *(vFloat32_8 *)&c4;
	Float32 c4s = c4[0];
	Float32 c5 [8] = {1.498030302f, 1.498030302f, 1.498030302f, 1.498030302f,
		1.498030302f, 1.498030302f, 1.498030302f, 1.498030302f};
	vFloat32_8 c5v = *(vFloat32_8 *)&c5;
	Float32 c5s = c5[0];
	Float32 c6 [8] = {1.72587999f, 1.72587999f, 1.72587999f, 1.72587999f,
		1.72587999f, 1.72587999f, 1.72587999f, 1.72587999f};
	vFloat32_8 c6v = *(vFloat32_8 *)&c6;
	Float32 c6s = c6[0];
	Float32 c7 [8] = {0.3520887068f, 0.3520887068f, 0.3520887068f, 0.3520887068f,
		0.3520887068f, 0.3520887068f, 0.3520887068f, 0.3520887068f};
	vFloat32_8 c7v = *(vFloat32_8 *)&c7;
	Float32 c7s = c7[0];
	
	
	union { float f; uint32_t i; } vx = { x };
	union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
	union { vFloat32_8 f; vUint32_8 i; } vX = { *(vFloat32_8 *)X };
	union { vUint32_8 i; vFloat32_8 f; } mX = { (vX.i & c1v) | c2v };
	
	Float32 ys = vx.i;
	vFloat32_8 y = __builtin_convertvector(vX.i,vFloat32_8);
	//        vFloat32_8 y = (vFloat32_8){(Float32)vX.i[0],(Float32)vX.i[1],
	//                                    (Float32)vX.i[2],(Float32)vX.i[3],
	//                                    (Float32)vX.i[4],(Float32)vX.i[5],
	//                                    (Float32)vX.i[6],(Float32)vX.i[7]};
	
	ys *= 1.1920928955078125e-7f;
	y *= c3v;
	
	log2x = ys - 124.22551499f
	- 1.498030302f * mx.f
	- 1.72587999f / (0.3520887068f + mx.f);
	*(vFloat32_8 *)log2X = y - c4v
	- c5v * mX.f
	- c6v / (c7v + mX.f);
	
	Float32* pff = (Float32*)&log2X;
	for(size_t i=0; i<8; i++){
		printf("%f, ",*pff++);
	}
	printf("\n");
	
	printf("scalar log output: %f\n",log2x);
}


void testNEONConditionals(){
	vFloat a = {0.0f,0.0f,0.0f,0.0f};
	vFloat b = {-1.0f,0.0f,1.0f,0.0f};
	vSInt32 c;
	c = a < b;
	for(size_t i=0; i<4; i++) printf("%i,",c[i]);
}

void testInRangef(){
	float testInput [32];
	float ll = 1;
	float ul = 2;
	float x=0;
	for(size_t i=0; i<32; i++) {
		testInput[i]=x;
		x += 0.25;
	}
	float testOutput [32];
	BMvInRangef(testInput,ll,ul,testOutput,32);
	
	for(size_t i=0; i<32; i++) printf("%f,",testInput[i]);
	printf("\n");
	for(size_t i=0; i<32; i++) printf("%f,",testOutput[i]);
}


/*
  *This function behaves like the function g(x) = x for x near zero but
 * its response tapers off asymptotically to +- limit as abs(x) increases.
 *
  *This FUNCTION DOES NOT WORK IN PLACE. (requires input != output)
 *
 * Function definition:
 *   f(x,limit) = x / (1 + |x/limit|)
 *
 * Function properties:
 *   f(0) = 0; f'(0) = 1
 *   f(+-inf) = +-limit; f'(+-inf) = 0
 *
 * Mathematica prototype code:
 *
 *  asymptoticLimit[input_, limit_] := input/(1 + Abs[input/limit])
 */
void asymptoticLimitTest1(float limit,
						  float* input,
						  float* output,
						  size_t numSamples){
	
	vFloat32_8 limitv = limit;
	vFloat32_8 onev = 1.0f;
	
	size_t n = numSamples;
	// process in chunks of 8 samples
	while(n >= 8){
		// vFloat32_8 d = vfabsf_8(*(vFloat32_8 *)input / limitv) + onev;
		// *(vFloat32_8 *)output = *(vFloat32_8 *)input / d;
		*(vFloat32_8 *)output = vfabsf_8(*(vFloat32_8 *)input);
		
		input += 8;
		output += 8;
		n -= 8;
	}
	
	// process the remaining samples individually
	while(n > 0){
		float d = fabsf(*input / limit) + 1.0f;
		*output = *input / d;
		
		input++;
		output++;
		n--;
	}
}


void asymptoticLimitTest2(float limit,
						  const float* input,
						  float* output,
						  size_t numSamples){
	
	// input / limit => output
	vDSP_vsdiv(input, 1, &limit, output, 1, numSamples);
	
	// abs(output) => output
	vDSP_vabs(output, 1, output, 1, numSamples);
	
	// 1 + output => output
	float one = 1.0;
	vDSP_vsadd(output, 1, &one, output, 1, numSamples);
	
	// input / output => output
	vDSP_vdiv(output, 1, input, 1, output, 1, numSamples);
}


void asymptoticLimitTest3(float limit,
						  const float* input,
						  float* output,
						  size_t numSamples){
	
	vFloat limitv = {limit, limit, limit, limit};
	vFloat onev = {1.0,1.0,1.0,1.0};
	
	size_t n = numSamples;
	// process in chunks of 4 samples
	while(n >= 4){
		vFloat d = vfabsf(*(vFloat *)input / limitv) + onev;
		*(vFloat *)output = *(vFloat *)input / d;
		
		input += 4;
		output += 4;
		n -= 4;
	}
	
	// process the remaining samples individually
	while(n > 0){
		float d = fabsf(*input / limit) + 1.0;
		*output = *input / d;
		
		input++;
		output++;
		n--;
	}
}


/*
 * On Intel, this is the fastest
 */
void asymptoticLimitTest4(float limit,
						  const float* input,
						  float* output,
						  size_t numSamples){
	
	//size_t n = numSamples;
	
	// process the remaining samples individually
	while(numSamples > 0){
		float d = fabsf(*input / limit) + 1.0;
		*output = *input / d;
		
		input++;
		output++;
		numSamples--;
	}
}



void asymptoticLimitTest5(float limit,
						  const float* input,
						  float* output,
						  size_t numSamples){
	
	for(size_t i=0; i<numSamples; i++)
		output[i] = input[i] / (fabsf(input[i] / limit) + 1.0);
}


/*
  *This is fastest on ARM processors
 */
void asymptoticLimitTest6(float limit,
						  float* input,
						  float* output,
						  size_t numSamples){
	
	
	size_t n = numSamples;
	// process in chunks of 32 samples
	while(n >= 32){
		vFloat32_32 d = vfabsf_32(*(vFloat32_32 *)input / limit) + 1.0f;
		*(vFloat32_32 *)output = *(vFloat32_32 *)input / d;
		
		input += 32;
		output += 32;
		n -= 32;
	}
	
	// process the remaining samples individually
	while(n > 0){
		float d = fabsf(*input / limit) + 1.0f;
		*output = *input / d;
		
		input++;
		output++;
		n--;
	}
}


void testAsymptoticLimit(){
	float t [256];
	float o [256];
	for(size_t i = 0; i<256; i++) t[i] = (float)i/64.0;
	
	// start a timer
	clock_t begin, end;
	double time_spent_1,time_spent_2, time_spent_3, time_spent_4, time_spent_5, time_spent_6;
	time_spent_1 = time_spent_2 = time_spent_3 = time_spent_4 = time_spent_5 = time_spent_6 = 0;
	
	for(size_t i=0; i < 1000000; i++){
		begin = clock();
		asymptoticLimitTest1(1.0, t, o, 256);
		end = clock ();
		time_spent_1 += (double)(end - begin);
		
		begin = clock();
		asymptoticLimitTest2(1.0, t, o, 256);
		end = clock ();
		time_spent_2 += (double)(end - begin);
		
		begin = clock();
		asymptoticLimitTest3(1.0, t, o, 256);
		end = clock ();
		time_spent_3 += (double)(end - begin);
		
		begin = clock();
		asymptoticLimitTest4(1.0, t, o, 256);
		end = clock ();
		time_spent_4 += (double)(end - begin);
		
		begin = clock();
		asymptoticLimitTest5(1.0, t, o, 256);
		end = clock ();
		time_spent_5 += (double)(end - begin);
		
		begin = clock();
		asymptoticLimitTest6(1.0, t, o, 256);
		end = clock ();
		time_spent_6 += (double)(end - begin);
	}
	
	// print the time taken to process
	end = clock();
	time_spent_1 /= CLOCKS_PER_SEC;
	time_spent_2 /= CLOCKS_PER_SEC;
	time_spent_3 /= CLOCKS_PER_SEC;
	time_spent_4 /= CLOCKS_PER_SEC;
	time_spent_5 /= CLOCKS_PER_SEC;
	time_spent_6 /= CLOCKS_PER_SEC;
	printf("time1: %f\ntime2: %f\ntime3: %f\ntime4: %f\ntime5: %f\ntime6: %f\n",
		   time_spent_1, time_spent_2, time_spent_3, time_spent_4, time_spent_5,time_spent_6);
	
}

/*
 * Vectorised softKnee function 1
 */
static inline void softKnee1(const float* overshoot,
							 const float* inTransition,
							 float* rect,
							 float halfKneeFraction,
							 float halfKneeWidthInDB,
							 size_t length) {
	//rect = inTransition*(halfKneeFraction * fastpow(overshoot + halfKneeWidthInDB, 2)) + ((1 - inTransition)*fmaxf(overshoot, 0));
	size_t n = length;
	
	// process in chunks of 32 samples
	while (n >= 32) {
		vFloat32_32 tv = *(vFloat32_32 *)overshoot + halfKneeWidthInDB;
		*(vFloat32_32 *)rect = *(vFloat32_32 *)inTransition
		* (halfKneeFraction*tv*tv)
		+ ((1.0 - *(vFloat32_32 *)inTransition)
		   * vfClipNeg32(*(vFloat32_32 *)overshoot));
		
		// advance pointers
		rect += 32;
		overshoot += 32;
		inTransition += 32;
		n -= 32;
	}
	
	// finish up the last few values
	while (n > 0) {
		// compute the soft knee
		float t = *overshoot + halfKneeWidthInDB;
		*rect = *inTransition * (halfKneeFraction * t * t)
		+ ((1.0 - *inTransition)*fmaxf(*overshoot, 0.0));
		
		// advance pointers
		rect++;
		overshoot++;
		inTransition++;
		n--;
	}
}


/*
 * Vectorised softKnee function 2
 */
static inline void softKnee2(const float* overshoot,
							 const float* inTransition,
							 float* rect,
							 float halfKneeFraction,
							 float halfKneeWidthInDB,
							 size_t length) {
	//rect = inTransition*(halfKneeFraction * fastpow(overshoot + halfKneeWidthInDB, 2)) + ((1 - inTransition)*fmaxf(overshoot, 0));
	size_t n = length;
	
	while (n > 0) {
		// compute the soft knee
		float t = *overshoot + halfKneeWidthInDB;
		*rect = *inTransition * (halfKneeFraction * t * t)
		+ ((1.0 - *inTransition)*fmaxf(*overshoot, 0.0));
		
		// advance pointers
		rect++;
		overshoot++;
		inTransition++;
		n--;
	}
}

/*
 * Vectorised softKnee function 3
 */
static inline void softKnee3(const float* overshoot,
							 const float* inTransition,
							 float* rect,
							 float halfKneeFraction,
							 float halfKneeWidthInDB,
							 size_t length) {
	//rect = inTransition*(halfKneeFraction * fastpow(overshoot + halfKneeWidthInDB, 2)) + ((1 - inTransition)*fmaxf(overshoot, 0));
	
	for(size_t i=0; i< length; i++){
		float t = overshoot[i] + halfKneeWidthInDB;
		rect[i] = inTransition[i] * (halfKneeFraction * t * t)
		+ ((1.0 - inTransition[i])*fmaxf(overshoot[i], 0.0));
	}
}


/*
 * Vectorised softKnee function 4
 *
 * temp1 may point to the same data as rect but it may not be the same as
 * overshoot or inTransition
 */
static inline void softKnee4(const float* overshoot,
							 const float* inTransition,
							 float* rect,
							 float* temp1,
							 float halfKneeFraction,
							 float halfKneeWidthInDB,
							 size_t length) {
	//rect = inTransition*(halfKneeFraction * fastpow(overshoot + halfKneeWidthInDB, 2)) + ((1 - inTransition)*fmaxf(overshoot, 0));
	
	// compute the soft knee
	// float t = *overshoot + halfKneeWidthInDB;
	// *rect = *inTransition * (halfKneeFraction * t * t)
	// + ((1.0 - *inTransition)*fmaxf(*overshoot, 0.0));
	
	// temp1 = fmaxf(overshoot,0.0);
	float zero = 0.0f;
	vDSP_vthres(overshoot, 1, &zero, temp1, 1, length);
	
	// temp2 = ((1.0 - *inTransition)
	float* temp2 = rect;
	vDSP_vneg(inTransition, 1, temp2, 1, length);
	float one = 1.0f;
	vDSP_vsadd(temp2, 1, &one, temp2, 1, length);
	
	// temp1 = ((1.0 - *inTransition)*fmaxf(*overshoot, 0.0))
	vDSP_vmul(temp1, 1, temp2, 1, temp1, 1, length);
	
	// the values in temp2 will not be used again;
	// we rename it to avoid confusion when reusing the array
	float* temp3 = temp2;
	
	// temp3 = *overshoot + halfKneeWidthInDB;
	vDSP_vsadd(overshoot, 1, &halfKneeWidthInDB, temp3, 1, length);
	
	// temp3 = t ^ 2;
	vDSP_vsq(temp3, 1, temp3, 1, length);
	
	// temp3 = halfKneeFraction * t * t
	vDSP_vsmul(temp3,1,&halfKneeFraction,temp3,1,length);
	
	// temp3 = *inTransition * (halfKneeFraction * t * t)
	vDSP_vmul(inTransition,1,temp3,1,temp3,1,length);
	
	// *rect = *inTransition * (halfKneeFraction * t * t)
	//       + ((1.0 - *inTransition)*fmaxf(*overshoot, 0.0));
	vDSP_vadd(temp1,1,temp2,1,rect,1,length);
}


void testSoftKnee(){
	float t [256];
	float o [256];
	float temp1 [256];
	float temp2 [256];
	float temp3 [256];
	float temp4 [256];
	for(size_t i = 0; i<256; i++) t[i] = 3.0 * sinf(i/64.0);
	size_t length = 256;
	float thresholdInDB = 0.0f;
	float halfKneeWidthInDB = 1.0f;
	
	
	/*
	 * Compute the volume in dB
	 */
	float* offsetInput = t;
	float* volumeLinear = temp1;
	vvfabsf(volumeLinear,offsetInput,(int*)&length);
	float* volumeDb = temp1;
	float one = 1.0f;
	vDSP_vdbcon(volumeLinear, 1, &one, volumeDb, 1, length, 1);
	
	
	/*
	 * Compute the overshoot in dB
	 * overshoot = volume - threshold
	 */
	float* overshoot = temp1;
	float negativeThreshold = -thresholdInDB;
	vDSP_vsadd(volumeDb, 1, &negativeThreshold, overshoot, 1, length);
	
	
	float negHalfKneeWidthInDB = - halfKneeWidthInDB;
	float* inTransition = temp2;
	
	// check which of the overshoot values are in the soft
	// knee section
	BMvInRangef(overshoot,
				negHalfKneeWidthInDB,
				halfKneeWidthInDB,
				inTransition,
				length);
	
	
	// start a timer
	clock_t begin, end;
	double time_spent_1,time_spent_2, time_spent_3, time_spent_4, time_spent_5, time_spent_6;
	time_spent_1 = time_spent_2 = time_spent_3 = time_spent_4 = time_spent_5 = time_spent_6 = 0;
	
	float* rect = temp3;
	float halfKneeFraction = 0.5;
	
	for(size_t i=0; i < 1000000; i++){
		begin = clock();
		softKnee1(overshoot, inTransition, rect, halfKneeFraction, halfKneeWidthInDB, length);
		end = clock ();
		time_spent_1 += (double)(end - begin);
		
		begin = clock();
		softKnee2(overshoot, inTransition, rect, halfKneeFraction, halfKneeWidthInDB, length);
		end = clock ();
		time_spent_2 += (double)(end - begin);
		
		begin = clock();
		softKnee3(overshoot, inTransition, rect, halfKneeFraction, halfKneeWidthInDB, length);
		end = clock ();
		time_spent_3 += (double)(end - begin);
		
		begin = clock();
		softKnee4(overshoot, inTransition, rect, temp4, halfKneeFraction, halfKneeWidthInDB, length);
		end = clock ();
		time_spent_4 += (double)(end - begin);
		
	}
	
	// print the time taken to process
	end = clock();
	time_spent_1 /= CLOCKS_PER_SEC;
	time_spent_2 /= CLOCKS_PER_SEC;
	time_spent_3 /= CLOCKS_PER_SEC;
	time_spent_4 /= CLOCKS_PER_SEC;
	time_spent_5 /= CLOCKS_PER_SEC;
	time_spent_6 /= CLOCKS_PER_SEC;
	printf("time1: %f\ntime2: %f\ntime3: %f\ntime4: %f\ntime5: %f\ntime6: %f\n",
		   time_spent_1, time_spent_2, time_spent_3, time_spent_4, time_spent_5,time_spent_6);
	
}


void testAsymptoticLimitOutputRange(){
	float t [256];
	float o [256];
	for(size_t i = 0; i<256; i++) t[i] = (float)i*64.0;
	
	asymptoticLimitTest1(1.0f, t, o, 256);
	
	printf("{");
	for(size_t i = 0; i<255; i++) printf("%f,",o[i]);
	printf("%f}\n",o[255]);
	
}





void arrayToFile(float* array, size_t length){
	char* name = "./arrayOut.csv";
	arrayToFileWithName(array, name, length);
}




void arrayXYToFile(float* arrayX,float* arrayY, size_t length){
	// open a file for writing
	FILE* audioFile;
	audioFile = fopen("./arrayOut.csv", "w+");
	
	// print out the entire frame in .csv format
	//fprintf(audioFile,"{");
	for (size_t i=0; i<length; i++) {
		fprintf(audioFile, "%f,%f\n",arrayX[i],arrayY[i]);
		
	}
	
	fclose(audioFile);
	
	system("pwd");
	printf("/arrayOut.csv\n");
}

#define BM_GAINSTAGE_TEST_OVERSAMPLE 8
#define BM_GAINSTAGE_TEST_SAMPLE_RATE 48000 * BM_GAINSTAGE_TEST_OVERSAMPLE
#define BM_GAINSTAGE_TEST_LENGTH BM_GAINSTAGE_TEST_SAMPLE_RATE / 7


#define BMCOMPRESSOR_TEST_SAMPLERATE 48000
#define BMCOMPRESSOR_TEST_LENGTH BMCOMPRESSOR_TEST_SAMPLERATE

void testCompressor(){
	BMCompressor c;
	
	BMCompressor_initWithSettings(&c, BMCOMPRESSOR_TEST_SAMPLERATE, -10.0f, 5.0f, 40.0f, .001f, 0.05f);
	
	// create sine wave for input test
	float frequency = 300.0f;
	float t [BMCOMPRESSOR_TEST_LENGTH];
	float o [BMCOMPRESSOR_TEST_LENGTH];
	float sr = BMCOMPRESSOR_TEST_SAMPLERATE;
	for (size_t i=0; i<BMCOMPRESSOR_TEST_LENGTH; i++) {
		t[i] = sinf(frequency * 2.0 * M_PI * (float)i / sr);
		
		if (i>8000 && i < 20000) t[i] *= BM_DB_TO_GAIN(-20.0f);
	}
	
	float minGain;
	BMCompressor_ProcessBufferMono(&c, t, o, &minGain, BMCOMPRESSOR_TEST_LENGTH);
	
	// print the output
	arrayToFile(o,BMCOMPRESSOR_TEST_LENGTH);
}


float asymThresholdTestFunction(float x, float threshold){
	float invThreshold = 1.0 / threshold;
	float xIfNeg = (x - fabsf(x)) * 0.5f;
	return x / (1.0f + (xIfNeg * invThreshold));
}

void testAsymptoticThreshold(){
	float t [256];
	for(size_t i=0; i<256; i++) {
		t[i] = asymThresholdTestFunction(1.0f - (float)i * 0.1f,-10.0);
	}
	
	arrayToFile(t,256);
}


void testCriticallyDampedFilterStepResponse(){
	float sampleRate = 48000.0f;
	
	float attackTimeInSeconds = 0.05;
	float alpha = 1.5; // an emperically determined magic number
	float filterFrequency = 1.0 / (alpha * attackTimeInSeconds);
	printf("filter frequency: %f\n", filterFrequency);
	
	BMMultiLevelBiquad This;
	BMMultiLevelBiquad_init(&This, 4, sampleRate, false, true, false);
	BMMultiLevelBiquad_setCriticallyDampedLP(&This, filterFrequency, 0, 4);
	
	float timeBeforeStart = 0.1;
	float stepLength = 1.0;
	size_t beforeStartNumSamples = sampleRate * timeBeforeStart;
	size_t stepLengthNumSamples = sampleRate * stepLength;
	size_t totalLength = beforeStartNumSamples + stepLengthNumSamples;
	
	float* input = malloc(sizeof(float) * totalLength);
	
	size_t i = 0;
	for(; i<beforeStartNumSamples; i++){
		input[i] = 0.0f;
	}
	for(; i<totalLength; i++){
		input[i] = 1.0f;
	}
	
	BMMultiLevelBiquad_processBufferMono(&This, input, input, totalLength);
	
	arrayToFile(input+beforeStartNumSamples, stepLengthNumSamples);
	
	// find time till -0.3db
	float pointThreeDb = BM_DB_TO_GAIN(-0.3f);
	
	bool done = false;
	for(i=beforeStartNumSamples; i<totalLength && !done; i++){
		if (input[i] >= pointThreeDb){
			printf("time to -0.3db: %f\n", ((float)i - (float)beforeStartNumSamples)/ (float)stepLengthNumSamples);
			done = true;
		}
	}
	
	free(input);
	input = NULL;
}



void testCriticallyDampedFilterFrequencyResponse(){
	float sampleRate = 48000.0f;
	
	float attackTimeInSeconds = 0.05;
	float alpha = 1.5; // an emperically determined magic number
	float filterFrequency = 1.0 / (alpha * attackTimeInSeconds);
	printf("filter frequency: %f\n", filterFrequency);
	
	BMMultiLevelBiquad This;
	size_t numLevels = 4;
	BMMultiLevelBiquad_init(&This, numLevels, sampleRate, false, true, false);
	BMMultiLevelBiquad_setCriticallyDampedLP(&This, filterFrequency, 0, numLevels);
	
	
	float frequencies [1024];
	float magnitudes [1024];
	
	for(size_t i=0; i < 1024; i++){
		float maxF = 300.0f;
		float alphaMax = log2f(maxF);
		float alpha = ((float)i / 1024.0f) * alphaMax;
		frequencies[i] = powf(2.0f,alpha);
	}
	
	BMMultiLevelBiquad_tfMagVector(&This, frequencies, magnitudes, 1024);
	
	for(size_t i=0; i < 1024; i++){
		magnitudes[i] = BM_GAIN_TO_DB(magnitudes[i]);
	}
	
	arrayToFile(magnitudes, 1024);
}


void testEnvelopeFollower(){
	BMEnvelopeFollower env;
	
	BMEnvelopeFollower_init(&env, 48000.0f);
	
	// create sine wave for input test
	float frequency = 80.0f;
	float t [BMCOMPRESSOR_TEST_LENGTH];
	float o [BMCOMPRESSOR_TEST_LENGTH];
	float sr = BMCOMPRESSOR_TEST_SAMPLERATE;
	for (size_t i=0; i<BMCOMPRESSOR_TEST_LENGTH; i++) {
		t[i] = sinf(frequency * 2.0 * M_PI * (float)i / sr);
		
		if (i>8000 && i < 20000) t[i] *= BM_DB_TO_GAIN(-20.0f);
	}
	
	// rectify the sine wave
	vDSP_vabs(t, 1, t, 1, BMCOMPRESSOR_TEST_LENGTH);
	
	BMEnvelopeFollower_processBuffer(&env, t, o, BMCOMPRESSOR_TEST_LENGTH);
	
	// print the output
	arrayToFile(o,BMCOMPRESSOR_TEST_LENGTH);
}

void testEnvAttackTime(){
	float attackTimeInSeconds = 0.100;
	float sampleRate = 48000.0f;
	
	// float correction = 1.0f / 1.0f; // one stage
	// float correction = 1.0f / 1.81f; // three stage
	// float correction = 1.0f / 1.35375f; // two stage
	// float correction = 1.0f / 0.85250f; // one stage
	
	BMEnvelopeFollower env;
	BMEnvelopeFollower_init(&env, sampleRate);
	BMEnvelopeFollower_setAttackTime(&env, attackTimeInSeconds);
	
	float t [BMCOMPRESSOR_TEST_LENGTH];
	float o [BMCOMPRESSOR_TEST_LENGTH];
	
	float one = 1.0f;
	vDSP_vfill(&one, t, 1, BMCOMPRESSOR_TEST_LENGTH);
	
	BMEnvelopeFollower_processBuffer(&env, t, o, BMCOMPRESSOR_TEST_LENGTH);
	
	bool done = false;
	for(size_t i=0; i<BMCOMPRESSOR_TEST_LENGTH && !done; i++){
		if (o[i] > 0.97){
			done = true;
			printf("reached 97 percent at time: %f\n", (float)i / sampleRate);
		}
	}
	
	// print the output
	arrayToFile(o,BMCOMPRESSOR_TEST_LENGTH);
}

void testEnvReleaseTime(){
	float attackTimeInSeconds = 0.001;
	float releaseTimeInSeconds = 0.050;
	float sampleRate = 48000.0f;
	
	BMEnvelopeFollower env;
	BMEnvelopeFollower_init(&env, sampleRate);
	BMEnvelopeFollower_setAttackTime(&env, attackTimeInSeconds);
	BMEnvelopeFollower_setReleaseTime(&env, releaseTimeInSeconds);
	
	float t [2*BMCOMPRESSOR_TEST_LENGTH];
	float o [2*BMCOMPRESSOR_TEST_LENGTH];
	
	float one = 1.0f; float zero = 0.0f;
	vDSP_vfill(&one, t, 1, BMCOMPRESSOR_TEST_LENGTH);
	vDSP_vfill(&zero, t+BMCOMPRESSOR_TEST_LENGTH,1,BMCOMPRESSOR_TEST_LENGTH);
	
	BMEnvelopeFollower_processBuffer(&env, t, o, 2*BMCOMPRESSOR_TEST_LENGTH);
	
	bool done = false;
	for(size_t i=BMCOMPRESSOR_TEST_LENGTH; i<2*BMCOMPRESSOR_TEST_LENGTH && !done; i++){
		if (o[i] < 0.0309){
			done = true;
			printf("reached -97 percent at time: %f\n", (float)(i-BMCOMPRESSOR_TEST_LENGTH) / sampleRate);
		}
	}
	
	// print the output
	arrayToFile(o,2*BMCOMPRESSOR_TEST_LENGTH);
}



void logSpeedTest(){
	// fill an array with random numbers
	float t [BMCOMPRESSOR_TEST_LENGTH];
	float s [BMCOMPRESSOR_TEST_LENGTH];
	for(size_t i=0; i<BMCOMPRESSOR_TEST_LENGTH; i++){
		t[i] = (double)arc4random() / (double) UINT32_MAX;
		s[i] = (double)arc4random() / (double) UINT32_MAX;
	}
	
	// test timing
	size_t numTests = 10000;
	clock_t start, end;
	clock_t vdspDB,vfdbtogain,vmul;
	vdspDB = vfdbtogain = vmul = 0;
	float one = 1.0f;
	for(size_t i =0; i<numTests; i++){
		start = clock();
		vDSP_vdbcon(t, 1, &one, t, 1, BMCOMPRESSOR_TEST_LENGTH, 1);
		end = clock();
		vdspDB += end - start;
		
		start = clock();
		vector_fastDbToGain(t, t, BMCOMPRESSOR_TEST_LENGTH);
		end = clock();
		vfdbtogain += end - start;
		
		start = clock();
		vDSP_vmul(t,1, s, 1, t, 1, BMCOMPRESSOR_TEST_LENGTH);
		end = clock();
		vmul += end - start;
	}
	
	printf("linearToDB: %f, dbToLinear: %f, mul: %f\n",(double)vdspDB/CLOCKS_PER_SEC,(double)vfdbtogain/CLOCKS_PER_SEC,(double)vmul/CLOCKS_PER_SEC);
}

void quadraticThresholdT(const float* input, float* output, float threshold, float softKneeWidth, size_t numFrames){
	assert(input != output);
	
	float l = threshold;
	float w = softKneeWidth;
	
	/*
	 * we want to express (- l x + (l + w + x)^2/4)/w
	 * as a polynomial A*x^2 * Bx + C
	 *
	 * Mathematica: (l^2 + 2 l w + w^2 + (-2 l + 2 w) x + x^2) / (4 w)
	 */
	float lPlusW = l + w;
	float lMinusW = l - w;
	float fourWinv = 1.0 / (4.0f * w);
	float C = lPlusW * lPlusW * fourWinv;
	float B = 2.0f * (-lMinusW) * fourWinv;
	float A = fourWinv;
	float coefficients [3] = {A, B, C};
	
	// clip values to [l-w,l+w];
	vDSP_vclip(input, 1, &lMinusW, &lPlusW, output, 1, numFrames);
	
	// process the polynomial
	vDSP_vpoly(coefficients, 1, output, 1, output, 1, numFrames, 2);
	
	// where the input is greater than the polynomial output, return the input
	vDSP_vmax(input,1,output,1,output,1,numFrames);
}

void quadraticThresholdTest(){
	float t [2*BMCOMPRESSOR_TEST_LENGTH];
	float o [2*BMCOMPRESSOR_TEST_LENGTH];
	
	float lengthInv = 1.0f / (float)BMCOMPRESSOR_TEST_LENGTH;
	for(size_t i=0; i<BMCOMPRESSOR_TEST_LENGTH; i++)
		t[i] = -10.0f + (10.0f * ((float)i * lengthInv));
	
	quadraticThresholdT(t,o,-4.0f,1.0f,BMCOMPRESSOR_TEST_LENGTH);
	
	// print the output
	arrayToFile(o,BMCOMPRESSOR_TEST_LENGTH);
}

//#define LT_SampleRate 48000.
//#define LT_Frequency 10.
//#define LT_TestRange 50
//void testBMLagTime(){
//	BMStereoLagTime stereoLagTime;
//	BMStereoLagTime_init(&stereoLagTime, 500, LT_SampleRate);
//
//
//	size_t length = LT_SampleRate *4;
//	float* audioX = malloc(sizeof(float)*length);
//	float* inLAudioY = malloc(sizeof(float)*length);
//	float* outLAudioY = malloc(sizeof(float)*length);
//	float* inRAudioY = malloc(sizeof(float)*length);
//	float* outRAudioY = malloc(sizeof(float)*length);
//
//	float scale = 1.;
//	float xFactor = 20;
//	for(int i=0;i<length;i++){
//		audioX[i] = (i*xFactor)/length;
//		inLAudioY[i] = sinf(LT_Frequency * 2.0 * M_PI * (float)i / LT_SampleRate) * scale;
//		inRAudioY[i] = sinf(LT_Frequency * 2.0 * M_PI * (float)i / LT_SampleRate) * scale;
//	}
//
//	float delayTime[LT_TestRange];
//	for(int i=0;i<30;i++){
//		delayTime[i] = (LT_TestRange - i)*1;
//	}
//
//	for(int i=10;i<LT_TestRange;i++){
//		delayTime[i] = i*2;
//	}
//	size_t processedSample = 0;
//	size_t processingSample = 1024;
//	for(int i=0;i<LT_TestRange;i++){
//		BMStereoLagTime_setDelayLeft(&stereoLagTime, delayTime[i]);
//		BMStereoLagTime_process(&stereoLagTime, inLAudioY+processedSample, inRAudioY+processedSample, outLAudioY+processedSample, outRAudioY+processedSample, processingSample);
//		processedSample += processingSample;
//	}
//
//	arrayXYToFile(audioX,outLAudioY, length);
//}

//#define BS_Range 10
//#define BS_Freq_Range 1000.0f
//#define BS_Freq_Start 1.0f
//#define BS_Freq_End 20000.0f
//void testBinauralSynthesis(){
//	float* input = malloc(sizeof(float) * BS_Range);
//	float* output = malloc(sizeof(float) * BS_Range);
//	float* frequency = malloc(sizeof(float)*BS_Freq_Range);
//	float* magnitude = malloc(sizeof(float)*BS_Freq_Range);
//	
//	float alpha = powf(BS_Freq_End/BS_Freq_Start, 1.0f/BS_Freq_Range);
//	
//	float* xArray = malloc(sizeof(float)*BS_Freq_Range);
//	
//	for(int i = 0;i< BS_Freq_Range;i++){
//		frequency[i] = BS_Freq_Start * powf(alpha, i);
//		xArray[i] = i * 10;
//	}
//	
//	
//	
//	BMBinauralSynthesis bs;
//	BMBinauralSynthesis_init(&bs, false, LT_SampleRate);
//	
//	BMBinauralSynthesis_setAngleLeft(&bs, 0);
//	
//	BMBinauralSynthesis_processMono(&bs, input, output, BS_Range);
//	
//	BMBinauralSynthesis_getTFMagVectorData(&bs,frequency, magnitude, BS_Freq_Range, 0);
//	//Convert from gain to db
//	for(int i=0;i<BS_Freq_Range;i++){
//		magnitude[i] = BM_GAIN_TO_DB(magnitude[i]) * 1000;
//	}
//	
//	arrayXYToFile(xArray,magnitude, BS_Freq_Range);
//	
//}



void testWetDryMixer() {
	size_t bufferSize = 600;
	size_t testLength = bufferSize * round(48000.0 / bufferSize);
	float* inputDryR = malloc(sizeof(float)*testLength);
	float* inputDryL = malloc(sizeof(float)*testLength);
	float* inputWetR = malloc(sizeof(float)*testLength);
	float* inputWetL = malloc(sizeof(float)*testLength);
	float* outputR = malloc(sizeof(float)*testLength);
	float* outputL = malloc(sizeof(float)*testLength);
	
	float zero = 0.0f;
	float fourth = 0.25f;
	float threeFourths = 0.75f;
	float one = 1.0f;
	
	vDSP_vfill(&zero,inputDryL, 1, testLength);
	vDSP_vfill(&fourth,inputDryR, 1, testLength);
	vDSP_vfill(&threeFourths, inputWetL, 1, testLength);
	vDSP_vfill(&one, inputWetR, 1, testLength);
	
	BMWetDryMixer m;
	BMWetDryMixer_init(&m, 48000);
	BMWetDryMixer_setMix(&m, 0.0f);
	
	size_t i;
	for(i=0; i<testLength/4; i+= bufferSize){
		BMWetDryMixer_processBufferInPhase(&m,
										   inputWetL+i, inputWetR+i,
										   inputDryL+i, inputDryR+i,
										   outputL+i, outputR+i,
										   bufferSize);
	}
	
	BMWetDryMixer_setMix(&m, 0.5f);
	
	for(; i<testLength; i+= bufferSize){
		BMWetDryMixer_processBufferInPhase(&m,
										   inputWetL+i, inputWetR+i,
										   inputDryL+i, inputDryR+i,
										   outputL+i, outputR+i,
										   bufferSize);
	}
	
	// fill an array of x coordinates for plotting the output
	float* xArray = malloc(sizeof(float)*testLength);
	vDSP_vramp(&zero,&one,xArray,1,testLength);
	
	// plot one channel of output
	arrayXYToFile(xArray, outputL, testLength);
	
	free(inputDryR);
	free(inputDryL);
	free(inputWetR);
	free(inputWetL);
	free(outputR);
	free(outputL);
}





void testUpsampler2x(){
	BMIIRUpsampler2x us;
	size_t numCoefficients = BMIIRUpsampler2x_init(&us, 110.0, 0.025, false);
	printf("\nnumber of allpass coefficients: %zu\n",numCoefficients);
	float sampleRate = 48000.0f;
	size_t testLength = sampleRate * 5;
	
	float* sineSweep = malloc(sizeof(float)*testLength);
	float* output = malloc(sizeof(float)*testLength*2);
	
	generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
	
	BMIIRUpsampler2x_processBufferMono(&us, sineSweep, output, testLength);
	
	arrayToFile(output,2*testLength);
	
	BMIIRUpsampler2x_free(&us);
	free(sineSweep);
	free(output);
}




void testDownsampler2x(){
	BMIIRDownsampler2x ds;
	size_t numCoefficients = BMIIRDownsampler2x_init(&ds, 110.0, 0.025, false);
	printf("\nnumber of allpass coefficients: %zu\n",numCoefficients);
	float sampleRate = 48000.0f;
	size_t testLength = sampleRate * 10;
	
	float* sineSweep = malloc(sizeof(float)*testLength);
	float* output = malloc(sizeof(float)*testLength/2);
	
	generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
	
	BMIIRDownsampler2x_processBufferMono(&ds, sineSweep, output, testLength);
	
	arrayToFile(output,testLength/2);
	
	BMIIRDownsampler2x_free(&ds);
	free(sineSweep);
	free(output);
}






void testUpDownsampler2x(){
	BMIIRUpsampler2x us;
	BMIIRDownsampler2x ds;
	BMIIRUpsampler2x_init(&us, 110.0, 0.025, false);
	size_t numCoefficients = BMIIRDownsampler2x_init(&ds, 110.0, 0.025, false);
	printf("\nnumber of allpass coefficients: %zu\n",numCoefficients);
	float sampleRate = 48000.0f;
	size_t testLength = sampleRate * 5;
	
	float* sineSweep = malloc(sizeof(float)*testLength);
	float* upsampled = malloc(sizeof(float)*testLength*2);
	float* output = malloc(sizeof(float)*testLength);
	
	generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
	
	size_t samplesRemaining = testLength;
	float* inputP = sineSweep;
	float* outputP = output;
	while(samplesRemaining > 0){
		size_t samplesProcessing = BM_MIN(samplesRemaining, 80);
		
		// upsample
		BMIIRUpsampler2x_processBufferMono(&us, inputP, upsampled, samplesProcessing);
		
		// waveshape
		//BMWaveshaper_processBufferBidirectional(upsampled, upsampled, testLength*2);
		
		// downsample
		BMIIRDownsampler2x_processBufferMono(&ds, upsampled, outputP, samplesProcessing*2);
		outputP += samplesProcessing;
		inputP += samplesProcessing;
		samplesRemaining -= samplesProcessing;
	}
	
	char* filename = "./sineSweepResponse.csv";
	arrayToFileWithName(output, filename, testLength);
	
	BMIIRUpsampler2x_free(&us);
	BMIIRDownsampler2x_free(&ds);
	free(sineSweep);
	free(upsampled);
	free(output);
}



void testUpsampler(){
	size_t upsampleFactor = 8;
	
	BMUpsampler us;
	BMUpsampler_init(&us, false, upsampleFactor, BMRESAMPLER_FULL_SPECTRUM);
	
	float sampleRate = 48000.0f;
	size_t testLength = sampleRate * 5 / upsampleFactor;
	
	float* sineSweep = malloc(sizeof(float)*testLength);
	float* output = malloc(sizeof(float)*testLength*upsampleFactor);
	
	generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
	
	BMUpsampler_processBufferMono(&us, sineSweep, output, testLength);
	
	arrayToFile(output,upsampleFactor*testLength);
	
	BMUpsampler_free(&us);
	free(sineSweep);
	free(output);
}




void testDownsampler(void){
	size_t downsampleFactor = 2;
	
	BMDownsampler ds;
	BMDownsampler_init(&ds, true, downsampleFactor, BMRESAMPLER_FULL_SPECTRUM);
	
	float sampleRate = 48000.0f;
	size_t testLength = sampleRate * 5 * downsampleFactor;
	
	float* sineSweepL = malloc(sizeof(float)*testLength);
	float* outputL = malloc(sizeof(float)*testLength/downsampleFactor);
	float* sineSweepR = malloc(sizeof(float)*testLength);
	float* outputR = malloc(sizeof(float)*testLength/downsampleFactor);
	
	generateSineSweep(sineSweepL, 20.0f, 24000.0f, 48000.0f, testLength);
	generateSineSweep(sineSweepR, 20.0f, 24000.0f, 48000.0f, testLength);
	
	size_t ic = 0;
	size_t sp = 0;
	size_t oc = 0;
	while(ic < testLength - 1){
		size_t samplesProcessing = BM_MIN(sp,testLength - (ic+1));
		size_t samplesOut;
		BMDownsampler_processBufferStereoOddInputLength(&ds, sineSweepL + ic, sineSweepR + ic, outputL + oc, outputR + oc, samplesProcessing, &samplesOut);
		ic += samplesProcessing;
		oc += samplesOut;
		sp++;
	}
	
	arrayToFile(outputL,testLength/downsampleFactor);
	
	BMDownsampler_free(&ds);
	free(sineSweepL);
	free(outputL);
	free(sineSweepR);
	free(outputR);
}


void testStaticDelay(void){
	// configure delay settings
	float sampleRate = 48000.0f;
	float delayTimeSeconds = 0.333f;
	float decayTimeSeconds = 1.0f;
	float wetGainDb = -6.0f;
	float crossMixAmount = 0.3f;
	float lowpassFc = 2000.0f;
	float highpassFc = 300.0f;
	
	BMStaticDelay d;
	BMStaticDelay_init(&d, delayTimeSeconds, decayTimeSeconds, wetGainDb, crossMixAmount, lowpassFc, highpassFc, sampleRate);
	BMStaticDelay_setDryGain(&d, -100.0f);
	
	// allocate memory for 5 seconds sine sweep test
	size_t testLength = sampleRate * 5;
	float* sineSweepL = malloc(sizeof(float)*testLength);
	float* outputL = malloc(sizeof(float)*testLength);
	float* sineSweepR = malloc(sizeof(float)*testLength);
	float* outputR = malloc(sizeof(float)*testLength);
	
	generateSineSweep(sineSweepL, 20.0f, 24000.0f, 48000.0f, testLength);
	generateSineSweep(sineSweepR, 20.0f, 24000.0f, 48000.0f, testLength);
	
	// reduce the volume
	float gain = BM_DB_TO_GAIN(-10.0f);
	vDSP_vsmul(sineSweepL, 1, &gain, sineSweepL, 1, testLength);
	vDSP_vsmul(sineSweepR, 1, &gain, sineSweepR, 1, testLength);
	
	size_t ic = 0;
	size_t sp = 1000;
	size_t oc = 0;
	while(ic < testLength - 1){
		size_t samplesProcessing = BM_MIN(sp,testLength - (ic+1));
		BMStaticDelay_processBufferStereo(&d, sineSweepL + ic, sineSweepR + ic, outputL + oc, outputR + oc, (int)samplesProcessing);
		ic += samplesProcessing;
		oc += samplesProcessing;
		sp++;
	}
	
	arrayToFile(outputL,testLength);
	
	// free memory
	BMStaticDelay_free(&d);
	free(sineSweepL);
	free(outputL);
	free(sineSweepR);
	free(outputR);
}



void testUpDownsampler(){
	size_t upsampleFactor = 16;
	
	BMDownsampler ds;
	BMUpsampler us;
	BMUpsampler_init(&us, false, upsampleFactor, BMRESAMPLER_FULL_SPECTRUM);
	BMDownsampler_init(&ds, false, upsampleFactor, BMRESAMPLER_FULL_SPECTRUM);
	
	float sampleRate = 48000.0f;
	size_t testLength = sampleRate * 5;
	
	float* sineSweep = malloc(sizeof(float)*testLength);
	float* upsampled = malloc(sizeof(float)*testLength*upsampleFactor);
	float* output = malloc(sizeof(float)*testLength);
	
	generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
	
	// gain boost
	//float boost = BM_DB_TO_GAIN(100.0f);
	//vDSP_vsmul(sineSweep,1,&boost,sineSweep,1,testLength);
	
	size_t samplesRemaining = testLength;
	float *inputP = sineSweep;
	float *outputP = output;
	while(samplesRemaining > 0){
		size_t samplesProcessing = BM_MIN(samplesRemaining, 80);
		if(upsampleFactor > 1)
			BMUpsampler_processBufferMono(&us, inputP, upsampled, samplesProcessing);
		else
			memcpy(upsampled,inputP,sizeof(float)*samplesProcessing);
		
		// waveshape
		//int length = (int)testLength*(int)upsampleFactor;
		//vvtanhf(upsampled, upsampled, &length);
		
		if(upsampleFactor > 1)
			BMDownsampler_processBufferMono(&ds, upsampled, outputP, samplesProcessing*upsampleFactor);
		else
			memcpy(outputP,upsampled,sizeof(float)*samplesProcessing);
		
		inputP += samplesProcessing;
		outputP += samplesProcessing;
		samplesRemaining -= samplesProcessing;
	}
	
	char* filename = "./sineSweepResponse.csv";
	arrayToFileWithName(output, filename, testLength);
	
	BMDownsampler_free(&ds);
	BMUpsampler_free(&us);
	free(sineSweep);
	free(upsampled);
	free(output);
}




void testCrossover(){
	
	BMCrossover4way cvr;
	BMCrossover4way_init(&cvr,
						 300.0f, 1600.0f, 8000.0f,
						 48000.0f,
						 false,
						 true);
	
	float sampleRate = 48000.0f;
	size_t testLength = sampleRate * 5;
	
	float* sineSweep = malloc(sizeof(float)*testLength);
	float* outputLLow = malloc(sizeof(float)*testLength*8);
	float* outputRLow = outputLLow + testLength;
	float* outputLHigh = outputRLow + testLength;
	float* outputRHigh = outputLHigh + testLength;
	float* outputLMid = outputRHigh + testLength;
	float* outputRMid = outputLMid + testLength;
	float* outputLHighMid = outputRMid + testLength;
	float* outputRHighMid = outputLHighMid + testLength;
	
	generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
	
	//	BMCrossover_processStereo(&cvr,
	//	sineSweep, sineSweep,
	//	outputLLow, outputRLow,
	//	outputLHigh, outputRHigh,
	//	testLength);
	
	BMCrossover4way_processStereo(&cvr,
								  sineSweep, sineSweep,
								  outputLLow, outputRLow,
								  outputLMid, outputRMid,
								  outputLHighMid, outputRHighMid,
								  outputLHigh, outputRHigh,
								  testLength);
	
	BMCrossover4way_recombine(outputLLow, outputRLow,
							  outputLMid, outputRMid,
							  outputLHighMid, outputRHighMid,
							  outputLHigh, outputRHigh,
							  outputLLow, outputRLow,
							  testLength);
	
	//	vDSP_vadd(outputLLow,1,outputLHigh,1,outputLLow,1,testLength);
	//	vDSP_vadd(outputRLow,1,outputRHigh,1,outputRLow,1,testLength);
	
	char* filename = "./sineSweepResponse.csv";
	arrayToFileWithName(outputLLow, filename, testLength);
	
	BMCrossover4way_free(&cvr);
	free(sineSweep);
	free(outputLLow);
}




void testVASVF(){
	float sampleRate = 48000;
	float testLength = 5000;
	
	float* input = malloc(sizeof(float)*testLength);
	float* output = malloc(sizeof(float)*testLength);
	input[0] = 1;
	for(int i =1;i<testLength;i++){
		input[1] = 0;
	}
	
	BMVAStateVariableFilter filter;
	BMVAStateVariableFilter_init(&filter, true, sampleRate, SVFLowpass);
	
	BMVAStateVariableFilter_setFilter(&filter, SVFBandpass, 1000, 1.0f/sqrtf(2.0f), 0);
	
	BMVAStateVariableFilter_processBufferStereo(&filter, input, input, output, output, testLength);
	
	for(int i =0;i<testLength;i++){
		output[i] = output[i]*100;
	}
	
	arrayToFile(output, testLength);
}



void fadeOut(const float* input, float* output, float rt60DecaySeconds, float sampleRate, size_t numSamples){
	double decayTimeSamples = rt60DecaySeconds * sampleRate;
	double commonRatio = pow(BM_DB_TO_GAIN(-60.0),1.0/decayTimeSamples);
	double gain = 1.0;
	for(size_t i=0; i<numSamples; i++){
		output[i] = input[i] * gain;
		gain *= commonRatio;
	}
}




void fadeIn(const float* input, float* output, float fadeTimeSeconds, float sampleRate, size_t numSamples){
	double fadeTimeSamples = fadeTimeSeconds * sampleRate;
	double commonRatio = pow(BM_DB_TO_GAIN(100.0),1.0/fadeTimeSamples);
	double gain = BM_DB_TO_GAIN(-100.0);
	for(size_t i=0; i<numSamples; i++){
		output[i] = input[i] * gain;
		if(gain <= 1.0)
			gain *= commonRatio;
	}
}



void testNoiseGate() {
	float sampleRate = 48000;
	size_t testLength = sampleRate * 5;
	float decayTime = 1.0 / 10.0;
	float threshold = -10.0;
	float toneFrequency = 800.0;
	BMNoiseGate ng1;
	
	BMNoiseGate_init(&ng1, threshold, sampleRate);
	
	// generate a 20 Hz tone for a test signal
	float* testSignal = malloc(sizeof(float)*testLength);
	generateSineSweep(testSignal, toneFrequency, toneFrequency, sampleRate, testLength);
	
	//fade the volume of the test signal so that it goes below the threshold
	float rt60DecayTime = 10.0f;
	fadeOut(testSignal, testSignal, rt60DecayTime, sampleRate, testLength/2);
	
	// fade back in
	float fadeInTime = 1.0f;
	fadeIn(testSignal + testLength/2, testSignal + testLength/2, fadeInTime, sampleRate, testLength/2);
	
	BMNoiseGate_processMono(&ng1, testSignal, testSignal, testLength);
	
	arrayToFile(testSignal, testLength);
}





void testOversamplerTransientResponse(){
	float sampleRate = 48000;
	size_t testLength = sampleRate * 5;
	float toneFrequency = 80.0;
	
	
	// generate a toneFrequency Hz tone for a test signal
	float* testSignal = malloc(sizeof(float)*testLength);
	generateSineSweep(testSignal, toneFrequency, toneFrequency, sampleRate, testLength);
	
	// zero out some sections of the signal
	for(size_t i=0; i<testLength; i++){
		if((int)ceil((float)i/sampleRate) % 2 == 0)
			testSignal[i]=0;
	}
	
	size_t upsampleFactor = 32;
	
	BMDownsampler ds;
	BMUpsampler us;
	BMUpsampler_init(&us, false, upsampleFactor, BMRESAMPLER_FULL_SPECTRUM);
	BMDownsampler_init(&ds, false, upsampleFactor, BMRESAMPLER_FULL_SPECTRUM);
	
	float* upsampled = malloc(sizeof(float)*testLength*upsampleFactor);
	float* output = malloc(sizeof(float)*testLength);
	
	//    // gain boost
	//    float boost = BM_DB_TO_GAIN(100.0f);
	//    vDSP_vsmul(sineSweep,1,&boost,sineSweep,1,testLength);
	
	if(upsampleFactor > 1)
		BMUpsampler_processBufferMono(&us, testSignal, upsampled, testLength);
	else
		memcpy(upsampled,testSignal,sizeof(float)*testLength);
	
	//    // waveshape
	//    int length = (int)testLength*(int)upsampleFactor;
	//    vvtanhf(upsampled, upsampled, &length);
	
	if(upsampleFactor > 1)
		BMDownsampler_processBufferMono(&ds, upsampled, output, testLength*upsampleFactor);
	else
		memcpy(output,upsampled,sizeof(float)*testLength);
	
	
	arrayToFile(output, testLength);
	
	BMDownsampler_free(&ds);
	BMUpsampler_free(&us);
	free(testSignal);
	free(upsampled);
	free(output);
}





void testOversamplerImpulseResponse(){
	
	float sampleRate = 48000;
	size_t testLength = 512;
	
	
	// generate the impulse input
	float* testSignal = malloc(sizeof(float)*testLength);
	float* emptySignal = malloc(sizeof(float)*testLength);
	memset(emptySignal,0,sizeof(float)*testLength);
	memset(testSignal,0,sizeof(float)*testLength);
	testSignal[0] = 1.0f;
	
	size_t upsampleFactor = 16;
	
	BMDownsampler ds;
	BMUpsampler us;
	BMUpsampler_init(&us, false, upsampleFactor, BMRESAMPLER_FULL_SPECTRUM);
	BMDownsampler_init(&ds, false, upsampleFactor, BMRESAMPLER_FULL_SPECTRUM);
	
	float* upsampled = malloc(sizeof(float)*testLength*upsampleFactor);
	float* output = malloc(sizeof(float)*testLength);
	
	if(upsampleFactor > 1)
		BMUpsampler_processBufferMono(&us, testSignal, upsampled, testLength);
	else
		memcpy(upsampled,testSignal,sizeof(float)*testLength);
	
	
	if(upsampleFactor > 1)
		BMDownsampler_processBufferMono(&ds, upsampled, output, testLength*upsampleFactor);
	else
		memcpy(output,upsampled,sizeof(float)*testLength);
	
	
	char* filename = "./impulseResponse.csv";
	arrayToFileWithName(output, filename, testLength);
	
	BMDownsampler_free(&ds);
	BMUpsampler_free(&us);
	free(testSignal);
	free(upsampled);
	free(output);
}



void SFMtoWNormalised(float* SFMwithoutZeroes, double* W, int values){
	double normalisationConstant = sqrt((double)FFTSIZE/2)/sqrt(pow(M_PI,2)-1);
	//meanSFM has to be e^(-eulergamma) to be white
	for (int i = 0; i<values; i++){
		W[i] = ((double)-1 * log((double)SFMwithoutZeroes[i]) - 0.57721566490153286060651209008240243104215933593992)*normalisationConstant;
		//        printf("original SFM : %f, W: %f \n", SFMwithoutZeroes[i], W[i]);
	}
}


void getZValue(double* data, int values, double* z, int sampleSize){
	
	//find mean of data
	double sample_mean = 0;
	
	for (int i = 0; i<values; i++){
		sample_mean += data[i]/(double)values;
	}
	
	double mu_0 = 0.0;
	double pop_var = 1;
	z[0] = (sample_mean - mu_0)/(pop_var/sqrt(values));
	
	//    printf("Sample mean is: %f \n", sample_mean);
	//    printf("Z-value is : %lf \n", z[0]);
	
}


void hypothesisTestAlphaBeta(float* data, int values, bool* result){
	float alpha = FFTSIZE*0.339989 + 0.495;
	float beta = FFTSIZE*0.265556 - 0.602;
	float mean_data = 0;
	
	for (int i = 0; i<values; i++){
		mean_data += data[i]/(float)values;
	}
	
	float p_value = BM_betaCDF(mean_data, alpha, beta);
	
	// two-tailed at 10% significance level
	if (p_value <= 0.05f || p_value >= 0.95f){
		result[0] = false; //fail the test
	}
	else{
		result[0] = true; //doesn't fail the test
	}
}


//p-value is 0.05
void hypothesisTestNormal(float* data, int values, bool* result, int sampleSize){
	
	//convert all SFM data into W
	double* W = malloc(sizeof(double)*values);
	SFMtoWNormalised(data, W, values);
	double z = 0;
	getZValue(W, values, &z, sampleSize);
	
	//two-tailed test, at 10% significance level
	if (z <= -1.645 || z>= 1.645){
		result[0] = false; //reject
		printf("H0 is rejected, z value is: %f \n", z);
	}
	else{
		result[0] = true; //does not reject
		//        printf("H0 is not rejected \n");
	}
	
	free(W);
	
}


void SFMStatsPerIR(float* IR, size_t irLength, float* SFMValue){
	//set up the FFT
	BMFFT fft;
	size_t fftSize = FFTSIZE;
	size_t outputLengthfft = fftSize / 2;
	
	BMFFT_init(&fft, fftSize);
	
	// compute the FFT on windows from the impulse response.
	size_t numWindows = irLength / fftSize;
	size_t numWindowsToSkip = 1; // skip the beginning of the IR
	
	// prepare arrays for result
	float* absFFTResult = malloc(sizeof(float)*outputLengthfft);
	float* tempFFTResult = malloc(sizeof(float)*fftSize);
    float* windowedFFTResult = malloc(sizeof(float)*fftSize);
	
	memset(absFFTResult, 0, outputLengthfft);
	
	for(size_t i=numWindowsToSkip; i<numWindows; i ++){
        //Do hamming window
        BMFFT_hammingWindow(&fft, IR + (fftSize*i), windowedFFTResult, fftSize);
		//Take FFT per window
		BMFFT_absFFTCombinedDCNQ(&fft, windowedFFTResult, tempFFTResult, fftSize);
		//store it to absFFTResult
		for (int j = 0; j<outputLengthfft; j++){
			absFFTResult[j] += tempFFTResult[j] / (float) outputLengthfft;
		}
	}
	
	
	// square the result
	vDSP_vsq(absFFTResult, 1, absFFTResult, 1, outputLengthfft);
	
	// find the sum of the logs of all elements
	float logSum = 0.0f;
	for(size_t i=0; i<outputLengthfft; i++){
		logSum += log2f(absFFTResult[i]);
		if(absFFTResult[i]==0)
			printf("zero");
	}
	
	float geometricMean = powf(2.0f,logSum/(float)outputLengthfft);
	
	float arithmeticMean;
	vDSP_meanv(absFFTResult, 1, &arithmeticMean, outputLengthfft);
	
	SFMValue[0] = geometricMean / arithmeticMean;
    
    free(windowedFFTResult);
    free(tempFFTResult);
    free(absFFTResult);
}


void SFMStatsMultiplePerIR(float* IR, size_t irLength, float* SFMValue, float* stats, bool write){
    // set up the SFM
    BMSFM sfm;
    
    size_t fftSize = FFTSIZE;
    
    BMSFM_init(&sfm, fftSize);
    
    // compute the SFM on windows from the impulse response.
    size_t numWindows = irLength / fftSize;
    size_t numWindowsToSkip = 0; // skip the beginning of the IR
    
    float* SFMResults = malloc(sizeof(float)*(numWindows-numWindowsToSkip));
    float* SFMResultsWithoutZeroes = malloc(sizeof(float)*(numWindows-numWindowsToSkip));
    
    for(size_t i=numWindowsToSkip; i<numWindows; i ++){
        //take SFM per window
        SFMResults[i-numWindowsToSkip] = BMSFM_process(&sfm, IR + (fftSize * i), fftSize);
        SFMValue[i-numWindowsToSkip] = SFMResults[i-numWindowsToSkip];
    }
    
    // calculate the min, max and mean
    float meanSFM = 0.0f;
    float maxSFM = FLT_MIN;
    float minSFM = FLT_MAX;
    float sumSFM = 0;
    int nonZeroSFM = 0;
    
    for(size_t i=0; i<numWindows-numWindowsToSkip; i ++){
        //skip if 0
        if (SFMResults[i] > 0.00000001){
            sumSFM += SFMResults[i];
            if(SFMResults[i] > maxSFM) maxSFM = SFMResults[i]; // max
            if(SFMResults[i] < minSFM) minSFM = SFMResults[i]; // min
            SFMResultsWithoutZeroes[nonZeroSFM] = SFMResults[i];
            nonZeroSFM++;
        }
    }
    
    meanSFM = sumSFM / (float)(nonZeroSFM);
    
    // calulate the std. dev.
    float variance = 0.0f;
    for(size_t i=0; i<numWindows-numWindowsToSkip; i ++){
        //skip if 0
        if (SFMResults[i] > 0.00000001){
            float diff = SFMResults[i] - meanSFM;
            variance += diff*diff;
        }
    }
    variance /= (float)(nonZeroSFM-1);
    float stdDev = sqrtf(variance);
    
//    printf("*** SFM Stats ***\n");
//    printf("mean: %f\n", meanSFM);
//    printf("std. dev.: %f\n", stdDev);
//    printf("min: %f\n", minSFM);
//    printf("max: %f\n", maxSFM);
    
    stats[0] = meanSFM;
    stats[1] = stdDev;
    stats[2] = minSFM;
    stats[3] = maxSFM;
    
    if (write){
        // output the impulse response to a file
        char *filename = malloc(sizeof(char)*128);
        sprintf(filename, "./SFMValuesPerIR.csv");

        arrayToFileWithName(SFMResultsWithoutZeroes, filename, nonZeroSFM);
    }
 
    free(SFMResultsWithoutZeroes);
    free(SFMResults);
}

//void SFMStatsMultiplePerIR(float* IR, size_t irLength, float* SFMValue, float* stats){
//	// set up the SFM
//	BMSFM sfm;
//	
//	size_t fftSize = FFTSIZE;
//	
//	BMSFM_init(&sfm, fftSize);
//	
//	// compute the SFM on windows from the impulse response.
//	size_t numWindows = irLength / fftSize;
//	size_t numWindowsToSkip = 0; // skip the beginning of the IR
//	
//	float* SFMResults = malloc(sizeof(float)*(numWindows-numWindowsToSkip));
//	float* SFMResultsWithoutZeroes = malloc(sizeof(float)*(numWindows-numWindowsToSkip));
//	
//	for(size_t i=numWindowsToSkip; i<numWindows; i ++){
//		//take SFM per window
//		SFMResults[i-numWindowsToSkip] = BMSFM_process(&sfm, IR + (fftSize * i), fftSize);
//		SFMValue[i-numWindowsToSkip] = SFMResults[i-numWindowsToSkip];
//	}
//	
//	// calculate the min, max and mean
//	float meanSFM = 0.0f;
//	float maxSFM = FLT_MIN;
//	float minSFM = FLT_MAX;
//	float sumSFM = 0;
//	int nonZeroSFM = 0;
//	
//	for(size_t i=0; i<numWindows-numWindowsToSkip; i ++){
//		//skip if 0
//		if (SFMResults[i] > 0.00000001){
//			sumSFM += SFMResults[i];
//			if(SFMResults[i] > maxSFM) maxSFM = SFMResults[i]; // max
//			if(SFMResults[i] < minSFM) minSFM = SFMResults[i]; // min
//			SFMResultsWithoutZeroes[nonZeroSFM] = SFMResults[i];
//			nonZeroSFM++;
//		}
//	}
//	
//	meanSFM = sumSFM / (float)(nonZeroSFM);
//	
//	// calulate the std. dev.
//	float variance = 0.0f;
//	for(size_t i=0; i<numWindows-numWindowsToSkip; i ++){
//		//skip if 0
//		if (SFMResults[i] > 0.00000001){
//			float diff = SFMResults[i] - meanSFM;
//			variance += diff*diff;
//		}
//	}
//	variance /= (float)(nonZeroSFM-1);
//	float stdDev = sqrtf(variance);
//	
//	//    printf("*** SFM Stats ***\n");
//	//    printf("mean: %f\n", meanSFM);
//	//    printf("std. dev.: %f\n", stdDev);
//	//    printf("min: %f\n", minSFM);
//	//    printf("max: %f\n", maxSFM);
//	
//	stats[0] = meanSFM;
//	stats[1] = stdDev;
//	stats[2] = minSFM;
//	stats[3] = maxSFM;
//	
//	free(SFMResultsWithoutZeroes);
//	free(SFMResults);
//	
//}




typedef struct StatsValue{
	float mean;
	float stdev;
	float min;
	float max;
}SFM_stats;

void print_SFMStats(SFM_stats* row1, SFM_stats* row2, SFM_stats* row3, SFM_stats* row4){
	printf("\t\t mean \t stdev\t min\t max\n");
	printf("mean \t %f\t %f \t%f \t%f \n", row1->mean, row1->stdev, row1->min, row1->max);
	printf("stdev \t %f\t %f \t%f \t%f \n", row2->mean, row2->stdev, row2->min, row2->max);
	printf("min \t %f\t %f \t%f\t %f \n", row3->mean, row3->stdev, row3->min, row3->max);
	printf("max \t %f\t %f \t%f \t%f \n", row4->mean, row4->stdev, row4->min, row4->max);
	
}

void computeStats(float* data, int repeat, SFM_stats* result){
	float mean_value = 0;
	
	float minVal = FLT_MAX_EXP;
	float maxVal = FLT_MIN_10_EXP -37;
	
	float variance = 0;
	float stdev = 0;
	
	for (int i = 0; i<repeat; i++){
		mean_value += data[i] / (float) repeat;
		if (data[i] > maxVal){
			maxVal = data[i];
		}
		if (data[i] < minVal){
			minVal = data[i];
		}
	}
	
	for (int i = 0; i<repeat; i++){
		variance += (data[i] - mean_value) * (data[i] - mean_value);
	}
	variance /= (float) (repeat-1);
	stdev = sqrtf(variance);
	
	result->mean = mean_value;
	result->stdev = stdev;
	result->min = minVal;
	result->max = maxVal;
}


void createAllStatsArray(float* SFMdata, int stride, int repeat, float* output){
	for (int i = 0; i<repeat; i++){
		output[i] = SFMdata[i*4 + stride];
	}
}



void computeStatsMatrix(float* SFMdata, int values){
	//create the 4x4 matrix
	//         mean stdev min max
	// mean
	// stdev
	// min
	// max
	
	int repeat = values/4;
	float* allMean = malloc(sizeof(float)*repeat);
	float* allStdev = malloc(sizeof(float)*repeat);
	float* allMin = malloc(sizeof(float)*repeat);
	float* allMax = malloc(sizeof(float)*repeat);
	
	createAllStatsArray(SFMdata, 0, repeat, allMean);
	createAllStatsArray(SFMdata, 1, repeat, allStdev);
	createAllStatsArray(SFMdata, 2, repeat, allMin);
	createAllStatsArray(SFMdata, 3, repeat, allMax);
	
	//each row is the data of which mean/stdev/min/max is computed
	SFM_stats* row1 = malloc(sizeof(SFM_stats));
	SFM_stats* row2 = malloc(sizeof(SFM_stats));
	SFM_stats* row3 = malloc(sizeof(SFM_stats));
	SFM_stats* row4 = malloc(sizeof(SFM_stats));
	
	computeStats(allMean, repeat, row1);
	computeStats(allStdev, repeat, row2);
	computeStats(allMin, repeat, row3);
	computeStats(allMax, repeat, row4);
	
	print_SFMStats(row1, row2, row3, row4);
	
	free(allMean);
	free(allStdev);
	free(allMin);
	free(allMax);
	
	free(row1);
	free(row2);
	free(row3);
	free(row4);
	
}

void generateWhiteNoise(float* input, int size){
	for (int i = 0; i<size; i++){
		input[i] = 2*(((float)rand()/(float)RAND_MAX)-0.5);
	}
}

void generateStdNormal(float* input, int size){
	for (int i = 0; i<size; i++){
		float S = 1;
		float U=0;
		float V=0;
		while (S>=1){
			U= 2*(((float)rand()/(float)RAND_MAX)-0.5);
			V= 2*(((float)rand()/(float)RAND_MAX)-0.5);
			S= U*U + V*V;
		}
		input[i] = U*sqrtf((-2*log(S))/S);
	}
}


/*@param repeat : how many IRs to produce
 *@param write : enable all IRs produce to file
 *@param method:
 - 0 -> take one SFM per IR, average absFFT per window
 - 1 -> take one SFM per IR, the last window only
 - 2 -> take multiple SFMs per IR
 */
void testFDN(int repeat, bool write, int method){
    // variables to set up the reverb
    float sampleRate = 48000;
    size_t numDelays = 16;
    float minDelayTime = 0.0070f;
    float maxDelayTime = 3.f * minDelayTime;
    size_t IRLength = (size_t)sampleRate * 3;
    float* IR = malloc(sizeof(float)*IRLength);
    float* SFMStats = malloc(sizeof(float)*repeat*4);
    float* SFMPerIR = malloc(sizeof(float)*repeat);
    
    assert(IRLength > maxDelayTime * sampleRate);
    
    float* sampledIR = malloc(sizeof(float)*FFTSIZE);
    
    for (int i = 0; i<repeat; i++){
        // initialise the reverberator
        BMSimpleFDN fdn;
        BMSimpleFDN_init(&fdn,
                         sampleRate,
                         numDelays,
                         DTM_RANDOM,
                         minDelayTime,
                         maxDelayTime,
                         FLT_MAX);

        
        //compute the impulse response
        BMSimpleFDN_impulseResponse(&fdn, IR, IRLength);
        
        // generate white noise for testing
     //   generateWhiteNoise(IR, (int) IRLength);

        //sample the last FFTSIZE samples from the IR
        for(int k = 0; k<FFTSIZE; k++){
            sampledIR[k] = IR[IRLength - FFTSIZE + k];
        }
        
        // compute SFM for the overall IR
        if (method == 0) SFMStatsPerIR(IR, IRLength, &SFMPerIR[i]);
        
        // compute 1 SFM per IR by taking the last window
        if (method == 1) SFMStatsMultiplePerIR(sampledIR, FFTSIZE, &SFMPerIR[i], &SFMStats[i*4], false);
        
        // compute multiple SFMs per IR
        if (method == 2) SFMStatsMultiplePerIR(IR, IRLength, &SFMPerIR[i], &SFMStats[i*4], false);
    
        if (write){
            // output the impulse response to a file
            char *filename = malloc(sizeof(char)*128);
            sprintf(filename, "./impulseResponse%d.csv", i);

            arrayToFileWithName(IR, filename, IRLength);
        }
        
        if (i%50 == 0)
        {
            printf("Completed : %f%% \n", (float)i/(float)repeat * 100);
        }
        
        BMSimpleFDN_free(&fdn);
    }
    
    //write the SFMs to file
    char* filename = "SFMs.csv";
    FILE* audioFile;
    audioFile = fopen(filename, "w+");
    
    for (size_t i=0; i<repeat; i++) {
        fprintf(audioFile, "%f,", SFMPerIR[i]);
    }
    
    fclose(audioFile);
    system("pwd");
    printf("%s",filename);
    printf("\n");
    
    
    if (method == 0 || method == 1){
        bool normalTestResult = false;
        hypothesisTestNormal(SFMPerIR, repeat, &normalTestResult, 1);
        printf("Based on Z-test, batch is white: %d \n", normalTestResult);
        SFM_stats* stats = malloc(sizeof(SFM_stats));
        computeStats(SFMPerIR, repeat, stats);
        
        printf("SFM Mean: %f\n SFM stdev: %f \n SFM min: %f \n SFM max: %f \n", stats->mean, stats->stdev, stats->min, stats->max);
        
    }

    
    //this method doesn't need a hypothesis test because mean of mean of SFM follows normal distribution
    if (method == 2){
        char* filename = "OverallStats.csv";
        
        FILE* audioFile;
        audioFile = fopen(filename, "w+");
        
        // print out all the stats in .csv format
        fprintf(audioFile,"mean, stdev, min, max\n");
        for (size_t i=0; i<repeat*4; i = i+4) {
            fprintf(audioFile, "%f,%f,%f,%f\n",SFMStats[i], SFMStats[i+1], SFMStats[i+2], SFMStats[i+3]);
        }
        
        fclose(audioFile);
        system("pwd");
        printf("%s",filename);
        printf("\n");

        //compute the statistic matrix
        computeStatsMatrix(SFMStats, repeat*4);
    }
    
    free(SFMStats);
    free(IR);
 
}


void testClipperSpeed(){
    size_t testLength = 1024;
    float *testData = malloc(sizeof(float)*testLength);
    
    // start a timer
    clock_t begin, end;
    double time_spent;
    begin = clock();
    
    for(size_t i=0; i<1000000; i++)
    BMAsymptoticLimit(testData, testData, 48000.0f, 1.0f, testLength);
    
    // print the time taken to process
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("asymptotic limiter time: %f\n", time_spent);
    
    
    begin = clock();
    
    int testLengthI = (int)testLength;
    for(size_t i=0; i<1000000; i++)
    vvtanhf(testData,testData,&testLengthI);
    
    // print the time taken to process reverb
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("tanh limiter time: %f\n", time_spent);
}


void testFormatConverter(){
    BMAudioStreamConverter converter;
    AudioStreamBasicDescription input, output;
    input.mSampleRate = 16000.0f;
    output.mSampleRate = 48000.0f;
    input.mChannelsPerFrame = 1;
    output.mChannelsPerFrame = 1;
    input.mFormatID = output.mFormatID = kAudioFormatLinearPCM;
    input.mFormatFlags = kAudioFormatFlagIsPacked | kAudioFormatFlagIsNonInterleaved;
    output.mFormatFlags = input.mFormatFlags | kAudioFormatFlagIsFloat;
    input.mFramesPerPacket = output.mFramesPerPacket = 1;
    input.mBytesPerFrame = sizeof(short) * input.mFramesPerPacket;
    output.mBytesPerFrame = sizeof(float) * output.mFramesPerPacket;

    BMAudioStreamConverter_init(&converter, input, output);
    
    size_t inputLength = 16000 * 5;
    
    // generate short int sine sweep
    float *testSignal_f = malloc(sizeof(float)*inputLength);
    generateSineSweep(testSignal_f, 20.0, 8000.0, input.mSampleRate, inputLength);
    float toShortScale = (float)INT16_MAX;
    vDSP_vsmul(testSignal_f, 1, &toShortScale, testSignal_f, 1, inputLength);
    short *testSignal_si = malloc(sizeof(short)*inputLength);
    vDSP_vfix16(testSignal_f, 1, testSignal_si, 1, inputLength);
    
    // convert the audio stream format
    size_t outputLength = 1 + (double)inputLength * output.mSampleRate / input.mSampleRate;
    float *outputSignal = malloc(8*inputLength*sizeof(float));
    float *outPointer = outputSignal;
    short *inPointer = testSignal_si;
    
    size_t numSamples = inputLength;
    outputLength = 0;
    while(numSamples > 0){
		size_t nextBufferSize = 1 + arc4random() % 100;
        size_t samplesIn = BM_MIN(numSamples, nextBufferSize);
        size_t samplesOut;
        samplesOut = BMAudioStreamConverter_convert(&converter, (const void**)&inPointer, (void**)&outPointer, samplesIn);
        numSamples -= samplesIn;
        outPointer += samplesOut;
        inPointer += samplesIn;
        printf("samplesOut: %zu\n",samplesOut);
        outputLength += samplesOut;
    }
	
//	// convert without changing sample rates
//	vDSP_vflt16(testSignal_si, 1, outputSignal, 1, inputLength);
//	float toFloatScale = 1.0f / toShortScale;
//	vDSP_vsmul(outputSignal, 1, &toFloatScale, outputSignal, 1, inputLength);
	
	char* filename = "./sineSweepResponse.csv";
//    arrayToFileWithName(outputSignal, filename, inputLength);
    arrayToFileWithName(outputSignal, filename, outputLength);
}


void testGroupDelay(){
	BMMultiLevelBiquad bqf;
	BMMultiLevelBiquad_init(&bqf, 1, 48000.0f, false, true, false);
	BMMultiLevelBiquad_setHighPass6db(&bqf, 3000.0f, 0);
	float groupDelay;
	float frequency;
	for(size_t i=0; i<100; i++){
		frequency = 0.01 + (float)i * 10000.0f / 100.0f;
		groupDelay = BMMultiLevelBiquad_groupDelay(&bqf, frequency);
		printf("group delay at %f Hz: %f\n", frequency, groupDelay);
	}
}


void testExtremeCompressor(){
    
}

void testExport(){
    uint32_t sr = 48000;
    uint32_t length = sr * 10;
    
    BMExportWavFile exportWavFile;
    BMExportWavFile_init(&exportWavFile,sr);
    
    float* dataL = malloc(sizeof(float)*length);
    memset(dataL, 0, sizeof(float)*length);
    dataL[0] = 1;
    
    char* filePath = "./sawtooth_test.wav";
    BMExportWavFile_exportAudioFloatToInt16(&exportWavFile,filePath, dataL, dataL, length);
}

void testRandomsInRange(){
    size_t randoms [20];
    BMReverbRandomsInRange(0, 100, randoms, 20);
    
    printf("{");
    for(size_t i=0; i<19; i++)
        printf("%zu, ",randoms[i]);
    printf("%zu}\n",randoms[19]);
}


void testCloudReverb(){
    uint32_t sr = 48000;
    uint32_t length = sr * 10;
    
    BMCloudReverb reverb;
    BMCloudReverb_init(&reverb, sr);
    
    float* outputL = malloc(sizeof(float)*length);
    float* outputR = malloc(sizeof(float)*length);
    //BMCloudReverb_impulseResponse(&reverb, outputL, outputR, length);
    
    //Export wav file
    BMExportWavFile exportWavFile;
    BMExportWavFile_init(&exportWavFile,sr);
    char* filePath = "./sawtooth_test.wav";
    BMExportWavFile_exportAudioFloatToInt16(&exportWavFile,filePath, outputL, outputR, length);
}


void testLongLoopFDN(){
	uint32_t sr = 48000;
    uint32_t length = sr * 10;
    
	BMLongLoopFDN llRv;
	size_t numTaps = 32;
	float maxDelay = 0.50;
	float minDelay = maxDelay / (float) numTaps;
//	BMLongLoopFDN_init(&llRv, numTaps, minDelay, maxDelay, true, 8, 1, sr);
	//BMLongLoopFDN_setInputPan(&llRv, 0.8);
	BMLongLoopFDN_setRT60Decay(&llRv, 100.0f);
    
	float* inputL = calloc(length,sizeof(float));
    float* inputR = calloc(length,sizeof(float));
	inputL[10] = 1.0f;
	inputR[10] = 1.0f;
    float* outputL = calloc(length,sizeof(float));
    float* outputR = calloc(length,sizeof(float));
	BMLongLoopFDN_process(&llRv, inputL, inputR, outputL, outputR, length);
    
    //Export wav file
    BMExportWavFile exportWavFile;
    BMExportWavFile_init(&exportWavFile,sr);
    char* filePath = "./longLoop_test.wav";
	float attenuation = 0.25;
	vDSP_vsmul(outputL, 1, &attenuation, outputL, 1, length);
	vDSP_vsmul(outputR, 1, &attenuation, outputR, 1, length);
    BMExportWavFile_exportAudioFloatToInt16(&exportWavFile,filePath, outputL, outputR, length);
}


void resampleToTarget(size_t sampleWidth, size_t targetSampleWidth, size_t *upsampleFactor, size_t *downsampleFactor){
    
    // more than 2x extra samples
    if (sampleWidth > 2 * targetSampleWidth){
        *upsampleFactor = 1;
        *downsampleFactor = floor((float)sampleWidth/(float)targetSampleWidth);
        return;
    }
    
    // more than 3/2 extra samples
    if (2*sampleWidth > 3*targetSampleWidth){
        *upsampleFactor = 2;
        *downsampleFactor = 3;
        return;
    }
    
    // near target samples (between 2/3 and 3/2 of target)
    if (3*sampleWidth > 2*targetSampleWidth){
        *upsampleFactor = 1;
        *downsampleFactor = 1;
        return;
    }
    
    // less than 1/2 of target samples
    if(2*sampleWidth < targetSampleWidth){
        *upsampleFactor = floor((float)targetSampleWidth/(float)sampleWidth);
        *downsampleFactor = 1;
        return;
    }
    
    // between 1/2 and 2/3 of target
    // else
    *upsampleFactor = 3;
    *downsampleFactor = 2;
}


void testBMSincUpsampler(){
	size_t length = 100;
	float *sineSweep = calloc(length,sizeof(float));
	generateSineSweep(sineSweep, 3000, 3000, 48000, length);
    //memset(sineSweep+10,0,sizeof(float)*20);
	
	BMSincUpsampler us;
	size_t upsampleFactor = 10;
	size_t numInterpolationPoints = 8;
	BMSincUpsampler_init(&us, numInterpolationPoints, upsampleFactor);
	
	float *upsampledSineSweep = calloc(length*upsampleFactor,sizeof(float));
	size_t outputLength = BMSincUpsampler_process(&us, sineSweep, upsampledSineSweep, length);
    
    printf("outputLength: %zu",outputLength);
    
    size_t leftPadding = BMSincUpsampler_inputPaddingBefore(&us);
    size_t rightPadding = BMSincUpsampler_inputPaddingAfter(&us);
	
	printf("{%f, ",upsampledSineSweep[0]);
	for(size_t i=1; i<outputLength-1; i++)
		printf("%f, ",upsampledSineSweep[i]);
	printf("%f}\n\n",upsampledSineSweep[outputLength-1]);
    
	
	free(upsampledSineSweep);
	free(sineSweep);
//    printf("{%f, ", sineSweep[leftPadding]);
//    for(size_t i=leftPadding+1; i<999-rightPadding; i++)
//        printf("%f, ", sineSweep[i]);
//    printf("%f}\n", sineSweep[999-rightPadding]);
}



void testBMSpectrogram(){
	BMSpectrogram sg;
	
	size_t maxFFTSize = 4096;
	size_t maxImageHeight = 1000;
	float sampleRate = 48000.0f;
	BMSpectrogram_init(&sg, maxFFTSize, maxImageHeight, sampleRate);
	
	size_t inputLength = sampleRate;
	float *testInput = malloc(sizeof(float)*inputLength);
	float *buffer = malloc(sizeof(float)*inputLength);
	
	generateSineSweep(testInput, 20.0f, 10000.0f, sampleRate, inputLength);
	generateSineSweep(buffer, 40.0f, 20000.0f, sampleRate, inputLength);
	vDSP_vadd(testInput,1,buffer,1,testInput,1,inputLength);
	
	// add noise
	float noiseLevel = BM_DB_TO_GAIN(-60.0);
	for(size_t i=0; i<inputLength; i++){
		testInput[i] += noiseLevel * (-1.0 + 2.0 * (double)arc4random() / (double)UINT32_MAX);
	}
	
	
	int startSampleIndex = 2000;
	int endSampleIndex = 40000;
	int fftSize = 2048;
	int pixelWidth = 300;
	int pixelHeight = 300;
	float minFrequency = 50.0f;
	float maxFrequency = sampleRate / 2;
    size_t bytesPerPixel = 4;
	uint8_t *imageOutput = malloc(sizeof(uint8_t)*bytesPerPixel*pixelWidth*pixelHeight);
	BMSpectrogram_process(&sg, testInput, (SInt32)inputLength, startSampleIndex, endSampleIndex, fftSize, imageOutput, pixelWidth, pixelHeight, minFrequency, maxFrequency);
	
	char* filename = "./spectrogramImage.txt";
	imageToFileWithName(imageOutput, filename, pixelWidth, pixelHeight);
	
    free(imageOutput);
}






void testDPWOscillator(){
    BMDPWOscillator osc;
    size_t length = 48000;
    size_t oversample = 1;
    BMDPWOscillator_init(&osc, BMDPWO_SAW, 10, oversample, length);

    float *output = malloc(sizeof(float)*length);
    float *frequencies = malloc(sizeof(float)*length);
    
    for(size_t i=0; i<length; i++)
        frequencies[i] = 300.0f;
    
    BMDPWOscillator_process(&osc, frequencies, output, length);
    
    //Export wav file
    BMExportWavFile exportWavFile;
    BMExportWavFile_init(&exportWavFile,48000);
    char* filePath = "./osc_test.wav";
    BMExportWavFile_exportAudioFloatToInt16(&exportWavFile,filePath, output, output, length);
}




void testCDBlepOscillator(){
	BMCDBlepOscillator osc;
	size_t length = 48000;
	size_t numBleps = 20;
	size_t filterOrder = 2;
	size_t oversampleFactor = 2;
	BMCDBlepOscillator_init(&osc, numBleps, filterOrder, oversampleFactor, length);

	float *output = malloc(sizeof(float)*length);
	float *frequencies = malloc(sizeof(float)*length);
	
	for(size_t i=0; i<length; i++)
		frequencies[i] = 1040.324f;
	
	BMCDBlepOscillator_process(&osc, frequencies, output, length);
	
	// drop the volume to avoid clipping
	float half = 0.5f;
	vDSP_vsmul(output, 1, &half, output, 1, length);
	
	//Export wav file
	BMExportWavFile exportWavFile;
	BMExportWavFile_init(&exportWavFile,48000);
	char* filePath = "./osc_test.wav";
	system("pwd");
	BMExportWavFile_exportAudioFloatToInt16(&exportWavFile,filePath, output, output, length);
}


void testBlipOscillator(){
    BMBlipOscillator osc;
    size_t length = 48000;
    size_t numBleps = 32;
    size_t filterOrder = 8;
    size_t oversampleFactor = 4;
    float sampleRate = 48000;
    BMBlipOscillator_init(&osc, sampleRate, oversampleFactor, filterOrder, numBleps);
    BMBlipOscilaltor_setLowpassFc(&osc, 4800.0f);

    float *output = malloc(sizeof(float)*length);
    float *frequencies = malloc(sizeof(float)*length);
    
    for(size_t i=0; i<length; i++)
        frequencies[i] = 251.123f;
    
    int length_i = (int)length;
    vvlog2f(frequencies, frequencies, &length_i);
    
    BMBlipOscillator_process(&osc, frequencies, output, length);
    
    // drop the volume to avoid clipping
    float half = 0.5f;
    vDSP_vsmul(output, 1, &half, output, 1, length);
    
    //Export wav file
    BMExportWavFile exportWavFile;
    BMExportWavFile_init(&exportWavFile,48000);
    char* filePath = "./osc_test.wav";
    system("pwd");
    BMExportWavFile_exportAudioFloatToInt16(&exportWavFile,filePath, output, output, length);
}




void testBlip(){
	BMBlip osc;
	size_t length = 48000;
	size_t numBleps = 20;
	size_t filterOrder = 16;
	size_t oversampleFactor = 1;
	float sampleRate = 48000;
	BMBlip_init(&osc, 2, 400.0f, sampleRate);
	BMBlip_restart(&osc, 0.5f);

	float *output = calloc(length, sizeof(float));
	
	BMBlip_process(&osc, output, length);
	
	// drop the volume to avoid clipping
	float half = 0.5f;
	vDSP_vsmul(output, 1, &half, output, 1, length);
	
	//Export wav file
	BMExportWavFile exportWavFile;
	BMExportWavFile_init(&exportWavFile,48000);
	char* filePath = "./osc_test.wav";
	system("pwd");
	BMExportWavFile_exportAudioFloatToInt16(&exportWavFile,filePath, output, output, length);
}


void randomVector(float *v, size_t length){
	// seed the random number generator
	srand48(time(0));
	
	// fill the array with random numbers
	for(size_t i=0; i<length; i++)
		v[i] = drand48();
	
	// normalise the result
	BMVectorNormalise(v, length);
}


void testIFFT(){
	size_t length = 128;
	BMFFT fft;
	BMFFT_init(&fft, length);
	
	float *input = malloc(sizeof(float)*length);
	float zero = 0.0f;
	float one = 1.0f;
	float unitaryMagnitude = 1.0f / sqrt(length);
	float increment = 1.0f / (6.493687f * length);
	//vDSP_vramp(&zero, &increment, input, 1, length);
	//vDSP_vfill(&unitaryMagnitude, input, 1, length);
	randomVector(input, length);
	
	// print the vector magnitude of the input
	float norm = BMVectorNorm(input,length);
	printf("input norm: %f\n", norm);
	
	// do the fft
	DSPSplitComplex zDomain;
	zDomain.imagp = calloc(length/2,sizeof(float));
	zDomain.realp = calloc(length/2,sizeof(float));
	BMFFT_FFTComplexOutput(&fft, input, &zDomain, length);
	
	// scale the zDomain data
//	float fftScale = 1.0 / (2.0 * sqrt(length));
//	vDSP_vsmul(zDomain.imagp, 1, &fftScale, zDomain.imagp, 1, length/2);
//	vDSP_vsmul(zDomain.realp, 1, &fftScale, zDomain.realp, 1, length/2);
	
	// print the vector magnitude of the zDomain data
	norm = BMVectorNormSplitComplex(zDomain,length/2);
	printf("zDomain norm: %f\n", norm);
	
	// print the vector magnitude of the real and imaginary components
	norm = BMVectorNorm(zDomain.realp,length/2);
	printf("real norm: %f\n", norm);
	norm = BMVectorNorm(zDomain.imagp,length/2);
	printf("imag norm: %f\n", norm);
	
	// do the ifft
	float *output = malloc(sizeof(float)*length);
	BMFFT_IFFT(&fft, &zDomain, output, length/2);
	
	// scale the output data
//	float ifftScale = 1.0 / sqrt(length);
//	float combinedScale = 1.0 / (2 * length);
//	vDSP_vsmul(output, 1, &combinedScale, output, 1, length);
	
	// print the vector magnitude of the output
	norm = BMVectorNorm(output, length);
	printf("output norm: %f\n", norm);
	
	// print the input and output
	printf("********** INPUT ************\n{");
	for(size_t i=0; i<length-1; i++) printf("%f, ", input[i]);
	printf("%f}\n",input[length-1]);
	printf("********** OUTPUT ************\n{");
	for(size_t i=0; i<length-1; i++) printf("%f, ", output[i]);
	printf("%f}\n",input[length-1]);
	
	// print the vector distance from input to output
	float distanceSq;
	vDSP_distancesq(input, 1, output, 1, &distanceSq, length);
	printf("distance: %f\n", sqrtf(distanceSq));
	
	free(input);
	free(output);
	free(zDomain.imagp);
	free(zDomain.realp);
}

//
//void testTherapyFilter (void) {
//	size_t length = 1024;
//	
//	BMVagusNerveTherapyFilter filter;
//	float sampleRate = 48000;
//	BMVagusNerveTherapyFilter_init(&filter, sampleRate);
//
//	size_t bufferSize = length;
//	float *inputL, *inputR;
//	float *outputL, *outputR;
//	float *buffer = (float*)malloc(sizeof(float) * bufferSize * 4);
//	memset(buffer, 0, sizeof(float) * bufferSize * 4);
//	inputL  = buffer + (bufferSize * 0);
//	inputR  = buffer + (bufferSize * 1);
//	outputL = buffer + (bufferSize * 2);
//	outputR = buffer + (bufferSize * 3);
//
//	// offset the impulse by 100 samples to prevent wierd things that happen
//	// in the first few samples of processing dynamic filters. However,
//	// if the length of the IR is short, start from 0 instead.
//	size_t impulseIndex = length < 200 ? 0 : 100;
//	inputL[impulseIndex] = 1.0f;
//	inputR[impulseIndex] = 1.0f;
//
//	BMVagusNerveTherapyFilter_process(&filter, inputL, inputR, outputL, outputR, bufferSize);
//
//	// print Impulse Response to console
//	printf("{");
//	size_t i=0;
//	for(; i<length-1; i++)
//		printf("%f, ", outputL[i]);
//	printf("%f}", outputL[i]);
//
//	free(buffer);
//	buffer = NULL;
//}


int main(int argc, const char * argv[]) {
	testStaticDelay();
    return 0;
}



