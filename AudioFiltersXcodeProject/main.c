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
#include "BMHIIRUpsampler2x.h"
#include "BMUpsampler.h"
#include "BMDownsampler.h"
#include "BMFIRFilter.h"
#include "BMGainStage.h"
#include "BMMultiLevelBiquad.h"
#include "BMVelvetNoiseDecorrelator.h"
#include "BMRoundRobin.h"
#include "BMBiquadArray.h"
#include "BMMultiTapDelay.h"
#include "BMReverb.h"
#include "BMStereoModulator.h"
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
#include "BMHIIRUpsampler2xFPU.h"
#include "BMHIIRDownsampler2xFPU.h"
#include "BMWaveshaping.h"

#define TESTBUFFERLENGTH 128


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
    BMVelvetNoiseDecorrelator_init(&vndFilter, true, 80,80, 100,100, 48000,80);
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

void testStereoModulator(){
    float sampleRate = 44100;
    BMStereoModulator stereoModFilter;
    
    BMStereoModulator_init(&stereoModFilter, sampleRate,S_Tightness_DV,S_Tightness_Max,S_Density_DV,S_Density_Max,S_Quad_DV);
    
    size_t frameCount = 256;
    float* input = malloc(sizeof(float)*frameCount);
    float* output = malloc(sizeof(float)*frameCount);
    int zero = 0;
    input[0] = 1;
    memset(input+1, zero, sizeof(float)*(frameCount-1));
    int random = 0;
    float factor = 0;
    while (true) {
        BMStereoModulator_processStereo(&stereoModFilter, input, input, output, output, frameCount);
        //Change parameter
        random = arc4random()%5;
        factor = (arc4random()%100)/100.;
        if(random==0){
            //Change numtap
            int newNumtap = roundf((S_Density_Max - S_Density_Min)*factor) + S_Density_Min;
            BMStereoModulator_setNumTap(&stereoModFilter, newNumtap);
        }else if(random==1){
            //Change delaytime
            int tightness = roundf((S_Tightness_Max-S_Tightness_Min)*factor) + S_Tightness_Min;
            BMStereoModulator_setDelayTime(&stereoModFilter, tightness);
        }else if(random==2){
            //Re-init struct
            BMStereoModulator_destroy(&stereoModFilter);
            BMStereoModulator_init(&stereoModFilter, sampleRate,S_Tightness_DV,S_Tightness_Max,S_Density_DV,S_Density_Max,S_Quad_DV);
        }
    }
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
    BMMultiLevelSVF_init(&svf, 1, 44100, false);
    
    BMMultiLevelSVF_setBell(&svf, 2000, -40, 2, 0);
    
    BMMultiLevelSVF_impulseResponse(&svf, 1000);
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
 * This function behaves like the function g(x) = x for x near zero but
 * its response tapers off asymptotically to +- limit as abs(x) increases.
 *
 * THIS FUNCTION DOES NOT WORK IN PLACE. (requires input != output)
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
 * This is fastest on ARM processors
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
    // open a file for writing
    FILE* audioFile;
    audioFile = fopen("./arrayOut.csv", "w+");
    
    // print out the entire frame in .csv format
    //fprintf(audioFile,"{");
    for (size_t i=0; i<length; i++) {
        fprintf(audioFile, "%f,%f\n",(float)i/(float)length,array[i]);
        
    }
    
    fclose(audioFile);
    
    system("pwd");
    printf("/arrayOut.csv\n");
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

void testGainStage(){
    BMGainStage gs;
    BMGainStage_init(&gs, BM_GAINSTAGE_TEST_SAMPLE_RATE, 10000.0f);
    
    // create sine wave for input test
    float frequency = 300.0f;
    float t [BM_GAINSTAGE_TEST_LENGTH];
    float sr = BM_GAINSTAGE_TEST_SAMPLE_RATE;
    for (size_t i=0; i<BM_GAINSTAGE_TEST_LENGTH; i++) {
        t[i] = sinf(frequency * 2.0 * M_PI * (float)i / sr);
        
        if (i>20000 && i < 30000) t[i] = 0.0f;
    }
    
    // set up output array
    float o [BM_GAINSTAGE_TEST_LENGTH];
    
    // set the tightness
    BMGainStage_setTightness(&gs, 0.6f);
    
    // set the bias
    BMGainStage_setBiasRatio(&gs, 0.0f);
    
    // set the input gain
    BMGainStage_setGain(&gs, BM_DB_TO_GAIN(15.0f));
    
    // process the test signal
    BMGainStage_processBufferMono(&gs, t, o, BM_GAINSTAGE_TEST_LENGTH);
    
    // print the output
    arrayToFile(o,BM_GAINSTAGE_TEST_LENGTH);
    
}

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
    
    BMMultiLevelBiquad bqf;
    BMMultiLevelBiquad_init(&bqf, 4, sampleRate, false, true, false);
    BMMultiLevelBiquad_setCriticallyDampedLP(&bqf, filterFrequency, 0, 4);
    
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
    
    BMMultiLevelBiquad_processBufferMono(&bqf, input, input, totalLength);
    
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
    
    BMMultiLevelBiquad bqf;
    size_t numLevels = 4;
    BMMultiLevelBiquad_init(&bqf, numLevels, sampleRate, false, true, false);
    BMMultiLevelBiquad_setCriticallyDampedLP(&bqf, filterFrequency, 0, numLevels);
    
    
    float frequencies [1024];
    float magnitudes [1024];
    
    for(size_t i=0; i < 1024; i++){
        float maxF = 300.0f;
        float alphaMax = log2f(maxF);
        float alpha = ((float)i / 1024.0f) * alphaMax;
        frequencies[i] = powf(2.0f,alpha);
    }
    
    BMMultiLevelBiquad_tfMagVector(&bqf, frequencies, magnitudes, 1024);
    
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

#define LT_SampleRate 48000.
#define LT_Frequency 10.
#define LT_TestRange 50
void testBMLagTime(){
    BMStereoLagTime stereoLagTime;
    BMStereoLagTime_init(&stereoLagTime, 500, LT_SampleRate);
    
    
    size_t length = LT_SampleRate *4;
    float* audioX = malloc(sizeof(float)*length);
    float* inLAudioY = malloc(sizeof(float)*length);
    float* outLAudioY = malloc(sizeof(float)*length);
    float* inRAudioY = malloc(sizeof(float)*length);
    float* outRAudioY = malloc(sizeof(float)*length);
    
    float scale = 1.;
    float xFactor = 20;
    for(int i=0;i<length;i++){
        audioX[i] = (i*xFactor)/length;
        inLAudioY[i] = sinf(LT_Frequency * 2.0 * M_PI * (float)i / LT_SampleRate) * scale;
        inRAudioY[i] = sinf(LT_Frequency * 2.0 * M_PI * (float)i / LT_SampleRate) * scale;
    }
    
    float delayTime[LT_TestRange];
    for(int i=0;i<30;i++){
        delayTime[i] = (LT_TestRange - i)*1;
    }
    
    for(int i=10;i<LT_TestRange;i++){
        delayTime[i] = i*2;
    }
    size_t processedSample = 0;
    size_t processingSample = 1024;
    for(int i=0;i<LT_TestRange;i++){
        BMStereoLagTime_setDelayLeft(&stereoLagTime, delayTime[i]);
        BMStereoLagTime_process(&stereoLagTime, inLAudioY+processedSample, inRAudioY+processedSample, outLAudioY+processedSample, outRAudioY+processedSample, processingSample);
        processedSample += processingSample;
    }
    
    arrayXYToFile(audioX,outLAudioY, length);
}

#define BS_Range 10
#define BS_Freq_Range 1000.0f
#define BS_Freq_Start 1.0f
#define BS_Freq_End 20000.0f
void testBinauralSynthesis(){
    float* input = malloc(sizeof(float) * BS_Range);
    float* output = malloc(sizeof(float) * BS_Range);
    float* frequency = malloc(sizeof(float)*BS_Freq_Range);
    float* magnitude = malloc(sizeof(float)*BS_Freq_Range);
    
    float alpha = powf(BS_Freq_End/BS_Freq_Start, 1.0f/BS_Freq_Range);
    
    float* xArray = malloc(sizeof(float)*BS_Freq_Range);
    
    for(int i = 0;i< BS_Freq_Range;i++){
        frequency[i] = BS_Freq_Start * powf(alpha, i);
        xArray[i] = i * 10;
    }
    
    
    
    BMBinauralSynthesis bs;
    BMBinauralSynthesis_init(&bs, false, LT_SampleRate);
    
    BMBinauralSynthesis_setAngleLeft(&bs, 0);
    
    BMBinauralSynthesis_processMono(&bs, input, output, BS_Range);
    
    BMBinauralSynthesis_getTFMagVectorData(&bs,frequency, magnitude, BS_Freq_Range, 0);
    //Convert from gain to db
    for(int i=0;i<BS_Freq_Range;i++){
        magnitude[i] = BM_GAIN_TO_DB(magnitude[i]) * 1000;
    }
    
    arrayXYToFile(xArray,magnitude, BS_Freq_Range);
    
}



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


double frequencyToPhaseIncrement(double frequency, double sampleRate){
    return (frequency * 2.0 * M_PI) / sampleRate;
}




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

/*
 * Linear sweep, matrix method
 */
//    BMVector2D state = {0.0, 1.0};
//    BM2x2MatrixD frequencyMatrix = BM2x2MatrixD_rotationMatrix(frequencyToPhaseIncrement(startFrequency, sampleRate));
//    BM2x2MatrixD frequencyIncreaseMatrix = BM2x2MatrixD_rotationMatrix(M_PI/(double)(numSamples-1));
//
//    for(size_t i=0; i<numSamples; i++){
//        output[i] = state.x;
//        state = BM2x2MatrixD_mvmul(frequencyMatrix, state);
//        frequencyMatrix = BM2x2MatrixD_mmul(frequencyMatrix, frequencyIncreaseMatrix);
//    }

}


void testUpsampler2x(){
    BMHIIRUpsampler2x us;
    float attenuation = BMHIIRUpsampler2x_init(&us, 7, 0.075);
    printf("\nstopband attenuation: %f\n",attenuation);
    float sampleRate = 48000.0f;
    size_t testLength = sampleRate * 5;
    
    float* sineSweep = malloc(sizeof(float)*testLength);
//    float* sineSweep4 = malloc(sizeof(float)*testLength*4);
    float* output = malloc(sizeof(float)*testLength*2);
    
    generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
    
    // duplicate each sample 4 times to get sineSweep4
//    for(size_t i=0; i<testLength; i++){
//        for(size_t j=0; j<4; j++)
//            sineSweep4[4*i + j] = sineSweep[i];
//    }
    
    BMHIIRUpsampler2x_processBufferMono(&us, sineSweep, output, testLength);
    
    arrayToFile(sineSweep,testLength);
    
    free(sineSweep);
    free(output);
}

void testSimdMulAdd(){
    simd_float4 x, y, z, result;
    x = 100.0f;
    y = 1.0f;
    z = 3.0f;
    
    // result = x + y*z;
    result = simd_muladd(x, y, z);
    
    printf("\nresult: %f\n",result[0]);
}

void testUpsampler2xFPU(){
    BMHIIRUpsampler2xFPU us;
    float attenuation = BMHIIRUpsampler2xFPU_init(&us, 7, 0.075);
    printf("\nstopband attenuation: %f\n",attenuation);
    float sampleRate = 48000.0f;
    size_t testLength = sampleRate * 5;
    
    float* sineSweep = malloc(sizeof(float)*testLength);
    //    float* sineSweep4 = malloc(sizeof(float)*testLength*4);
    float* output = malloc(sizeof(float)*testLength*2);
    
    generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
    
    BMHIIRUpsampler2xFPU_processBufferMono(&us, sineSweep, output, testLength);
    
    arrayToFile(output,2*testLength);
    
    free(sineSweep);
    free(output);
}



void testDownsampler2xFPU(){
    BMHIIRDownsampler2xFPU ds;
    float attenuation = BMHIIRDownsampler2xFPU_init(&ds, 7, 0.075);
    printf("\nstopband attenuation: %f\n",attenuation);
    float sampleRate = 48000.0f;
    size_t testLength = sampleRate * 10;
    
    float* sineSweep = malloc(sizeof(float)*testLength);
    //    float* sineSweep4 = malloc(sizeof(float)*testLength*4);
    float* output = malloc(sizeof(float)*testLength/2);
    
    generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
    
    BMHIIRDownsampler2xFPU_processBufferMono(&ds, sineSweep, output, testLength);
    
    arrayToFile(output,testLength/2);
    
    free(sineSweep);
    free(output);
}





void testUpDownsampler2xFPU(){
    BMHIIRDownsampler2xFPU ds;
    BMHIIRUpsampler2xFPU us;
    BMHIIRUpsampler2xFPU_init(&us, 100.0, 0.02);
    float attenuation = BMHIIRDownsampler2xFPU_init(&ds, 100.0, 0.02);

    printf("\nstopband attenuation: %f\n",attenuation);
    float sampleRate = 48000.0f;
    size_t testLength = sampleRate * 10;
    
    float* sineSweep = malloc(sizeof(float)*testLength);
    float* upsampled = malloc(sizeof(float)*testLength*2);
    float* output = malloc(sizeof(float)*testLength);
    
    generateSineSweep(sineSweep, 20.0f, 24000.0f, 48000.0f, testLength);
    
    // upsample
    BMHIIRUpsampler2xFPU_processBufferMono(&us, sineSweep, upsampled, testLength);
    
//    // amplify
//    float gainBoost = 2.0f;
//    vDSP_vsmul(upsampled, 1, &gainBoost, upsampled, 1, testLength*2);
    
    // waveshape
    BMWaveshaper_processBufferBidirectional(upsampled, upsampled, testLength*2);
                                            
    // downsample
    BMHIIRDownsampler2xFPU_processBufferMono(&ds, upsampled, output, testLength*2);
    
    arrayToFile(output,testLength);
    
    free(sineSweep);
    free(output);
}



int main(int argc, const char * argv[]) {
    // testGainStage();
    // testAsymptoticLimitOutputRange();
    // testCompressor();
    // testAsymptoticThreshold();
    // testCriticallyDampedFilterStepResponse();
    // testCriticallyDampedFilterFrequencyResponse();
    // testEnvelopeFollower();
    // testEnvReleaseTime();
    // logSpeedTest();
//    quadraticThresholdTest();
//    testBMLagTime();
    // testBinauralSynthesis();
//    testWetDryMixer();
//    testSimdMulAdd();
//    testUpsampler2xFPU();
    testUpDownsampler2xFPU();
    return 0;
}



