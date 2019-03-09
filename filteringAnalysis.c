/*****************************************************************************
*
* File Name: filteringAnalysis.c
*
* Description: This program takes a timedomain signal in from a file, intepolates by a harcoded amount,
*				filters the signal via coefficients that are hardcoded
*				into the program and then decimates the filter output again by a hardcoded amount to provide the required output sample rate. 
* 
* 
* 				It will output one file with columns for 
*				the input signal, 
*				the output signal, 
*				frequency scale, 
*				frequency spectra for the input and the output signal and 
*				the frequency response of the filter.
*
* Programmer: Sean Mullery, Shane Ormonde
*
* Date: 09/03/19
*
* Version: 0.1
*				
*
******************************************************************************/

/******************************include files*********************************/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"wave.h"

/****************************defined constants*******************************/
#define NUMCOEFF 7 /*the number of coefficients in the filter.*/
#define SIGSIZE 110000 /*Size of the signals, 110013 is the exact size of the 
						* PianoC6 file, which has be rounded here.
						* In future versions this will be replaced
						* with a variable, and malloc used to allocate the 
						* memory for all the arrays of this size.*/
#define SIGSEGSIZE 4000 /*The size of the segment to be analysed by the DFT */
#define DFTSIZE 2001 /*Length of DFT */
#define NYQUIST 22050.0 /* Input file Nyquist Frequency in Hz, must be floating point*/
#define OUTPUT_NYQUIST 8000.00 /*output file nyquist frequency in Hz, must be floating point*/
#define INTERPOLATE_AMOUNT 160 /*interpolation factor*/
#define DECIMATE_AMOUNT    441 /*decimation factor*/
#define INTER_ARRAY_SIZE SIGSIZE*INTERPOLATE_AMOUNT /*array length required to store interpolated signal*/
#define DECI_ARRAY_SIZE  INTER_ARRAY_SIZE/DECIMATE_AMOUNT /*array length required to store decimated signal*/

#define TRUE 1
#define FALSE 0


/****************************Function Prototypes******************************/
void dft(const float* signal, float* mag, float* pha, const int sigSize);

void filter(const float* sigIn, float* sigOut, const float* aCoeff, 
									const float* bCoeff, const int sigSize); 
									
void interpolate(float *signal, float *interSignalOut , const int oldSigSize);

void decimate(float *signal, float *deciSignalOut, const int newSigSize);

void takeInFrom16BitWav(float* signal, const int sigSize);

void sendOutTo16BitWav(float* signal, const int sigSize);

void takeInFromCsv(float* signal, const int sigSize);


struct HEADER header;

/*****************************main function**********************************/
int main(void) {
	/*Note that if SIGSIZE and SIGSEGSIZE are very large they may cause stack 
	memory to run out and cause program failure. These will be replaced by 
	dynamic allocation of the arrays in the next version of this program. */
	float signalIn[SIGSIZE] = {0}; 
	float sigInFreqMag[SIGSEGSIZE]={0};
	float sigInFreqPha[SIGSEGSIZE]={0};
	
	float *interpolateSignalOut = malloc(INTER_ARRAY_SIZE * sizeof(float)); /*allocate memory for interpolated signal*/
	float *decimateSignalOut    = malloc(DECI_ARRAY_SIZE  * sizeof(float)); /*allocate memory for decimated signal*/
	
	
	float *filterSignalOut = malloc(INTER_ARRAY_SIZE * sizeof(float)); /*allocate memry for output from filter*/
	
	float sigOutFreqMag[SIGSEGSIZE]={0};
	float sigOutFreqPha[SIGSEGSIZE]={0};
		
	float impulse[SIGSEGSIZE] = {1};
	
	float impulseResp[SIGSEGSIZE] = {0};
	float freqRespMag[SIGSEGSIZE] = {0};
	float freqRespPha[SIGSEGSIZE] = {0};
	
	/******************Filter Coefficients*************************************/	
	/*The current values are for a six pole Lowpass chebyshev filter with 
	* a cutoff of 0.1 of the sampling Frequency 
	* see http://www.dspguide.com/ch20/2.htm for lists of filter coefficients*/
	
	
	// filter coefficients for 0.00099206 of sample rate (44.1*0.00099206) = 7Khz
	float aCoeff[NUMCOEFF] = {1.8344e-18 ,  1.1006e-17  , 2.7516e-17  , 3.6688e-17 ,  2.7516e-17  , 1.1006e-17  , 1.8344e-18};
	float bCoeff[NUMCOEFF] = {  -5.99564 ,  14.97823 , -19.95650  , 14.95655 ,  -5.97830  ,  0.99566};
	
	int i;	
	
	/*********************Initialise File IO**********************************/
	FILE * outFile;
	
	outFile = fopen("outputsAnalysis.csv", "w");
			
	if(outFile == NULL){
		printf("problem Opening output file\n");
		return 1;
	}
	
	
	/***********Take in the input signal**************************************/	
	/*uncomment the following line (and comment out the one after if you wish to take
	* the signal in fromm a csv file instead of a wav file.*/
	//takeInFromCsv(signalIn,SIGSIZE);
	takeInFrom16BitWav(signalIn,SIGSIZE);

	/**Get the inpulse response of the filter and store in impulseResp*********/	
	filter(impulse,impulseResp,aCoeff,bCoeff,SIGSEGSIZE);
	
	/**take the input signal and interpolate it by the defined amount**/
	interpolate(signalIn,interpolateSignalOut,SIGSIZE);
	
	
	/*****Filter the input signal signalIn and put the filtered version in 
		the output signal signalOut********************************************/
	filter(interpolateSignalOut,filterSignalOut,aCoeff,bCoeff,INTER_ARRAY_SIZE);
	
	
	/****decimate the output from the filter *****/
	
	decimate(filterSignalOut, decimateSignalOut, DECI_ARRAY_SIZE);
	
	
	/*****Get the Frequency Response of the filter****************************/
	dft(impulseResp, freqRespMag, freqRespPha, SIGSEGSIZE);
	
	/****Get the input and output signal's frequency spectra*****************/
	dft(signalIn, sigInFreqMag, sigInFreqPha, SIGSEGSIZE);
	dft(decimateSignalOut, sigOutFreqMag, sigOutFreqPha, SIGSEGSIZE);
	

	/*********Dump the outputs to a csv file **************************************/
	fprintf(outFile,"Input Signal, Output Signal, Freq Scale, Mag Spec Input,"
					" Phase Spec Input, Impulse Resp, Freq Resp Mag,"
					" Freq Resp Phase, Mag Spec Out, Mag Phase Out\n");
					
	/*Only going to print to DFTSIZE so that all columnns are equal length *******
		and we don't go off the end of the freq arrays****************************/						
	for(i=0; i<DFTSIZE; i++){
				
		
		fprintf(outFile,"%f,", signalIn[i]); /*Input Signal*/
		fprintf(outFile,"%f,", decimateSignalOut[i]);/*Output Signal*/ 
		
		fprintf(outFile,"%f,", i*NYQUIST/DFTSIZE);
		fprintf(outFile,"%f,", sigInFreqMag[i]);/*Mag Freq Spect Input Signal*/
		fprintf(outFile,"%f,", sigInFreqPha[i]);/*Pha Freq Spect Input Signal*/
		
		fprintf(outFile,"%f,", impulseResp[i]);/*Impulse Response of Filter*/
		fprintf(outFile,"%f,", freqRespMag[i]);/*Mag Freq Response of Filter*/
		fprintf(outFile,"%f,", freqRespPha[i]);/*Pha Freq Response of Filter*/
				
		fprintf(outFile,"%f,", sigOutFreqMag[i]);/*Mag Freq Spect Output Signal*/
		fprintf(outFile,"%f,", sigOutFreqPha[i]);/*Pha Freq Spect Output Signal*/
		
		fprintf(outFile,"\n");
	}
	
	
	
	sendOutTo16BitWav(decimateSignalOut, DECI_ARRAY_SIZE);

	/***********************Clean up and close files*****************************/


	fclose(outFile);
	free(interpolateSignalOut);
	free(decimateSignalOut);
	system("pause");
	return 0;	
	
}

/******************************************************************************
*
* Function Name: interpolate
*
* Input Parameters: pointer input signal array, pointer to output signal array an int containing the size of the input array
* 							
*
* Output Parameters: float* interSignalOut is the output
*							
*
* Returns: void
*
* Purpose of Function: take in an array, interpolate it be a set number, and return a pointer to a new array of interpolated data
*
******************************************************************************/

void interpolate(float *signal, float *interSignalOut, const int oldSigSize) {
	
	int i;
	int h;
	int x;
	
	
	for (i=0; i<oldSigSize; i++) {                 /*iterate through each element of input array*/
		interSignalOut[h] = signal[i];             /*copy input array element to output array*/
		for (x=1; x<=INTERPOLATE_AMOUNT; x++){     /*fill interSignalOut with amount of zeroes specified by (INTERPOLATE_AMOUNT)*/
			interSignalOut[h+x] = 0;			   
		}
		h = h + (INTERPOLATE_AMOUNT-1);            /*now that zeros have been added , update h to point to latest element after newest zeroes*/
			
	   
	}
	
	
	
	
	
}

/******************************************************************************
*
* Function Name: decimate
*
* Input Parameters: pointer input signal array, pointer to output signal array an int containing the size of the input array
* 						
*
* Output Parameters: float* deciSignalOut is the output
*
* Returns: void
*
* Purpose of Function: take in an array, decimate it by a set number, and return a pointer to a new array of decimated data
*
******************************************************************************/
void decimate(float *signal, float *deciSignalOut, const int newSigSize){
	int i,h=0;
	
	for (i=0; i<newSigSize; i++) {         /*iterate through each element of input array*/
		deciSignalOut[i] = signal[h];      /*copy input signal element to deciSignalOut array*/
		h = h + (DECIMATE_AMOUNT-1);       /*skip ahead by DECIMATE_AMOUNT-1 to perfrom decimation action*/
	}
	
	
}


/******************************************************************************
*
* Function Name: takeInFromCsvFile
*
* Input Parameters: const sigSize is the length of the signal to taken in
*							
*
* Output Parameters: float* signal is where the signal from the csv file will 
*							be stored for later use
*
* Returns: void
*
* Purpose of Function: Take in data from a csv file and store in a array
*						It is assumed the data is floating point data
*
******************************************************************************/
void takeInFromCsv(float* signal, const int sigSize){
	FILE * inFile;	
	inFile = fopen("signal.csv", "r");
	int i;

	if(inFile == NULL){
		printf("problem Opening input csv file\n");
		system("pause");
		exit(1);
	}
	/*take in the data*/
	for(i=0; i<SIGSIZE; i++){		
		fscanf(inFile, "%f\n", &signal[i]);
	}

	fclose(inFile);
	return;	
}


/******************************************************************************
*
* Function Name: takeInFrom16BitWav
*
* Input Parameters: const sigSize is the length of the signal to taken in
*							
*
* Output Parameters: float* signal is where the signal from the Wav file will 
*							be stored for later use
*
* Returns: void
*
* Purpose of Function: Take in data from a Wav file and store in a array
*						It is assumed the data is 16 bits per sample stereo
*						Stereo means 2 channels.
*						Only one of the channels is stored in the array
* 						The 16 bits are modified to be floating point numbers
*						between -1 and +1.  
*
******************************************************************************/
void takeInFrom16BitWav(float* signal, const int sigSize){
	unsigned int buffer4;
	unsigned short int buffer2; 
	FILE *inFile;
	unsigned int  xChannels = 2;
	long numSamples;
	long int  i =0;
    short dataBuffer;
	FILE * outFile;

	
	printf("Opening Wav file..\n");
	/* following file will have to be in the same directory as the executible
	* You can change this for the name of another file you wish to try */
 	inFile = fopen("PianoC6.wav", "rb");

 	if (inFile == NULL) {
    	printf("Error opening input wav file\n");
		system("pause");
    	exit(1);
 	}	
 
 // read header parts
 	fread(header.riff, sizeof(header.riff), 1, inFile);
	fread((unsigned char*)&header.overallSize, 4, 1, inFile);
	printf("Overall Size of input WAV file is %ld\n", header.overallSize);
	fread(header.wave, sizeof(header.wave), 1, inFile);
	fread(header.fmtChunkMarker, sizeof(header.fmtChunkMarker), 1, inFile);
	fread((unsigned char*)&header.lengthOfFmt, sizeof(buffer4), 1, inFile);
 	fread((unsigned char*)&header.formatType, 2, 1, inFile); 
	fread((unsigned char*)&header.channels, 2, 1, inFile);
	printf("Number of Channels in input Wav file is %d\n", header.channels);
	fread((unsigned char*)&header.sampleRate, 4, 1, inFile);
	printf("Sample Rate used in input WAV file is %d Hz\n", header.sampleRate);
	fread((unsigned char*)&header.byterate, 4, 1, inFile);
	fread((unsigned char*)&header.blockAlign, 2, 1, inFile);
	fread((unsigned char*)&header.bitsPerSample, 2, 1, inFile);
	fread(header.dataChunkHeader, sizeof(header.dataChunkHeader), 1, inFile);
	fread((unsigned char*)&header.dataSize, 4, 1, inFile);
	printf("Data Size of input WAV file is %ld\n", header.dataSize);
	
 	/*calculate the number of samples in the file*/	
	numSamples = (8 * header.dataSize) / (header.channels * header.bitsPerSample);
	/*Report the number of samples in the file, note only sigSize samples will be 
	stored in the signal array */

	printf("Number of Samples in input WAV file: %ld bytes\n\n", numSamples);
   
    for (i=0; i<sigSize; i++) {
                
        for (xChannels = 0; xChannels < header.channels; xChannels ++ ) {    
			fread((unsigned char*)&dataBuffer, sizeof(dataBuffer), 1, inFile);
            if(xChannels ==0){/*we will only take a single channel*/			   
				signal[i]= dataBuffer/32768.0;
        	}
        }
           
	}
	fclose(inFile);
	//fclose(outFile);
}


/******************************************************************************
*
* Function Name: sendOutTo16BitWav
*
* Input Parameters: const sigSize is the length of the signal that will be 
*					written to the WAV file
*							
*
* Output Parameters: float* signal is that is to be written to the WAV file. 
*
* Returns: void
*
* Purpose of Function: Take the data stored in signal and write it to a WAV file
*						It is assumed the data is 16 bits per sample stereo
*						Stereo means 2 channels.
*						Only one channel is stored in the array
*						This one channel will be duplicated on both channels
*						The array contains floating point numbers between -1 
*						and +1.
* 						These will be changed to be 16 bit signed integer 
*						numbers for output to WAV file  
*
******************************************************************************/
void sendOutTo16BitWav(float* signal, const int sigSize){
	unsigned int buffer4;
	unsigned short int buffer2; 
	
 	FILE *outfile;
	 
	long numSamples;
	long sizeOfEachSample;	
	long i = 0;
	short int dataBuffer;
    
	printf("Opening output WAV file...\n");
	outfile = fopen("outputWav.wav", "wb");
	
	if(outfile == NULL){
		printf("Error opening output WAV file\n");
		system("pause");
		exit(1);		
	}


	printf("Signal Size to be output is %d Samples\n", sigSize);	
 	strcpy(header.riff, "RIFF");
	header.overallSize = sigSize*4 + 36;
	strcpy(header.wave, "WAVE");
	strcpy(header.fmtChunkMarker, "fmt ");
	header.lengthOfFmt = 16;
	header.formatType = 1;
	header.channels = 2;
	header.sampleRate = OUTPUT_NYQUIST*2;
	header.byterate = OUTPUT_NYQUIST*2*4;
	header.blockAlign = 4;
	header.bitsPerSample = 16;
	strcpy(header.dataChunkHeader,"data");
	header.dataSize = sigSize*2;
	
	
	/* writing WAV file header parts*/
 	fwrite((unsigned char*)&header.riff, sizeof(header.riff), 1, outfile);
 	fwrite((unsigned char*)&header.overallSize, sizeof(buffer4), 1, outfile);
 	fwrite((unsigned char*)&header.wave, sizeof(header.wave), 1, outfile);
 	fwrite((unsigned char*)&header.fmtChunkMarker, sizeof(header.fmtChunkMarker), 1, outfile);
 	fwrite((unsigned char*)&header.lengthOfFmt, sizeof(buffer4), 1, outfile);
 	fwrite((unsigned char*)&header.formatType, sizeof(buffer2), 1, outfile); 
	fwrite((unsigned char*)&header.channels, sizeof(buffer2), 1, outfile);
 	fwrite((unsigned char*)&header.sampleRate, sizeof(buffer4), 1, outfile);
 	fwrite((unsigned char*)&header.byterate, sizeof(buffer4), 1, outfile);
 	fwrite((unsigned char*)&header.blockAlign, sizeof(buffer2), 1, outfile);
 	fwrite((unsigned char*)&header.bitsPerSample, sizeof(buffer2), 1, outfile);
 	fwrite((unsigned char*)&header.dataChunkHeader, sizeof(header.dataChunkHeader), 1, outfile);
 	fwrite((unsigned char*)&header.dataSize, sizeof(buffer4), 1, outfile);
 	

	/*Writing out the data, in this case the same signal will be duplicated on both channels*/
	for (i =0; i < (sigSize); i++) {
			dataBuffer = (short int)(((*(signal+i)))*32768 );	
			fwrite((unsigned char*)(&dataBuffer), 2, 1, outfile); /*Data for Left Channel*/
			fwrite((unsigned char*)(&dataBuffer), 2, 1, outfile);/*Data for Right Channel*/	
	}
}

/******************************************************************************
*
* Function Name: filter
*
* Input Parameters: const float* sigIn is the input signal to be filtered
*	
* 					const float* aCoeff are the a coefficients of the filter
*
*					const float* bCoeff are the b coefficients of the filter
*
*					const sigSize is the length of the input and output signal
*							Note: that as this is an IIR filter the output 
*								signal could be theoretically infinite size
*								but it is truncated to the same size as the 
*								input signal.
*
* Output Parameters: float* sigOut is the signal after it is filtered
*
* Returns: void
*
* Purpose of Function: Filter the input signal with the coefficients provided.
*
******************************************************************************/
void filter(const float* sigIn, float* sigOut, const float* aCoeff, 
								const float* bCoeff, const int sigSize) { 
	int i;
	int j;
	
	/*********Filtering Algorithm************************************************/
	for(i=0;i<sigSize;i++){
		sigOut[i]=0;
		for(j=0; j< NUMCOEFF; j++){
			if((i-j)>=0){			
				sigOut[i] += sigIn[i-j]*aCoeff[j]+sigOut[i-j]*bCoeff[j];
			}
		}			
	}
	
	return; /*void*/
}



/******************************************************************************
*
* Function Name: dft
*
* Input Parameters: const float* signal is the input signal to by analysed
*	
*					const sigSize is the length of the input signal
*							Note: The lenght of the dft is as defined by DFTSIZE
*
* Output Parameters: float* mag is the frequency magnitude analysis of the 
*								input signal.
*
*					float* pha is the frequency phase analysis of the 
*								input signal.	
*
* Returns: void
*
* Purpose of Function: Perform the discrete fourier transform of the input
*						signal producing the Magnitude and Phase frequency 
*						analysis produced by the DFT.
*
******************************************************************************/
void dft(const float* signal, float* mag, float* pha, const int sigSize) {
	
	int i;
	int j;
	int k;
	
	const float PI = 3.14159;	
	float reX[DFTSIZE] = {0};
	float imX[DFTSIZE] = {0};
	
	for(k=0;k<DFTSIZE;k++) { 
        for(i=0;i<sigSize;i++) {                    
            reX[k]=reX[k]+signal[i]*cos(2*PI*k*i/sigSize);
            imX[k]=imX[k]-signal[i]*sin(2*PI*k*i/sigSize);       
        }       
    }
    
    /* This is where we convert from Rectangular to Polar form and print 
      the result out to a file*/ 
    for(k=0;k<DFTSIZE;k++) {
        mag[k]= sqrt((reX[k]*reX[k])+(imX[k]*imX[k]));
        
		/* for negligably small magnitudes, set the phase to 0*/
        if (mag[k] < 0.03){
            pha[k] = 0;
        }
        else {
            pha[k] = atan(imX[k]/reX[k]); /*polar nuicances not dealt with see 
										 * http://www.dspguide.com/ch8/9.htm */
        }
		/*Change Magnitude to Log range*/
		mag[k]=20*log10(mag[k]);
       
    }
    return; /*void*/
}
