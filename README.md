# DFT-and-Filtering-of-WAV-file

This program analyses and filters signals taken in from wav files or csv files.

It is currently set up for wav files but can be changed to csv input by modifying the code under the section 
/***********Take in the input signal**************************************/ in the filteringAnalysis.c file.

The filtering requires hard coding of filter coefficients into the file. It is currently set up as an IIR Low pass Chebyshev filter with a cutoff frequency of 0.1 of the sampling frequency.

The program will output a wav file with the filtered signal named output.wav.
The input in this program is a Piano C6 note in the PianoC.wav file (source:http://freewavesamples.com/alesis-fusion-bright-acoustic-piano-c6)
This must be in the same directory as the executable file (or c file if you are compiling and running).
You can change the input wav file in the takeInFrom16BitWav(float* signal, const int sigSize); function.

The program will also output an analysis csv file called outputAnalysis.csv.
You can open this in most spreadsheet programs such as Microsoft Excel and then graph the data.
Note if you open this file and then try to re-run the program while the file is open you will get an error. So make sure to close it (or save it as another name) before re-running the program.

The columns in the output analysis are as follows.
Column one (A) is the input signal (but only points up to DFTSIZE are shown).
Column two (B) is the output signal (but only points up to DFTSIZE are shown).
Column three (C) is a frequency scale which should be used as the x-axis in any graph of the frequency spectra or response.
Columns four (D) and five (E) are the Magnitude and Phase specturm respectively of the input signal.
Column six (F) is the Impulse Response (time domain and only points up to DFTSIZE are shown) of the filter.
Columns seven (G) and eight (H) are the Magnitude and Phase Frequency Response respectively of the filter.
Columns eight (I) and nine (J) are the Magnitude and Phase spectrum respectively of the output signal.
The DFT algorithm is as given in http://www.dspguide.com/ch8/6.htm.

The size of the DFT can be changed to show increased (or decreased) frequency resolution by changing the DFTSIZE value and SIGSEGSIZE.
SIGSEGSIZE should be DFTSIZE/2+1 although if you do not have enough signal you can still increase the DFTSIZE by padding the signal with zeros up to the required SIGSEGSIZE.
NOTE: This is a DFT not an FFT and as such increasing the size fo the DFT can have a dramatic effect on processing time even for a powerful computer.
Only the first SIGSEGSIZE samples are analysed in the program. If you wish to analyse other parts you will have to make that change yourself.

The entire signal (up to SIGSIZE) is filtered and output to the outputWav.wav file.
At present, you have to have some indication of the size of the input files and set this with SIGSIZE.
In future implementations, it is planned to detemine this at runtime and allocate all the necessary arrays with malloc.
